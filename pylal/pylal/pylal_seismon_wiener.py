#!/usr/bin/python

import os
import numpy as np
import scipy.linalg

import pylal.pylal_seismon_utils

import gwpy.time, gwpy.timeseries
import gwpy.spectrum, gwpy.spectrogram
import gwpy.plotter

__author__ = "Michael Coughlin <michael.coughlin@ligo.org>"
__date__ = "2012/8/26"
__version__ = "0.1"

# =============================================================================
#
#                               DEFINITIONS
#
# =============================================================================

def wiener(params, target_channel, segment):
    """@calculates wiener filter for given channel and segment.

    @param params
        seismon params dictionary
    @param target_channel
        seismon channel structure
    @param segment
        [start,end] gps
    """

    ifo = pylal.pylal_seismon_utils.getIfo(params)

    gpsStart = segment[0]
    gpsEnd = segment[1]

    # set the times
    duration = np.ceil(gpsEnd-gpsStart)

    dataAll = []

    for channel in params["channels"]:
        # make timeseries
        dataFull = pylal.pylal_seismon_utils.retrieve_timeseries(params, channel, segment)
        if dataFull == []:
            continue

        dataFull = dataFull / channel.calibration
        indexes = np.where(np.isnan(dataFull.data))[0]
        meanSamples = np.mean(np.ma.masked_array(dataFull.data,np.isnan(dataFull.data)))
        for index in indexes:
            dataFull[index] = meanSamples
        dataFull -= np.mean(dataFull.data)

        if np.mean(dataFull.data) == 0.0:
            print "data only zeroes... continuing\n"
            continue
        if len(dataFull.data) < 2*channel.samplef:
            print "timeseries too short for analysis... continuing\n"
            continue
    
        dataAll.append(dataFull)

    X = []
    y = []
    for dataFull in dataAll:
        if dataFull.channel.name == target_channel.station:
            tt = np.array(dataFull.times)
            y = dataFull.data
        else:
            if X == []:
                X = dataFull.data
            else:
                try:
                    X = np.vstack([X,dataFull.data])
                except:
                    continue

    if len(y) == 0:
        print "No data for target channel... continuing"
        return

    originalASD = []
    residualASD = []
    FFASD = []

    N = 1000
    gpss = np.arange(gpsStart,gpsEnd,2*params["fftDuration"])
    create_filter = True
    for i in xrange(len(gpss)-1):
        tt = np.array(dataFull.times)
        indexes = np.intersect1d(np.where(tt >= gpss[i])[0],np.where(tt <= gpss[i+1])[0])

        if len(indexes) == 0:
            continue

        indexMin = np.min(indexes)
        indexMax = np.max(indexes)

        yCut = y[indexMin:indexMax]
        XCut = X[:,indexMin:indexMax]

        XCut = XCut.T
        if create_filter:
            W,R,P = miso_firwiener(N,XCut,yCut)
            create_filter = False
            continue
            
        residual, FF = subtractFF(W,XCut,yCut,target_channel)

        gpsStart = tt[indexMin]
        dataOriginal = gwpy.timeseries.TimeSeries(yCut, epoch=gpsStart, sample_rate=target_channel.samplef)
        dataResidual = gwpy.timeseries.TimeSeries(residual, epoch=gpsStart, sample_rate=target_channel.samplef)
        dataFF = gwpy.timeseries.TimeSeries(FF, epoch=gpsStart, sample_rate=target_channel.samplef)

        #cutoff = 1.0
        #dataOriginal = dataOriginal.lowpass(cutoff, amplitude=0.9, order=3, method='scipy')
        #dataResidual = dataResidual.lowpass(cutoff, amplitude=0.9, order=3, method='scipy')
        #dataFF = dataFF.lowpass(cutoff, amplitude=0.9, order=3, method='scipy')

        # calculate spectrum
        NFFT = params["fftDuration"]
        #window = None
        dataOriginalASD = dataOriginal.asd(NFFT,NFFT,'welch')
        dataResidualASD = dataResidual.asd(NFFT,NFFT,'welch')
        dataFFASD = dataFF.asd(NFFT,NFFT,'welch')

        freq = np.array(dataOriginalASD.frequencies)
        indexes = np.where((freq >= params["fmin"]) & (freq <= params["fmax"]))[0]
        freq = freq[indexes]

        dataOriginalASD = np.array(dataOriginalASD.data)
        dataOriginalASD = dataOriginalASD[indexes]
        dataOriginalASD = gwpy.spectrum.Spectrum(dataOriginalASD, f0=np.min(freq), df=(freq[1]-freq[0]))

        dataResidualASD = np.array(dataResidualASD.data)
        dataResidualASD = dataResidualASD[indexes]
        dataResidualASD = gwpy.spectrum.Spectrum(dataResidualASD, f0=np.min(freq), df=(freq[1]-freq[0]))

        dataFFASD = np.array(dataFFASD.data)
        dataFFASD = dataFFASD[indexes]
        dataFFASD = gwpy.spectrum.Spectrum(dataFFASD, f0=np.min(freq), df=(freq[1]-freq[0]))

        originalASD.append(dataOriginalASD)
        residualASD.append(dataResidualASD)
        FFASD.append(dataFFASD)

    dt = gpss[1] - gpss[0]
    epoch = gwpy.time.Time(gpss[0], format='gps')
    originalSpecgram = gwpy.spectrogram.Spectrogram.from_spectra(*originalASD, dt=dt,epoch=epoch)
    residualSpecgram = gwpy.spectrogram.Spectrogram.from_spectra(*residualASD, dt=dt,epoch=epoch)
    FFSpecgram = gwpy.spectrogram.Spectrogram.from_spectra(*FFASD, dt=dt,epoch=epoch)

    freq = np.array(originalSpecgram.frequencies)
    bins,originalSpecvar = pylal.pylal_seismon_utils.spectral_histogram(originalSpecgram)
    original_spectral_variation_50per = pylal.pylal_seismon_utils.spectral_percentiles(originalSpecvar,bins,50)
    bins,residualSpecvar = pylal.pylal_seismon_utils.spectral_histogram(residualSpecgram)
    residual_spectral_variation_50per = pylal.pylal_seismon_utils.spectral_percentiles(residualSpecvar,bins,50)
    bins,FFSpecvar = pylal.pylal_seismon_utils.spectral_histogram(FFSpecgram)
    FF_spectral_variation_50per = pylal.pylal_seismon_utils.spectral_percentiles(FFSpecvar,bins,50)

    if params["doPlots"]:

        plotDirectory = params["path"] + "/Wiener/" + target_channel.station_underscore
        pylal.pylal_seismon_utils.mkdir(plotDirectory)

        fl, low, fh, high = pylal.pylal_seismon_NLNM.NLNM(2)

        pngFile = os.path.join(plotDirectory,"psd.png")

        plot = gwpy.plotter.Plot(figsize=[14,8])
        kwargs = {"linestyle":"-","color":"b","label":"Original"}
        plot.add_line(freq, original_spectral_variation_50per, **kwargs)
        kwargs = {"linestyle":"-","color":"g","label":"Residual"}
        plot.add_line(freq, residual_spectral_variation_50per, **kwargs)
        kwargs = {"linestyle":"-.","color":"k"}
        plot.add_line(fl, low, **kwargs)
        plot.add_line(fh, high, **kwargs)
        plot.axes.set_yscale("log")
        plot.axes.set_xscale("log")
        plot.xlim = [params["fmin"],params["fmax"]]
        #plot.ylim = [np.min(bins), np.max(bins)]
        plot.xlabel = "Frequency [Hz]"
        plot.ylabel = "Amplitude Spectrum [(m/s)/rtHz]"
        plot.add_legend(loc=1,prop={'size':10})
        plot.save(pngFile,dpi=200)
        plot.close()

def xcorr(x, y, normed=True, maxlags=None):
    """
    Call signature::

    xcorr(x, y, normed=True, detrend=mlab.detrend_none,
    usevlines=True, maxlags=10, **kwargs)

    """

    Nx = len(x)
    if Nx!=len(y):
        raise ValueError('x and y must be equal length')

    c = np.correlate(x, y, mode='full')

    if normed: c/= np.sqrt(np.dot(x,x) * np.dot(y,y))

    if maxlags is None: maxlags = Nx - 1

    if maxlags >= Nx or maxlags < 1:
        raise ValueError('maglags must be None or strictly positive < %d'%Nx)

    lags = np.arange(-maxlags,maxlags+1)
    c = c[Nx-1-maxlags:Nx+maxlags]

    return c,lags

def miso_firwiener(N,X,y):

    # MISO_FIRWIENER Optimal FIR Wiener filter for multiple inputs.
    # MISO_FIRWIENER(N,X,Y) computes the optimal FIR Wiener filter of order
    # N, given any number of (stationary) random input signals as the columns
    # of matrix X, and one output signal in column vector Y.
    # Author: Keenan Pepper
    # Last modified: 2007/08/02
    # References:
    # [1] Y. Huang, J. Benesty, and J. Chen, Acoustic MIMO Signal
    # Processing, SpringerVerlag, 2006, page 48

    # Number of input channels.
    try:
        junk, M = X.shape
    except:
        M = 1

    # Input covariance matrix.
    R = np.zeros([M*(N+1),M*(N+1)])
    for m in xrange(M):
        for i in xrange(m,M):
            rmi,lags = xcorr(X[:,m],X[:,i],maxlags=N,normed=False)
            Rmi = scipy.linalg.toeplitz(np.flipud(rmi[range(N+1)]),r=rmi[range(N,2*N+1)])
            top = m*(N+1)
            bottom = (m+1)*(N+1)
            left = i*(N+1)
            right = (i+1)*(N+1)
            #R[range(top,bottom),range(left,right)] = Rmi

            for j in xrange(top,bottom):
                for k in xrange(left,right):
                    R[j,k] = Rmi[j-top,k-left]
         
            if not i == m:
                #R[range(left,right),range(top,bottom)] = Rmi  # Take advantage of hermiticity.

                #RmiT = Rmi.T
                for j in xrange(left,right):
                    for k in xrange(top,bottom):
                        R[j,k] = Rmi[j-left,k-top]

    # Crosscorrelation vector.
    P = np.zeros([M*(N+1),])
    for i in xrange(M):
        top = i*(N+1)
        bottom = (i+1)*(N+1)
        p, lags = xcorr(y,X[:,i],maxlags=N,normed=False)

        P[range(top,bottom)] = p[range(N,2*N+1)]

    # The following step is very inefficient because it fails to exploit the
    # block Toeplitz structure of R. Its done the same way in the builtin
    # function "firwiener".
    # P / R
    #W = 

    W = np.linalg.lstsq(R.T, P.T)[0].T
    W = np.reshape(W,(N+1,M))

    return W,R,P

def subtractFF(W,SS,S,channel):

    # Subtracts the filtered SS from S using FIR filter coefficients W.
    # Routine written by Jan Harms. Routine modified by Michael Coughlin.
    # Modified: August 17, 2012
    # Contact: michael.coughlin@ligo.org

    N = len(W)-1
    ns = len(S)

    FF = np.zeros([ns-N,])

    for k in xrange(N+1,ns):
        tmp = SS[range(k-N,k+1),:] * W
        FF[k-N] = np.sum(tmp)

    cutoff = 1.0
    dataFF = gwpy.timeseries.TimeSeries(FF, sample_rate=channel.samplef)
    dataFFLowpass = dataFF.lowpass(cutoff, amplitude=0.9, order=3, method='scipy')
    FF = np.array(dataFFLowpass)

    residual = S[range(ns-N)]-FF

    return residual, FF