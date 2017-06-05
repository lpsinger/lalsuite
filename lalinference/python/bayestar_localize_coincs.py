#
# Copyright (C) 2013-2017  Leo Singer
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
"""
Produce GW sky maps for all coincidences in a LIGO-LW XML file.

The filename of the (optionally gzip-compressed) LIGO-LW XML input is an
optional argument; if omitted, input is read from stdin.

The distance prior is controlled by the --prior-distance-power argument.
If you set --prior-distance-power=k, then the distance prior is
proportional to r^k. The default is 2, uniform in volume.

If the --min-distance argument is omitted, it defaults to zero. If the
--max-distance argument is omitted, it defaults to the SNR=4 horizon
distance of the most sensitive detector.

A FITS file is created for each sky map, having a filename of the form

  "X.toa_phoa_snr.fits"
  "X.toa_snr_mcmc.fits"
  "X.toa_phoa_snr_mcmc.fits"

where X is the LIGO-LW row id of the coinc and "toa" or "toa_phoa_snr"
identifies whether the sky map accounts for times of arrival (TOA),
PHases on arrival (PHOA), and amplitudes on arrival (SNR).

OpenMP and MPI Parallelism
--------------------------

bayestar_localize_coincs supports two kinds of parallelism: OpenMP and MPI.
With OpenMP, BAYESTAR uses multiple cores on the current machine to speed up
numerically intensive calculations for each sky map. With MPI, BAYESTAR uses
multiple machines to evaluate many sky maps at once. It is possible to use
both OpenMP and MPI at the same time, but usually it is best to use only
OpenMP if you are using a single workstation or only MPI if you are submitting
jobs to a computing cluster.

The number of OpenMP threads is set by the OMP_NUM_THREADS environment
variable. If not set, then the default is the number of cores. Here is
an example of running BAYESTAR with OpenMP and 8 threads:

    $ export OMP_NUM_THREADS=8
    $ bayestar_localize_coincs coinc.xml

or equivalently:

    $ OMP_NUM_THREADS=8 bayestar_localize_coincs coinc.xml

Running BAYESTAR with MPI requires launching bayestar_localize_coincs
using mpiexec, which is a standard tool that is included with any MPI
installation. It also requires the latest development version of mpi4py
(pip install --user git+https://github.com/mpi4py/mpi4py). Here is an
example of running BAYESTAR using 8 processes and no OpenMP threads:

    $ OMP_NUM_THREADS=1 mpiexec -n 8 bayestar_localize_coincs coinc.xml

Of course you can use both OpenMP and MPI:

    $ OMP_NUM_THREADS=8 mpiexec -n 8 bayestar_localize_coincs coinc.xml

Submission to an HTCondor cluster
---------------------------------
Here is an example HTCondor submit file, bayestar.sub (doesn't work yet).

    universe = parallel
    getenv = true
    machine_count = 32
    request_memory = 1000 MB
    executable = /usr/bin/env
    arguments = "mpiexec -n $(machine_count) bayestar_localize_coincs coinc.xml"
    error = $(Node).err
    log = $(Node).log
    queue

Submit with:

    $ condor_submit bayestar.sub

Submission to a PBS cluster
---------------------------
Here is an example PBS submit script, bayestar.pbs, for NASA's Pleiades
supercomputer (https://www.nas.nasa.gov/hecc/resources/pleiades.html).

    #PBS -V -S /bin/sh -l select=4:ncpus=16:mpiprocs=16:ompthreads=1:model=san
    mpiexec bayestar_localize_coincs coinc.xml

Submit with:

    $ qsub bayestar.pbs
"""


# Command line interface.
import argparse
from lalinference.bayestar import command

methods = '''
    toa_phoa_snr
    toa_phoa_snr_mcmc
    '''.split()
default_method = 'toa_phoa_snr'
command.skymap_parser.add_argument(
    '--method', choices=methods, default=[default_method], nargs='*',
    help='Sky localization methods [default: %(default)s]')
parser = command.ArgumentParser(
    parents=[
        command.waveform_parser, command.prior_parser, command.skymap_parser])
parser.add_argument(
    'input', metavar='INPUT.xml[.gz]', default='-', nargs='+',
    type=argparse.FileType('rb'),
    help='Input LIGO-LW XML file [default: stdin] or PyCBC HDF5 files. '
         'For PyCBC, you must supply the coincidence file '
         '(e.g. "H1L1-HDFINJFIND.hdf" or "H1L1-STATMAP.hdf"), '
         'the template bank file (e.g. H1L1-BANK2HDF.hdf), '
         'the single-detector merged PSD files '
         '(e.g. "H1-MERGE_PSDS.hdf" and "L1-MERGE_PSDS.hdf"), '
         'and the single-detector merged trigger files '
         '(e.g. "H1-HDF_TRIGGER_MERGE.hdf" and "L1-HDF_TRIGGER_MERGE.hdf"), '
         'in any order.')
parser.add_argument(
    '--pycbc-sample', default='foreground',
    help='sample population [PyCBC only; default: %(default)s]')
parser.add_argument(
    '--output', '-o', default='.',
    help='output directory [default: current directory]')
opts = parser.parse_args()

#
# Logging
#

import logging
log = logging.getLogger('BAYESTAR')

# BAYESTAR imports.
from lalinference.io import fits, events
from lalinference.bayestar.sky_map import localize

# Other imports.
import os
import six
from functools import partial

# Squelch annoying and uniformative LAL log messages.
import lal
lal.ClobberDebugLevel(lal.LALNDEBUG)

# Read coinc file.
log.info('%s:reading input files', ','.join(file.name for file in opts.input))
event_source = events.open(*opts.input, sample=opts.pycbc_sample)

command.mkpath(opts.output)


def process(opts, int_coinc_event_id):
    event = event_source[int_coinc_event_id]

    coinc_event_id = 'coinc_event:coinc_event_id:{}'.format(int_coinc_event_id)

    # Loop over sky localization methods
    for method in opts.method:
        log.info("%s:method '%s':computing sky map", coinc_event_id, method)
        if opts.chain_dump:
            chain_dump = '%s.chain.npy' % int_coinc_event_id
        else:
            chain_dump = None
        sky_map = localize(
            event, opts.waveform, opts.f_low, opts.min_distance,
            opts.max_distance, opts.prior_distance_power, opts.cosmology,
            method=method, nside=opts.nside, chain_dump=chain_dump,
            enable_snr_series=opts.enable_snr_series,
            f_high_truncate=opts.f_high_truncate)
        sky_map.meta['objid'] = coinc_event_id
        log.info(
            "%s:method '%s':saving sky map",
            coinc_event_id, method)
        filename = '%d.%s.fits' % (int_coinc_event_id, method)
        fits.write_sky_map(
            os.path.join(opts.output, filename), sky_map, nest=True)


# Determine if we are running under MPI
try:
    from emcee.mpi_pool import MPIPool
    pool = MPIPool(loadbalance=True)
except (ImportError, ValueError):
    log.info('Not using MPI')
    map(partial(process, opts), six.iterkeys(event_source))
else:
    try:
        log.info('Using %d MPI processes', pool.size)
        pool.map(partial(process, opts), event_source.keys())
    finally:
        pool.close()
