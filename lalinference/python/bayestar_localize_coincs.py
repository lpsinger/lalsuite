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
"""
from __future__ import division
from __future__ import print_function


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
    '--keep-going', '-k', default=False, action='store_true',
    help='Keep processing events if a sky map fails to converge [default: no]')
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
group = parser.add_argument_group(
    'parallelization options', 'options for parallel or remote execution')
submit_arg = group.add_argument(
    '--submit', default='openmp', const='auto', nargs='?',
    choices=('auto', 'openmp', 'condor', 'pbs'),
    help='parallelize locally or submit to a cluster using a job scheduler')
group.add_argument(
    '-i', '--i', type=int,
    help='job number (internal, used by --submit=condor or --submit=pbs)')
group.add_argument(
    '-j', '--jobs', nargs='?', default=1, type=int,
    help='number of jobs/processes')
opts = parser.parse_args()

#
# Logging
#

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('BAYESTAR')

# BAYESTAR imports.
from lalinference.io import fits, events
from lalinference.bayestar.sky_map import (
    localize, get_num_threads, set_num_threads)

# Other imports.
import os
from distutils.spawn import find_executable
import sys
import tempfile

# Read coinc file.
log.info('%s:reading input files', ','.join(file.name for file in opts.input))
event_source = events.open(*opts.input, sample=opts.pycbc_sample)

# Decide on the type of parallelism.
if opts.submit == 'auto':
    if find_executable('qsub'):
        opts.submit = 'pbs'
    elif find_executable('condor'):
        opts.submit = 'condor'
if opts.submit in ('auto', 'openmp'):
    try:
        if opts.jobs is None:
            opts.jobs = get_num_threads()
        else:
            set_num_threads(opts.jobs)
    except NotImplementedError:
        if opts.submit == 'openmp':
            raise NotImplementedError(
                'BAYESTAR was not built with OpenMP support')
        else:
            opts.submit = None
    else:
        opts.submit = 'openmp'

# Prepare basic arguments for job submission
if opts.submit in ('condor', 'pbs'):
    if opts.i is not None:
        raise ValueError(
            'must not set -i with --submit=condor or --submit=pbs')

    clean_args = list(sys.argv)
    for choice in submit_arg.choices:
        try:
            clean_args.remove('--submit={}'.format(choice))
        except ValueError:
            pass
        try:
            i = clean_args.index('--submit')
        except ValueError:
            pass
        else:
            try:
                next_arg = clean_args[i + 1]
            except IndexError:
                j = i + 1
            else:
                j = i + 2
            del clean_args[i:j]

    if opts.jobs is None or opts.jobs == 1:
        raise RuntimeError(
            '"--jobs j" with j > 1 is required with "--submit=condor" or '
            '"--submit=pbs"')

    if opts.submit == 'condor':
        clean_args += ['-i', '$(Process)']
        cmd = ['condor_submit', 'accounting_group=ligo.dev.o3.cbc.pe.bayestar',
               'on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)',
               'on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)',
               'on_exit_hold_reason = (ExitBySignal == True ? '
               'strcat("The job exited with signal ", ExitSignal) : '
               'strcat("The job exited with signal ", ExitCode))',
               'request_memory = 1000 MB', 'universe=vanilla', 'getenv=true',
               'executable=' + sys.executable, 'JobBatchName=BAYESTAR',
               'error=' + os.path.join(opts.output, '$(Cluster).$(Process).err'),
               'log=' + os.path.join(opts.output, '$(Cluster).$(Process).log'),
               'arguments="-B ' + ' '.join(arg for arg in clean_args) + '""',
               '-append', 'queue {}'.format(len(event_source)), '/dev/null']
        os.execvp('condor_submit', cmd)
    elif opts.submit == 'pbs':
        clean_args += ['-i', '${PBS_ARRAY_INDEX}']
        with tempfile.NamedTemporaryFile(mode='w', prefix='bayestar',
                                         suffix='.pbs', dir='.',
                                         delete=False) as f:
            filename = f.name
            print('#PBS', '-l', 'select=1:ncpus=1:model=san', file=f)
            print('#PBS', '-S', '/bin/sh', file=f)
            print('#PBS', '-V', file=f)
            print('#PBS', '-J', '0-{}'.format(opts.jobs - 1), file=f)
            print(sys.executable, '-B', *clean_args, file=f)
        os.execvp('qsub', ['qsub', filename])

command.mkpath(opts.output)
count_sky_maps_failed = 0

if opts.i is not None:
    from lalinference.io.events.base import EventSource

    class RangeEventSource(EventSource):
        def __init__(self, source, i, j):
            self._parent = source
            keys = list(source.keys())
            n = len(keys)
            m = -(-n // j)
            self._keys = keys[(i * m):((i + 1) * m)]

        def __len__(self):
            return len(self._keys)

        def __iter__(self):
            return iter(self._keys)

        def __getitem__(self, key):
            if key not in self._keys:
                raise KeyError
            return self._parent[key]

    event_source = RangeEventSource(event_source, opts.i, opts.jobs)

for int_coinc_event_id, event in event_source.items():
    coinc_event_id = 'coinc_event:coinc_event_id:{}'.format(int_coinc_event_id)

    # Loop over sky localization methods
    for method in opts.method:
        log.info("%s:method '%s':computing sky map", coinc_event_id, method)
        if opts.chain_dump:
            chain_dump = '%s.chain.npy' % int_coinc_event_id
        else:
            chain_dump = None
        try:
            sky_map = localize(
                event, opts.waveform, opts.f_low, opts.min_distance,
                opts.max_distance, opts.prior_distance_power, opts.cosmology,
                method=method, nside=opts.nside, chain_dump=chain_dump,
                enable_snr_series=opts.enable_snr_series,
                f_high_truncate=opts.f_high_truncate)
            sky_map.meta['objid'] = coinc_event_id
        except (ArithmeticError, ValueError):
            log.exception(
                "%s:method '%s':sky localization failed",
                coinc_event_id, method)
            count_sky_maps_failed += 1
            if not opts.keep_going:
                raise
        else:
            log.info(
                "%s:method '%s':saving sky map",
                coinc_event_id, method)
            filename = '%d.%s.fits' % (int_coinc_event_id, method)
            fits.write_sky_map(
                os.path.join(opts.output, filename), sky_map, nest=True)

if count_sky_maps_failed > 0:
    raise RuntimeError("{0} sky map{1} did not converge".format(
        count_sky_maps_failed, 's' if count_sky_maps_failed > 1 else ''))
