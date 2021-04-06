#!/usr/bin/env python3

"""
Download TENDL ENDF 2019 data from official site
and convert it to a HDF5 library for use with OpenMC.
"""

import argparse
import tarfile
from multiprocessing import Pool
from pathlib import Path
from shutil import rmtree
from urllib.parse import urljoin

import openmc.data
from utils import download, process_neutron


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=CustomFormatter
)
parser.add_argument('-d', '--destination', type=Path, default=None,
                    help='Directory to create new library in')
parser.add_argument('--download', action='store_true',
                    help='Download files from PSI')
parser.add_argument('--no-download', dest='download', action='store_false',
                    help='Do not download files from PSI')
parser.add_argument('--extract', action='store_true',
                    help='Extract tar files')
parser.add_argument('--no-extract', dest='extract', action='store_false',
                    help='Do not extract tar files')
parser.add_argument('--libver', choices=['earliest', 'latest'],
                    default='latest', help="Output HDF5 versioning. Use "
                    "'earliest' for backwards compatibility or 'latest' for "
                    "performance")
parser.add_argument('-r', '--release', choices=['2019'],
                    default='2019')
parser.add_argument('--cleanup', action='store_true',
                    help="Remove download directories when data has "
                    "been processed")
parser.add_argument('--no-cleanup', dest='cleanup', action='store_false',
                    help="Do not remove download directories when data has "
                    "been processed")
parser.add_argument('--temperatures', type=float,
                    default=[250.0, 293.6, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 
                             900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 2500.0],
                    help="Temperatures in Kelvin", nargs='+')
parser.set_defaults(download=True, extract=True, cleanup=False)
args = parser.parse_args()


library_name = 'tendl'

cwd = Path.cwd()

endf_files_dir = cwd.joinpath('-'.join([library_name, args.release, 'endf']))
download_path = cwd.joinpath('-'.join([library_name, args.release, 'download']))
# the destination is decided after the release is known
# to avoid putting the release in a folder with a misleading name
if args.destination is None:
    args.destination = Path('-'.join([library_name, args.release, 'hdf5']))

# This dictionary contains all the unique information about each release.
# This can be extended to accommodate new releases
release_details = {
    '2019': {
        'base_url': 'https://tendl.web.psi.ch/tendl_2019/tar_files/',
        'compressed_files': ['TENDL-n.tgz'],
        'neutron_files': endf_files_dir.glob('n-*.tendl'),
        'compressed_file_size': '2.78 GB',
        'uncompressed_file_size': '12 GB'
    }
}

download_warning = """
WARNING: This script will download {} of data.
Extracting and processing the data requires {} of additional free disk space.
""".format(release_details[args.release]['compressed_file_size'],
           release_details[args.release]['uncompressed_file_size'])

# ==============================================================================
# DOWNLOAD FILES FROM WEBSITE

if args.download:
    print(download_warning)
    for f in release_details[args.release]['compressed_files']:
        # Establish connection to URL
        download(urljoin(release_details[args.release]['base_url'], f),
                 output_path=download_path)


# ==============================================================================
# EXTRACT FILES FROM TGZ
if args.extract:
    for f in release_details[args.release]['compressed_files']:
        with tarfile.open(download_path / f) as zf:
            print('Extracting {0}...'.format(f))
            zf.extractall(path=endf_files_dir)

    if args.cleanup and download_path.exists():
        rmtree(download_path)

# ==============================================================================
# GENERATE HDF5 LIBRARY -- NEUTRON FILES

# Get a list of all ENDF files
neutron_files = release_details[args.release]['neutron_files']

# Create output directory if it doesn't exist
args.destination.mkdir(parents=True, exist_ok=True)

library = openmc.data.DataLibrary()

with Pool() as pool:
    results = []
    for filename in sorted(neutron_files):

        func_args = (filename, args.destination, args.libver, args.temperatures)
        r = pool.apply_async(process_neutron, func_args)
        results.append(r)

    for r in results:
        r.wait()

# Register with library
for p in sorted((args.destination).glob('*.h5')):
    library.register_file(p)

# Write cross_sections.xml
library.export_to_xml(args.destination / 'cross_sections.xml')
