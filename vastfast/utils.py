import logging

import os
import glob
import tarfile

logger = logging.getLogger(__name__)

def remove_fits(dir):
    """remove fits files"""
    fits_files = glob.glob(dir+"/*.fits")
    if len(fits_files) > 0:
        for item in fits_files:
            print(item)
            os.remove(item)

def tar_file(tar_name, dir):
    """tar files"""
    with tarfile.open(tar_name, "w:gz") as tar:
        tar.add(dir, arcname=os.path.basename(dir))
