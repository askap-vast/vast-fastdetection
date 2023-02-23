# default prefix for output files
import configparser

OUT_PREFIX = "output"

# list of maps 
KTYPELIST = ['chisquare', 'peak', 'std', 'gaussian']


# Gaussian width
G_WIDTH = 4

# chunksize used in dask.array when generating Gaussian map
# chunksize = (CHUNK_0, CHUNK_1, time_dim), where time_dim (the number of time slices) is determined at run time
# the default value is (100, 100, time_dim), which has been verified during the course of devlopment at ADACS
# NOTICE: if the chunksize is too big, it will result in out-of-memory issue or kill the process
CHUNK_0 = 100
CHUNK_1 = 100


try:
    config = configparser.ConfigParser()
    config.read("/input/app_settings.ini")
except Exception as e:
    print(e)
else:
    try:
        KTYPELIST = config["RUN_SETTINGS"]["KTYPELIST"]
    except:
        pass

import logging
logger = logging.getLogger()

logger.info("============")
logger.info("KTYPELIST")
print(KTYPELIST)
