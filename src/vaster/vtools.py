import os
import logging
logger = logging.getLogger(__name__)

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>"



def measure_running_time(start_time, end_time, nround=2):
    total_time = end_time - start_time
    if total_time <= 60:
        logger.info('Total running time %s seconds', round(total_time, nround))
    elif total_time <= 60*60:
        logger.info('Total running time %s minutes', round(total_time/60, nround))
    elif total_time <= 60*60*24:
        logger.info('Total running time %s hours', round(total_time/60/60, nround))
    else:
        logger.info('Total running time %s days', round(total_time/60/60/24, nround))


def process_txt(args, txt):
    if args.dry_run:
        logger.warning('Dry run: skip "%s"', txt)
    else:
        logger.info('Executing "%s"', txt)
        exit_code = os.system(txt)
        logger.debug('Exit with %s: %s', exit_code, txt)
