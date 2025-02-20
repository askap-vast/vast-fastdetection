#!/usr/bin/env python
"""
Copyright (C) Swinburne 2025

Given an input measurement sets and config files, produce output measurement sets 
basic info and print out the "intervels-out" parameter - i.e. the number of short images
based on TIMESTEP and total duration of the observation 
"""
import os
import sys
import yaml
import argparse
import numpy as np
from collections import Counter
from casacore.tables import table

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>"


def _main():

    parser = argparse.ArgumentParser(
        prog='VOevent', 
        description='VO Event trigger', 
        epilog='Example usage: python ~/scripts/notebooks/notes/template.py -h', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        )
    parser.add_argument('msfile', type=str, help='path for the ms file')
    parser.add_argument('--config', type=str, default=None, help='path for the config file')
    parser.add_argument('--savename', type=str, default='measurements.txt', help='path to save measurementsets info')
    args = parser.parse_args()

    # ===============
    # main program 
    # ===============

    if args.config is not None:
        config = read_config(args.config)
        timestep = config['TIMESTEP']
        intervals_out = read_ms(args, timestep=timestep)
    else:
        intervals_out = read_ms(args, timestep='int')
    
    print(intervals_out)


    
def read_ms(args, timestep='int'):
    '''
    timestep: str "int" (integration time) or any float (unit of seconds) 
    '''
    tb = table(args.msfile)
    # observation time list 
    times = np.array(list(Counter(tb.getcol('TIME')).keys()))
    # observation duration 
    duration = times.max() - times.min()
    # sampling/integration time 
    interval = np.unique(tb.getcol('INTERVAL'))[0]

    tb.close()
    
    if timestep == 'int':
        intervals_out = times.shape[0]
    else:
        intervals_out = int(round(duration / float(timestep)))

    with open(args.savename, 'w') as fw:
        fw.write('FILENAME:          ' + args.msfile + '\n')
        fw.write('INTEGRATION TIME:  ' + str(interval) + ' seconds' + '\n')
        fw.write('DURATION:          ' + str(duration) + ' seconds' + '\n')
        fw.write('                   ' + str(duration/60) + ' minutes' + '\n') 
        fw.write('                   ' + str(duration/3600) + ' hours' + '\n') 
        fw.write('NUM INTEGRATIONS:  ' + str(times.shape[0]) + '\n' )
        fw.write('TIMESTEP:          ' + str(timestep) + ' seconds' + '\n')
        fw.write('INTERVALS OUT:     ' + str(intervals_out) + '\n')
        if timestep == 'int':
            fw.write('TIMESTEP (actual): ' + str(interval) + ' seconds' + '\n')
        else:
            fw.write('TIMESTEP (actual): ' + str(duration/intervals_out) + ' seconds' + '\n')

    return intervals_out


def read_config(fname):
    # read configuration file 
    if os.path.isfile(fname):
        with open(fname, 'r') as yaml_file:
            config = yaml.safe_load(yaml_file)
        return config
    
    else:
        print('Cannot find configuration file')
        print('%s does not exists' % fname)
        sys.exit()


if __name__ == '__main__':
    _main()
