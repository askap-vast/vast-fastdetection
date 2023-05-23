#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 20:41:50 2023
@author: ywan3191

Get ready for various scripts. 

For each beam, we want: 
1. Query CASDA to get visibility download url + selavy catalogue download url. Generate **bash_GETDATA_beam??.sh** and run the script to download visibilities and catalogues. (Double check the completion of visibilities downloading and untar visibilities) 
2. Run slurm_FIXDATA_beam??.sh (including rescale and fix dir)
3. Run slurm_MODDEEP_beam??.sh (deep modeling)
4. Run slurm_IMAFAST_beam??.sh (fast imaging, in either time interval)
5. Run slurm_SELCAND_beam??.sh (select candidates)
"""





