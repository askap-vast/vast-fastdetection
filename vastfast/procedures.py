import logging

import os
import glob
import tarfile
import  shutil

from .cube import Cube, Filter
from .plot import Candidates, Products, combine_csv
from .exceptions import *
from .setting import OUT_PREFIX, KTYPELIST, G_WIDTH

logger = logging.getLogger(__name__)

def create_out_beam(outdir, beam):
    outdir_beam = outdir + "/output_beam{:02}".format(beam)
    if not os.path.exists(outdir_beam):
        os.mkdir(outdir_beam)
    return outdir_beam

class Procedures():
    def __init__(self, imagelist, beam, catalogue, outdir, outdir_beam, deepimage=None, nprocess=4):
        self.imagelist = imagelist
        self.beam = beam
        self.catalogue = catalogue
        self.outdir = outdir
        self.outdir_beam = outdir_beam
        self.deepimage = deepimage
        self.nprocess = nprocess
        self._create_outfile_names()

    def _create_outfile_names(self):
        # output file prefix
        self.name = OUT_PREFIX + "_beam{:02}".format(self.beam)
        # output tar file
        self.tar_name = self.outdir + "/output_beam{:02}.tar.gz".format(self.beam)

    def run_all_steps(self):
        sigcube = self.get_sigcube()
        f = self.get_filter(sigcube)
        self.get_map_all(f)
        self.read_map_list()
        self.get_candidates()
        self.combine_catalogue()
        self.plot_final_candidates()
        self.finalize_output()        

    def _get_sigcube(self):
        """get sigcube"""
        logger.info("============")
        logger.info("Building the cube...")
        cube = Cube(self.imagelist)
        cube.icube()
        logger.info("Finish buidling the cube.")
        logger.info("============")
        logger.info("Removing bad images...")
        cube.remove_bad_images()
        logger.info("cube shape: {}".format(cube.sigcube.shape))
        return cube.sigcube

    def get_sigcube(self):
        try:
            sigcube = self._get_sigcube()
            return sigcube

        except Exception:
            logger.exception("Fail to generate sigcube; skip the beam{:02}".format(self.beam))
            return -1

    def get_filter(self, sigcube):
        try:
            f = Filter(sigcube)
            return f
        except Exception:
            logger.exception("Fail to generate rms cube; skip the beam{:02}".format(self.beam))
            return -1

    def _get_map(self, f, ktype):
        logger.info("===== Matched Filter =====")
        logger.info("Kernel match filter '{}'...".format(ktype))
        f.fmap(ktype, width=G_WIDTH, nprocess=self.nprocess)
        logger.info("Kernel match Done")
        f.tofits(fitsname="{}/{}_{}.fits".format(self.outdir_beam, self.name, ktype), imagename=self.imagelist[0])
        logger.info("Save the results to {}_{}.fits".format(self.name, ktype))

    def get_map_all(self, f):
        """get maps of given ktype"""
        for ktype in KTYPELIST:
            try:
                self._get_map(f, ktype)
            except Exception:
                logger.exception("Fail to generate {} map".format(ktype))

    def read_map_list(self):
        maps = {}
        for ktype in KTYPELIST:
            map_name = ktype+ "_map"
            maps[map_name] = "{}/{}_{}.fits".format(self.outdir_beam, self.name, ktype)

        self.map_dict = maps

    def _get_candidates_single_map(self, maptype):
        c = Candidates(self.map_dict)

        # find local maximum
        logger.info("Finding local maximum for {}...".format(maptype))
        c.local_max(maptype=maptype, min_distance=30, sigma=5)
        logger.info("Finish finding local maximum. ")
        
        
        # plot a map with all of candidates above the threshold 
        c.plot_fits(fitsname=self.map_dict[maptype], 
                imagename="{}/{}_{}1".format(self.outdir_beam, self.name, maptype))
        
        logger.info("Deep image catalogue {}".format(self.catalogue))
        c.select_candidates(deepcatalogue=self.catalogue)
        
        # save the table
        c.save_csvtable(tablename="{}/{}_{}_cand".format(self.outdir_beam, self.name, maptype), savevot=True)
        logger.info("Save candidates in table.")

        # plot a final map with promising candidates 
        c.plot_fits(fitsname=self.map_dict[maptype], 
                imagename="{}/{}_{}2".format(self.outdir_beam, self.name, maptype))

    def get_candidates(self):
        self.valid_map = Candidates(self.map_dict).valid_map
        for maptype in self.valid_map:
            self._get_candidates_single_map(maptype)
            
    def combine_catalogue(self):
        logger.info("=========Combine catalogue==========")
        namelist = ['{}/{}_{}_cand.csv'.format(self.outdir_beam, self.name, maptype) 
                    for maptype in self.valid_map ]
        combine_csv(namelist, tablename="{}/{}_final".format(self.outdir_beam, self.name), 
                        savevot=True)

    def plot_final_candidates(self):
        logger.info("========= Plotting =============")
        final_csv = "{}/{}_final.csv".format(self.outdir_beam, self.name)

        if os.path.exists(final_csv):
            p = Products(final_csv)
            p.generate_slices(imagelist=self.imagelist, 
                            savename='{}/{}_slices'.format(self.outdir_beam, self.name))

            if self.deepimage is None:
                logger.warning("No deepcutout or lightcurve files are generated without deep image.")
            else:
                p.generate_cutout(fitsname=self.deepimage, 
                                savename='{}/{}_deepcutout'.format(self.outdir_beam, self.name))
                p.generate_lightcurve(imagelist=self.imagelist, 
                                    deepname=self.deepimage, 
                                    savename='{}/{}_lightcurve'.format(self.outdir_beam, self.name))
    
    def _remove_fits(self):
        """remove existing fits files"""
        fits_files = glob.glob(self.outdir_beam+"/*.fits")
        if len(fits_files) > 0:
            for item in fits_files:
                os.remove(item)
                logger.info("Remove file {}.".format(item))

    def _tar_file(self):
        """create tar file"""
        nfile = len(glob.glob(self.outdir_beam+"/*"))
        if nfile == 0:
            logger.warning("No output files available for beam{:02}".format(self.beam))
        else:
            with tarfile.open(self.tar_name, "w:gz") as tar:
                tar.add(self.outdir_beam, arcname=os.path.basename(self.outdir_beam))

    def finalize_output(self):
        logger.info("Removing fits files...")
        self._remove_fits()

        logger.info("Creating tarfile...")
        self._tar_file()

        logger.info("Cleaning up output directory...")
        shutil.rmtree(self.outdir_beam)


 # # plot!
        # for i, candname in enumerate(c.cand_name):
        #     logger.info("Plot slices {}/{}: {}".format(i, len(c.cand_name), candname))
        #     plot.plot_slices(src_name=candname, 
        #                      imagelist=imagelist, 
        #                      name="{}/{}_{}".format(outdir_beam, name, candname))
        
