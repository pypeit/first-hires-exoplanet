#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 10:18:22 2021

@author: pheeren
"""

# Import packages
import pyodine
from pyodine.timeseries.misc import robust_mean, robust_std

import os
import matplotlib.pyplot as plt
import time
import logging

import argparse
import importlib


def combine_velocity_results(Pars, res_files=None, comb_res_in=None, 
                             plot_dir=None, comb_res_out=None, vels_out=None, 
                             reject_files=None, bary_dict=None, temp_vel=None,
                             ref_vel=None, error_log=None, info_log=None,
                             compact=False, quiet=False):
    """Weight and combine chunk velocities from modelling results
    
    :param Pars: The parameters to use in the routine.
    :type Pars: :class:`Timeseries_Parameters`
    :param res_files: A pathname to a text-file with pathnames of individual 
        results to load for the combination, or, alternatively, a list of 
        pathnames to individual results. If this is None, hand an existing 
        saved :class:`pyodine.timeseries.base.CombinedResults` object to 
        'comb_res_in'!
    :type res_files: str, list, tuple, or None
    :param comb_res_in: A pathname to a saved 
        :class:`pyodine.timeseries.base.CombinedResults` object to load. 
        If this is None, hand individual results to 'res_files'!
    :type comb_res_in: str, or None
    :type plot_dir: str, or None
    :param comb_res_out: The pathname where to save the final 
        :class:`pyodine.timeseries.base.CombinedResults` object into. If None, 
        the results are not saved.
    :type comb_res_out: str, or None
    :param vels_out: The pathname of a text-file to write chosen timeseries 
        results into. If None, no results are written.
    :type vels_out: str, or None
    :param reject_files: A pathname to a text-file with pathnames of individual 
        results to reject from the combination, or, alternatively, a list of 
        pathnames to individual results. If None, all results are used in the 
        combination algorithm.
    :type reject_files: str, list, tuple, or None
    :param bary_dict: A dictionary with stellar (and possible observatory) 
        information that should be used in the computation of barycentric 
        velocity corrections. If None, the info from the model results is 
        used.
        Possible entries:
        'star_ra' and 'star_dec' (in deg), 'star_pmra' and 
        'star_pmdec' (in mas/yr), 'star_rv0' (in m/s), 'star_name' (e.g. 
        HIPXXX), 'instrument_lat' and 'instrument_long' (in deg), 
        'instrument_alt' (in m).
    :type bary_dict: dict, or None
    :param temp_vel: Velocity offset of the template used in the computation of
        the observation velocities (required if precise barycentric correction 
        is performed). If not supplied here, the template velocity from the
        results file is used.
    :type temp_vel: float, int, or None
    :param ref_vel: Velocity offset of the reference spectrum used for the 
        template velocity estimate (required if precise barycentric correction 
        is performed). If not supplied here, it will be 0.
    :type ref_vel: float, int, or None
    :param error_log: A pathname of a log-file used for error messages. If 
        None, no errors are logged.
    :type error_log: str, or None
    :param info_log: A pathname of a log-file used for info messages. If 
        None, no info is logged.
    :type info_log: str, or None
    :param compact: If True, use a compact version of the CombinedResults
        (only when loading individual results), where only the bare minimum of
        parameters is loaded (to prevent memory crashed for very large time
        series). Defaults to False.
    :type compact: bool
    :param quiet: Whether or not to print info messages to terminal. Defaults 
        to False (messages are printed).
    :type quiet: bool
    
    :return: The final :class:`pyodine.timeseries.base.CombinedResults` object, 
        containing the timeseries results.
    :rtype: :class:`pyodine.timeseries.base.CombinedResults`
    """
    
    # Check whether a logger is already setup. If no, setup a new one
    #if not logging.getLogger().hasHandlers():
    pyodine.lib.misc.setup_logging(
            config_file=Pars.log_config_file, level=Pars.log_level,
            error_log=error_log, info_log=info_log, quiet=quiet)
    
    try:
        # Start timer
        start_t = time.time()
        
        # Get git information
        branch_name = pyodine.lib.git_check.get_git_branch_name()
        branch_hash = pyodine.lib.git_check.get_git_revision_short_hash()
        
        logging.info('')
        logging.info('Branch: {}'.format(branch_name))
        logging.info('Hash: {}'.format(branch_hash))
        logging.info('---------------------------')
        logging.info('Weighting and combining velocities')
        
        ###########################################################################
        ## Set up the environment, and load all neccessary data and parameters
        ###########################################################################
        
        # Set up the CombinedResults object, and load from file
        # Either load a list of individual fit results, handed directly through
        # res_files or in a text-file, or load a previously saved CombinedResults 
        # object if a filename has been supplied through comb_res_in
        Results = pyodine.timeseries.base.CombinedResults()
        
        if isinstance(res_files, str):
            with open(res_files, 'r') as f:
                res_names = [l.strip() for l in f.readlines()]
            Results.load_individual_results(res_names, compact=compact)
            
        elif isinstance(res_files, (list,tuple)):
            res_names = res_files
            Results.load_individual_results(res_names, compact=compact)
            
        elif isinstance(comb_res_in, str):
            Results.load_combined(comb_res_in)
            
        else:
            raise ValueError('Either hand individual fit results through "res_files"' +
                             'as list or tuple or in a text-file, or an existing' +
                             'CombinedResults object through "comb_res_in"!')
        
        # Output directory for plots (setup the directory structure if non-existent)
        if isinstance(plot_dir, str):
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
        
        # Final output name for the CombinedResults object (setup the directory 
        # structure if non-existent)
        if isinstance(comb_res_out, str):
            comb_res_dir = os.path.dirname(comb_res_out)
            if not os.path.exists(comb_res_dir) and comb_res_dir != '':
                os.makedirs(comb_res_dir)
        
        # Output name for the RV text file (in .vels format or any other defined in
        # the parameter input file) (setup the directory structure if non-existent)
        if isinstance(vels_out, str):
            vels_out_dir = os.path.dirname(vels_out)
            if not os.path.exists(vels_out_dir) and vels_out_dir != '':
                os.makedirs(vels_out_dir)
        
        # Load a list of files that should be rejected in the timeseries
        if isinstance(reject_files, (list,tuple)):
            reject_names = reject_files
        elif isinstance(reject_files, str):
            with open(reject_files, 'r') as f:
                reject_names = [l.strip() for l in f.readlines()]
        else:
            reject_names = None
        
        
        ###########################################################################
        ## Now do the velocity weighting and combination, as prescribed in the 
        ## parameter input file
        ###########################################################################
        
        logging.info('')
        logging.info('Star: {}'.format(Results.info['star_name']))
        
        # Possibly first throw out bad observations
        if isinstance(reject_names, (list,tuple)):
            logging.info('')
            logging.info('Rejecting files...')
            if Pars.reject_type == 'obs_files':
                Results.remove_observations(obs_names=reject_names)
            else:
                Results.remove_observations(res_names=reject_names)
        
        logging.info('')
        logging.info('Velocity weighting and combination...')
        
        if Pars.weighting_algorithm == 'song':
            Results.create_timeseries(weighting_pars=Pars.weighting_pars_song, 
                                      do_crx=Pars.do_crx, crx_pars=Pars.crx_pars)
        elif Pars.weighting_algorithm == 'lick':
            Results.create_timeseries_dop(weighting_pars=Pars.weighting_pars_lick, 
                                          do_crx=Pars.do_crx, 
                                          crx_pars=Pars.crx_pars)
        
        # Do precise barycentric velocity correction (using the multiplicative
        # method taking the absolute measured Doppler shift into account)
        if Pars.compute_bvc == 'precise':
            logging.info('')
            logging.info('Doing the precise barycentric velocity correction...')
            Results.compute_bvcs(use_hip=Pars.use_hip_for_bvc,
                                bary_dict=bary_dict, precise=True, 
                                temp_vel=temp_vel, ref_vel=ref_vel,
                                solar=Pars.solar_bvc)
            
        # Or estimate predictive barycentric velocities using barycorrpy, and
        # correct RVs using the additive formula? (less precise)
        elif Pars.compute_bvc == 'predictive':
            logging.info('')
            logging.info('Predictive barycentric velocity computation...')
            Results.compute_bvcs(use_hip=Pars.use_hip_for_bvc, 
                                 bary_dict=bary_dict, precise=False,
                                 solar=Pars.solar_bvc)
        
        
        if isinstance(vels_out, str):
            logging.info('')
            Results.results_to_txt(vels_out, outkeys=Pars.txt_outkeys, 
                                   delimiter=Pars.txt_delimiter, 
                                   header=Pars.txt_header,
                                   outformat=Pars.txt_outformat,
                                   detailed=Pars.txt_detailed,
                                   flux_chunk=Pars.txt_flux_chunk)
        
        ###########################################################################
        ## Possibly save the CombinedResults object and create analysis plots
        ###########################################################################
        
        if Pars.save_comb_res and isinstance(comb_res_out, str):
            logging.info('')
            Results.save_combined(comb_res_out)
        
        if Pars.plot_analysis and isinstance(plot_dir, str):
            logging.info('')
            logging.info('Creating and saving analysis plots to:')
            logging.info('\t{}'.format(plot_dir))
            
            # Plot velocity results
            fig = plt.figure(figsize=(10,6))
            plt.errorbar(Results.bary_date, Results.rv_bc, yerr=Results.rv_err, 
                         fmt='.', alpha=0.7, label='Weighted velocities:\n{:.2f}+-{:.2f} m/s'.format(
                                 robust_mean(Results.rv_bc),
                                 robust_std(Results.rv_bc)))
            plt.plot(Results.bary_date, Results.mdvel+Results.bary_vel_corr+robust_mean(Results.rv_bc), 
                     '.', alpha=0.5, label='Median velocities:\n{:.2f}+-{:.2f} m/s'.format(
                             robust_mean(Results.mdvel+Results.bary_vel_corr),
                             robust_std(Results.mdvel+Results.bary_vel_corr)))
            plt.legend()
            plt.xlabel('JD')
            plt.ylabel('RV [m/s]')
            plt.title('{}, RV time series'.format(Results.info['star_name']))
            plt.savefig(os.path.join(plot_dir, 'RV_timeseries.png'), format='png', dpi=300)
            plt.close()
            
            # Same as above, but outliers rejected
            mask_rvs, good_rvs, bad_rvs = pyodine.lib.misc.chauvenet_criterion(Results.rv_bc)
            rv_good = Results.rv_bc[good_rvs]
            bjd_good = Results.bary_date[good_rvs]
            rv_err_good = Results.rv_err[good_rvs]
            mdvel_good = Results.mdvel[good_rvs]
            bvc_good = Results.bary_vel_corr[good_rvs]
            
            fig = plt.figure(figsize=(10,6))
            plt.errorbar(bjd_good, rv_good, yerr=rv_err_good, fmt='.', alpha=0.7,
                         label='Weighted velocities:\n{:.2f}+-{:.2f} m/s'.format(
                                 robust_mean(rv_good),
                                 robust_std(rv_good)))
            plt.plot(bjd_good, mdvel_good+bvc_good, '.', alpha=0.5,
                     label='Median velocities:\n{:.2f}+-{:.2f} m/s'.format(
                             robust_mean(mdvel_good+bvc_good),
                             robust_std(mdvel_good+bvc_good)))
            plt.legend()
            plt.xlabel('JD')
            plt.ylabel('RV [m/s]')
            plt.title('{}, RV time series, without {} outliers'.format(
                    Results.info['star_name'], len(bad_rvs[0])))
            plt.savefig(os.path.join(plot_dir, 'RV_timeseries_goodobs.png'), format='png', dpi=300)
            plt.close()
            
            # Print the outliers to file, if desired
            if Pars.print_outliers:
                logging.info('')
                logging.info('Observations with outlier RVs:')
                if len(bad_rvs[0]) > 0:
                    for i in bad_rvs[0]:
                        logging.info(Results.res_filename[i])
            
            # Plot chunk-to-chunk scatter of observations
            fig = plt.figure(figsize=(10,6))
            plt.plot(Results.bary_date, Results.c2c_scatter, '.', alpha=0.7, 
                     label='Mean: {:.2f}+-{:.2f} m/s'.format(
                             robust_mean(Results.c2c_scatter),
                             robust_std(Results.c2c_scatter)))
            plt.legend()
            plt.xlabel('JD')
            plt.ylabel('Chunk scatter [m/s]')
            plt.title('{}, chunk scatter of observations'.format(Results.info['star_name']))
            plt.savefig(os.path.join(plot_dir, 'c2c_scatter.png'), format='png', dpi=300)
            plt.close()
            
            # Plot of chunk sigmas
            fig = plt.figure(figsize=(10,6))
            plt.plot(Results.auxiliary['chunk_sigma'], '.', alpha=0.7, 
                     label='Mean: {:.2f}+-{:.2f} m/s'.format(
                             robust_mean(Results.auxiliary['chunk_sigma']), 
                             robust_std(Results.auxiliary['chunk_sigma'])))
            plt.legend()
            plt.xlabel('Chunks')
            plt.ylabel('Chunk sigmas [m/s]')
            plt.title('{}, chunk sigmas'.format(Results.info['star_name']))
            plt.savefig(os.path.join(plot_dir, 'chunk_sigma.png'), format='png', dpi=300)
            plt.close()
            
            # Plot of chunk-to-chunk offsets
            fig = plt.figure(figsize=(10,6))
            plt.plot(Results.auxiliary['chunk_offsets'], '.', alpha=0.7,
                     label='Mean: {:.2f}+-{:.2f}'.format(
                             robust_mean(Results.auxiliary['chunk_offsets']), 
                             robust_std(Results.auxiliary['chunk_offsets'])))
            plt.legend()
            plt.xlabel('Chunks')
            plt.ylabel('Chunk offsets [m/s]')
            plt.title('{}, chunk offsets from observation means'.format(Results.info['star_name']))
            plt.savefig(os.path.join(plot_dir, 'chunk_offsets.png'), format='png', dpi=300)
            plt.close()
            
            # 3D plot of velocities (corrected by chunk offsets & barycentric velocities)
            vel_corrected = Results.params['velocity'] - Results.auxiliary['chunk_offsets']
            vel_corrected = vel_corrected.T + Results.bary_vel_corr
            
            fig = plt.figure(figsize=(12,10))
            plt.imshow(vel_corrected, aspect='auto')
            plt.colorbar()
            plt.xlabel('Observations')
            plt.ylabel('Chunks')
            plt.title('{}, offset-BV-corrected chunk velocities'.format(Results.info['star_name']))
            plt.savefig(os.path.join(plot_dir, 'chunk_vels_corr.png'), format='png', dpi=300)
            plt.close()
            
            # 3D plot of chunk deviations
            fig = plt.figure(figsize=(12,10))
            plt.imshow(Results.auxiliary['chunk_dev'].T, aspect='auto')
            plt.colorbar()
            plt.xlabel('Observations')
            plt.ylabel('Chunks')
            plt.title('{}, chunk deviations'.format(Results.info['star_name']))
            plt.savefig(os.path.join(plot_dir, 'chunk_devs.png'), format='png', dpi=300)
            plt.close()
            
            # Histogram of the final velocity weights
            fig = plt.figure(figsize=(10,6))
            plt.hist(Results.auxiliary['chunk_weights'].flatten(), bins=100, alpha=0.7,
                     label='Mean: {:.3e}+-{:.3e} (m/s)$^{{-2}}$'.format(
                             robust_mean(Results.auxiliary['chunk_weights']), 
                             robust_std(Results.auxiliary['chunk_weights'])))
            plt.legend()
            plt.xlabel(r'Weights [(m/s)$^{-2}$]')
            plt.title('{}, chunk weights'.format(Results.info['star_name']))
            plt.savefig(os.path.join(plot_dir, 'chunk_weights_hist.png'), format='png', dpi=300)
            plt.close()
            
            # Plot chromatic indices (if any)
            if Pars.do_crx:
                fig = plt.figure(figsize=(10,6))
                plt.errorbar(Results.bary_date, Results.crx, yerr=Results.crx_err, 
                             fmt='.', alpha=0.7, label='Mean: {:.2f}+-{:.2f} (m/s)/Np'.format(
                                     robust_mean(Results.crx), 
                                     robust_std(Results.crx)))
                plt.legend()
                plt.xlabel('JD')
                plt.ylabel('CRX [(m/s)/Np]')
                plt.title('{}, CRX time series'.format(Results.info['star_name']))
                plt.savefig(os.path.join(plot_dir, 'CRX_timeseries.png'), format='png', dpi=300)
                plt.close()
            
        ###########################################################################
        ## Everything's done now, return the CombinedResults object
        ###########################################################################
        
        work_time = time.time() - start_t
        logging.info('')
        logging.info('All done! Full work time: {}'.format(work_time))
        
        return Results
    
    except Exception as e:
        logging.error('Something went wrong!', exc_info=True)


if __name__ == '__main__':
    
    # Set up the parser for input arguments
    parser = argparse.ArgumentParser(
            description='Weight and combine velocities from observation modelling')
    
    # Required input arguments:
    # utilities_dir, ostar_files, temp_files, temp_outname, (plot_dir=None, par_file=None)
    parser.add_argument('par_file', type=str, help='The pathname to the timeseries parameters file to use.')
    parser.add_argument('--res_files', type=str, help='A pathname to a text-file with the pathnames of modelling results.')
    parser.add_argument('--comb_res_in', type=str, help='The pathname to a saved CombinedResults object.')
    parser.add_argument('--plot_dir', type=str, help='The pathname to a directory where to save analysis plots.')
    parser.add_argument('--comb_res_out', type=str, help='The pathname where to save the CombinedResults object.')
    parser.add_argument('--vels_out', type=str, help='The pathname of a text-file where to save chosen timeseries results.')
    parser.add_argument('--reject_files', type=str, help='A pathname of a text-file with the pathnames of results to reject.')
    parser.add_argument('--temp_vel', type=float, help='Optional template velocity offset to use in barycentric correction.')
    parser.add_argument('--ref_vel', type=float, help='Optional reference velocity offset to use in barycentric correction.')
    parser.add_argument('--error_file', type=str, help='The pathname to the error log file.')
    parser.add_argument('--info_file', type=str, help='The pathname to the info log file.')
    parser.add_argument('-c', '--compact', action='store_true', dest='compact', help='Load compact version of combined results.')
    parser.add_argument('-q', '--quiet', action='store_true', dest='quiet', help='Do not print messages to the console.')
    
    # Parse the input arguments
    args = parser.parse_args()
    
    par_file     = args.par_file
    res_files    = args.res_files
    comb_res_in  = args.comb_res_in
    plot_dir     = args.plot_dir
    comb_res_out = args.comb_res_out
    vels_out     = args.vels_out
    reject_files = args.reject_files
    temp_vel     = args.temp_vel
    ref_vel      = args.ref_vel
    error_file   = args.error_file
    info_file    = args.info_file
    compact      = args.compact
    quiet        = args.quiet
    
    # Import and load the timeseries parameters
    par_file = os.path.splitext(par_file)[0].replace('/', '.')
    timeseries_parameters = importlib.import_module(par_file)
    Pars = timeseries_parameters.Timeseries_Parameters()
    
    # And run the velocity weighting routine
    combine_velocity_results(Pars, res_files=res_files, comb_res_in=comb_res_in, 
                             plot_dir=plot_dir, comb_res_out=comb_res_out, 
                             vels_out=vels_out, reject_files=reject_files, 
                             temp_vel=temp_vel, ref_vel=ref_vel,
                             error_log=error_file, info_log=info_file, 
                             compact=compact, quiet=quiet)
