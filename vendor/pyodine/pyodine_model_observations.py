#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 17:21:06 2020

@author: pheeren
"""

# Import packages
import pyodine
import pipe_lib

import os
import sys
import time
import numpy as np
import logging
from pathos.multiprocessing import Pool
import traceback

import argparse
import importlib


def model_single_observation(utilities, Pars, obs_file, temp_file, 
                             iod=None, orders=None, order_correction=None, 
                             normalizer=None, tellurics=None, plot_dir=None, 
                             res_names=None, error_log=None, info_log=None, 
                             quiet=False, live=False):
    """Model a single observation
    
    This routine models a stellar observation spectrum with I2, using a stellar
    template spectrum, I2 template and LSF model. The modelling can be
    performed in multiple runs, in order to allow a better determination of the
    fit parameters. The results and analysis plots can be saved to file.
    
    :param utilities: The utilities module for the instrument used in this 
        analysis.
    :type utilities: library
    :param Pars: The parameter input object to use.
    :type Pars: :class:`Parameters`
    :param obs_file: The pathname of the stellar observation to model.
    :type obs_file: str
    :param temp_file: The pathname of the deconvolved stellar template to use 
        in the modelling.
    :type temp_file: str
    :param iod: The I2 template to use in the modelling. If None, it is loaded 
        as specified in the parameter input object.
    :type iod: :class:`IodineTemplate`, or None
    :param orders: The orders of the observation to work on. If None, they are 
        defined as specified by the template and the parameter input object.
    :type orders: list, ndarray, or None
    :param order_correction: Possible order shift between template and 
        observation. If not given, the order shift is estimated in the code.
    :type order_correction: int, or None
    :param normalizer: The normalizer to use. If None, it is created as 
        specified in the parameter input object.
    :type normalizer: :class:`SimpleNormalizer`, or None
    :param tellurics: The tellurics to use. If None, it is created as specified 
        in the parameter input object.
    :type tellurics: :class:`SimpleTellurics`, or None
    :param plot_dir: The directory name where to save plots. If the directory 
        structure does not exist yet, it will be created in the process. If 
        None is given, no plots will be saved (default).
    :type plot_dir: str, or None
    :param res_names: The pathname under which to save the results file. If you 
        want to save results from multiple runs, you should supply a list with 
        pathnames for each run. If the directory structure does not exist yet, 
        it will be created in the process. If None is given, no results will be 
        saved (default).
    :type res_names: str, list, or None
    :param error_log: A pathname of a log-file used for error messages. If 
        None, no errors are logged.
    :type error_log: str, or None
    :param info_log: A pathname of a log-file used for info messages. If 
        None, no info is logged.
    :type info_log: str, or None
    :param quiet: Whether or not to print info messages to terminal. Defaults 
        to False (messages are printed).
    :type quiet: bool
    :param live: If True, then the modelling is performed in live-mode, i.e.
        each modelled chunk is plotted and the best-fit parameters printed to
        terminal. Defaults to False.
    :type live: bool
    """
    
    # Check whether a logger is already setup. If no, setup a new one
    #if not logging.getLogger().hasHandlers():
    pyodine.lib.misc.setup_logging(
            config_file=Pars.log_config_file, level=Pars.log_level,
            error_log=error_log, info_log=info_log, quiet=quiet)
    
    # I put everything in try - except, so that when it runs in parallel 
    # processes and something goes wrong, the overall routine does not crash.
    # (Is this elegant though?)
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
        logging.info('Working on: {}'.format(obs_file))
        
        ###########################################################################
        ## Set up the environment, and load all neccessary data and parameters
        ###########################################################################
        
        # Load observation
        obs = utilities.load_pyodine.ObservationWrapper(obs_file)
        
        # And log an idea of the flux level of the observation
        logging.info('')
        logging.info('Median flux of the observation: {:.0f}'.format(
                np.median([np.median(obs[o].flux) for o in obs.orders])))
        
        # Load the deconvolved stellar template
        template = pyodine.template.base.StellarTemplate_Chunked(temp_file)
        
        # If the I2 template is not given, load it
        if not iod:
            iod = utilities.load_pyodine.IodineTemplate(Pars.i2_to_use)
        
        # Output directory for plots (setup the directory structure if non-existent)
        if isinstance(plot_dir, str):
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
        
        # Pathnames for the result files
        if isinstance(res_names, list) and isinstance(res_names[0], str):
            res_names = res_names
        elif isinstance(res_names, str):
            res_names = [res_names]
        else:
            res_names = None
        
        
        ###########################################################################
        ## Now prepare the modelling: Choose the orders, compute weights,
        ## maybe tellurics...
        ###########################################################################
        
        # If no order range is given, define it
        if not isinstance(orders, (list, tuple, np.ndarray)):
            orders = template.orders_unique
            if isinstance(Pars.order_range, (list, tuple)) and Pars.order_range[0] is not None:
                orders = orders[Pars.order_range[0]:Pars.order_range[1]+1]
        
        # Compute possible order shifts between template and observation,
        # by searching for the best coverage of first template order in
        # observation (if not supplied in the top)
        if not isinstance(order_correction, int):
            template_ind = template.get_order_indices(template[0].order)
            obs_order_min, min_coverage = obs.check_wavelength_range(
                    template[template_ind[0]].w0, template[template_ind[-1]].w0)
            order_correction = obs_order_min - template[template_ind[0]].order
        logging.info('')
        logging.info('Order correction: {}'.format(order_correction))
        
        # Compute weights array
        weight = obs.compute_weight(weight_type=Pars.weight_type, rel_noise=Pars.rel_noise)
        # Pixels 1001 - 1004 in Lick orders 31 - 60 are corrupted -> set weights to zero there
        #for i in range(len(orders)):
        #    weight[1001:1005] = 0.
        
        # Potentially compute a bad pixel mask for the observation
        if Pars.bad_pixel_mask is True:
            mask = pyodine.bad_pixels.BadPixelMask(obs, cutoff=Pars.bad_pixel_cutoff)
            weight[np.where(mask.mask == 1.)] = 0.
        
        # Correct spectrum at weight == 0 if desired
        if Pars.correct_obs is True:
            obs = pyodine.bad_pixels.correct_spec.correct_spectrum(obs, weight, orders)
        
        # Initialize Normalizer with reference spectrum, then cross-correlate 
        # to obtain the velocity guess of the observation relative to the
        # template velocity
        if not normalizer:
            normalizer = pyodine.template.normalize.SimpleNormalizer(reference=Pars.ref_spectrum)
        
        ref_velocity = normalizer.guess_velocity(
                obs[Pars.velgues_order_range[0]:Pars.velgues_order_range[1]],
                delta_v=Pars.delta_v, maxlag=Pars.maxlag)
        obs_velocity = ref_velocity - template.velocity_offset
        logging.info('')
        logging.info('Measured velocity rel. to reference spectrum: {0:.3f} km/s'.format(ref_velocity*1e-3))
        logging.info('Template velocity: {0:.3f} km/s'.format(template.velocity_offset*1e-3))
        logging.info('Velocity guess: {0:.3f} km/s (relative to template)\n'.format(obs_velocity*1e-3))
        logging.info('(Barycentric velocity of observation: {0:.3f} km/s)'.format(obs.bary_vel_corr*1e-3))
        logging.info('(Barycentric velocity of template: {0:.3f} km/s)\n'.format(template.bary_vel_corr*1e-3))
        
        # Load the tellurics (if desired and not given)
        if not tellurics and Pars.telluric_mask is not None:
            tellurics = pyodine.tellurics.SimpleTellurics(tell_type=Pars.telluric_mask, wave_min=Pars.tell_wave_range[0],
                                                          wave_max=Pars.tell_wave_range[1], disp=Pars.tell_dispersion)
        
        # If the observation spectrum should be normalized prior to fitting, 
        # this is done here
        if Pars.normalize_chunks is True:
            for o in obs.orders:
                obs._flux[o] = (obs[o].flux / obs[o].cont)
        
        # Now create the chunks, using the algorithm (and corresponding parameters) 
        # as defined in the parameter input file
        if Pars.chunking_algorithm == 'auto_wave_comoving':
            obs_chunks = pyodine.chunks.auto_wave_comoving(
                    obs, template, orders=orders, padding=Pars.chunk_padding, 
                    order_correction=order_correction, delta_v=Pars.chunk_delta_v)
        else:
            raise KeyError('Algorithm {} not known! (Only option right now: auto_wave_comoving)'.format(
                    Pars.chunking_algorithm))
            
        nr_chunks_total  = len(obs_chunks)
        nr_chunks_order0 = len(obs_chunks.get_order(obs_chunks.orders[0]))
        nr_orders_chunks = len(obs_chunks.orders)
        
        logging.info('')
        logging.info('Total number of chunks: {}'.format(nr_chunks_total))
        logging.info('Nr. chunks in order 0: {}'.format(nr_chunks_order0))
        logging.info('First and last covered pixel of chunks in order 0: {}, {}'.format(
                obs_chunks[0].abspix[0], obs_chunks[nr_chunks_order0-1].abspix[-1]))
        logging.info('Orders: {} - {} ({} in total)'.format(
                obs_chunks[0].order, obs_chunks[-1].order, nr_orders_chunks))
        
        # Produce the chunk weight array
        chunk_weight = []
        for chunk in obs_chunks:
            chunk_weight.append(np.array(weight[chunk.order, chunk.abspix[0]:chunk.abspix[-1]+1]))
        #chunk_weight = np.array(chunk_weight)
        
        ###########################################################################
        ## Build and fit the model in as many runs as specified in the 
        ## input parameter file.
        ## Other things specified there:
        ## - which parameters to vary, and possibly limits and starting values
        ## - should output plots be created
        ## - smoothing of wavelength parameters
        ###########################################################################
        
        # Set up the dictionary holding all results from the runs
        run_results = {}
        for run_id in range(len(Pars.model_runs)):
            run_results[run_id] = {
                    'velocity_guess': obs_velocity,
                    'obs_bvc': obs.bary_vel_corr,
                    'results': [],
                    'fitting_failed': [],
                    'red_chi_sq': np.zeros((nr_chunks_total)),
                    'chunk_w': [], #np.ones((chunk_weight.shape)),
                    'median_pars': pyodine.models.base.ParameterSet(),
                    'std_pars': pyodine.models.base.ParameterSet()
                    }
        
        # Run for run now...
        for run_id, run_dict in Pars.model_runs.items():
            logging.info('')
            logging.info('----------------------')
            logging.info('RUN {}'.format(run_id))
            logging.info('----------------------')
            
            # Build the desired LSF, wavelength and continuum models
            lsf_model = run_dict['lsf_model']
            # For the LSF model: Potentially adapt it to the instrument
            if 'lsf_setup_dict' in run_dict.keys() and isinstance(run_dict['lsf_setup_dict'], dict):
                lsf_model.adapt_LSF(run_dict['lsf_setup_dict'])
            
            wave_model = run_dict['wave_model'] #pyodine.models.wave.LinearWaveModel
            cont_model = run_dict['cont_model'] #pyodine.models.cont.LinearContinuumModel
            
            # If the LSF model is a fixed LSF, try and smooth LSF results from
            # an earlier run
            if lsf_model.name() == 'FixedLSF':
                if 'smooth_lsf_run' not in run_dict.keys():
                    smooth_lsf_run = run_id - 1
                else:
                    smooth_lsf_run = run_dict['smooth_lsf_run']
                # Check if the run LSFs to smooth over actually exist
                if smooth_lsf_run in run_results.keys() and len(run_results[smooth_lsf_run]['results']) == nr_chunks_total:
                    # Smooth the desired LSFs over adjacent chunks and orders
                    logging.info('')
                    logging.info('Smooth LSF...')
                    manual_redchi2 = None
                    if 'smooth_manual_redchi' in run_dict.keys() and run_dict['smooth_manual_redchi']:
                        manual_redchi2 = run_results[smooth_lsf_run]['red_chi_sq']
                    smooth_osample = None
                    if 'smooth_osample' in run_dict.keys() and run_dict['smooth_osample'] > 0:
                        smooth_osample = int(run_dict['smooth_osample'])
                    
                    lsf_smoothed = pyodine.lib.misc.smooth_lsf(
                            obs_chunks, run_dict['smooth_pixels'], run_dict['smooth_orders'], 
                            run_dict['order_separation'], run_results[smooth_lsf_run]['results'],
                            redchi2=manual_redchi2, osample=smooth_osample)
                    logging.info('')
                    logging.info('LSFs with nans: {}'.format(len(np.unique(np.argwhere(np.isnan(lsf_smoothed))[:,0]))))
                    
                    LSFarr = pyodine.models.lsf.LSF_Array(lsf_smoothed, np.array([ch.order for ch in obs_chunks]),
                                                          np.array([ch.abspix[0] for ch in obs_chunks]))
                    
                    # Now create the full model, including the fixed LSF
                    model = pyodine.models.spectrum.SimpleModel(
                            lsf_model, wave_model, cont_model, iod, 
                            stellar_template=template, lsf_array=LSFarr,
                            osample_factor=Pars.osample_obs, conv_width=Pars.lsf_conv_width)
                        
                else:
                    raise KeyError('smooth_lsf_run not in run_results or no results for that run yet!')
            
            else:
                lsf_smoothed = None
                # If no fixed LSF model, set up a "normal" full model
                model = pyodine.models.spectrum.SimpleModel(
                        lsf_model, wave_model, cont_model, iod, stellar_template=template, 
                        osample_factor=Pars.osample_obs, conv_width=Pars.lsf_conv_width)
            
            # Set up the fitter
            fitter = pyodine.fitters.lmfit_wrapper.LmfitWrapper(model)
            
            # Loop through all chunks to pre-guess the parameters
            starting_pars = []
            for i, chunk in enumerate(obs_chunks):
                starting_pars.append(model.guess_params(chunk))
            
            # Fit wavelength guess within each order if desired (smoother input)
            # (this really only makes sense for the very first run, otherwise
            # use smoothed wavelength results from previous runs)
            if 'pre_wave_slope_deg' in run_dict.keys() and run_dict['pre_wave_slope_deg'] > 0:
                poly_pars = pyodine.lib.misc.smooth_parameters_over_orders(
                        starting_pars, 'wave_slope', obs_chunks, deg=run_dict['pre_wave_slope_deg'])
                # Write the smoothed dispersion into starting_pars
                for i in range(nr_chunks_total):
                    starting_pars[i]['wave_slope'] = poly_pars[i]
                    
            if 'pre_wave_intercept_deg' in run_dict.keys() and run_dict['pre_wave_intercept_deg'] > 0:
                poly_pars = pyodine.lib.misc.smooth_parameters_over_orders(
                        starting_pars, 'wave_intercept', obs_chunks, deg=run_dict['pre_wave_intercept_deg'])
                # Write the smoothed intercepts into starting_pars
                for i in range(nr_chunks_total):
                    starting_pars[i]['wave_intercept'] = poly_pars[i]
            
            # Convert starting parameters to lmfit-parameters 
            lmfit_params = []
            for i in range(nr_chunks_total):
                lmfit_params.append(fitter.convert_params(starting_pars[i], to_lmfit=True))
            
            # Constrain the parameters, using the definitions as supplied in
            # the parameter input file
            lmfit_params = Pars.constrain_parameters(
                    lmfit_params, run_id, run_results, fitter
                    )
            
            ###########################################################################
            ## Finally loop over the chunks and model each of them.
            ###########################################################################
            
            use_chauvenet = False
            if 'use_chauvenet_pixels' in run_dict.keys() and run_dict['use_chauvenet_pixels']:
                use_chauvenet = True
            
            modelling_return = pipe_lib.model_all_chunks(
                    obs_chunks, chunk_weight, fitter, lmfit_params, 
                    tellurics, use_chauvenet=use_chauvenet, compute_redchi2=True, 
                    use_progressbar=Pars.use_progressbar, live=live)
            
            (run_results[run_id]['results'], run_results[run_id]['chunk_w'], 
             run_results[run_id]['fitting_failed'], chauvenet_outliers, 
             run_results[run_id]['red_chi_sq']) = modelling_return
            
            
            ###########################################################################
            ## Now we do some run diagnostics and analysis of the model results.
            ## Plotting and results saving is done only if a res_dir was supplied.
            ###########################################################################
            
            # For which chunks could uncertainties not be estimated?! (Why?)
            uncertainties_failed = [i for i in range(nr_chunks_total) if \
                                    'uncertainties could not be estimated' in run_results[run_id]['results'][i].report]
            # For which chunks is the red. Chi2 nan?!
            nan_rchi_fit = [i for i in range(nr_chunks_total) if np.isnan(run_results[run_id]['results'][i].redchi)]
            logging.info('')
            logging.info('Number of chunks with no uncertainties: {}'.format(len(uncertainties_failed)))
            logging.info('Number of chunks with outliers: {}'.format(len(chauvenet_outliers)))
            logging.info('Number of chunks with nan fitted red. Chi2: {}'.format(len(nan_rchi_fit)))
            
            
            # Save the fitting result to file
            # Check in the list for the correct results name
            if isinstance(res_names, list) and 'save_result' in run_dict.keys() and run_dict['save_result'] is True:
                if (len(res_names) - 1) >= run_id:
                    res_save_name = res_names[run_id]
                else:
                    res_save_name = res_names[0]
                # Create the directory structure if it does not exist yet
                res_save_name_dir = os.path.dirname(res_save_name)
                if not os.path.exists(res_save_name_dir) and res_save_name_dir != '':
                    os.makedirs(res_save_name_dir)
                # Save it under the correct file type
                if 'save_filetype' in run_dict.keys() and run_dict['save_filetype'] == 'dill':
                    pyodine.fitters.results_io.save_results(
                            res_save_name, run_results[run_id]['results'], filetype='dill')
                else:
                    pyodine.fitters.results_io.save_results(
                            res_save_name, run_results[run_id]['results'], filetype='h5py')
            
            
            wave_slope_fit = None
            wave_intercept_fit = None
            
            # Fit best-fit wavelength slopes within each order with a polynomial
            if 'wave_slope_deg' in run_dict.keys() and run_dict['wave_slope_deg'] > 0:
                wave_slope_fit = pyodine.lib.misc.smooth_fitresult_over_orders(
                        run_results[run_id]['results'], 'wave_slope', deg=run_dict['wave_slope_deg'])
                run_results[run_id]['wave_slope_fit'] = wave_slope_fit
            
            # Fit best-fit wavelength intercepts within each order with a polynomial
            if 'wave_intercept_deg' in run_dict.keys() and run_dict['wave_intercept_deg'] > 0:
                wave_intercept_fit = pyodine.lib.misc.smooth_fitresult_over_orders(
                        run_results[run_id]['results'], 'wave_intercept', deg=run_dict['wave_intercept_deg'])
                run_results[run_id]['wave_intercept_fit'] = wave_intercept_fit
            
            
            # Calculate median parameter results from this run
            param_names = [k for k in run_results[run_id]['results'][0].params.keys()]
            params = pyodine.models.base.ParameterSet() #{k: np.zeros(nr_chunks_total) for k in param_names}
            errors = pyodine.models.base.ParameterSet() #{k: np.zeros(nr_chunks_total) for k in param_names}
            for p in param_names:
                params[p] = np.array([r.params[p] for r in run_results[run_id]['results']])
                errors[p] = np.array([r.errors[p] for r in run_results[run_id]['results']])
                run_results[run_id]['median_pars'][p] = np.nanmedian(params[p])
                run_results[run_id]['std_pars'][p] = np.nanstd(params[p])
            
            # Write the median parameter results to file
            if plot_dir and 'save_median_pars' in run_dict.keys() and run_dict['save_median_pars']:
                # Print median parameters into file
                pars_file = os.path.join(plot_dir, 'r{}_median_pars.txt'.format(run_id))
                
                with open(pars_file, 'a') as f:
                    for p in param_names:
                        f.write('\t' + p + ':\t' + str(run_results[run_id]['median_pars'][p]) + \
                                '\t +/- \t' + str(run_results[run_id]['std_pars'][p]) + '\n')
                    f.write('\n')
            
            
            # More plots
            if plot_dir and 'plot_analysis' in run_dict.keys() and run_dict['plot_analysis'] is True:
                # Plot residuals scatter and histogram, red. Chi**2, possibly evaluated chunks,
                # the wavelength slopes and intercepts along with their fits (if accessible),
                # the LSF parameters (if desired), the fitting success (if desired),
                # and some exemplary smoothed LSFs (if accessible)
                plot_chunks = None
                plot_lsf_pars = False
                nr_chunks_order, nr_orders = None, None
                if 'plot_chunks' in run_dict.keys() and isinstance(run_dict['plot_chunks'], (list,tuple)):
                    plot_chunks = run_dict['plot_chunks']
                if 'plot_lsf_pars' in run_dict.keys() and run_dict['plot_lsf_pars']:
                    plot_lsf_pars = True
                if 'plot_success' not in run_dict.keys() or not run_dict['plot_success']:
                    uncertainties_failed = None
                    nan_rchi_fit = None
                    chauvenet_outliers = None
                # If all orders are split into equal number of chunks
                if len(np.unique([len(obs_chunks.get_order(o)) for o in obs_chunks.orders])) == 1:
                    nr_chunks_order = nr_chunks_order0
                    nr_orders = nr_orders_chunks
                
                logging.info('')
                logging.info('Creating analysis plots...')
                
                pipe_lib.create_analysis_plots(
                        run_results[run_id]['results'], plot_dir, run_id=run_id, 
                        tellurics=tellurics, red_chi_sq=run_results[run_id]['red_chi_sq'], 
                        nr_chunks_order=nr_chunks_order, nr_orders=nr_orders, 
                        chunk_weight=run_results[run_id]['chunk_w'], plot_chunks=plot_chunks, 
                        chunks=obs_chunks, 
                        wave_intercept_fit=wave_intercept_fit, wave_slope_fit=wave_slope_fit,
                        plot_lsf_pars=plot_lsf_pars,
                        uncertainties_failed=uncertainties_failed,
                        nan_rchi_fit=nan_rchi_fit, chauvenet_outliers=chauvenet_outliers,
                        lsf_array=lsf_smoothed, live=live)
            
            ###########################################################################
            ## Run finished, proceeding to next run (unless all through)
            ###########################################################################
    
        
        ###########################################################################
        ## Now all runs are done, do some very basic analysis of the 
        ## velocity results to provide instant feedback (only if plot_dir
        ## exists).
        ## Then exit.
        ###########################################################################
        
        if plot_dir and isinstance(Pars.vel_analysis_plots, int):
            if Pars.vel_analysis_plots in list(run_results.keys()) \
            or abs(Pars.vel_analysis_plots) <= len(list(run_results.keys())):
                
                run_id = list(run_results.keys())[Pars.vel_analysis_plots]
                
                nr_chunks_order, nr_orders = None, None
                # If all orders are split into equal number of chunks
                if len(np.unique([len(obs_chunks.get_order(o)) for o in obs_chunks.orders])) == 1:
                    nr_chunks_order = nr_chunks_order0
                    nr_orders = nr_orders_chunks
                    
                pipe_lib.velocity_results_analysis(
                        run_results[run_id], plot_dir, nr_chunks_order, nr_orders, obs.orig_filename)
            else:
                logging.warning('')
                logging.warning('Desired run id for velocity analysis plots {} not existent!'.format(
                        Pars.vel_analysis_plots))
        
        modelling_time = time.time() - start_t
        logging.info('')
        logging.info('Time to model this observation: {}'.format(modelling_time))
    
    except Exception as e:
        """
        with open(error_file, 'a') as f:
            f.write(parameters[5]+'\n')
        """
        traceback.print_exc()
        logging.info('Something went wrong with input file {}'.format(obs_file), 
                     exc_info=True)
        print()
        print(e)
        


def model_multi_observations(utilities, Pars, obs_files, temp_files, 
                             order_corrections=None, plot_dirs=None, 
                             res_files=None, error_files=None, info_files=None, 
                             quiet=False):
    """Model multiple observations at the same time
    
    This function can parallelize the modelling of multiple observations,
    taking advantage of Python's :class:`pathos.multiprocessing.Pool` capabilities.
    The number of parallel processes is defined in the parameter input object.
    
    :param utilities: The utilities module for the instrument used in this 
        analysis.
    :type utilities: library
    :param Pars: The parameter input object for the used instrument.
    :type Pars: :class:`Parameters`
    :param obs_files: A pathname to a text-file with pathnames of stellar 
        observations for the modelling, or, alternatively, a list with the 
        pathnames.
    :type obs_files: str or list
    :param temp_files: A pathname to a text-file with pathnames of the 
        deconvolved stellar templates to use for each of the observations, or, 
        alternatively, a list with the pathnames.
    :type temp_files: str or list
    :param order_corrections: A pathname to a text-file with possible order 
        shifts between templates and observations, or, alternatively, a list
        with the order shifts. If not given, the order shifts are estimated in 
        the code.
    :type order_corrections: str, list, or None
    :param plot_dirs: A pathname to a text-file with directory names for each 
        observation where to save plots, or, alternatively, a list with the 
        directory names. If the directory structure does not exist yet, it will 
        be created in the process. If None is given, no results/plots will be 
        saved (default).
    :type plot_dirs: str, list, or None
    :param res_files: A pathname to a text-file with pathnames under which to 
        save the results file(s), or, alternatively, a list with the 
        pathname(s). If you want to save results from multiple runs, you should 
        supply pathnames for each run. If the directory structure does not 
        exist yet, it will be created in the process. If None is given, no 
        results will be saved (default).
    :type res_files: str, list, or None
    :param error_files: A pathname to a text-file with pathnames to log-files 
        used for error messages, or, alternatively, a list with the 
        pathname(s). If None, no errors are logged (default).
    :type error_files: str, list, or None
    :param info_files: A pathname to a text-file with pathnames to log-files 
        used for info messages, or, alternatively, a list with the pathname(s). 
        If None, no info is logged (default).
    :type info_files: str, or None
    :param quiet: Whether or not to print info messages to terminal. Defaults 
        to False (messages are printed).
    :type quiet: bool
    """
    
    ###########################################################################
    ## Some modules and data are loaded here already and then passed to the
    ## individual parallel modelling sessions, so that they do not need to
    ## be loaded in each one individually:
    ## - the stellar template spectrum
    ## - the I2 template spectrum
    ## - the orders to work on
    ## - the normalizer
    ## - potentially tellurics
    ## - the plot directories
    ###########################################################################
    
    # Start timer
    fulltime_start = time.time()
    
    # Load the pathnames of the observations
    if isinstance(obs_files, list):
        obs_names = obs_files
    elif isinstance(obs_files, str):
        with open(obs_files, 'r') as f:
            obs_names = [l.strip() for l in f.readlines()]
    
    # Load the pathnames of the deconvolved stellar templates
    if isinstance(temp_files, list):
        temp_names = temp_files
    elif isinstance(temp_files, str):
        with open(temp_files, 'r') as f:
            temp_names = [l.strip() for l in f.readlines()]
    #template = pyodine.template.base.StellarTemplate_Chunked(temp_file)
    
    # Load the iodine atlas from file
    iod = utilities.load_pyodine.IodineTemplate(Pars.i2_to_use)
    
    # Initialize Normalizer with reference spectrum
    normalizer = pyodine.template.normalize.SimpleNormalizer(reference=Pars.ref_spectrum)
    
    # Load the tellurics (if desired)
    if Pars.telluric_mask is not None:
        tellurics = pyodine.tellurics.SimpleTellurics(tell_type=Pars.telluric_mask, wave_min=Pars.tell_wave_range[0],
                                                      wave_max=Pars.tell_wave_range[1], disp=Pars.tell_dispersion)
    else:
        tellurics = None
    
    # Load the plot directory names for each observation
    if isinstance(plot_dirs, list) and isinstance(plot_dirs[0], str):
        plot_dir_names = plot_dirs
    elif isinstance(plot_dirs, str):
        with open(plot_dirs, 'r') as f:
            plot_dir_names = [l.strip() for l in f.readlines()]
    else:
        plot_dir_names = [None] * len(obs_names)
    
    # Pathnames for the result files
    # Several result names for each run of each observation can be supplied
    # either as a list of lists of names, or a file with several pathnames
    # in each line (for each observation)
    if isinstance(res_files, list):
        res_names = res_files
    elif isinstance(res_files, str):
        res_names = []
        with open(res_files, 'r') as f:
            for l in f.readlines():
                names = l.split() #.strip()
                res_names.append(names)
    else:
        res_names = [None] * len(obs_names)
    
    # Pathnames for the error log files
    if isinstance(error_files, list):
        error_logs = error_files
    elif isinstance(error_files, str):
        error_logs = []
        with open(error_files, 'r') as f:
            for l in f.readlines():
                error_logs.append(l.strip())
    else:
        error_logs = [None] * len(obs_names)
    
    # Pathnames for the info log files
    if isinstance(info_files, list):
        info_logs = info_files
    elif isinstance(info_files, str):
        info_logs = []
        with open(info_files, 'r') as f:
            for l in f.readlines():
                info_logs.append(l.strip())
    else:
        info_logs = [None] * len(obs_names)
    
    # Order corrections
    if isinstance(order_corrections, list):
        order_shifts = order_corrections
    elif isinstance(order_corrections, str):
        order_shifts = []
        with open(order_corrections, 'r') as f:
            for l in f.readlines():
                order_shifts.append(int(l.strip()))
    else:
        order_shifts = [None] * len(obs_names)
    
    
    ###########################################################################
    ## Now parallelize the modelling of the observations, by initializing the
    ## Pool workers and distribute the observations.
    ## Changed this (similar to https://www.py4u.net/discuss/237878)
    ###########################################################################
    
    # Prepare the input arguments list for all the jobs (corresponding
    # to the arguments of the function model_single_observation)
    input_arguments = [
            (utilities, Pars, obs_name, temp_name,
             ) for obs_name, temp_name in zip(obs_names, temp_names)]
    # Prepare the keyword arguments list for all the jobs (corresponding
    # to the keywords of the function model_single_observation)
    input_keywords  = [
            {'iod': iod, 'order_correction': order_shift,
             'normalizer': normalizer, 'tellurics': tellurics, 
             'plot_dir': plot_dir_name, 'res_names': res_name,
             'error_log': error_log, 'info_log': info_log, 'quiet': quiet
             } for plot_dir_name, res_name, error_log, info_log, order_shift in zip(
             plot_dir_names, res_names, error_logs, info_logs, order_shifts)]
    
    # Setup the Pool object, distribute the arguments and start the jobs
    with Pool(Pars.number_cores) as p:
        jobs = [
            p.apply_async(model_single_observation,
                          args=input_arg, kwds=input_kwd
                          ) for input_arg, input_kwd in zip(input_arguments, input_keywords)]
        
        for job in jobs:
            # Wait for the modelling in all workers to finish
            job.wait()
            #job.get()
            # regularly you'd use `job.get()`, but it would `raise` the exception,
            # which is not suitable for this example, so we dig in deeper and just use
            # the `._value` it'd return or raise:
            #print(params, type(job._value), job._value)
    
    # And done
    
    full_modelling_time = time.time() - fulltime_start
    print('\nDone, full working time: ', full_modelling_time)
    


if __name__ == '__main__':
    
    # Set up the parser for input arguments
    parser = argparse.ArgumentParser(
            description='Model a number of observations')
    
    # Required input arguments:
    # utilities_dir, obs_files, temp_file, (plot_dirs=None, res_files=None, par_file=None)
    parser.add_argument('utilities_dir', type=str, help='The pathname to the utilities directory for this instrument.')
    parser.add_argument('obs_files', type=str, help='A pathname to a text-file with pathnames of stellar observations for the modelling.')
    parser.add_argument('temp_files', type=str, help='A pathname to a text-file with pathnames of deconvolved stellar templates to use.')
    parser.add_argument('--order_corrections', type=str, help='A pathname to a text-file with order corrections between observations and templates.')
    parser.add_argument('--plot_dirs', type=str, help='A pathname to a text-file with directory names for each observation where to save analysis plots.')
    parser.add_argument('--res_files', type=str, help='A pathname to a text-file with pathnames under which to save modelling results.')
    parser.add_argument('--par_file', type=str, help='The pathname of the parameter input file to use.')
    parser.add_argument('--error_files', type=str, help='The pathname to a text-file with pathnames of error log files.')
    parser.add_argument('--info_files', type=str, help='The pathname to a text-file with pathnames of info log files.')
    parser.add_argument('-q', '--quiet', action='store_true', dest='quiet', help='Do not print messages to the console.')
    
    # Parse the input arguments
    args = parser.parse_args()
    
    utilities_dir = args.utilities_dir
    obs_files     = args.obs_files
    temp_files    = args.temp_files
    order_corrections = args.order_corrections
    plot_dirs     = args.plot_dirs
    res_files     = args.res_files
    par_file      = args.par_file
    error_files    = args.error_files
    info_files     = args.info_files
    quiet         = args.quiet
    
    # Import and load the utilities
    sys.path.append(os.path.abspath(utilities_dir))
    utilities_dir = utilities_dir.strip('/').split('/')[-1]
    utilities = importlib.import_module(utilities_dir)
    
    # Import and load the reduction parameters
    if par_file == None:
        module = utilities_dir + '.pyodine_parameters'
        pyodine_parameters = importlib.import_module(module)
        Pars = pyodine_parameters.Parameters()
    else:
        par_file = os.path.splitext(par_file)[0].replace('/', '.')
        pyodine_parameters = importlib.import_module(par_file)
        Pars = pyodine_parameters.Parameters()\
    
    # And run the multiprocessing observation modelling routine
    model_multi_observations(utilities, Pars, obs_files, temp_files, 
                             order_corrections=order_corrections,
                             plot_dirs=plot_dirs, res_files=res_files,
                             error_files=error_files, info_files=info_files, 
                             quiet=quiet)
