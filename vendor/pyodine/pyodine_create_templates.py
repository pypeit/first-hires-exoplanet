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

import argparse
import importlib


def create_template(utilities, Pars, ostar_files, temp_files, temp_outname, 
                    plot_dir=None, res_files=None, obs_sum_outname=None,
                    error_log=None, info_log=None, quiet=False, live=False):
    """Create a deconvolved stellar template
    
    This routine takes a list of hot star observations, which are modelled
    with the I2 template and LSF. The best-fit parameters (mainly LSF) are then
    used to deconvolve the stellar template observation(s), resulting in a
    high-S/N, oversampled stellar template for later observation modelling of
    the same star.
    
    :param utilities: The utilities module for the instrument used in this 
        analysis.
    :type utilities: library
    :param Pars: The parameter input object to use.
    :type Pars: :class:`Template_Parameters`
    :param ostar_files: A pathname to a text-file with pathnames of hot star 
        observations for the modelling, or, alternatively, a list with the 
        pathnames.
    :type ostar_files: str or list 
    :param temp_files: A pathname to a text-file with pathnames of stellar 
        template observations to use, or, alternatively, a list with the 
        pathnames.
    :type temp_files: str or list 
    :param temp_outname: The pathname where to save the deconvolved stellar 
        template. If the directory structure does not exist yet, it will be 
        created in the process.
    :type temp_outname: str
    :param plot_dir: The directory name where to save plots and modelling 
        results. If the directory structure does not exist yet, it will be 
        created in the process. If None is given, no results/plots will be 
        saved (default).
    :type plot_dir: str or None
    :param res_files: A pathname to a text-file with pathnames under which to 
        save the results file(s), or, alternatively, a list with the 
        pathname(s). If you want to save results from multiple runs, you should 
        supply pathnames for each run. If the directory structure does not 
        exist yet, it will be created in the process. If None is given, no 
        results will be saved (default).
    :type res_files: str or list or None
    :param obs_sum_outname: A pathname under which to save the sum of the 
        template observations, both normalized and unnormalized, as fits file.
        If None is given, this data product is not saved.
    :type obs_sum_outname: str or None
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
    #print('Create logger')
    #print(logging.getLogger().hasHandlers())
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
        logging.info('Aiming to create template: {}'.format(temp_outname))
        
        ###########################################################################
        ## Set up the environment, and load all neccessary data and parameters
        ###########################################################################
        
        # Load the hot star observation(s)
        if isinstance(ostar_files, list):
            ostar_names = ostar_files
        elif isinstance(ostar_files, str):
            with open(ostar_files, 'r') as f:
                ostar_names = [l.strip() for l in f.readlines()]
        
        all_ostar_obs = []
        for filename in ostar_names:
            all_ostar_obs.append(utilities.load_pyodine.ObservationWrapper(filename))
        
        # Take the sum of the O-star spectra
        ostar = pyodine.components.SummedObservation(*all_ostar_obs)
        
        # And log the number of hot star observation(s)
        logging.info('')
        logging.info('Hot star observations ({}):'.format(len(ostar_names)))
        for ostar_name in ostar_names:
            logging.info(ostar_name)
        
        # Load the stellar spectra obtained without I2 (stellar template observations)
        if isinstance(temp_files, list):
            temp_names = temp_files
        elif isinstance(temp_files, str):
            with open(temp_files, 'r') as f:
                temp_names = [l.strip() for l in f.readlines()]
        
        all_temp_obs = []
        for filename in temp_names:
            all_temp_obs.append(utilities.load_pyodine.ObservationWrapper(filename))
        
        # Take the sum of the stellar template spectra
        temp_obs = pyodine.components.SummedObservation(*all_temp_obs)
        
        # And log the number of stellar template observation(s)
        logging.info('')
        logging.info('Stellar template observations ({}):'.format(len(temp_names)))
        for temp_name in temp_names:
            logging.info(temp_name)
        
        # And log an idea of the flux level of the template observation
        logging.info('')
        logging.info('Median flux of the summed stellar template observation: {:.0f}'.format(
                np.median([np.median(temp_obs[o].flux) for o in temp_obs.orders])))
        
        # Load the iodine atlas from file
        iod = utilities.load_pyodine.IodineTemplate(Pars.i2_to_use)
        
        # Final template output name (setup the directory structure if non-existent)
        temp_outname_dir = os.path.dirname(temp_outname)
        if not os.path.exists(temp_outname_dir) and temp_outname_dir != '':
            os.makedirs(temp_outname_dir)
        
        # Output directory for plots (setup the directory structure if non-existent)
        if isinstance(plot_dir, str):
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
        
        # Pathnames for the result files
        if isinstance(res_files, list):
            res_names = res_files
        elif isinstance(res_files, str):
            with open(res_files, 'r') as f:
                res_names = [l.strip() for l in f.readlines()]
        else:
            res_names = None
        
        
        ###########################################################################
        ## Now prepare the modelling: Choose the orders, compute weights,
        ## maybe tellurics...
        ###########################################################################
        
        # Choose orders for template generation
        if isinstance(Pars.temp_order_range, (list, tuple)) and Pars.temp_order_range[0] is not None:
            orders = np.arange(Pars.temp_order_range[0], Pars.temp_order_range[1]+1, dtype='int')
        else:
            orders = np.array(ostar.orders) 
        logging.info('')
        logging.info('Orders to use: {}'.format(orders))
        
        # Compute weights array for the combined hot star observation
        weight = ostar.compute_weight(weight_type=Pars.weight_type, rel_noise=Pars.rel_noise)
        
        # Potentially compute a bad pixel mask, using the stellar template observation
        if Pars.bad_pixel_mask is True:
            mask = pyodine.bad_pixels.BadPixelMask(temp_obs, cutoff=Pars.bad_pixel_cutoff)
            weight[np.where(mask.mask == 1.)] = 0.
        
        # Correct stellar template observation at weight == 0 if desired
        # (e.g. to correct for bad pixels)
        if Pars.correct_obs is True:
            temp_obs = pyodine.bad_pixels.correct_spec.correct_spectrum(temp_obs, weight, orders)
        
        # Initialize Normalizer with reference spectrum, then cross-correlate 
        # to obtain the velocity guess of the stellar template observation
        # (observed velocity, no barycentric correction)
        normalizer = pyodine.template.normalize.SimpleNormalizer(reference=Pars.ref_spectrum)
        
        temp_velocity = normalizer.guess_velocity(
                temp_obs[Pars.velgues_order_range[0]:Pars.velgues_order_range[1]],
                delta_v=Pars.delta_v, maxlag=Pars.maxlag)
        bary_v = temp_obs.observations[0].bary_vel_corr
        logging.info('')
        logging.info('Velocity guess: {0:.3f} km/s'.format(temp_velocity*1e-3))
        logging.info('Barycentric velocity: {0:.3f} km/s'.format(bary_v*1e-3))
        
        # Normalize the template observation (this is used as input to the
        # deconvolver later)
        norm_temp_obs = normalizer.normalize_obs(temp_obs, temp_velocity, orders=orders)
        
        # If the summed template observation should be saved, do this here
        # (set up the directory structure if non-existent yet)
        if isinstance(obs_sum_outname, str):
            obs_sum_outname_dir = os.path.dirname(obs_sum_outname)
            if not os.path.exists(obs_sum_outname_dir) and obs_sum_outname_dir != '':
                os.makedirs(obs_sum_outname_dir)
            norm_temp_obs.save_norm(obs_sum_outname)
            logging.info('')
            logging.info('Saved summed, normalized template observations to:')
            logging.info(obs_sum_outname)
        
        # Load the tellurics (if desired)
        if Pars.telluric_mask is not None:
            tellurics = pyodine.tellurics.SimpleTellurics(tell_type=Pars.telluric_mask, wave_min=Pars.tell_wave_range[0],
                                                          wave_max=Pars.tell_wave_range[1], disp=Pars.tell_dispersion)
        else:
            tellurics = None
        
        # If the hot star spectra should be normalized prior to fitting, this is
        # done here
        if Pars.normalize_chunks is True:
            for o in ostar.orders:
                ostar._flux[o] = (ostar[o].flux / ostar[o].cont / len(all_ostar_obs))
                #ostar._flux[o] = (ostar[o].flux / pyodine.template.normalize.top(ostar[o].flux, deg=3))
        
        # Now create the chunks, using the algorithm (and corresponding parameters) 
        # as defined in the parameter input file
        if Pars.chunking_algorithm == 'auto_equal_width':
            ostar_chunks = pyodine.chunks.auto_equal_width(
                    ostar, width=Pars.chunk_width, orders=orders, padding=Pars.chunk_padding, 
                    chunks_per_order=Pars.chunks_per_order, pix_offset0=Pars.pix_offset0
                    )
        elif Pars.chunking_algorithm == 'wavelength_defined':
            ostar_chunks = pyodine.chunks.wavelength_defined(
                    ostar, Pars.wavelength_dict, Pars.chunk_padding
                    )
        else:
            raise KeyError('Algorithm {} not known! (Must be one of auto_equal_width, wavelength_defined)'.format(
                    Pars.chunking_algorithm))
        
        nr_chunks_total  = len(ostar_chunks)
        nr_chunks_order0 = len(ostar_chunks.get_order(ostar_chunks.orders[0]))
        nr_orders_chunks = len(ostar_chunks.orders)
        
        logging.info('')
        logging.info('Total number of chunks: {}'.format(nr_chunks_total))
        logging.info('Nr. chunks in order 0: {}'.format(nr_chunks_order0))
        logging.info('First and last covered pixel of chunks in order 0: {}, {}'.format(
                ostar_chunks[0].abspix[0], ostar_chunks[nr_chunks_order0-1].abspix[-1]))
        logging.info('Orders: {} - {} ({} in total)'.format(
                ostar_chunks[0].order, ostar_chunks[-1].order, nr_orders_chunks))
        
        # Produce the chunk weight array
        chunk_weight = []
        for chunk in ostar_chunks:
            chunk_weight.append(weight[chunk.order, chunk.abspix[0]:chunk.abspix[-1]+1])
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
                    logging.info('Smoothing LSF...')
                    manual_redchi2 = None
                    if 'smooth_manual_redchi' in run_dict.keys() and run_dict['smooth_manual_redchi']:
                        manual_redchi2 = run_results[smooth_lsf_run]['red_chi_sq']
                    smooth_osample = None
                    if 'smooth_osample' in run_dict.keys() and run_dict['smooth_osample'] > 0:
                        smooth_osample = int(run_dict['smooth_osample'])
                    
                    lsf_smoothed = pyodine.lib.misc.smooth_lsf(
                            ostar_chunks, run_dict['smooth_pixels'], run_dict['smooth_orders'], 
                            run_dict['order_separation'], run_results[smooth_lsf_run]['results'],
                            redchi2=manual_redchi2, osample=smooth_osample)
                    logging.info('LSFs with nans: ', len(np.unique(np.argwhere(np.isnan(lsf_smoothed))[:,0])))
                    
                    LSFarr = pyodine.models.lsf.LSF_Array(lsf_smoothed, np.array([ch.order for ch in ostar_chunks]),
                                                          np.array([ch.abspix[0] for ch in ostar_chunks]))
                    
                    # Now create the full model, including the fixed LSF
                    model = pyodine.models.spectrum.SimpleModel(
                            lsf_model, wave_model, cont_model, iod, lsf_array=LSFarr,
                            osample_factor=Pars.osample_obs, conv_width=Pars.lsf_conv_width)
                        
                else:
                    raise KeyError('smooth_lsf_run not in run_results or no results for that run yet!')
            
            else:
                lsf_smoothed = None
                # If no fixed LSF model, set up a "normal" full model
                model = pyodine.models.spectrum.SimpleModel(
                        lsf_model, wave_model, cont_model, iod,
                        osample_factor=Pars.osample_obs, conv_width=Pars.lsf_conv_width)
            
            # Set up the fitter
            fitter = pyodine.fitters.lmfit_wrapper.LmfitWrapper(model)
            
            # Loop through all chunks to pre-guess the parameters
            starting_pars = []
            for i, chunk in enumerate(ostar_chunks):
                starting_pars.append(model.guess_params(chunk))
            
            # Fit wavelength guess within each order if desired (smoother input)
            # (this really only makes sense for the very first run, otherwise
            # use smoothed wavelength results from previous runs)
            if 'pre_wave_slope_deg' in run_dict.keys() and run_dict['pre_wave_slope_deg'] > 0:
                poly_pars = pyodine.lib.misc.smooth_parameters_over_orders(
                        starting_pars, 'wave_slope', ostar_chunks, deg=run_dict['pre_wave_slope_deg'])
                # Write the smoothed dispersion into starting_pars
                for i in range(nr_chunks_total):
                    starting_pars[i]['wave_slope'] = poly_pars[i]
                    
            if 'pre_wave_intercept_deg' in run_dict.keys() and run_dict['pre_wave_intercept_deg'] > 0:
                poly_pars = pyodine.lib.misc.smooth_parameters_over_orders(
                        starting_pars, 'wave_intercept', ostar_chunks, deg=run_dict['pre_wave_intercept_deg'])
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
                    ostar_chunks, chunk_weight, fitter, lmfit_params, 
                    tellurics, use_chauvenet=use_chauvenet, compute_redchi2=True, 
                    use_progressbar=Pars.use_progressbar, live=live)
            
            (run_results[run_id]['results'], run_results[run_id]['chunk_w'], 
             run_results[run_id]['fitting_failed'], chauvenet_outliers, 
             run_results[run_id]['red_chi_sq']) = modelling_return
            
            
            ###########################################################################
            ## Now we do some run diagnostics and analysis of the model results.
            ## Plotting and results saving is done only if a plot_dir was supplied.
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
            
            # Write them to file
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
                if Pars.chunking_algorithm == 'auto_equal_width':
                    nr_chunks_order = nr_chunks_order0
                    nr_orders = nr_orders_chunks
                
                logging.info('')
                logging.info('Creating analysis plots...')
                
                pipe_lib.create_analysis_plots(
                        run_results[run_id]['results'], plot_dir, run_id=run_id, 
                        tellurics=tellurics, red_chi_sq=run_results[run_id]['red_chi_sq'], 
                        nr_chunks_order=nr_chunks_order, nr_orders=nr_orders, 
                        chunk_weight=run_results[run_id]['chunk_w'], plot_chunks=plot_chunks, 
                        chunks=ostar_chunks, 
                        wave_intercept_fit=wave_intercept_fit, wave_slope_fit=wave_slope_fit,
                        plot_lsf_pars=plot_lsf_pars,
                        uncertainties_failed=uncertainties_failed,
                        nan_rchi_fit=nan_rchi_fit, chauvenet_outliers=chauvenet_outliers,
                        lsf_array=lsf_smoothed)
            
            ###########################################################################
            ## Run finished, proceeding to next run (unless all through)
            ###########################################################################
        
        
        ###########################################################################
        ## Jansson deconvolution
        ## Deconvolve orders in 'norm_temp_obs' using the O-star LSFs 
        ## and produce a StellarTemplate_Chunked object.
        ###########################################################################
        
        logging.info('')
        logging.info('----------------------')
        logging.info('All runs done, continuing to template deconvolution...')
        
        # Parameter results to use
        ostar_params = []
        for i in range(nr_chunks_total):
            ostar_params.append(run_results[Pars.jansson_run_model]['results'][i].params)
        
        # You can supply the ChunkedDeconvolver with weights for each chunk,
        # e.g. from the fitting red. Chi**2
        chunk_weights = None
        if Pars.chunk_weights_redchi is True:
            chunk_weights = 1./run_results[Pars.jansson_run_model]['red_chi_sq']
        
        # If a smoothed LSF should be used for deconvolution
        if Pars.jansson_lsf_smoothing['do_smoothing'] is True:
            logging.info('')
            logging.info('Smoothing LSF...')
            # Compute smoothed LSFs of desired run, with template oversampling
            manual_redchi2 = None
            if 'smooth_manual_redchi' in Pars.jansson_lsf_smoothing.keys() and \
                Pars.jansson_lsf_smoothing['smooth_manual_redchi']:
                manual_redchi2 = run_results[Pars.jansson_lsf_smoothing['smooth_lsf_run']]['red_chi_sq']
            
            lsf_smoothed_deconv = pyodine.lib.misc.smooth_lsf(
                    ostar_chunks, Pars.jansson_lsf_smoothing['smooth_pixels'], 
                    Pars.jansson_lsf_smoothing['smooth_orders'], Pars.jansson_lsf_smoothing['order_separation'], 
                    run_results[Pars.jansson_lsf_smoothing['smooth_lsf_run']]['results'], 
                    redchi2=manual_redchi2, osample=Pars.deconvolution_pars['osample_temp'])
            logging.info('')
            logging.info('LSFs with nans: {}'.format(len(np.unique(np.argwhere(np.isnan(lsf_smoothed_deconv))[:,0]))))
            
            LSFarr = pyodine.models.lsf.LSF_Array(
                    lsf_smoothed_deconv, np.array([ch.order for ch in ostar_chunks]), 
                    np.array([ch.abspix[0] for ch in ostar_chunks]))
            
            # Build the desired model for deconvolution
            lsf_model = Pars.model_runs[Pars.jansson_run_model]['lsf_model']
            wave_model = Pars.model_runs[Pars.jansson_run_model]['wave_model'] #pyodine.models.wave.LinearWaveModel
            cont_model = Pars.model_runs[Pars.jansson_run_model]['cont_model'] #pyodine.models.cont.LinearContinuumModel
            
            ostar_model = pyodine.models.spectrum.SimpleModel(
                    lsf_model, wave_model, cont_model, iod, lsf_array=LSFarr,
                    osample_factor=Pars.osample_obs, conv_width=Pars.lsf_conv_width)
            
            # Set up the deconvolver
            deconvolver = pyodine.template.deconvolve.ChunkedDeconvolver(
                    ostar_chunks, ostar_model, ostar_params)
            
            # Deconvolve
            template = deconvolver.deconvolve_obs(
                    norm_temp_obs, temp_velocity, bary_v, 
                    weights=chunk_weights, lsf_fixed=lsf_smoothed_deconv,
                    deconv_pars=Pars.deconvolution_pars)
        
        # No smoothed LSF
        else:
            # Build the desired model for deconvolution
            lsf_model = Pars.model_runs[Pars.jansson_run_model]['lsf_model']
            wave_model = Pars.model_runs[Pars.jansson_run_model]['wave_model'] #pyodine.models.wave.LinearWaveModel
            cont_model = Pars.model_runs[Pars.jansson_run_model]['cont_model'] #pyodine.models.cont.LinearContinuumModel
            
            ostar_model = pyodine.models.spectrum.SimpleModel(
                    lsf_model, wave_model, cont_model, iod,
                    osample_factor=Pars.osample_obs, conv_width=Pars.lsf_conv_width)
            
            # Set up the deconvolver
            deconvolver = pyodine.template.deconvolve.ChunkedDeconvolver(
                    ostar_chunks, ostar_model, ostar_params)
            
            # Deconvolve
            template = deconvolver.deconvolve_obs(
                    norm_temp_obs, temp_velocity, bary_v,
                    weights=chunk_weights, deconv_pars=Pars.deconvolution_pars)
        
        ###########################################################################
        ## Now it's done, save the template and exit
        ###########################################################################
        
        template.save(temp_outname)
        
        modelling_time = time.time() - start_t
        logging.info('')
        logging.info('Done, full working time: {}'.format(modelling_time))
    
    except Exception as e:
        logging.error('Something went wrong!', exc_info=True)



if __name__ == '__main__':
    
    # Set up the parser for input arguments
    parser = argparse.ArgumentParser(
            description='Create a deconvolved stellar template')
    
    # Required input arguments:
    # utilities_dir, ostar_files, temp_files, temp_outname, (plot_dir=None, par_file=None)
    parser.add_argument('utilities_dir', type=str, help='The pathname to the utilities directory for this instrument.')
    parser.add_argument('ostar_files', type=str, help='A pathname to a text-file with pathnames of hot star observations for the modelling.')
    parser.add_argument('temp_files', type=str, help='A pathname to a text-file with pathnames of stellar template observations to use.')
    parser.add_argument('temp_outname', type=str, help='The pathname where to save the deconvolved stellar template.')
    parser.add_argument('--plot_dir', type=str, help='The directory name where to save plots.')
    parser.add_argument('--res_files', type=str, help='A pathname to a text-file with pathnames under which to save modelling results.')
    parser.add_argument('--par_file', type=str, help='The pathname of the parameter input file to use.')
    parser.add_argument('--obs_sum_outname', type=str, help='The pathname where to save the summed, normalized template observation (not deconvolved yet).')
    parser.add_argument('--error_file', type=str, help='The pathname to the error log file.')
    parser.add_argument('--info_file', type=str, help='The pathname to the info log file.')
    parser.add_argument('-q', '--quiet', action='store_true', dest='quiet', help='Do not print messages to the console.')
    parser.add_argument('-l', '--live', action='store_true', dest='live', help='Run the code in live mode.')
    
    # Parse the input arguments
    args = parser.parse_args()
    
    utilities_dir   = args.utilities_dir
    ostar_files     = args.ostar_files
    temp_files      = args.temp_files
    temp_outname    = args.temp_outname
    plot_dir        = args.plot_dir
    res_files       = args.res_files
    par_file        = args.par_file
    obs_sum_outname = args.obs_sum_outname
    error_file      = args.error_file
    info_file       = args.info_file
    quiet           = args.quiet
    live            = args.live
    
    # Import and load the reduction parameters
    if par_file == None:
        module = utilities_dir + '.pyodine_parameters'
        pyodine_parameters = importlib.import_module(module)
        Pars = pyodine_parameters.Template_Parameters()
    else:
        par_file = os.path.splitext(par_file)[0].replace('/', '.')
        pyodine_parameters = importlib.import_module(par_file)
        Pars = pyodine_parameters.Template_Parameters()
    
    # Import and load the utilities
    sys.path.append(os.path.abspath(utilities_dir))
    utilities_dir = utilities_dir.strip('/').split('/')[-1]
    utilities = importlib.import_module(utilities_dir)
    
    # And run the template creation routine
    create_template(utilities, Pars, ostar_files, temp_files, temp_outname, 
                    plot_dir=plot_dir, res_files=res_files, 
                    obs_sum_outname=obs_sum_outname, error_log=error_file, 
                    info_log=info_file, quiet=quiet, live=live)
