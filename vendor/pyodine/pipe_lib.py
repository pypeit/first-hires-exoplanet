#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 16:25:21 2021

@author: pheeren
"""

import numpy as np
import pyodine
import os
import matplotlib.pyplot as plt

from progressbar import ProgressBar
import logging
import sys


def model_all_chunks(chunks, chunk_weight, fitter, lmfit_params, 
                     tellurics=None, use_chauvenet=True, compute_redchi2=True, 
                     use_progressbar=True, live=False):
    """Loop over all chunks and model them
    
    :params chunks: The chunks of the observation to model.
    :type chunks: :class:`ChunkArray`
    :param chunk_weight: The weights array for the chunks.
    :type chunk_weight: list #ndarray[nr_chunks,nr_pix]
    :param fitter: The fitter instance to use.
    :type fitter: :class:`LmfitWrapper`
    :param lmfit_params: A list of :class:`lmfit.Parameters` objects for the 
        chunks.
    :type lmfit_params: list[:class:`lmfit.Parameters`]
    :param tellurics: The tellurics to use. If None, they are not used.
    :type tellurics: :class:`SimpleTellurics`, or None
    :param use_chauvenet: Whether to use Chauvenet criterion in the modelling. 
        Defaults to True.
    :type use_chauvenet: bool
    :param compute_redchi2: Whether to manually compute red. Chi**2 values for 
        the chunks. Defaults to True.
    :type compute_redchi2: bool
    :param use_progressbar: Whether to show a progressbar during the modelling. 
        Defaults to True.
    :type use_progressbar: bool
    :param live: If True, then the modelling is performed in live-mode, i.e.
        each modelled chunk is plotted.
    :type live: bool
    
    :return: The best-fit results of the modelled chunks.
    :rtype: list[:class:`LmfitResult`]
    :return: An array with updated chunk weights.
    :rtype: list #ndarray[nr_chunks,nr_pix]
    :return: A list of chunk indices where the fitting failed.
    :rtype: list
    :return: A list of chunk indices where the Chauvenet criterion took effect.
    :rtype: list
    :return: An array of manually computed red. Chi**2 values (or zeros if 
        keyword **compute_redchi2** == False).
    :rtype: ndarray[nr_chunks]
    
    """
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')
    
    chauvenet_outliers = []
    fitting_failed = []
    results = []
    if compute_redchi2:
        red_chi_sq = np.zeros((len(chunks)))
    else:
        red_chi_sq = None
    chunk_w = [] #np.ones((chunk_weight.shape))
    
    # Use a progressbar?
    if use_progressbar:
        bar = ProgressBar(max_value=len(chunks), redirect_stdout=True)
        bar.update(0)
    
    for i, chunk in enumerate(chunks):
        
        # These are the nominal pixel weights of this chunk
        ch_w = chunk_weight[i].copy()
        # Check which chunk pixels are potentially affected by tellurics
        # Weights of these pixels are set to zero
        if tellurics:
            waves = fitter.model.wave_model.eval(
                    chunk.pix, fitter.convert_params(lmfit_params[i], from_lmfit=True).filter(prefix='wave'))
            ind = np.where(tellurics.is_affected(waves))
            ch_w[ind[0]] = 0.
        
        # Find bad pixels, set weights there to 0.
        # (compare to dop_pre)
        bad_pix = np.where(ch_w <= 0.) #& (ch_w > 5. * np.median(ch_w)))  
        if len(bad_pix[0]) > 0:
            ch_w[bad_pix[0]] = 0.
        
        # Now model the chunk
        try:
            result = fitter.fit(chunk, lmfit_params[i], weight=ch_w, chunk_ind=i)
            if result.lmfit_result is not None and use_chauvenet is True:
                mask, mask_true, mask_false = pyodine.lib.misc.chauvenet_criterion(result.residuals)
                if any(mask)==False:
                    logging.info('Fit again, reject outliers... (chunk {})'.format(i))
                    chauvenet_outliers.append([i, chunk])
                    ch_w[mask_false] = 0.
                    # Fit again
                    result = fitter.fit(chunk, lmfit_params[i], weight=ch_w, chunk_ind=i)
        except Exception as e:
            logging.warning('Chunk {}:'.format(i))
            logging.warning(e)
            fitting_failed.append(i)
            result = fitter.LmfitResult(chunk, fitter.model, None, chunk_ind=i)
        results.append(result)
        
        # Compute manual red. Chi2
        if compute_redchi2:
            if result.lmfit_result is not None:
                rchi2 = np.sum(ch_w * result.residuals**2) / \
                        (len(chunk) - len(bad_pix[0]) - result.lmfit_result.nvarys) # divided by degrees of freedom
                if rchi2 < 0.:
                    rchi2 = 0.
            else:
                rchi2 = 0.
            red_chi_sq[i] = rchi2**0.5
        chunk_w.append(ch_w) #[i] = ch_w
        
        # If live, plot the chunk and print fit report
        if live:
            logging.info(result.report)
            if i==0:
                live_fig, live_ax = None, None
            try:
                live_fig, live_ax = pyodine.plot_lib.live_chunkmodel(
                        result, chunks, i, tellurics=tellurics, weight=ch_w, 
                        fig=live_fig, ax=live_ax)
            except Exception as e:
                logging.warning(e)
        
        # Update the progressbar
        if use_progressbar:
            bar.update(i+1)
    
    if use_progressbar:
        bar.finish()
    
    return results, chunk_w, fitting_failed, chauvenet_outliers, red_chi_sq


def create_analysis_plots(fit_results, save_dir, run_id=None, tellurics=None, 
                          red_chi_sq=None, nr_chunks_order=None, nr_orders=None, 
                          chunk_weight=None, plot_chunks=None, chunks=None, 
                          wave_intercept_fit=None, wave_slope_fit=None, 
                          plot_lsf_pars=False, uncertainties_failed=None,
                          nan_rchi_fit=None, chauvenet_outliers=None,
                          lsf_array=None, live=False):
    """Create analysis plots for a modelling run
    
    :param fit_results: A list containing instances of :class:`LmfitResult`
        for each chunk.
    :type fit_results: list[:class:`LmfitResult`]
    :param save_dir: The directory name where to save plots.
    :type save_dir: str
    :param run_id: The run ID (to include in plot titles and savenames, 
        if supplied).
    :type run_id: int, or None
    :param tellurics: An instance of tellurics. If None, they are not included 
        (default).
    :type tellurics: :class:`SimpleTellurics` or None
    :param red_chi_sq: An array of red. Chi**2 values to plot along with the 
        fit-results values. If None, they are not plotted.
    :type red_chi_sq: ndarray[nr_chunks], or None
    :param nr_chunks_order: Number of chunks per order. If this and nr_orders 
        is given, the order borders are plotted.
    :type nr_chunks_order: int, or None
    :param nr_orders: Number of orders. If this and nr_chunks_orders is given, 
        the order borders are plotted.
    :type nr_orders: int, or None
    :param chunk_weight: The weights array for the chunks.
    :type chunk_weight: ndarray[nr_chunks,nr_pix]
    :param plot_chunks: Which chunks should be plotted. Defaults to None.
    :type plot_chunks: int, list, tuple, ndarray, or None
    :param chunks: The chunks of the observation to model.
    :type chunks: :class:`ChunkArray`
    :param wave_intercept_fit: Fitted wave intercepts for each chunk. If None, 
        the corresponding plots are not created.
    :type wave_intercept_fit: ndarray[nr_chunks], or None
    :param wave_slope_fit: Fitted wave slopes for each chunk. If None, the 
        corresponding plot is not created.
    :type wave_slope_fit: ndarray[nr_chunks], or None
    :param plot_lsf_pars: Whether to plot the LSF parameters. Defaults to False.
    :type plot_lsf_pars: bool
    :param uncertainties_failed: Chunk indices where the uncertainties were not 
        computed. If None, no fit success plot is created.
    :type uncertainties_failed: list, tuple, ndarray, or None
    :param nan_rchi_fit: Chunk indices where the fitting failed. If None, no 
        fit success plot is created.
    :type nan_rchi_fit: list, tuple, ndarray, or None
    :param chauvenet_outliers: Chunk indices where the Chauvernet criterion 
        took effect. If None, no fit success plot is created.
    :type chauvenet_outliers: list, tuple, ndarray, or None
    :param lsf_array: An array of evaluated LSFs to plot some exemplary ones. 
        Defaults to None.
    :type lsf_array: ndarray[nr_chunks, nr_pix], or None
    :param live: If True, then the modelling is performed in live-mode, i.e.
        each modelled chunk is plotted.
    :type live: bool
    """
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')
    
    # I put everything in try - except, so that it does not destroy 
    #everything else
    try:
        # Extract the parameters
        param_names = [k for k in fit_results[0].params.keys()]
        params = pyodine.models.base.ParameterSet()
        for p in param_names:
            params[p] = np.array([r.params[p] for r in fit_results])
        
        # Check if the run_id is given (for the titles and savenames)
        run_title_str = ''
        run_save_str = ''
        if isinstance(run_id, (str, int)):
            run_title_str = ', run {}'.format(run_id)
            run_save_str = 'r{}_'.format(run_id)
        
        # Plot residuals histogramm
        hist_result = pyodine.plot_lib.plot_residual_hist(
                fit_results, residual_arr=None, tellurics=tellurics, 
                title='Histogram of residuals' + run_title_str, 
                savename=os.path.join(save_dir, run_save_str + 'residuals_hist.png'), 
                dpi=300, show_plot=live)
        
        # Plot chunk residuals
        if tellurics:
            (all_res, sub_res) = hist_result
            
            pyodine.plot_lib.plot_chunk_scatter(
                    scatter=sub_res, scatter_alpha=0.7,
                    nr_chunks_order=nr_chunks_order, nr_orders=nr_orders, 
                    ylabel='Chunk residuals [%]', title='Chunk residuals (outside tellurics)' + run_title_str, 
                    savename=os.path.join(save_dir, run_save_str + 'residuals.png'), 
                    dpi=300, show_plot=live)
        
        else:
            all_res = hist_result
            
            pyodine.plot_lib.plot_chunk_scatter(
                    scatter=all_res, scatter_alpha=0.7,
                    nr_chunks_order=nr_chunks_order, nr_orders=nr_orders, 
                    ylabel='Chunk residuals [%]', title='Chunk residuals' + run_title_str, 
                    savename=os.path.join(save_dir, run_save_str + 'residuals.png'), 
                    dpi=300, show_plot=live)
        
        # Plot red. Chi**2
        fit_red_chi2 = [r.redchi for r in fit_results]
        scatter_label = 'Fit red. Chi**2 med={:.2f}'.format(np.nanmedian(fit_red_chi2))
        if isinstance(red_chi_sq, (list,tuple,np.ndarray)):
            fit_red_chi2 = [fit_red_chi2, red_chi_sq]
            scatter_label = [scatter_label,
                             'Manual red. Chi**2 med={:.2f}'.format(np.nanmedian(red_chi_sq))]
        pyodine.plot_lib.plot_chunk_scatter(
                scatter=fit_red_chi2, scatter_alpha=0.7, scatter_label=scatter_label,
                nr_chunks_order=nr_chunks_order, nr_orders=nr_orders, 
                ylabel='log(red. Chi**2)', ylog=True, title='Chunk red. Chi**2' + run_title_str, 
                savename=os.path.join(save_dir, run_save_str + 'redchi.png'), 
                dpi=300, show_plot=live)
        
        # Plot some evaluated chunks
        if isinstance(plot_chunks, (list,tuple,np.ndarray,int)) and chunks:
            if isinstance(plot_chunks, int):
                plot_chunks = [plot_chunks]
            if fit_results[0].model.stellar_template:
                template = True
            else:
                template = False
            for chunk_id in plot_chunks:
                pyodine.plot_lib.plot_chunkmodel(
                        fit_results, chunks, chunk_id, template=template, 
                        tellurics=tellurics, show_plot=live, 
                        savename=os.path.join(save_dir, run_save_str + 'chunk{}.png'.format(chunk_id)), 
                        weight=chunk_weight[chunk_id])
        
        # Plot the wavelength slopes and fits
        if isinstance(wave_slope_fit, (np.ndarray, list, tuple)):
            pyodine.plot_lib.plot_chunk_scatter(
                    scatter=params['wave_slope'], scatter_alpha=0.7, scatter_label='Fit results',
                    curve=wave_slope_fit, curve_alpha=0.7, curve_label='Polynomial fit',
                    ylabel=r'Wave slope [$\AA$/pix]', title='Wave slope results' + run_title_str, 
                    savename=os.path.join(save_dir, run_save_str + 'wave_slope.png'), 
                    dpi=300, show_plot=live)
        
        # Plot the wavelength intercepts and fits, and the residuals between both
        if isinstance(wave_intercept_fit, (np.ndarray, list, tuple)):
            pyodine.plot_lib.plot_chunk_scatter(
                    scatter=params['wave_intercept'], scatter_alpha=0.7, scatter_label='Fit results',
                    curve=wave_intercept_fit, curve_alpha=0.7, curve_label='Polynomial fit',
                    ylabel=r'Wave intercept [$\AA$]', title='Wave intercept results' + run_title_str, 
                    savename=os.path.join(save_dir, run_save_str + 'wave_intercept.png'), 
                    dpi=300, show_plot=live)
            
            pyodine.plot_lib.plot_chunk_scatter(
                    scatter=params['wave_intercept'] - wave_intercept_fit, scatter_alpha=0.7,
                    ylabel=r'Residuals [$\AA$]', title='Wave intercepts - fit' + run_title_str, 
                    savename=os.path.join(save_dir, run_save_str + 'wave_residuals.png'), 
                    dpi=300, show_plot=live)
        
        # Plot the lsf parameters
        if plot_lsf_pars:
            lsf_parnames = fit_results[0].params.filter(prefix='lsf')
            for p in lsf_parnames:
                k = 'lsf_' + p
                pyodine.plot_lib.plot_chunk_scatter(
                    scatter=params[k], scatter_fmt='.', scatter_alpha=0.7, grid=False,
                    nr_chunks_order=nr_chunks_order, nr_orders=nr_orders, 
                    title='{}'.format(k) + run_title_str, 
                    savename=os.path.join(save_dir, run_save_str + '{}.png'.format(k)), 
                    dpi=300, show_plot=live)
        
        # Plot where fitting went wrong for chunks
        if (isinstance(uncertainties_failed, (list,tuple,np.ndarray)) and
            isinstance(nan_rchi_fit, (list,tuple,np.ndarray)) and
            isinstance(chauvenet_outliers, (list,tuple,np.ndarray))):
            
            fig = plt.figure(figsize=(12,6))
            plt.vlines(uncertainties_failed, ymin=0.05, ymax=0.3, 
                       color='r', label='No fit uncertainties')
            plt.vlines(nan_rchi_fit, ymin=0.35, ymax=0.6, 
                       color='g', label='NaN fitted red. Chi**2')
            if len(chauvenet_outliers) != 0:
                plt.vlines(chauvenet_outliers[:][0], ymin=0.65, ymax=0.9, 
                           color='k', label='Chauvenet outliers')
            if nr_chunks_order and nr_orders:
                plt.vlines([nr_chunks_order*o for o in range(nr_orders)],
                            ymin=0, ymax=1, color='k', alpha=0.5, linestyle=':')
            plt.legend()
            plt.xlabel('Chunk #')
            plt.title('Fitting success' + run_title_str)
            
            plt.savefig(os.path.join(save_dir, run_save_str + 'fitting_success.png'),
                        format='png', dpi=300)
            if live:
                plt.show()
            plt.close()
        
        # Plot exemplary smoothed LSFs in 9 sectors
        if isinstance(lsf_array, (np.ndarray,list,tuple)):
            pyodine.plot_lib.plot_lsfs_grid(
                    lsf_array, chunks, x_nr=3, y_nr=3, 
                    alpha=0.7, xlim=(12,36), grid=True, 
                    savename=os.path.join(save_dir, run_save_str + 'smoothed_lsfs.png'), 
                    dpi=300, show_plot=live)
    
    except Exception as e:
        logging.error('Run results analysis failed...', exc_info=True)


def velocity_results_analysis(run_result, save_dir, nr_chunks_order, 
                              nr_orders, obs_filename):
    """Perform a short analysis of velocity results
    
    :param run_result: Then results dictionary of the final modelling run.
    :type run_result: dict
    :param save_dir: The directory name where to save plots.
    :type save_dir: str
    :param nr_chunks_order: Number of chunks per order. If this and nr_orders 
        is given, the order borders are plotted.
    :type nr_chunks_order: int, or None
    :param nr_orders: Number of orders. If this and nr_chunks_orders is given, 
        the order borders are plotted.
    :type nr_orders: int, or None
    :param obs_filename: The filename of the modelled observation (to use in
        the plots' savenames.)
    :type obs_filename: str
    
    """
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')
    
    # I put everything in try - except, so that it does not destroy 
    #everything else
    try:
        # Process and print overview of results
        obs_params = []
        obs_errors = []
        for r in run_result['results']:
            obs_params.append(r.params)
            obs_errors.append(r.errors)
        
        velocities = np.array([p['velocity'] for p in obs_params], dtype=np.float64)
        vel_errors = np.array([e['velocity'] for e in obs_errors], dtype=np.float64)
        
        velocity_robust_mean = pyodine.timeseries.misc.robust_mean(velocities)
        velocity_robust_std = pyodine.timeseries.misc.robust_std(velocities)
        vel_error_robust_mean = pyodine.timeseries.misc.robust_mean(vel_errors)
        vel_error_robust_std = pyodine.timeseries.misc.robust_std(vel_errors)
        
        logging.info('Velocity robust mean: {}'.format(velocity_robust_mean))
        logging.info('Velocity robust std: {}'.format(velocity_robust_std))
        logging.info('')
        
        logging.info('Velocity error robust mean: {}'.format(vel_error_robust_mean))
        logging.info('Velocity error robust std: {}'.format(vel_error_robust_std))
        logging.info('')
        
        try:
            logging.info('Minimum error: {} in chunk {}'.format(np.nanmin(vel_errors), np.nanargmin(vel_errors)))
            logging.info('Maximum error: {} in chunk {}'.format(np.nanmax(vel_errors), np.nanargmax(vel_errors)))
        except Exception as e:
            logging.error(e)
        logging.info('')
        
        logging.info('Chunks with nan velocities: {}'.format(np.where(np.isnan(velocities))[0]))
        logging.info('Chunks with nan velocity errors: {}'.format(np.where(np.isnan(vel_errors))[0]))
        
        #velocities_finite = velocities[np.where(np.isfinite(vel_errors))]
        vel_errors_finite = vel_errors[np.where(np.isfinite(vel_errors))]
        
        # Plot velocities with errors
        errorbar_label = 'Rob. mean = {:.3f} km/s\nRob. std = {:.3f} km/s'.format(
                velocity_robust_mean*1e-3, velocity_robust_std*1e-3)
        savename = os.path.join(save_dir, 'vel_' + os.path.splitext(os.path.basename(obs_filename))[0] + '.png')
        pyodine.plot_lib.plot_chunk_scatter(
                errorbar=velocities*1e-3, errorbar_yerr=vel_errors*1e-3, 
                errorbar_fmt='o', errorbar_alpha=0.7, errorbar_label=errorbar_label,
                hlines=[velocity_robust_mean*1e-3,
                        (velocity_robust_mean+velocity_robust_std)*1e-3,
                        (velocity_robust_mean-velocity_robust_std)*1e-3],
                hlines_fmt=':', hlines_alpha=0.7, hlines_color=['r','g','g'],
                yrange=((velocity_robust_mean-5*velocity_robust_std)*1e-3, 
                        (velocity_robust_mean+5*velocity_robust_std)*1e-3),
                ylabel='Velocity [km/s]', title='Velocities', grid=True,
                savename=savename, dpi=300, show_plot=False
                )
        
        # Plot histogram of velocity errors
        percent = 80
        if len(vel_errors_finite) != 0:
            percentiles = np.percentile(vel_errors_finite, percent)
            logging.info('Velocity error at percentile {}: {}'.format(percent, percentiles))
            
            maxbin = 1000.
            
            fig = plt.figure(figsize=(12,6))
            plt.hist(vel_errors_finite, bins=np.arange(0,maxbin,20), alpha=0.75)#, range=(0,2000));
            plt.xlabel('Velocity error [m/s]')
            #plt.xlim((0,400.))
            plt.axvline(x=percentiles, color='r')
            plt.text(x=percentiles+10, y=60, s='{}%: {:.1f} m/s'.format(percent, percentiles), color='r', fontsize=14)
            plt.legend(['Rob. mean = {:5.1f} m/s\nRob. std = {:5.1f} m/s'.format(vel_error_robust_mean, 
                                                                      vel_error_robust_std)])
            plt.savefig(os.path.join(save_dir, 'velerr_hist_' + os.path.splitext(os.path.basename(obs_filename))[0] + '.png'),
                        format='png', dpi=300)
            plt.close()
        
        # Plot velocities phased over orders
        if isinstance(nr_orders, int):
            fig = plt.figure(figsize=(12,6))
            for i in range(nr_orders):
                plt.plot(np.arange(nr_chunks_order), velocities[i*nr_chunks_order:(i+1)*nr_chunks_order]*1e-3, 
                         #yerr=vel_errors[i*nr_chunks_order:(i+1)*nr_chunks_order]*1e-3, fmt='.-', alpha=0.5)
                         '-', alpha=0.5)
            
            plt.xlabel('Chunk # within order')
            plt.ylabel('Velocity [km/s]')
            plt.minorticks_on()
            plt.ylim((velocity_robust_mean-5.*velocity_robust_std)*1e-3, 
                     (velocity_robust_mean+5.*velocity_robust_std)*1e-3)
            plt.title('Velocities within orders')
            plt.savefig(os.path.join(save_dir, 'vel_orders_{}.png'.format(os.path.splitext(os.path.basename(obs_filename))[0])),
                                     format='png', dpi=300)
            plt.close()
        
        # Plot distribution of velocities
        fig = plt.figure(figsize=(12,6))
        plt.hist(velocities*1e-3, bins=200, alpha=0.75)#, range=(0,2000));
        plt.xlabel('Velocity [km/s]')
        #plt.xlim((0,400.))
        plt.axvline(x=velocity_robust_mean*1e-3, color='r', alpha=0.5)
        plt.axvline(x=(velocity_robust_mean-velocity_robust_std)*1e-3, color='g', alpha=0.5)
        plt.axvline(x=(velocity_robust_mean+velocity_robust_std)*1e-3, color='g', alpha=0.5)
        plt.xlim((velocity_robust_mean-5.*velocity_robust_std)*1e-3, 
                 (velocity_robust_mean+5.*velocity_robust_std)*1e-3)
        plt.legend(['Rob. mean = {:.3f} km/s\nRob. std = {:.3f} km/s'.format(velocity_robust_mean*1e-3, velocity_robust_std*1e-3)])
        plt.savefig(os.path.join(save_dir, 'vel_hist_{}.png'.format(os.path.splitext(os.path.basename(obs_filename))[0])),
                                 format='png', dpi=300)
        plt.close()
        
        # Wavelengths (residuals to fit) phase-folded over orders 
        if isinstance(nr_orders, int):
            fig = plt.figure(figsize=(12,6))
            for i in range(nr_orders):
                plt.plot(np.arange(nr_chunks_order), 
                         [obs_params[j]['wave_intercept'] for j in range(i*nr_chunks_order,(i+1)*nr_chunks_order)] - \
                         run_result['wave_intercept_fit'][i*nr_chunks_order:(i+1)*nr_chunks_order],
                         '-', alpha=0.5)
            plt.ylim(-0.02,0.02)
            plt.xlabel('Chunk # within order')
            plt.ylabel(r'Wavelength residuals [$\AA$]')
            plt.minorticks_on()
            plt.title('Wavelength residuals within orders')
            plt.savefig(os.path.join(save_dir, 'wave_res_orders_{}.png'.format(os.path.splitext(os.path.basename(obs_filename))[0])),
                                     format='png', dpi=300)
            plt.close()
    
    except Exception as e:
        logging.error('Velocity results analysis failed...', exc_info=True)
    