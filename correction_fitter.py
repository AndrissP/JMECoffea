#!/usr/bin/env python
    # coding: utf-8
### correction_fitter.py
### File automatically converted using ConvertJupyterToPy.ipynb from correction_fitter.ipynb
### Formatting or commets may not be preserved by the conversion.


# #### Imports

import matplotlib.pyplot as plt
import numpy as np
import os
from JetEtaBins import JetEtaBins, PtBins
from helpers import read_data
from fileNames.available_datasets import legend_labels
# from correction_fitter_helpers import save_correction_txt_file, init_vals_2014, init_two_gaus, fit_corrections

#### some newer versions of pyplot and mplhep, aren't good friends with jupyter
#### To make the plots be formatted directly well, we need to make a dummy plot and rerun the import
### (a very silly solution)
plt.figure(num=None, figsize=(2, 2), dpi=80)
plt.plot([1,2,3],[1,3,3])
import matplotlib.pyplot as plt

from data_tools import read_or_recreate_data, read_or_recreate_data_txt
from collections.abc import Iterable

from plotters.pltStyle import pltStyle
pltStyle(style='hep', font_frac=1.1)
plt.rcParams['figure.dpi'] = 110

# ### Fitting the inverse median responses

def my_mapping(flav):
    return 'a' if flav=='all' else flav

def correction_fitter(saveplots = True, do_simfit = False, do_Mikkofit = False, correction_for = 'Py', eta_binning="HCalPart", combine_antiflavour=True, figdir="fig/median_correction_fits"  ):
    '''
    correction_for = 'Py' #'Py', 'Her' or 'Py_MGQCD' , 'Her_MGQCD', standPy
    eta_binning  = "CoarseCalo"  ### HCalPart, JERC, CoarseCalo, CaloTowers, onebin, Summer20Flavor
    combine_antiflavour = True for flavor-antiflavor uncertainties
    '''
    
    if not os.path.exists(figdir):
        os.mkdir(figdir)
    
    eta_binning_str = '_'+eta_binning if eta_binning != "HCalPart" else ''    
    combine_antiflavour_txt = '_split_antiflav' if not combine_antiflavour else ''
    
    jeteta_bins = JetEtaBins(eta_binning, absolute=True)
    pt_bins = PtBins("MC_truth")
    
    # saveplots = True

    # do_simfit = False
    # do_Mikkofit = False
    
    ttbarlab = legend_labels['ttbar']['lab']
    
    pltStyle(style='hep')
    # plt.rcParams['figure.subplot.left'] = 0.162
    plt.rcParams['figure.dpi'] = 110
    
    plt.rcParams['figure.subplot.top'] = 0.93
    plt.rcParams['figure.subplot.right'] = 0.96
    # pltStyle(style='hep')
    # plt.rcParams['figure.subplot.bottom'] = 0.37
    plt.rcParams['figure.subplot.left'] = 0.175
    plt.rcParams['font.size'] = plt.rcParams['font.size']/1.04
    
    out_txt_path = 'out_txt'
    def read_data4plot(tag, closure=1, path=out_txt_path):
        '''Read the Mean, MeanStd, Median, MedianStd and RecoPt values of the data with tag `tag`.
        If closure==1, there is no clusure, otherwise it has to be of the same shape as the data read
        '''
    #     file_path = f'../out_txt/fit_results_L5_{tag}.json'
    #     with open(file_path, 'r') as json_file:
    #         json_data = json.load(json_file)
        
        data = read_or_recreate_data(tag, out_txt_path)['data']
    
        if not isinstance(closure, Iterable):
            closure_tmp = np.array([closure])
        else:
            closure_tmp = np.array(closure).copy()
            closure_tmp[closure_tmp==0] = np.nan
        
        close = ["Median", "Mean"]
        for flav in data:
            for mean_name in close:
                data[flav][mean_name] = data[flav][mean_name]/closure_tmp #[2:]
            for typeii in ["MedianStd", "MeanStd", "MeanRecoPt"]:
                data[flav][typeii] = np.array(data[flav][typeii])
    
        return data
        
    
    # from collections.abc import Iterable
    # out_txt_path = 'out_txt'
    
    # def read_data4plot(flav, tag, closure=1, path=out_txt_path):
    #     '''Read the Mean, MeanStd and RecoPt values of the data with tag `tag` and flavor `flav`.
    #     If closure==1, there is no clusure, otherwise it has to be of the same shape as the data read '''
    #     mean_name = "Median" #or 'Mean'
    #     mean_name_std = mean_name+'Std'
    #     if not isinstance(closure, Iterable):
    #         closure_tmp = np.array([closure])
    #     else:
    #         closure_tmp = np.array(closure).copy()
    #         closure_tmp[closure_tmp==0] = np.nan
    #     median = read_or_recreate_data(mean_name, flav, tag, path)/closure_tmp #[2:]
    #     median[(median==0) | (np.abs(median)==np.inf)] = np.nan
    #     medianstd = read_data(mean_name_std, flav, tag, path) #[2:]
    #     reco = read_data("MeanRecoPt", flav, tag, path)
    # #     median = 1/median
    # #     medianstd = median**2*medianstd
    #     return [median, medianstd, reco]
    
    import correction_fitter_helpers as fits
    
    # ### Individual fits
    
    def do_fits(tags, names, flavors, pion_fit=False):
        ### Put the minimum limit on the relative uncertainty to min_rel_uncert
        # min_rel_uncert = 0.0005
        min_rel_uncert_relative = 0.05 if combine_antiflavour else False
        inverse = True if combine_antiflavour else False
        
        if combine_antiflavour:
            fits2plot = { ##name: [function, initial values, # parameters]
        #              "MC truth extended": [fits.two_gaus_fnc, [0]*9, 9],
        #              "MC truth simp": [fits.response_fnc_simp, [0]*4, 4],
                     "Poly, n=4": [fits.poly4, [ 1, 1, 1, 1, 1], 5],
                     "Poly, n=3": [fits.poly3, [ 1, 1, 1, 1], 4],
        #              "Poly, n=3 + 1/x": [fits.poly3_1_x, [ 1, 1, 1, 1, 1], 5],
                     "Poly, n=3 + a/x^b": [fits.poly3_1_xn, [ 1, 1, 1, 1, 0.1, 0.1], 6],    
                       }
        else:
            fits2plot = { ##name: [function, initial values, # parameters]
             "Poly, n=3": [fits.poly3, [ 1, 1, 1, 1], 4],
             "Poly, n=4": [fits.poly4, [ 1, 1, 1, 1, 1], 5],
               }
    
            
    
        closure = read_data("Median", "all", '_L5_QCD-Py'+eta_binning_str)
    
        fit_res_all_tags = {}
        for data_tag, data_name in np.transpose([tags, names]):
            fit_res_all = {}
    
            for flav in flavors:
                fit_res = []
                ## if fit flavor-antiflavor, closure by the L5 correction
                if not combine_antiflavour:
                    flav2 = flav[:-3] if 'bar' in flav else flav
                    closure = read_data4plot(data_tag[:-len(combine_antiflavour_txt)])[flav2]["Median"]
    #                 closure = read_data("Median", flav2, data_tag[:-len(combine_antiflavour_txt)])
                    if pion_fit:
                        closure = read_data("Median", flav, '_L5_not_scaled_pion'+'_split_antiflav'+eta_binning_str)
                    
                data = {data_tag: read_data4plot(data_tag, closure=closure )[flav]}
    
                for etaidx in range(jeteta_bins.nbins):
    #             for etaidx in range(1):
                    print('Fitting subsample: ', flav, f'Eta idx = {etaidx}, eta bin = [{jeteta_bins.edges[etaidx]}; {jeteta_bins.edges[etaidx+1]} ]', )
                    fit_res_new = fits.fit_corrections(etaidx, data, flav=flav, data_tags=[data_name],
                                                         fits2plot=fits2plot, main_fit="Poly, n=4",
                                      figdir2=figdir,
                                      jeteta_bins=jeteta_bins, pt_bins=pt_bins,
                                      plot_initial_val=[0],
                                      use_recopt=True,
                                      maxlimit_static_pnt=True,
                                      min_pt_val=15,
        #                               max_ptval=4000,
                                      min_rel_uncert_relative=min_rel_uncert_relative,
    #                                   min_rel_uncert_relative=False,
                                      show_original_errorbars=True,
                                      saveplots=saveplots,
                                                       inverse=inverse,
                                                       combine_antiflavour=combine_antiflavour,
                                                      ncoefs_out=7)
    
                    fit_res = np.concatenate([fit_res, fit_res_new])
    
                num = int(5 + fit_res[2])
                fit_res = np.reshape(fit_res,((len(fit_res)//num), num))
                correction_name = my_mapping(flav)
                correction_name += 'T' if 'TTBAR' in data_tag else 'J' 
                fit_res_all[correction_name] = fit_res
            fit_res_all_tags[data_tag] = fit_res_all
            
        return fit_res_all_tags
    
    from collections import OrderedDict
    def merge_DY_QCD_TTBAR(data, tags, pt_merge=60):
        ''' Merge DY up to 60 GeV and QCD from 60 GeV. A method recommended by Mikko.
        '''
        QCD_tag = tags[np.where(['QCD' in tag for tag in tags])[0][0]]
        DY_tag = tags[np.where(['DY' in tag for tag in tags])[0][0]]
        TTBAR_tag = tags[np.where(['TTBAR' in tag for tag in tags])[0][0]]
        dataQCD = data[QCD_tag]
        dataDY = data[DY_tag]
        dataTTBAR = data[TTBAR_tag]
        merge_idx = pt_bins.get_bin_idx(pt_merge)
    
    #     dataQCD2 = {}
        for key in dataQCD:
            dataQCD[key][:merge_idx,:] = np.nan
        return OrderedDict({QCD_tag: dataQCD, TTBAR_tag: dataTTBAR, DY_tag: dataDY})
    
    def do_fits_simfit(tags, names, colors, flavors):
        ### Put the minimum limit on the relative uncertainty to min_rel_uncert
        min_rel_uncert = 0.0005
        min_rel_uncert_relative = 0.01
    
        fits2plot = { ##name: [function, initial values, # parameters]
                     "Poly, n=4": [fits.poly4, [ 1, 1, 1, 1, 1], 5],
    #                  "Poly, n=3 + a/x^b": [fits.poly3_1_xn, [ 1, 1, 1, 1, 0.1, 0.1], 6],
                     "Poly, n=3": [fits.poly3, [ 1, 1, 1, 1], 4],
                     }
    
        closure = read_data("Median", "all", '_L5_QCD-Py'+eta_binning_str)
    
        fit_res_all = {}
    
        for flav in flavors:
            fit_res = []
            if not combine_antiflavour:
                assert("Not defined in the flavor antiflavor fit scheme")
            data = {tag: read_data4plot(tag, closure=closure )[flav] for tag in tags}
            data = merge_DY_QCD_TTBAR(data, tags, pt_merge=60)
            for etaidx in range(jeteta_bins.nbins):
                print('Fitting flavor: ', flav, f'Eta idx = {etaidx}, eta bin = [{jeteta_bins.edges[etaidx]}; {jeteta_bins.edges[etaidx+1]} ]', )
                fit_res_new = fits.fit_corrections(etaidx, data, flav=flav, data_tags=names,
                                      fits2plot=fits2plot, main_fit="Poly, n=4",
                                      figdir2=figdir,
                                      jeteta_bins=jeteta_bins, pt_bins=pt_bins,
    #                                   plot_initial_val=["Mikko fun ud"],
                                      use_recopt=True,
                                      min_pt_val=17,
                                      maxlimit_static_pnt=True,
                                      min_rel_uncert_relative=min_rel_uncert_relative,
                                      show_original_errorbars=True,
                                      inflate_smallest_std_bool=True,
                                      saveplots=saveplots,
                                      colors=colors)
    
                fit_res = np.concatenate([fit_res, fit_res_new])
    
            num = int(5 + fit_res[2])
            fit_res = np.reshape(fit_res,((len(fit_res)//num), num))
            correction_name = my_mapping(flav)+'S'
            fit_res_all[correction_name] = fit_res
    
        return fit_res_all
    
    def normalize_QCD(data, tags, norm_QCD):
        QCD_tag = tags[np.where(['QCD' in tag for tag in tags])[0][0]]
        data[QCD_tag]['Median'] = data[QCD_tag]['Median']/norm_QCD
        
        return data
    
    def do_fits_simfit(tags, names, colors, flavors):
        ### Put the minimum limit on the relative uncertainty to min_rel_uncert
        min_rel_uncert = 0.0005
        min_rel_uncert_relative = 0.01
    
        fits2plot = { ##name: [function, initial values, # parameters]
                     "Poly, n=4": [fits.poly4, [ 1, 1, 1, 1, 1], 5],
    #                  "Poly, n=3 + a/x^b": [fits.poly3_1_xn, [ 1, 1, 1, 1, 0.1, 0.1], 6],
                     "Poly, n=3": [fits.poly3, [ 1, 1, 1, 1], 4],
                     }
    
        closure = read_data("Median", "all", '_L5_QCD-Py'+eta_binning_str)
    
        fit_res_all = {}
    
        for flav in ['ud']:
            fit_res = []
            if not combine_antiflavour:
                assert("Not defined in the flavor antiflavor fit scheme")
            data = {tag: read_data4plot(tag, closure=closure )[flav] for tag in tags}
            QCD_tag = tags[np.where(['QCD' in tag for tag in tags])[0][0]]
            norm_QCD = data[QCD_tag]['Median'].copy()
            norm_QCD_pt = data[QCD_tag]['MeanRecoPt'].copy()
            data = merge_DY_QCD_TTBAR(data, tags, pt_merge=60)
            for etaidx in range(jeteta_bins.nbins):
    #         for etaidx in range(2):
                print('Fitting flavor: ', flav, f'Eta idx = {etaidx}, eta bin = [{jeteta_bins.edges[etaidx]}; {jeteta_bins.edges[etaidx+1]} ]', )
                fit_res_new = fits.fit_corrections(etaidx, data, flav=flav, data_tags=names,
                                      fits2plot=fits2plot, main_fit="Poly, n=4",
                                      figdir2=figdir,
                                      jeteta_bins=jeteta_bins, pt_bins=pt_bins,
    #                                   plot_initial_val=["Mikko fun ud"],
                                      use_recopt=True,
                                      ncoefs_out=7,
                                      min_pt_val=17,
                                      maxlimit_static_pnt=True,
                                      min_rel_uncert_relative=min_rel_uncert_relative,
                                      show_original_errorbars=True,
                                      inflate_smallest_std_bool=True,
                                      saveplots=saveplots,
                                      colors=colors)
    
                fit_res = np.concatenate([fit_res, fit_res_new])
    #             breakpoint()
                norm_QCD[:-3, etaidx] = norm_QCD[:-3, etaidx]*fits.poly4(norm_QCD_pt[:-3, etaidx],*fit_res_new[-7:-2])
    #             breakpoint()
            num = int(5 + fit_res[2])
            fit_res = np.reshape(fit_res,((len(fit_res)//num), num))
            correction_name = my_mapping(flav)+'S'
            fit_res_all[correction_name] = fit_res
    
        flavors = [flav for flav in flavors if 'ud' not in flav]
        for flav in flavors:
            fit_res = []
            if not combine_antiflavour:
                assert("Not defined in the flavor antiflavor fit scheme")
            data = {tag: read_data4plot(tag, closure=closure )[flav] for tag in tags}
            data = normalize_QCD(data, tags, norm_QCD)
            for etaidx in range(jeteta_bins.nbins):
    #         for etaidx in range(2):
                print('Fitting flavor: ', flav, f'Eta idx = {etaidx}, eta bin = [{jeteta_bins.edges[etaidx]}; {jeteta_bins.edges[etaidx+1]} ]', )
                fit_res_new = fits.fit_corrections(etaidx, data, flav=flav, data_tags=names,
                                      fits2plot=fits2plot, main_fit="Poly, n=4",
                                      figdir2=figdir,
                                      jeteta_bins=jeteta_bins, pt_bins=pt_bins,
    #                                   plot_initial_val=["Mikko fun ud"],
                                      use_recopt=True,
                                      min_pt_val=17,
                                      maxlimit_static_pnt=True,
                                      min_rel_uncert_relative=min_rel_uncert_relative,
                                      show_original_errorbars=True,
                                      inflate_smallest_std_bool=True,
                                      saveplots=saveplots,
                                      colors=colors,
                                      ncoefs_out=7)
    
                fit_res = np.concatenate([fit_res, fit_res_new])
    
            num = int(5 + fit_res[2])
            fit_res = np.reshape(fit_res,((len(fit_res)//num), num))
            correction_name = my_mapping(flav)+'S'
            fit_res_all[correction_name] = fit_res
            
        return fit_res_all
    
    def merge_DY_QCD(data, tags, pt_merge=60):
        ''' Merge DY up to 60 GeV and QCD from 60 GeV. A method recommended by Mikko.
        '''
        QCD_tag = tags[np.where(['QCD' in tag for tag in tags])[0][0]]
        DY_tag = tags[np.where(['DY' in tag for tag in tags])[0][0]]
        dataQCD = data[QCD_tag]
        dataDY = data[DY_tag]
        merge_idx = pt_bins.get_bin_idx(pt_merge)
    
        data_merged = {}
        for key in dataQCD:
            data_merged[key] = np.vstack([dataDY[key][:merge_idx,:], dataQCD[key][merge_idx:,:]])
        return {QCD_tag: dataQCD, DY_tag: dataDY, "merged": data_merged}
    
    def do_fits_simfit_Mikko(tags, names, colors, flavors):
        ### Put the minimum limit on the relative uncertainty to min_rel_uncert
        min_rel_uncert = 0.0006
        min_rel_uncert_relative = 0.022
    
        ##### ud fit
        fits2plot = { ##name: [function, initial values, # parameters]
                     "Poly, n=4": [fits.poly4, [ 1, 1, 1, 1, 1], 5],
                     "Mikko fun py": [fits.Mikko_fun_ud,[0.8657 ,-0.8604 ,-0.2572 ,0.2917 ,-0.1611 ,0.8742], 6],  #[0.8657 ,-0.8603 ,-0.2572 ,0.2917 ,-0.1611 ,0.8741], 6],
                     "Mikko fun ud": [None, [0.8657 ,-0.8604 ,-0.2572 ,0.2917 ,-0.1611 ,0.8742], 6, None, fits.Mikkofun_ud], 
                     }
    
        closure = read_data("Median", "all", '_L5_QCD-Py'+eta_binning_str)
    
        fit_res_all = {}
    
        for flav in ['ud']:
            fit_res = []
            if not combine_antiflavour:
                assert("Not defined in the flavor antiflavor fit scheme")
            data = {tag: read_data4plot(tag, closure=closure )[flav] for tag in tags}
            data_merged = merge_DY_QCD(data, tags)
    
            for etaidx in range(jeteta_bins.nbins):
    
                print('Fitting flavor: ', flav, f'Eta idx = {etaidx}, eta bin = [{jeteta_bins.edges[etaidx]}; {jeteta_bins.edges[etaidx+1]} ]', )
                fit_res_new = fits.fit_corrections(etaidx, data_merged, flav=flav, data_tags=names,
                                      fits2plot=fits2plot, main_fit="Mikko fun ud",
                                      figdir2=figdir,
                                      jeteta_bins=jeteta_bins, pt_bins=pt_bins,
                                      plot_initial_val=[],
                                      ncoefs_out=6,
                                      use_recopt=True,
                                      inverse=False,
                                      min_pt_val=17,
                                      maxlimit_static_pnt=False,
                                      min_rel_uncert_relative=min_rel_uncert_relative,
                                      show_original_errorbars=False,
                                      inflate_smallest_std_bool=True,
                                      saveplots=saveplots,
                                      colors=colors,
                                      fit_sample=2)
    
                fit_res = np.concatenate([fit_res, fit_res_new])
                
    
            num = int(5 + fit_res[2])
            fit_res = np.reshape(fit_res,((len(fit_res)//num), num))
            fit_res_ud = fit_res.copy()
            correction_name = my_mapping(flav)+'M'
            fit_res_all[correction_name] = fit_res
        
        flavors = [flav for flav in flavors if 'ud' not in flav]
        merge_pt = {flav: 60 for flav in flavors}
        merge_pt['b'] = 15
        merge_pt['c'] = 15
        merge_pt['all'] = 15
        closure_ud = data_merged.copy()
            
        ##### s fit
        fits2plot = { ##name: [function, initial values, # parameters, ROOT string]
                 "Poly, n=4": [fits.poly4, [ 1, 1, 1, 1, 1], 5],
                 "Mikko fun py": [fits.Mikko_fun2,[-1,-0.7, -0.025,np.log(150.),0.3], 5],  #[0.8657 ,-0.8603 ,-0.2572 ,0.2917 ,-0.1611 ,0.8741], 6],
                 "Mikko fun s": [None, [-0.35,-0.17, -0.025,np.log(150.),0.3], 5, [None, None,None, [np.log(2),np.log(200)],None], fits.Mikkofun_s], 
                 }
        
        for flav in 's':
            fit_res = []
            if not combine_antiflavour:
                assert("Not defined in the flavor antiflavor fit scheme")
            data = {tag: read_data4plot(tag, closure=closure*closure_ud['merged']['Median'] )[flav] for tag in tags}
            data = merge_DY_QCD(data, tags, merge_pt[flav])
    
            for etaidx in range(jeteta_bins.nbins):
    
                print('Fitting flavor: ', flav, f'Eta idx = {etaidx}, eta bin = [{jeteta_bins.edges[etaidx]}; {jeteta_bins.edges[etaidx+1]} ]', )
                fit_res_new = fits.fit_corrections(etaidx, data, flav=flav, data_tags=names,
                                      fits2plot=fits2plot, main_fit="Mikko fun s",
                                      figdir2=figdir,
                                      jeteta_bins=jeteta_bins, pt_bins=pt_bins,
    #                                   plot_initial_val=["Mikko fun s"],
                                      ncoefs_out=5,
                                      use_recopt=True,
                                      inverse=False,
                                      min_pt_val=17,
                                      maxlimit_static_pnt=False,
                                      min_rel_uncert=min_rel_uncert,
                                      min_rel_uncert_relative=min_rel_uncert_relative,
                                      show_original_errorbars=True,
                                      inflate_smallest_std_bool=True,
                                      saveplots=saveplots,
                                      colors=colors,
                                      fit_sample=2,
                                      custom_jet_legend=f's/ud jets')
    
                fit_res = np.concatenate([fit_res, fit_res_new, fit_res_ud[etaidx][5:]])
    
            num = int(5 + fit_res[2]+6)
            fit_res = np.reshape(fit_res,((len(fit_res)//num), num))
            fit_res[:,2] = fit_res[:,2]+6
            correction_name = my_mapping(flav)+'M'
            fit_res_all[correction_name] = fit_res
            
        flavors = [flav for flav in flavors if 's' not in flav]
        fits2plot = { ##name: [function, initial values, # parameters, ROOT string]
                 "Poly, n=4": [fits.poly4, [ 1, 1, 1, 1, 1], 5],
                 "Mikko fun py": [fits.Mikko_fun1,[-1,-0.3,0.01,0.0001], 4],  #[0.8657 ,-0.8603 ,-0.2572 ,0.2917 ,-0.1611 ,0.8741], 6],
                 "Mikko fun": [None, [-1,-0.3,0.01,0.0001], 4, None, fits.Mikkofun], 
                 }
        
        for flav in flavors:
            fit_res = []
            if not combine_antiflavour:
                assert("Not defined in the flavor antiflavor fit scheme")
            data = {tag: read_data4plot(tag, closure=closure*closure_ud['merged']['Median'] )[flav] for tag in tags}
            data = merge_DY_QCD(data, tags, merge_pt[flav])
    
            for etaidx in range(jeteta_bins.nbins):
    
                print('Fitting flavor: ', flav, f'Eta idx = {etaidx}, eta bin = [{jeteta_bins.edges[etaidx]}; {jeteta_bins.edges[etaidx+1]} ]', )
                fit_res_new = fits.fit_corrections(etaidx, data, flav=flav, data_tags=names,
                                      fits2plot=fits2plot, main_fit="Mikko fun",
                                      figdir2=figdir,
                                      jeteta_bins=jeteta_bins, pt_bins=pt_bins,
    #                                   plot_initial_val=["Mikko fun"],
                                      ncoefs_out=4,
                                      use_recopt=True,
                                      inverse=False,
                                      min_pt_val=17,
                                      maxlimit_static_pnt=False,
                                      min_rel_uncert_relative=False,
    #                                   min_rel_uncert_relative=min_rel_uncert_relative,
                                      show_original_errorbars=True,
                                      inflate_smallest_std_bool=True,
                                      saveplots=saveplots,
                                      colors=colors,
                                      fit_sample=2,
                                      custom_jet_legend=f'{flav}/ud jets')
    
                fit_res = np.concatenate([fit_res, fit_res_new, fit_res_ud[etaidx][5:]])
    
            num = int(5 + fit_res[2]+6)
            fit_res = np.reshape(fit_res,((len(fit_res)//num), num))
            fit_res[:,2] = fit_res[:,2]+6
            correction_name = my_mapping(flav)+'M'
            fit_res_all[correction_name] = fit_res
        
        return fit_res_all
    
    from uncertainty_plotters import color_scheme
    
    unc_eta_str = '_'+eta_binning if eta_binning != "Summer20Flavor" else ''
    if correction_for == 'MGQCD_Py': #used only for the flavor uncertainties
        tags, names = ['_L5_QCD-MG-Py', '_L5_Pythia-TTBAR'], ['QCD, MG+Py8', f"{ttbarlab} Pow+Py8"]
        tags_simfit, names_simfit = ['_L5_QCD-MG-Py', '_L5_Pythia-TTBAR', '_L5_DY-MG-Py'], ['QCD, MG+Py8', f"{ttbarlab} Pow+Py8", 'DY, MG+Py8']
        colors = {names_simfit[0]: color_scheme['QCD']['color'], names_simfit[1]: color_scheme['TTBAR']['color'], names_simfit[2]: color_scheme['DY']['color']}
        tags_simfitMikko, names_simfitMikko = ['_L5_QCD-MG-Py', '_L5_DY-MG-Py'], ['QCD, MG+Py8', 'DY, MG+Py8']
        txtfile_outname = f'Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs{combine_antiflavour_txt}_MGQCD{unc_eta_str}.txt'
    elif correction_for == 'MGQCD_Her': #used only for the flavor uncertainties
        tags, names = ['_L5_QCD-MG-Her', '_L5_Herwig-TTBAR'], ['QCD, MG+Her7', f"{ttbarlab} Pow+Her7"]
        tags_simfit = ['_L5_QCD-MG-Her', '_L5_Herwig-TTBAR', '_L5_DY-MG-Her']
        names_simfit = ['QCD, MG+Her7', f"{ttbarlab} Pow+Her7", 'DY, MG+Her7']
        tags_simfitMikko = [tags_simfit[0],tags_simfit[2]]
        names_simfitMikko = [names_simfit[0],names_simfit[2]]
        colors = {names[0]: color_scheme['QCD']['color'], names[1]: color_scheme['TTBAR']['color'], 'DY, MG+Her7':color_scheme['DY']['color']}
        txtfile_outname = f'Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs{combine_antiflavour_txt}_MGQCD_Her{unc_eta_str}.txt'
    elif correction_for == 'Her':
        txtfile_outname = f'Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs{combine_antiflavour_txt}_Her{unc_eta_str}.txt'
        tags, names = ['_L5_Herwig-TTBAR', '_L5_DY-MG-Py'], [f"{ttbarlab} Pow+Her7", 'DY, MG+Her7']
        tags_simfit, names_simfit = ['_L5_Herwig-TTBAR', '_L5_DY-MG-Her'], [f"{ttbarlab} Pow+Her7", 'DY, MG+Her7']
        tags_simfitMikko, names_simfitMikko = ['_L5_QCD-MG-Her', '_L5_DY-MG-Her'], ['QCD, MG+Her7', 'DY, MG+Her7']
        colors = {names[0]: color_scheme['TTBAR']['color'], names[1]:color_scheme['DY']['color']}
    elif correction_for == 'Py': 
        ''' used only for the central flavor correction but not for the flavor uncertainties
        Uses the standalone QCD-Py sample which has no corresponging Herwig sample
        '''
        do_Mikkofit = False
        txtfile_outname = f'Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs{combine_antiflavour_txt}{unc_eta_str}.txt'
        tags, names =  ['_L5_QCD-Py', '_L5_Pythia-TTBAR'], ['QCD, Py8', f"{ttbarlab}, Pow+Py8"]
        tags_simfit, names_simfit = ['_L5_QCD-Py', '_L5_Pythia-TTBAR', '_L5_DY-MG-Py'], ['QCD, Py8', f"{ttbarlab}, Pow+Py8", 'DY, MG+Py8']
    #     tags_simfitMikko, names_simfitMikko = ['_L5_QCD-MG-Py', '_L5_DY-MG-Py'], ['QCD, MG+Py8', 'DY, MG+Py8']
        colors = {names_simfit[0]: color_scheme['QCD']['color'], names_simfit[1]: color_scheme['TTBAR']['color'], names_simfit[2]: color_scheme['DY']['color']}
    elif correction_for == 'standPy':
        do_simfit = False
        do_Mikkofit = False
        txtfile_outname = f'Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs{combine_antiflavour_txt}_standPy8{unc_eta_str}.txt'
        tags, names =  ['_L5_QCD-Py'], ['QCD, Py8']
    #     tags_simfit, names_simfit = ['_L5_QCD-Py', '_L5_Pythia-TTBAR', '_L5_DY-MG-Py'], ['QCD, Py8', f"{ttbarlab}, Pow+Py8", 'DY, MG+Py8']
    #     tags_simfitMikko, names_simfitMikko = ['_L5_QCD-MG-Py', '_L5_DY-MG-Py'], ['QCD, MG+Py8', 'DY, MG+Py8']
    #     colors = {names[0]: color_scheme['QCD']['color'], names[1]: color_scheme['TTBAR']['color'], 'DY, MG+Py8':color_scheme['DY']['color']}
    else:
        raise ValueError("Please provide correction_for from ['Py', 'Her', 'Her_MGQCD' or 'Py_MGQCD']")
        
    if combine_antiflavour == False:
        do_simfit=False
        do_Mikkofit=False
        if correction_for == 'Her':
            txtfile_outname = f'Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs{combine_antiflavour_txt}_Her{unc_eta_str}.txt'
            tags, names =  ['_L5_Herwig-TTBAR'], [f"{ttbarlab} Pow+Her7"]
        elif correction_for == 'Py':
            txtfile_outname = f'Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs{combine_antiflavour_txt}{unc_eta_str}.txt'
            tags, names =  ['_L5_Pythia-TTBAR'], [f"{ttbarlab} Pow+Py8"]
            
    tags = [tag+eta_binning_str+combine_antiflavour_txt for tag in tags]
    if do_simfit:
        tags_simfit = [tag+eta_binning_str+combine_antiflavour_txt for tag in tags_simfit]
    if do_Mikkofit:
        tags_simfitMikko = [tag+eta_binning_str+combine_antiflavour_txt for tag in tags_simfitMikko]
    
    if not combine_antiflavour:
        flavors = ['bbar', 'b', 'c', 's', 'ud', 'q', 'cbar', 'sbar', 'udbar', 'qbar']
    else:
        flavors = ['b', 'c', 'g', 's', 'ud', 'q', 'u', 'd', 'all']

    fit_res_all_tags = do_fits(tags, names, flavors) #:-1 to exclude DY
    if do_simfit:
        fit_res_all_tags['simfit'] = do_fits_simfit(tags_simfit, names_simfit, colors, flavors)
    if do_Mikkofit:
        fit_res_all_tags['Mikko'] = do_fits_simfit_Mikko(tags_simfitMikko, names_simfit+["merged"], colors=None, flavors=flavors)

    fits.save_correction_txt_file_Mikko(txtfile_outname, fit_res_all_tags)
    
    print('-----'*10)
    print("All done. Congrats!")
    
if __name__ == "__main__":
    correction_fitter()