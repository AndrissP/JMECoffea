#!/usr/bin/env python
    # coding: utf-8
### plotters/Plotting_comparison.py
### File automatically converted using ConvertJupyterToPy.ipynb from plotters/Plotting_comparison.ipynb
### No comments or formatting is preserved by the transfer.
def main():
    
    # import sys
    # top_path = '../'
    # if top_path not in sys.path:
    #     sys.path.append(top_path)
        
    # coffea_path = '/afs/cern.ch/user/a/anpotreb/top/JERC/coffea/'
    # if coffea_path not in sys.path:
        # sys.path.insert(0,coffea_path)

    import sys
    top_path = '../'
    if top_path not in sys.path:
        sys.path.append(top_path)
    
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    # import matplotlib as mpl
    from coffea import util
    from make_comparison_plot import make_comparison_plot, make_double_ratio_plot
    from coffea.lookup_tools import extractor
    from JetEtaBins import JetEtaBins, PtBins
    from collections.abc import Iterable


    from data_tools import read_or_recreate_data, read_or_recreate_data_txt
    out_txt_path = '../out_txt'

    def read_data2(mean_name, samp, tag1):
        return read_or_recreate_data_txt(mean_name, samp, tag1, out_txt_path)

    from common_binning import JERC_Constants
    
    from fileNames.available_datasets import legend_labels
    ttbarlab = legend_labels['ttbar']['lab']

    from pltStyle import pltStyle
    # from scipy.interpolate import CubicSpline
    pltStyle(style='hep')

    #### some newer versions of pyplot and mplhep, aren't good friends with jupyter
    #### To make the plots be formatted directly well, we need to make a dummy plot and rerun the import
    ### (a very silly solution)
    plt.figure(num=None, figsize=(2, 2), dpi=80)
    plt.plot([1,2,3],[1,3,3])
    import matplotlib.pyplot as plt
    pltStyle('hep', size_frac=3.0)
    # plt.rcParams['figure.dpi'] = 100

    def read_data4plot_txt(flav, tag, closure=1, path=out_txt_path, mean_name="Median"):
        '''Read the Mean, MeanStd and RecoPt values of the data with tag `tag` and flavor `flav`.
        If closure==1, there is no clusure, otherwise it has to be of the same shape as the data read
        mean_name = "Median" or 'Mean'
        '''
        mean_name_std = mean_name+'Std'
        if not isinstance(closure, Iterable):
            closure_tmp = np.array([closure])
        else:
            closure_tmp = np.array(closure).copy()
            closure_tmp[closure_tmp==0] = np.nan
        median = read_or_recreate_data_txt(mean_name, flav, tag, path=path)/closure_tmp #[2:]
        medianstd = read_or_recreate_data_txt(mean_name_std, flav, tag, path=path) #[2:]
        reco = read_or_recreate_data_txt("MeanRecoPt", flav, tag, path=path)
        return [median, medianstd, reco]
        
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
        
    # load_fit_res=True
    # subsamples = ['all', 'b', 'c', 'd', 'u', 's', 'g', 'bbar', 'cbar', 'ubar', 'dbar', 'sbar']
    # flavors = ['all'] #, 'b', 'c', 'd', 'u', 's', 'g', 'ud', 'q', 'unmatched']
    # flavors = ['b', 'g', 'ud', 'c', 's', 'all', 'q', 'unmatched', 'd', 'u']
    # flavors = ['c_prompt', 'c_gluon_splitting', 'c', 'b_prompt', 'b_gluon_splitting', 'b', 'ud', 's', 'all', 'q', 'unmatched', 'd', 'u']
    # flavors = ['b', 'c', 'g', 'ud', 's', 'all', 'q', 'unmatched', 'd', 'u']
    plotvspt = True
    flavors = ['b', 'g']

    eta_binning  = "HCalPart"  ### HCalPart, JERC, CoarseCalo, CaloTowers, Summer20Flavor, onebin;
    eta_binning_str = '_'+eta_binning if eta_binning != "HCalPart" else ''
    etabins = JetEtaBins(eta_binning, absolute=True)
    ptbins = PtBins("MC_truth")

    iso_str = ''
    tag1 = '_L5_QCD-MG-Py-evtgen'
    tag2 = '_L5_QCD-MG-Py-2000toInf'
    tag3 = '_L5_QCD-MG-Her-2000toInf'
    # tag1 = '_L5_QCD-Py'+iso_str+eta_binning_str
    # tag1_gen_gen = '_L5_QCD-Py_gen_gen'+eta_binning_str
    # tag1_reco_reco = '_L5_QCD-Py_reco_reco'+eta_binning_str
    # tag1_gen_reco = '_L5_QCD-Py_recoforalpha'+eta_binning_str
    # tag1_genwt = '_L5_QCD-Py_genwt'+eta_binning_str
    # tag1_noiso2 = '_L5_QCD-Py_noiso2'+eta_binning_str
    # tag2 = '_L5_QCD-MG-Py'+iso_str+eta_binning_str
    # tag2Her = '_L5_QCD-MG-Her'+iso_str+eta_binning_str
    # tag3 = '_L5_Pythia-TTBAR'+iso_str+eta_binning_str
    # tag3 = '_L5_Pythia-semilep-TTBAR'+eta_binning_str
    # tag3_100files = '_L5_Pythia-TTBAR_100files'+eta_binning_str
    # tag3Her = '_L5_Herwig-TTBAR'+iso_str+eta_binning_str
    # tag4noiso = '_L5_DY-MG-Py_noiso'+eta_binning_str
    # tag4 = '_L5_DY-MG-Py'+iso_str+eta_binning_str
    # tag4Her = '_L5_DY-MG-Her'+iso_str+eta_binning_str
    # tag5 = '_L5_Pythia-TTBAR_ISR-FSR'+eta_binning_str

    mean_name = "Median"
    mean_name_std = mean_name+'Std'
    # import pdb; pdb.set_trace()
    closure_QCD = read_data4plot('_L5_QCD-Py'+eta_binning_str)['all'][mean_name] #read_data2(mean_name, 'all', '_L5_QCD-Py'+eta_binning_str)
    # closure_QCD = read_data2(mean_name, 'all', '_L5_QCD-Py_iso_cut'+eta_binning_str)
    # closure_QCD_iso = read_data2(mean_name, 'all', tag26)
    # closure_TTBAR = read_data2(mean_name, 'all', '_L5_QCD-Py'+eta_binning_str)
    # closure_tag29 = read_data2(mean_name, 'all', tag29)
    # closure_tag30 = read_data2(mean_name, 'all', tag30)
    # closure_Py2 = read_data2(mean_name, 'all', tag3)
    # closure_Py3 = read_data2(mean_name, 'all', tag5)
    # closure_Py = read_data(mean_name, 'all', tag3)

    data_to_read = {
        # label_on_plot: data_list,
        r"QCD MG+Py8": read_data4plot(tag2, closure_QCD), #, closure_Py3),
        r"QCD MG+Her7": read_data4plot(tag3, closure_QCD), #, closure_Py3),
        r"QCD MG+Py8, evtGen": read_data4plot(tag1, closure_QCD),
    #     f"{ttbarlab} Pow+Py8": read_data4plot(tag3, closure_QCD), #[:,:-1,:],
    #     f"ZJets MG+Py8": read_data4plot(tag4, closure_QCD), #, closure_QCD)), #[:,:-1,:],
    #     f"QCD MG+Py8": read_data4plot(tag2, closure_QCD),#         f"QCD Py8": np.array(read_data4plot(flav, '_L5_QCD-Py_leading_gen_jet')), #[:,:-1,:],
    #     f"QCD Py8": read_data4plot(tag1, closure_QCD), #[:,:-1,:],
    #     f"{ttbarlab} Pow+Her7": read_data4plot(tag3Her, closure_QCD), #[:,:-1,:],
    #     f"QCD Py8, gen-gen": read_data4plot(tag1_gen_gen), #[:,:-1,:],
    #     f"QCD Py8, reco-reco": read_data4plot(tag1_reco_reco), #[:,:-1,:],
    #     f"QCD Py8, reco-gen": read_data4plot(tag1_gen_reco), #[:,:-1,:],
    #     f"QCD Py8, gen-reco": read_data4plot(tag1_noiso2), #[:,:-1,:],
    #     f"{ttbarlab} Pow+Her7": read_data4plot(tag3Her, closure_QCD), #[:,:-1,:],
    #     f"ZJets MG+Py8, noiso": read_data4plot(tag4noiso, closure_QCD), #, closure_QCD)), #[:,:-1,:],
        # f"ZJets MG+Py8": read_data4plot(tag4, closure_QCD), #, closure_QCD)), #[:,:-1,:],
        # f"ZJets MG+Her7": read_data4plot(tag4Her, closure_QCD), #, closure_QCD)), #[:,:-1,:],
    #     f"QCD MG+Her7": read_data4plot(tag2Her, closure_QCD),#         f"QCD Py8": np.array(read_data4plot(flav, '_L5_QCD-Py_leading_gen_jet')), #[:,:-1,:],

        #     f"QCD MG+Her7": read_data4plot(tag2Her, closure_QCD),#         f"QCD Py8": np.array(read_data4plot(flav, '_L5_QCD-Py_leading_gen_jet')), #[:,:-1,:],

    }


    for flav in flavors:
        data = {tag:data_to_read[tag][flav] for tag in data_to_read}
        data = {key:np.array([data[key][mean_name], data[key][mean_name_std], data[key]["MeanRecoPt"]]) for key in data}
    #     data = {
    #         # label_on_plot: data_list,
    #         f"{ttbarlab} Pow+Py8": np.array(read_data4plot(flav, tag3, closure_QCD)), #[:,:-1,:],
    #         f"ZJets MG+Py8": np.array(read_data4plot(flav, tag4, closure_QCD)), #, closure_QCD)), #[:,:-1,:],
    #         f"QCD MG+Py8": np.array(read_data4plot(flav, tag2, closure_QCD)),#         f"QCD Py8": np.array(read_data4plot(flav, '_L5_QCD-Py_leading_gen_jet')), #[:,:-1,:],
    #         f"QCD Py8": np.array(read_data4plot(flav, tag1, closure_QCD)), #[:,:-1,:],
    # #         f"{ttbarlab} Pow+Her7": np.array(read_data4plot(flav, tag3Her, closure_QCD)), #[:,:-1,:],
    # #         f"ZJets MG+Her7": np.array(read_data4plot(flav, tag4Her, closure_QCD)), #, closure_QCD)), #[:,:-1,:],
    # #         f"QCD MG+Her7": np.array(read_data4plot(flav, tag2Her, closure_QCD)),#         f"QCD Py8": np.array(read_data4plot(flav, '_L5_QCD-Py_leading_gen_jet')), #[:,:-1,:],
    # #         f"{ttbarlab} Pow+Py8+ISR-FSR": np.array(read_data4plot(flav, tag5, closure_QCD)), #[:,:-1,:],
    #     }

    #     for k in range(len(etabins_abs)-1):
        for k in etabins.get_bin_idx([0, 1.305, 2.5]):
        # for k in etabins.get_bin_idx([0, 1.305, 2.5, 4]):
    #     for k in ptbins.get_bin_idx([20, 35, 150, 400]):
            print('Fitting subsample: ', flav)
            print('Eta: ', etabins.centres[k]) if plotvspt else print('pt: ', ptbins.centres[k])
            if plotvspt:
                if not np.any(data[list(data.keys())[0]][2][:,k]>-0.1):
                    continue
            
            make_comparison_plot(data, 
                            {},
                            etabins, ptbins,
                            binidx=k, flav=flav, pt_min=15, ratio_name='ratio,\nHer7/Py8', inverse=False, plotvspt=plotvspt, ratio_ylim=[0.985, 1.015], extra_text='')

    
    # etabins_abs = etabins[(len(etabins)-1)//2:]
    # etabins_c = (etabins_abs[:-1]+etabins_abs[1:])/2 #output['ptresponse'].axis('jeteta').centers()
    
    # import matplotlib as mpl
    
    # # (read_data2("Median", "b", '_L5_QCD-MG-Py')[0:,0])/read_data2("Median", "b", '_L5_QCD-MG-Her')[0:,0]
    
    # color_scheme = {key: cycler_vals
    #     for cycler_vals, key in zip(plt.rcParams['axes.prop_cycle'], ['g', 'q', 'c', 'b', 'QCD', 'DY', 'TTBAR', 'DY200'])}
    
    # leggend_dict = {'g': 'Gluons', 'q': 'Quarks', 'b': 'Bottom', 'c': 'Charm'}
    
    # # 1/closure
    
    # def read_data4plot(flav, tag, closure=1):
    #     '''Read the Mean, MeanStd and RecoPt values of the data with tag `tag` and flavor `flav`.
    #     If closure==1, there is no clusure, otherwise it has to be of the same shape as the data read '''
    #     mean_name = "Median" #or 'Mean'
    #     mean_name_std = mean_name+'Std'
    #     median = read_data2(mean_name, flav, tag)/closure #[2:]
    #     medianstd = read_data2(mean_name_std, flav, tag) #[2:]
    #     reco = read_data2("MeanRecoPt", flav, tag)
    #     return [median, medianstd, reco]
        
    
    # # np.array(read_data4plot(flav, tag1, closure_Py))[:,:-1,:].shape
    # # np.array(read_data4plot(flav, tag3, closure_Py2)).shape
    
    # # load_fit_res=True
    # # subsamples = ['all', 'b', 'c', 'd', 'u', 's', 'g', 'bbar', 'cbar', 'ubar', 'dbar', 'sbar']
    # flavors = ['b', 'g'] #, 'b', 'c', 'd', 'u', 's', 'g', 'ud', 'q', 'unmatched']
    # # flavors = ['all', 'all_unmatched', 'b', 'c', 'd', 'u', 's', 'g', 'ud', 'q', 'unmatched']
    
    # # subsamples = ['unmatched']
    # tag1 = '_L5_QCD-MG-Py-evtgen'
    # tag2 = '_L5_QCD-MG-Py-2000toInf'
    # tag3 = '_L5_QCD-MG-Her-2000toInf'
    # # tag4 = '_L5_QCD-MG-Her_alphacut_0p2'
    # # tag5 = '_L5_QCD-MG-Py_alphacut_0p2_gen15'
    # # tag6 = '_L5_QCD-MG-Her_alphacut_0p2_gen15'
    # # tag7 = '_L5_QCD-MG-Her_alphacut_0p2_promptlep_gen15'
    # # tag3 = '_L5_QCD-JME'
    # # tag1 = '_L5_Pythia-TTBAR'
    # # tag2 = '_L5_Herwig-TTBAR'
    # # tag1 = '_L5_DY-MG-Py'
    # # tag2 = '_L5_DY-MG-Her'
    # # tag3 = '_L5_DY-JME'
    
    # mean_name = "Median"
    # mean_name_std = mean_name+'Std'
    # # closure_Py = read_data2(mean_name, 'all', tag1)
    # # closure_Py2 = read_data2(mean_name, 'all', tag3)
    # # closure_Py3 = read_data2(mean_name, 'all', tag5)
    # # closure_Py = read_data(mean_name, 'all', tag3)
    
    # for flav in flavors:
    #     data = {
    #         # label_on_plot: data_list,
    # #         "QCD MG+Py8": np.array(read_data4plot(flav, tag1, closure_Py))[:,:-1,:],
    # #         "QCD MG+Her7": np.array(read_data4plot(flav, tag2))[:,:-1,:],
    # #         r"QCD MG+Py8, $\alpha<0.2$": read_data4plot(flav, tag3),
    #         r"QCD MG+Py8, evtGen": read_data4plot(flav, tag1),
    #         r"QCD MG+Py8": read_data4plot(flav, tag2), #, closure_Py3),
    #         r"QCD MG+Her7": read_data4plot(flav, tag3), #, closure_Py3),
    # #         r"QCD MG+Her7, alpha<0.2 ": read_data4plot(flav, tag7), #, closure_Py3),
    # #, $\alpha<0.4$, gen15
    #         #, $\alpha<0.2$, promtlep, gen15
    #         #         "TTBAR MG+Her7": [median_2, medianstd_2, reco_pt2],
    # #         "DY MGFxFx+Py8": [median_3, medianstd_2, reco_pt2],
    # #         "QCD MG+Her7": [median_2, medianstd_2, reco_pt2],
    
    #        }
    
    #     for k in range(len(etabins_abs)-1):
    #         print('Fitting subsample: ', flav)
    #         print('Eta: ', k)
    #         if not np.any(data[list(data.keys())[0]][2][:,k]>-0.1):
    #             continue
            
    #         make_comparison_plot(data, 
    #                           {},
    #                           etabins_abs, ptbins[:27],
    #                           etaidx=k, flav=flav, ratio_name='*/ \n Py8')
            
    # 3;
    
    # len(ptbins)
    
    # # load_fit_res=True
    # # subsamples = ['all', 'b', 'c', 'd', 'u', 's', 'g', 'bbar', 'cbar', 'ubar', 'dbar', 'sbar']
    # subsamples = ['all', 'b', 'c', 'd', 'u', 's', 'g', 'ud', 'q', 'untagged']
    # # subsamples = ['ud', 'all']
    # # subsamples = ['all', 'b']
    # # CoffeaJERCOutputs_L5_DY-MG-Py.coffea
    # tag1 = '_L5_QCD-JME'
    # tag2 = '_L5_Herwig-QCD'
    # tag3 = '_L5_QCD_MG_Py8'
    # # tag4 = '_L5_QCD-JME-weights-leading3jets'
    # # tag5 = '_L5_Herwig-QCD'
    # # tag3 = '_L5_QCD-divided'
    
    # mean_name = "Median"
    # mean_name_std = mean_name+'Std'
    
    # # closure_Py = read_data(mean_name, 'q', tag1)*0.8+read_data(mean_name, 'g', tag1)*0.2
    # # closure_Py[( read_data(mean_name, 'q', tag1)==0 ) | (read_data(mean_name, 'g', tag1)==0)] = np.nan
    # # closure_Her = read_data(mean_name, 'q', tag2)*0.8 + read_data(mean_name, 'g', tag2)*0.2
    # # closure_Her[( read_data(mean_name, 'q', tag2)==0 ) | (read_data(mean_name, 'g', tag2)==0)] = np.nan
    # closure_Py = read_data(mean_name, 'all', tag1)
    # # closure_3 = read_data(mean_name, 'all', tag3)
    # # closure_4 = read_data(mean_name, 'all', tag4)
    # # closure_Her = read_data(mean_name, 'all', tag5)
    
    # # k2 = np.where(etabins_mod<=0)[0][-1]
    # # k4 = np.where(etabins_mod<=1.3)[0][-1]
    # # k6 = np.where(etabins_mod<=2.5)[0][-1]
    # # k8 = np.where(etabins_mod<=3.0)[0][-1]
    # # ks = [k2, k4, k6, k8]
    # # closure_Aut18 = evaluator['Autumn18_V3_MC_Pythia8_all_L2Relative_AK4PFchs']
    # # closure_Aut18_Her = evaluator['Autumn18_V3_MC_Herwig7_all_L2Relative_AK4PFchs']
    
    # # def closure_Aut18(pt,eta):
    # #     return 0.8*evaluator['Autumn18_V3_MC_Pythia8_ud_L2Relative_AK4PFchs'](pt, eta) + 0.2*evaluator['Autumn18_V3_MC_Pythia8_g_L2Relative_AK4PFchs'](pt, eta)
    # # def closure_Aut18_Her(pt,eta):
    # #     return 0.8*evaluator['Autumn18_V3_MC_Herwig7_ud_L2Relative_AK4PFchs'](pt, eta) + 0.2*evaluator['Autumn18_V3_MC_Herwig7_g_L2Relative_AK4PFchs'](pt, eta)
    
    # # ks = [k2, k4] #, k6, k8]
    # for samp in subsamples:
    #     samp_Aut18 = samp
    # #     samp_Sum20 = '_'+samp
    #     samp_Aut18 = '_ud' if samp_Aut18=='u' or samp_Aut18=='d' else '_'+samp_Aut18
    #     median_base = read_data(mean_name, samp, tag1)/closure_Py #[2:]
    #     medianstd_base = read_data(mean_name_std, samp, tag1) #[2:]
    #     reco_pt = read_data("MeanRecoPt", samp, tag1)
    #     median_2 = read_data(mean_name, samp, tag2)/closure_Py
    #     medianstd_2 = read_data(mean_name_std, samp, tag2)
    #     reco_pt2 = read_data("MeanRecoPt", samp, tag2)
    #     median_3 = read_data(mean_name, samp, tag3)/closure_Py
    #     medianstd_3 = read_data(mean_name_std, samp, tag3)
    #     reco_pt3 = read_data("MeanRecoPt", samp, tag3)
    # #     median_3 = read_data(mean_name, samp, tag3)/closure_3
    # #     medianstd_3 = read_data(mean_name_std, samp, tag3)
    # #     reco_pt3 = read_data("MeanRecoPt", samp, tag3)
    # #     median_4 = read_data(mean_name, samp, tag4)/closure_4
    # #     medianstd_4 = read_data(mean_name_std, samp, tag4)
    # #     reco_pt4 = read_data("MeanRecoPt", samp, tag4)
    # #     median_5 = read_data(mean_name, samp, tag5)/closure_1
    # #     medianstd_5 = read_data(mean_name_std, samp, tag5)
    # #     reco_pt5 = read_data("MeanRecoPt", samp, tag5)
    #     data = {
    #         "Py8": [median_base, medianstd_base, reco_pt],
    #         "MG+Her7": [median_2, medianstd_2, reco_pt2],
    #         "MG+Py8": [median_3, medianstd_3, reco_pt3],
    # #         "divided": [median_3, medianstd_3, ptbins_c],
    #        }
        
    # #     evo_Her = evaluator[f'Autumn18_V3_MC_Herwig7{samp_Aut18}_L2Relative_AK4PFchs']
    # #     evo = evaluator[f'Autumn18_V3_MC_Pythia8{samp_Aut18}_L2Relative_AK4PFchs']
        
    # #     functions = {
    # #             "Autumn18_Py":    [evo, closure_Aut18],
    # #             "Autumn18_Her":   [evo_Her, closure_Aut18_Her],
    # # #             "Autumn18_Her_Her":   [evo_Her, closure_Aut18_Her],
                
    # #            }
    
    #     for k in range(len(etabins_abs)-1):
    #         print('Fitting subsample: ', samp)
    #         print('Eta: ', k)
    #         if not np.any(median_base[:,k]>-0.1):
    #             continue
            
    #         make_comparison_plot(data, 
    #                           {},
    #                           etaidx=k, flav=samp, ratio_name='*/ \n Py8')
            
    # 3;
    
    # corr_loc_Sum20_Py = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Py_fineeta.txt"]
    # corr_loc_Sum20_Her = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Her_fineeta.txt"]
    # unc = ["* * ../Summer19UL18_V5_MC/Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.junc.txt"]
    
    # ext = extractor()
    # ext.add_weight_sets(corr_loc_Sum20_Py+corr_loc_Sum20_Her+unc)
    # ext.finalize()
    # evaluator2 = ext.make_evaluator()
    
    # # corr_loc_Sum20_Py = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Py.txt"]
    # # corr_loc_Sum20_Her = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Her.txt"]
    
    # corr_loc_Sum20_Py = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Py_pt30to500.txt"]
    # corr_loc_Sum20_Her = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Her_pt30to500.txt"]
    # # corr_loc_Sum20_Py = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Py_fineeta.txt"]
    # # corr_loc_Sum20_Her = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Her_fineeta.txt"]
    # # unc = ["* * ../Summer19UL18_V5_MC/Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.junc.txt"]
    
    # ext = extractor()
    # ext.add_weight_sets(corr_loc_Sum20_Py+corr_loc_Sum20_Her)
    # ext.finalize()
    # evaluator = ext.make_evaluator()
    
    # # load_fit_res=True
    # # flavoes = ['all', 'b', 'c', 'd', 'u', 's', 'g', 'bbar', 'cbar', 'ubar', 'dbar', 'sbar']
    # flavoes = ['b', 'c', 'd', 'u', 's', 'g', 'ud', 'q']
    
    # tag1 = '_L5_QCD-MG-Py'
    # tag2 = '_L5_QCD-MG-Her'
    # # tag1 = '_L5_Pythia-TTBAR'
    # # tag2 = '_L5_Herwig-TTBAR'
    
    # mean_name = "Median"
    # mean_name_std = mean_name+'Std'
    
    # closure_Py = read_data2(mean_name, 'all', tag1)
    
    # # ks = [k2, k4] #, k6, k8]
    # for flav in flavoes:
    #     data = {
    #         "QCD, MG+Py8": read_data4plot(flav, tag1, closure_Py), #[median_base, medianstd_base, reco_pt],
    #         "QCD, MG+Her7": read_data4plot(flav, tag2, closure_Py),
    # #         "MG+Py8": [median_3, medianstd_3, reco_pt3],
    # #         "divided": [median_3, medianstd_3, ptbins_c],
    #        }
        
    #     flav2 = 'a' if flav=='all' else flav
    #     evo_Her = evaluator[f'Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Her_pt30to500_{flav2}J']
    #     evo_Py = evaluator[f'Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Py_pt30to500_{flav2}J']
        
    #     functions = {
    #             "Summer20_Py":    [evo_Py, None],
    #             "Summer20_Her":   [evo_Her, None],
    #                        }
    
    #     for k in range(len(etabins_abs)-1):
    #         print('Fitting subsample: ', flav)
    #         print('Eta: ', k)
    #         median_base = data[list(data.keys())[0]]
    #         if not np.any(median_base[2][:,k]>-0.1):
    #             continue
            
    #         make_comparison_plot(data, 
    #                           functions,
    #                              etabins_abs, ptbins,
    #                           etaidx=k, flav=flav, ratio_name='*/ \n Py8')
            
    # 3;
    
    # # load_fit_res=True
    # # subsamples = ['all', 'b', 'c', 'd', 'u', 's', 'g', 'bbar', 'cbar', 'ubar', 'dbar', 'sbar']
    # subsamples = ['all', 'b', 'c', 'd', 'u', 's', 'g', 'ud', 'q']
    
    # # tag1 = '_L5_QCD-MG-Py'
    # # tag2 = '_L5_QCD-MG-Her'
    # # tag1 = '_L5_Pythia-TTBAR'
    # # tag2 = '_L5_Herwig-TTBAR'
    # # tag1 = '_L5_QCD-JME'
    # # tag2 = '_L5_QCD-JME-leading3jets'
    # tag1 = '_L5_QCD-MG-Py'
    # tag2 = '_L5_QCD-MG-Py_leading_jets'
    # tag3 = '_L5_QCD-MG-Her'
    # tag4 = '_L5_QCD-MG-Her_leading_jets'
    
    # mean_name = "Median"
    # mean_name_std = mean_name+'Std'
    
    # closure_Py = read_data2(mean_name, 'all', tag1)
    
    # # ks = [k2, k4] #, k6, k8]
    # for samp in subsamples:
    # #     median_base = read_data2(mean_name, samp, tag1) #/closure_Py #[2:]
    # #     medianstd_base = read_data2(mean_name_std, samp, tag1) #[2:]
    # #     reco_pt = read_data2("MeanRecoPt", samp, tag1)
    # #     median_2 = read_data2(mean_name, samp, tag2) #/closure_Py
    # #     medianstd_2 = read_data2(mean_name_std, samp, tag2)
    # #     reco_pt2 = read_data2("MeanRecoPt", samp, tag2)
    
    #     data = {
    # #         "QCD, Py8": read_data4plot(samp, tag1), 
    # #         "QCD, Py8, 3lead jets": read_data4plot(samp, tag2),
    #         "QCD, MG+Py8": np.array(read_data4plot(samp, tag1, closure_Py))[:,:-1,:], 
    #         "QCD, MG+Py8, 3lead jets": read_data4plot(samp, tag2, closure_Py[:-1,:]),
    #         "QCD, MG+Her7": np.array(read_data4plot(samp, tag3, closure_Py))[:,:-1,:], 
    #         "QCD, MG+Her7, 3lead jets": read_data4plot(samp, tag4, closure_Py[:-1,:]),
    # #         "TTBAR, Pow+Py8": [median_base, medianstd_base, reco_pt],
    # #         "TTBAR, Pow+Her7": [median_2, medianstd_2, reco_pt2],
    # #         "MG+Py8": [median_3, medianstd_3, reco_pt3],
    # #         "divided": [median_3, medianstd_3, ptbins_c],
    #        }
        
    # #     samp2 = 'a' if samp=='all' else samp
    # #     evo_Her = evaluator[f'Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Her_{samp2}T']
    # #     evo_Py = evaluator[f'Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Py_{samp2}T']
        
    # #     functions = {
    # #             "fit Pow+Py8":    [evo_Py, None],
    # #             "fit Pow+Her7":   [evo_Her, None],
    # #                        }
    
    #     for k in range(len(etabins_abs)-1):
    #         print('Fitting subsample: ', samp)
    #         print('Eta: ', k)
    #         median_base = data[list(data.keys())[0]]
    #         if not np.any(median_base[2][:,k]>-0.1):
    #             continue
            
    #         make_comparison_plot(data, 
    #                           {}, #functions,
    #                           etabins_abs, ptbins[:-1],
    #                           etaidx=k, flav=samp, ratio_name='*/ \n Py8')
            
    # 3;
    
    # # ### Comparison with old corrections
    
    # Aut18_samples = ['all', 'b', 'c', 's', 'ud', 'g' ]
    # Sum16_samples = ['b', 'c', 's', 'ud', 'g' ]
    
    # list_corr_Sum16 = ["Summer16_07Aug2017_V15_Flavor_Pythia8_MC_"+samp+"_L5Flavor_AK4PFchs.txt" for samp in Sum16_samples]
    # list_corr_Sum16.append("Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L2Relative_AK4PFchs.txt")
    # list_corr_Sum16.append("Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L3Absolute_AK4PFchs.txt")
    # corr_loc_Sum16 = ["* * ../Summer16_07Aug2017_V15_Flavor_Pythia8_MC/"+corr for corr in list_corr_Sum16]
    # list_corr_Aut18 = ["Autumn18_V3_MC_Pythia8_"+samp+"_L2Relative_AK4PFchs.txt" for samp in Aut18_samples]
    # corr_loc_Aut18 = ["* * ../Autumn18_V3_MC_Pythia8/"+corr for corr in list_corr_Aut18]
    # corr_loc_Winter14 = ["* * ../Winter14_V8_MC_L5Flavor/Winter14_V8_MC_L5Flavor_AK5PFchs.txt"]
    # corr_loc_Sum20 = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Py-etaAut18.txt"]
    
    # list_corr_Aut18 = ["Autumn18_V3_MC_Herwig7_"+samp+"_L2Relative_AK4PFchs.txt" for samp in Aut18_samples]
    # corr_loc_Aut18_Her = ["* * ../Autumn18_V3_MC_Pythia8/"+corr for corr in list_corr_Aut18]
    # # corr_loc_Sum20_Her = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Her-etaAut18.txt"]
    
    # # corr_loc_Sum16+corr_loc_Aut18+corr_loc_Winter14+corr_loc_Sum20+corr_loc_Aut18_Her+corr_loc_Sum20_Her
    
    # ext = extractor()
    # ext.add_weight_sets(corr_loc_Sum16+corr_loc_Aut18+corr_loc_Winter14+corr_loc_Sum20+corr_loc_Aut18_Her) #+corr_loc_Sum20_Her)
    # # ext.add_weight_sets(corr_loc_Winter14)
    # ext.finalize()
    
    # # ext._names
    # evaluator = ext.make_evaluator()
    # # evo = evaluator['Autumn18_V3_MC_Pythia8_b_L2Relative_AK4PFchs']
    # # evaluator['Autumn18_V3_MC_Herwig7_b_L2Relative_AK4PFchs']
    # # evaluator['Autumn18_V3_MC_Herwig7_b_L2Relative_AK4PFchs']
    
    # # plt.rcParams['figure.subplot.top'] = 0.93
    # # plt.rcParams['figure.subplot.right'] = 0.99
    # # pltStyle(style='hep')
    # # plt.rcParams['figure.subplot.bottom'] = 0.37
    # # plt.rcParams['figure.subplot.left'] = 0.50
    # # plt.rcParams['figure.figsize'] = [10,20]
    # plt.rcParams['font.size'] = plt.rcParams['font.size']/0.98
    
    # # load_fit_res=True
    # # subsamples = ['all', 'b', 'c', 'd', 'u', 's', 'g']
    # flavors = ['b', 'c', 'd', 'u', 's', 'g', 'ud', 'q']
    # # flavors = ['all']
    # etabins = np.array(JERC_Constants.etaBinsEdges_Aut18_full())
    # etabins_abs = etabins[(len(etabins)-1)//2:]
    

    
if __name__ == "__main__":
    main()