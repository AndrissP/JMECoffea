''' The script makes the comparison plots for the datasets for flavour with gluon splitting and without for Herwig and Pythia.
Run using `python plotters/plot_all_flavors.py`
'''

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
out_txt_path = 'out_txt'

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



plotvspt = True

eta_binning  = "HCalPart"  ### HCalPart, JERC, CoarseCalo, CaloTowers, Summer20Flavor, onebin;
eta_binning_str = '_'+eta_binning if eta_binning != "HCalPart" else ''
etabins = JetEtaBins(eta_binning, absolute=True)
ptbins = PtBins("MC_truth")

tag1 = '_L5_QCD-Py'+eta_binning_str
tag1_gen_gen = '_L5_QCD-Py_gen_gen'+eta_binning_str
tag1_reco_reco = '_L5_QCD-Py_reco_reco'+eta_binning_str
tag1_gen_reco = '_L5_QCD-Py_recoforalpha'+eta_binning_str
tag1_genwt = '_L5_QCD-Py_genwt'+eta_binning_str
tag1_noiso2 = '_L5_QCD-Py_noiso2'+eta_binning_str
tag2 = '_L5_QCD-MG-Py'+eta_binning_str
tag2Her = '_L5_QCD-MG-Her'+eta_binning_str
tag3 = '_L5_Pythia-TTBAR'+eta_binning_str
tag3Her = '_L5_Herwig-TTBAR'+eta_binning_str
tag4 = '_L5_DY-MG-Py'+eta_binning_str
tag4Her = '_L5_DY-MG-Her'+eta_binning_str
tag5 = '_L5_Pythia-TTBAR_ISR-FSR'+eta_binning_str


mean_name = "Median"
mean_name_std = mean_name+'Std'
closure_QCD = read_data4plot('_L5_QCD-Py'+eta_binning_str)['all'][mean_name] #read_data2(mean_name, 'all', '_L5_QCD-Py'+eta_binning_str)

######### Draw Herwig

data_to_read = {
    # label_on_plot: data_list,
    f"{ttbarlab} Pow+Her7": read_data4plot(tag3Her, closure_QCD),
    f"ZJets MG+Her7": read_data4plot(tag4Her, closure_QCD),
    f"QCD MG+Her7": read_data4plot(tag2Her, closure_QCD),
}

flavors = ['c_prompt', 'c_gluon_splitting', 'c']
for flav in flavors:
    data = {tag:data_to_read[tag][flav] for tag in data_to_read}
    data = {key:np.array([data[key][mean_name], data[key][mean_name_std], data[key]["MeanRecoPt"]]) for key in data}
    for k in etabins.get_bin_idx([0, 1.305, 2.5, 4]):
        print('Fitting subsample: ', flav)
        print('Eta: ', etabins.centres[k]) if plotvspt else print('pt: ', ptbins.centres[k])
        if plotvspt:
            if not np.any(data[list(data.keys())[0]][2][:,k]>-0.1):
                continue
        
        make_comparison_plot(data, 
                          {},
                          etabins, ptbins,
                          binidx=k, flav=flav, pt_min=15, ratio_name=f'ratio, */{ttbarlab}', inverse=False, plotvspt=plotvspt, ratio_ylim=[0.98, 1.02], ylim=[0.955,1.066])
        
flavors = ['b_prompt', 'b_gluon_splitting', 'b']
for flav in flavors:
    data = {tag:data_to_read[tag][flav] for tag in data_to_read}
    data = {key:np.array([data[key][mean_name], data[key][mean_name_std], data[key]["MeanRecoPt"]]) for key in data}
    for k in etabins.get_bin_idx([0, 1.305, 2.5, 4]):
        print('Fitting subsample: ', flav)
        print('Eta: ', etabins.centres[k]) if plotvspt else print('pt: ', ptbins.centres[k])
        if plotvspt:
            if not np.any(data[list(data.keys())[0]][2][:,k]>-0.1):
                continue
        
        make_comparison_plot(data, 
                          {},
                          etabins, ptbins,
                          binidx=k, flav=flav, pt_min=15, ratio_name=f'ratio, */{ttbarlab}', inverse=False, plotvspt=plotvspt, ratio_ylim=[0.98, 1.02], ylim=[0.87,1.075])


######### Draw Pythia

data_to_read = {
    # label_on_plot: data_list,
    f"{ttbarlab} Pow+Py8": read_data4plot(tag3, closure_QCD), #[:,:-1,:],
    f"ZJets MG+Py8": read_data4plot(tag4, closure_QCD), #, closure_QCD)), #[:,:-1,:],
    f"QCD MG+Py8": read_data4plot(tag2, closure_QCD),#         f"QCD Py8": np.array(read_data4plot(flav, '_L5_QCD-Py_leading_gen_jet')), #[:,:-1,:],
    f"QCD Py8": read_data4plot(tag1, closure_QCD), #[:,:-1,:],
}

flavors = ['c_prompt', 'c_gluon_splitting', 'c']
for flav in flavors:
    data = {tag:data_to_read[tag][flav] for tag in data_to_read}
    data = {key:np.array([data[key][mean_name], data[key][mean_name_std], data[key]["MeanRecoPt"]]) for key in data}
    for k in etabins.get_bin_idx([0, 1.305, 2.5, 4]):
        print('Fitting subsample: ', flav)
        print('Eta: ', etabins.centres[k]) if plotvspt else print('pt: ', ptbins.centres[k])
        if plotvspt:
            if not np.any(data[list(data.keys())[0]][2][:,k]>-0.1):
                continue
        
        make_comparison_plot(data, 
                          {},
                          etabins, ptbins,
                          binidx=k, flav=flav, pt_min=15, ratio_name=f'ratio, */{ttbarlab}', inverse=False, plotvspt=plotvspt, ratio_ylim=[0.98, 1.02], ylim=[0.955,1.075])
        
flavors = ['b_prompt', 'b_gluon_splitting', 'b']
for flav in flavors:
    data = {tag:data_to_read[tag][flav] for tag in data_to_read}
    data = {key:np.array([data[key][mean_name], data[key][mean_name_std], data[key]["MeanRecoPt"]]) for key in data}
    for k in etabins.get_bin_idx([0, 1.305, 2.5, 4]):
        print('Fitting subsample: ', flav)
        print('Eta: ', etabins.centres[k]) if plotvspt else print('pt: ', ptbins.centres[k])
        if plotvspt:
            if not np.any(data[list(data.keys())[0]][2][:,k]>-0.1):
                continue
        
        make_comparison_plot(data, 
                          {},
                          etabins, ptbins,
                          binidx=k, flav=flav, pt_min=15, ratio_name=f'ratio, */{ttbarlab}', inverse=False, plotvspt=plotvspt, ratio_ylim=[0.98, 1.02], ylim=[0.915,1.055])