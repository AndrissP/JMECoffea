import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib as mpl
import os
from JetEtaBins import JetEtaBins, PtBins
from helpers import read_data
# from correction_fitter_helpers import save_correction_txt_file, init_vals_2014, init_two_gaus, fit_corrections

from plotters.pltStyle import pltStyle
import mplhep as hep
pltStyle(style='hep')
plt.rcParams['figure.dpi'] = 110
# plt.rcParams['image.cmap'] = 'viridis'

from data_tools import read_or_recreate_data, read_or_recreate_data_txt
from collections.abc import Iterable
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
    

flavors = ['b', 'c', 's', 'd', 'u', 'ud', 'q']
# flavors = ['b', 'bbar', 'c', 'cbar', 's', 'sbar', 'ud', 'udbar', 'q', 'qbar', 'unmatched', 'all']
plotvspt = True
eta_binning  = "HCalPart"  ### HCalPart, JERC, CoarseCalo, CaloTowers, Summer20Flavor, onebin;
eta_binning_str = '_'+eta_binning if eta_binning != "HCalPart" else ''
etabins = JetEtaBins(eta_binning, absolute=True)
ptbins = PtBins("MC_truth")

outdir = 'fig/pion_scaling'
os.makedirs(outdir, exist_ok=True)

tag1 = '_L5_not_scaled_pion'+'_split_antiflav'+eta_binning_str
tag2 = '_L5_scaled_pion'+'_split_antiflav'+eta_binning_str
# tag3 = '_L5_scaled_times2_pion'+'_split_antiflav'+eta_binning_str
tagx5 = '_L5_scaled_times5_pion'+'_split_antiflav'+eta_binning_str
tagx10 = '_L5_scaled_times10_pion'+'_split_antiflav'+eta_binning_str

for flav in flavors:
    data1 = read_data4plot(tag1)
    data2 = read_data4plot(tag2)
    # data3 = read_data4plot(tag3)
    datax5 = read_data4plot(tagx5)
    datax10 = read_data4plot(tagx10)

    data_div = data1[flav]['Median'].copy()
    data_div[data_div==0] = 1
    data_div_qbar = data1[flav+'bar']['Median'].copy()
    data_div_qbar[data_div_qbar==0] = 1
    mean_q = data2[flav]['Median']/data_div
    mean_qbar = data2[flav+'bar']['Median']/data_div_qbar

    mean_q_x5 = datax5[flav]['Median']/data_div
    mean_qbar_x5 = datax5[flav+'bar']['Median']/data_div_qbar
    mean_q_x10 = datax10[flav]['Median']/data_div
    mean_qbar_x10 = datax10[flav+'bar']['Median']/data_div_qbar
    # mean_q = data_div/data2[flav]['Median']
    # mean_qbar = data_div_qbar/data2[flav+'bar']['Median']

    uncertainty_q = np.sqrt((data2[flav]['MedianStd']/data_div)**2 + (data2[flav]['Median']*data1[flav]['MeanStd']/data_div**2)**2)
    uncertainty_qbar = np.sqrt((data2[flav+'bar']['MedianStd']/data_div_qbar)**2 + (data2[flav+'bar']['Median']*data1[flav+'bar']['MeanStd']/data_div_qbar**2)**2)
    uncertainty_q_x10 = np.sqrt((datax10[flav]['MeanStd']/data_div)**2 + (datax10[flav]['Median']*data1[flav]['MeanStd']/data_div**2)**2)
    uncertainty_qbar_x10 = np.sqrt((datax10[flav+'bar']['MeanStd']/data_div_qbar)**2 + (datax10[flav+'bar']['Median']*data1[flav+'bar']['MeanStd']/data_div_qbar**2)**2)
    uncertainty_q_x5 = np.sqrt((datax5[flav]['MeanStd']/data_div)**2 + (datax5[flav]['Median']*data1[flav]['MeanStd']/data_div**2)**2)
    uncertainty_qbar_x5 = np.sqrt((datax5[flav+'bar']['MeanStd']/data_div_qbar)**2 + (datax5[flav+'bar']['Median']*data1[flav+'bar']['MeanStd']/data_div_qbar**2)**2)
     

    uncertainty_q = uncertainty_q/100
    uncertainty_qbar = uncertainty_qbar/100
    uncertainty_q_x10 = uncertainty_q_x10/100
    uncertainty_qbar_x10 = uncertainty_qbar_x10/100
    uncertainty_q_x5 = uncertainty_q_x5/100
    uncertainty_qbar_x5 = uncertainty_qbar_x5/100

    difference = mean_q - mean_qbar
    uncertainty = np.sqrt(uncertainty_q**2 + uncertainty_qbar**2)
    difference_x10 = (mean_q_x10 - mean_qbar_x10)
    uncertainty_x10 = np.sqrt(uncertainty_q_x10**2 + uncertainty_qbar_x10**2)
    difference_x5 = (mean_q_x5 - mean_qbar_x5)
    uncertainty_x5 = np.sqrt(uncertainty_q_x5**2 + uncertainty_qbar_x5**2)

    for binidx in etabins.get_bin_idx([0, 1.305, 2.5, 4]):
        eta_string = etabins.idx2str(binidx)
        fig, ax = plt.subplots()
        # breakpoint()
        # print('flav: ', flav)
        # print(np.transpose([data1[flav]['MeanRecoPt'][:,binidx], mean_q[:, binidx], mean_qbar[:, binidx]]))
        # print('mean_q: ', mean_q[:, binidx])
        # print('mean_qbar: ', mean_qbar[:, binidx])
        ax.errorbar(data1[flav]['MeanRecoPt'][:,binidx], difference[:, binidx], yerr=uncertainty[:, binidx], fmt='o', label='x1')
        ax.errorbar(data1[flav]['MeanRecoPt'][:,binidx], difference_x10[:, binidx], yerr=uncertainty_x10[:, binidx], fmt='o', label='x10')
        ax.errorbar(data1[flav]['MeanRecoPt'][:,binidx], difference_x5[:, binidx], yerr=uncertainty_x5[:, binidx], fmt='o', label='x5')
        ax.set_xlabel('Reco Pt')
        ax.set_ylabel('$R(scaled, q)/R(central, q)$'+ f'\n'+ '$ - R(scaled, \overline{q})/R(central, \overline{q})$')
        ax.set_xscale('log')
        good_xlims = ax.get_xlim()
        ax.hlines(0,1, 10000, linestyles='--',color="black",
            linewidth=1,)
        ax.set_xticks([20, 50, 100, 500, 1000, 5000])
        ax.set_xlim(good_xlims)

        ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
        hep.label.exp_text(text=f'{etabins.idx2plot_str(binidx)}\n{flav} jets', loc=2, ax=ax)
        figname = f'{outdir}/difference_{flav}_pion_{eta_string}'
        hep.cms.label("Private work", loc=0, data=False, ax=ax, rlabel='')
        ax.legend(loc="upper right")
        [left_lim, right_lim] = ax.get_ylim()
        lim_pad = (right_lim - left_lim)/5
        ax.set_ylim(left_lim, right_lim+lim_pad)
        # inclrease the figure left margin so that the y-axis label is not cut off
        plt.subplots_adjust(left=0.23)

        plt.savefig(figname+'.png')
        plt.savefig(figname+'.pdf')
        print(f'Saved {figname}.png /.pdf')
        plt.close()
    # ax.set_title(f'{flav} quark')
