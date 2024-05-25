import matplotlib.pyplot as plt
from pltStyle import pltStyle
import mplhep as hep
import matplotlib as mpl
# import numpy as np

import os
pltStyle('hep')
plt.rcParams['font.size'] = plt.rcParams['font.size']/0.98
from coffea import util
from helpers import get_xsec_dict, hist_div, sum_subhist
from fileNames.available_datasets import dataset_dictionary



HT_Py = util.load('out/draw_HT_spectrum_L5_QCD-MG-Py.coffea')
HT_Her = util.load('out/draw_HT_spectrum_L5_QCD-MG-Her.coffea')

hists_HT_merged = {}
hists_HT_theory_merged = {}
legend_labels = []

for tag in ['QCD-MG-Her', 'QCD-MG-Py']:
    xsec_dict, legend_label = get_xsec_dict(tag, dataset_dictionary)
    hists = util.load('out/Processor_HT_spectrum_L5_'+tag+'.coffea')
    # h_HT = sum([hists[key]['h_HT'] for key in hists])
    # h_HT_theory = sum([hists[key]['h_HT_theory'] for key in hists])
    # hs_HT = {key: hists[key]['h_HT'] for key in hists}
    # hs_HT_theory = {key: hists[key]['h_HT_theory'] for key in hists}

    keys = hists.keys()
    Nev = {key: hists[key]['cutflow']['all_events'].value for key in keys}
    scale_factors_theory = hist_div(xsec_dict, Nev)
    scale_factors = hist_div({key: 1 for key in keys}, Nev)
    all_histo_keys = hists[next(iter(hists.keys()))].keys()
    h_HT_theory = sum_subhist(hists, 'h_HT_theory', scale_factors_theory)
    h_HT        = sum_subhist(hists, 'h_HT', scale_factors)

    hists_HT_merged[tag] = h_HT
    hists_HT_theory_merged[tag+'theory'] = h_HT_theory
    legend_labels.append(legend_label)


# assert False

def plot_spectra(histdict, legend_labels, saveplot=True):

    samples = list(histdict.keys())
    fig = plt.figure()
    gs = fig.add_gridspec(nrows=2, ncols=1, hspace=0, height_ratios=[3, 1])
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    ### normalize the histograms
    spectra = {samp: histdict[samp]/histdict[samp].sum().value for samp in samples}

    xbins_ed = histdict[samples[0]].axes[0].edges
    bin_widths = (xbins_ed[1:]-xbins_ed[:-1])
    # legend = {legend_labels}
    colors = iter(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    for key in samples:
#         mc = next(ax._get_lines.prop_cycler)
#         colors
#         print(mc['color'])
        artist = (spectra[key]/bin_widths).plot1d(ax=ax, label=key+'blahhhh', color=next(colors), linewidth=0.95) #, markersize=1.5) #color=plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(stack)])
        artist[0].errorbar[0].set_markersize(0)
#         artist[0].stairs.set_lw(0.95)
#     mc = next(ax._get_lines.prop_cycler)
#     (pt_spectrumHer/bin_widths).plot1d(ax=ax, color = mc['color'], label='QCD Py8')
    lims = [15,5000]
    # ax.get_xlim()
    # ax.set_xticks([10, 20, 50, 100, 200, 500, 1000, 2000, 5000])
    # ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_xlim(lims)
    # breakpoint()
    ax.legend(legend_labels) #, 'QCD, Py8 weighted', 'QCD, Py8'])
    ax.set_xlim([0,6])

    denom = spectra[samples[0]].values()
    colors = iter(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    artist = (spectra[samples[0]]/denom).plot1d(ax=ax2, color=next(colors), linewidth=0.95)
    artist[0].errorbar[0].set_markersize(0)
#     artist[0].stairs.set_lw(0.95)
#     assert False
    for key in samples[1:]:
        artist = (spectra[key]/denom).plot1d(ax=ax2, color=next(colors), linewidth=0.95)
        artist[0].errorbar[0].set_markersize(0)
#         artist[0].stairs.set_lw(0.95)

    ax.set_xscale('log')
    ax2.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_ylim((0.0,2))
    ax2.set_xticks([])
    ax2.set_xticks([10, 20, 50, 100, 200, 500, 1000, 2000, 5000])
    ax2.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
        
    ax.set_xticks(ax2.get_xticks())
    ax.set_xticklabels([])
    ax.set_xlim(lims)
    ax2.set_xlim(lims)
        
    ax2.hlines(1,-10, 10000, linestyles='--',color="black", 
       linewidth=1,)
    ax2.set_ylabel("ratio")
    ax2.set_xlabel('$H_T$ (GeV)')
    ax.set_ylabel("$dN/dH_{T}$ (GeV)")
#     hep.label.exp_text(text=f'{bins.idx2plot_str(eta_idx)}, {flav} jets', loc=2, ax=ax)
    hep.cms.label("Preliminary", loc=0, data=False, ax=ax)
    #To make the minor ticks appear and be without labels even if it gets busy with major ticks.
    ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10, numticks=20))  #numticks - the maximum number of ticks. 
    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    
    if saveplot:
        if not os.path.exists("fig/pt_spectra"):
            os.mkdir("fig/pt_spectra")

        fig_name = 'fig/pt_spectra/HT_spectra_'+"_".join(samples)
        print("Saving plot with the name = ", fig_name)
        plt.savefig(fig_name+'.pdf');
        plt.savefig(fig_name+'.png');

    plt.show()
1;

plot_spectra(hists_HT_merged, legend_labels)
plot_spectra(hists_HT_theory_merged, legend_labels)