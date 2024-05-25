'''
Run as `python plotters/draw_spectra.py`.
Requires the coffea output files with tags:
['QCD-Py_genwt_test', 'QCD-MG-Py', 'QCD-MG-Her', 'QCD-Py']
'''

from fileNames.available_datasets import legend_labels

include_unmatched = True
from uncertainty_helpers import get_output, sum_output, combine_flavors
# from uncertainty_plotters import plot_spectra
from helpers import get_xsec_dict, sum_neg_pos_eta, rebin_hist, hist_mult
from fileNames.available_datasets import dataset_dictionary
from JetEtaBins import JetEtaBins, PtBins
from plotters.pltStyle import pltStyle
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import mplhep as hep
pltStyle('hep')
plt.rcParams['figure.subplot.right'] = 0.955
plt.rcParams['figure.subplot.left'] = 0.18

combine_antiflavour = True
eta_binning  = "Summer20Flavor"  ### HCalPart, CoarseCalo, CaloTowers, one_bin, Summer20Flavor;       
                         ### HCalPart: bin in HCal sectors, CaloTowers: the standard JERC binning,
                         ### CoarseCalo: like 'CaloTowers' but many bins united;

jeteta_bins = JetEtaBins(eta_binning, absolute=True)
pt_bins = PtBins("MC_truth")
eta_binning_str = '_'+eta_binning if eta_binning != "HCalPart" else ''

eta_idx = jeteta_bins.get_bin_idx(0)

# Her_legends = ['QCD MG+Her7', 'DY MG+Her7', legend_labels["ttbar"]["lab"]+' Pow+Her7']
Py_legends = ['QCD MG+Py8', 'DY MG+Py8', legend_labels["ttbar"]["lab"]+' Pow+Py8']
# samples = ['QCD', 'DY', 'TTBAR']
# Her_samples = ['_QCD-MG-Py_3rd_jet', '_DY-MG-Her', '_Herwig-TTBAR']
# Py_samples = ['_QCD-Py_3rd_jet', '_DY-MG-Py', '_Pythia-TTBAR'] #_sel_67ac6c3

qfrac_dict = {}
qfrac_var_dict = {}
qfrac_spline_dict = {}
# qfrac_spline_dict2 = {}
qfrac_spline_dict2D = {}
legend_labels = []
flavors = ['g', 'c', 'b', 'ud', 's', 'unmatched'] if include_unmatched==True else ['g', 'c', 'b', 'ud', 's']
# if not combine_antiflavour:
    # flavors = get_flavor_antiflavor_list(flavors)
flavors_to_obtain = flavors+['all'] if include_unmatched==True else flavors
    
saveplot = include_unmatched if combine_antiflavour else True

sample_plot = {}
hists_rebinned_dict = {}
Neffs = {}

def plot_spectra(histdict, labels, flav, etaidx, jeteta_bins, ptbins, saveplot=True, plotvspt=True):
    samples = list(histdict.keys())
    xbins = ptbins if plotvspt else jeteta_bins
    bins = jeteta_bins if plotvspt else ptbins
    xbins_c = xbins.centres
    xbins_ed = xbins.edges
    
    fig = plt.figure()
    gs = fig.add_gridspec(nrows=2, ncols=1, hspace=0, height_ratios=[3, 1])
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    Neff = {samp: histdict[samp][flav].sum().value**2/(histdict[samp][flav].sum().variance) for samp in samples}
    if plotvspt:
        spectra = {samp: histdict[samp][flav][:,sum] for samp in samples}
    else:
        spectra = {samp: histdict[samp][flav][sum,:] for samp in samples}
    spectra['QCD-Py_weights'] = spectra['QCD-Py_genwt_test']
    # for key in ['QCD-Py_genwt', 'QCD-Py']:
#     pt_spectrumPy = histsPy[flav][:,sum]
#     pt_spectrumHer = histsHer[flav][:,sum]
#     ed = pt_spectrum.axes[0].edges
#     centres = pt_spectrum.axes[0].centers
    bin_widths = (xbins_ed[1:]-xbins_ed[:-1])
    colors = iter(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    for key in samples:
#         mc = next(ax._get_lines.prop_cycler)
#         colors
#         print(mc['color'])
        artist = (spectra[key]/bin_widths).plot1d(ax=ax, label=key, color=next(colors), linewidth=0.95) #, markersize=1.5) #color=plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(stack)])
        artist[0].errorbar[0].set_markersize(2.5)
#         artist[0].stairs.set_lw(0.95)
#     mc = next(ax._get_lines.prop_cycler)
#     (pt_spectrumHer/bin_widths).plot1d(ax=ax, color = mc['color'], label='QCD Py8')
    # lims = ax.get_xlim()
    # lims = [np.min(centres), np.max(centres)]
    if plotvspt:
        lims = [15,5000]
    else:
        lims = [-0.2,5.3]
    # ax.get_xlim()
    # ax.set_xticks([10, 20, 50, 100, 200, 500, 1000, 2000, 5000])
    # ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_xlim(lims)
    ax.legend(labels)
    ax.set_xlim([0,6])

    denom = spectra[samples[0]].values()
    colors = iter(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    artist = (spectra[samples[0]]/denom).plot1d(ax=ax2, color=next(colors), linewidth=0.95)
    artist[0].errorbar[0].set_markersize(2.5)
#     artist[0].stairs.set_lw(0.95)
#     assert False
    for key in samples[1:]:
        artist = (spectra[key]/denom).plot1d(ax=ax2, color=next(colors), linewidth=0.95)
        artist[0].errorbar[0].set_markersize(2.5)
#         artist[0].stairs.set_lw(0.95)
    if plotvspt:
        ax.set_xscale('log')
        ax2.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_ylim((0.0,2))
    if plotvspt:
        ax2.set_xticks([])
        ax2.set_xticks([10, 20, 50, 100, 200, 500, 1000, 2000, 5000])
        ax2.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
        
    ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10, numticks=30))  #numticks - the maximum number of ticks. 
    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_xticks(ax2.get_xticks())
    ax.set_xticklabels([])
    ax.set_xlim(lims)
    ax2.set_xlim(lims)
        
    ax2.hlines(1,-10, 10000, linestyles='--',color="black", 
       linewidth=1,)
    ax2.set_ylabel("ratio")
    ax2.set_xlabel('$p_{T,ptcl}$ (GeV)') if plotvspt else ax2.set_xlabel('$|\eta|$') 
    ax.set_ylabel("$dN/dp_{T,ptcl}$ (GeV)") if plotvspt else ax.set_ylabel("$dN/d|\eta|$")
#     hep.label.exp_text(text=f'{bins.idx2plot_str(eta_idx)}, {flav} jets', loc=2, ax=ax)
    hep.cms.label("Preliminary", loc=0, data=False, ax=ax)
    
    if saveplot:
        if not os.path.exists("fig/pt_spectra"):
            os.mkdir("fig/pt_spectra")

        fig_name = 'fig/pt_spectra/pt_spectra_'+"_".join(samples)
        print("Saving plot with the name = ", fig_name)
        plt.savefig(fig_name+'.pdf')
        plt.savefig(fig_name+'.png')

    plt.show()
    # fig.close()


for sample in ['QCD-Py_genwt_test', 'QCD-MG-Py', 'QCD-MG-Her', 'QCD-Py']:
#         print(sample)
    output = get_output(sample)
    file_dict, legend_label = get_xsec_dict(sample, dataset_dictionary)
    output = sum_output(output, sample, file_dict)
    hists = combine_flavors(output, flavors_to_obtain, sumeta=False, combine_antiflavour=combine_antiflavour)
    Neffs[sample] = output['sum_weights']['sum_weights'].value
    hists_rebinned = {flav: rebin_hist(sum_neg_pos_eta(hists[flav]), 'jeteta', jeteta_bins.edges) for flav in hists.keys() }
    if include_unmatched==False: #recalculate the all
        hists_rebinned['all'] = sum([hists_rebinned[flav] for flav in flavors])
    hists_rebinned_dict[sample] = hists_rebinned
    # hists_vals = {}
    # for flav in hists_rebinned.keys():
    #     vals = hists_rebinned[flav].values().copy()
    #     vals[vals==0] = np.nan
    #     hists_vals[flav] = vals
    # qfracs = {flav: hists_vals[flav]/hists_vals['all'] for flav in flavors}
#         qfracs_var = {flav: hists_rebinned[flav].variances()/hists_rebinned['all'].variances() for flav in flavors}
    # qfrac_var_all = hists_rebinned['all'].variances()/hists_vals['all']**2
    # qfrac_var = {}
    # for flav in flavors:
    #     qfrac_var[flav] = qfracs[flav]**2*(hists_rebinned[flav].variances()/hists_vals[flav]**2 + qfrac_var_all**2)
    # Efrac_splines = {key: 
    #                     np.array([get_spline(qfracs[key][:,eta_idx], pt_bins) 
    #                     for eta_idx in range(jeteta_bins.nbins)
    #                     ])
    #                 for key in qfracs.keys()
    #                 }
    # Efrac_2Dsplines = {key: RegularGridInterpolator((np.log10(pt_bins.centres), jeteta_bins.centres), qfracs[key], fill_value=None) 
                        # for key in qfracs.keys()}
    
#         yval = qfracs['g'][:,0]

#         Efrac_fit_dict[sample] = Efrac_fits
    # qfrac_dict[sample] = qfracs
    # qfrac_var_dict[sample] = qfrac_var
    # zero_spline = get_spline(np.array([1e-15]*pt_bins.nbins), pt_bins)
    # if add_composed_flav:
        # Efrac_splines['q'] = np.array([zero_spline]*jeteta_bins.nbins)
        # if not combine_antiflavour:
            # Efrac_splines['qbar'] = np.array([zero_spline]*jeteta_bins.nbins)
        # qfrac_spline_dict[sample] = FlavorFractions(Efrac_splines, eta_binning) #Efrac_splines #Efrac_2Dsplines
#         qfrac_spline_dict2[sample] = Efrac_splines
    
    # qfrac_spline_dict2D[sample] = Efrac_2Dsplines
    # sample_plot[sample] = [qfracs, qfrac_var, Efrac_splines, Efrac_2Dsplines]
    legend_labels.append(legend_label) if sample != 'QCD-Py_genwt_test' else legend_labels.append(legend_label+', weights')

### to make
hists_rebinned_dict['QCD-MG-Her'] = hist_mult(hists_rebinned_dict['QCD-MG-Her'], Neffs['QCD-MG-Py']/Neffs['QCD-MG-Her'])
plot_spectra(hists_rebinned_dict, legend_labels, 'all', eta_idx, jeteta_bins, pt_bins, saveplot=True, plotvspt=True)
