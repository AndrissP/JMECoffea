import matplotlib.pyplot as plt

import matplotlib as mpl
import os
import mplhep as hep
import numpy as np
from uncertainty_helpers import get_ratio, ptmin_global, ptmax_global
from coffea.lookup_tools import extractor
from scipy.optimize import curve_fit
from correction_fitter_helpers import inflate_smallest_std
from common_binning import JERC_Constants

hep_label = "Private work"

color_scheme = {key: cycler_vals
    for cycler_vals, key in zip(plt.rcParams['axes.prop_cycle'], ['g', 'ud', 'c', 'b', 'QCD', 'DY', 'TTBAR', 'DY200', 'unmatched', 's', 'q', 'cs', 'u', 'd'])}
color_scheme_antiflav = {key: cycler_vals
    for cycler_vals, key in zip(plt.rcParams['axes.prop_cycle'], ['g', 'udbar', 'cbar', 'bbar', 'QCD', 'DY', 'TTBAR', 'DY200', 'unmatched', 'sbar', 'qbar'])}
color_scheme = color_scheme | color_scheme_antiflav

color_scheme['b_prompt'] = color_scheme['b']
color_scheme['b_gluon_splitting'] = color_scheme['c']
color_scheme['c_prompt'] = color_scheme['b']
color_scheme['c_gluon_splitting'] = color_scheme['c']


legend_dict = {'g': 'Gluons', 'q': 'Quarks', 'ud':'UpDown', 'b': 'Bottom', 'c': 'Charm', 's': 'Strange', 'unmatched': 'Unmatched', 'u': 'Up', 'd': 'Down',}
from fileNames.available_datasets import legend_labels
legend_dict_short = {'g': 'g',
                     'ud': 'ud', 'q':'q', 'b': 'b', 'c': 'c', 's':'s', 'cs':'cs', 'd': 'd', 'u': 'u',
                     'unmatched': 'unmatched',
                     'udbar': '$\overline{ud}$', 'qbar':'$\overline{q}$', 'bbar': '$\overline{b}$', 'cbar': '$\overline{c}$', 'sbar':'$\overline{s}$',
                     'QCD': legend_labels["QCD"]["lab"], 'TTBAR': legend_labels["ttbar"]["lab"], 'DY': legend_labels["DY"]["lab"],
                      'b_gluon_splitting':'b gluon split', 'b_prompt':' b prompt',
                      'c_gluon_splitting':'c gluon split', 'c_prompt':' c prompt' }

legend_dict_q_vs_qbar = {'ud': 'ud/$(ud+\overline{ud}$)', 'q':'q/($q+\overline{q}$)', 'b': '$b/(b+\overline{b})$',
                         'c': '$c/(c+\overline{c})$', 's':'$s/(s+\overline{s})$', 'cs':'$cs/(cs+\overline{cs})$'}

def plot_Efractions(sampledict, etaidx, jeteta_bins, ptbins, legenddict=None, saveplot=False):
    samples = list(sampledict.keys())
    ptbins_c = ptbins.centres

#     ### Check that Herwig is the first sample and Pythia the second
#     if not ('Her' in samples[0] and 'Py' in samples[1]):
#         raise ValueError('key in the dictionary happened to get reversed')
    
    
    qfracs0, qfrac_var0, spline0, spline2D0 = sampledict[samples[0]]
    qfracs1, qfrac_var1, spline1, spline2D1 = sampledict[samples[1]]
    
    plot_range = range(0, np.searchsorted(ptbins_c,1250)) if 'DY' in "".join(samples) else range(0, np.searchsorted(ptbins_c,2750))
    ptbins_c_plot = ptbins_c[plot_range]
    
    fig, ax = plt.subplots()
    xplot = np.geomspace(ptbins_c_plot.min() - (1), ptbins_c_plot.max(),1000)
    xplot2 = np.geomspace(ptbins_c_plot.min(), ptbins_c_plot.max(),1000)
    points_ls = []
    for flav in qfracs0.keys():
        lab = legend_dict_short[flav]
#         mc = next(ax._get_lines.prop_cycler)

        points = ax.errorbar(ptbins_c_plot, qfracs0[flav][plot_range, etaidx],
                             yerr=np.sqrt(qfrac_var0[flav][plot_range, etaidx]),
                             linestyle='none', label=lab,  **color_scheme[flav], capsize=1.6, capthick=0.7, linewidth=1.0)
        points2 = ax.errorbar(ptbins_c_plot, qfracs1[flav][plot_range, etaidx],
                              yerr=np.sqrt(qfrac_var1[flav][plot_range, etaidx]),
                              linestyle='none', mfc='none', markeredgewidth=1.2, **color_scheme[flav], capsize=1.6, capthick=0.7, linewidth=1.0)

        valid_fit_val = ~(np.isnan(qfracs1[flav]) | np.isinf(qfracs1[flav]) | (qfracs1[flav]==0))
        
#         ax.plot(xplot, spline0[flav](np.log10(xplot)), '--', markersize=0, **mc, linewidth=1.0)
#         sp1 = ax.plot(xplot, spline1[flav](np.log10(xplot)), '--', markersize=0, **mc, linewidth=1.0)
        # ax.plot(xplot2, spline2D0[flav]((np.log10(xplot2), np.repeat([jeteta_bins.centres[etaidx]],len(xplot2)))),
        #               '-.', markersize=0, **color_scheme[flav], linewidth=1.0)
        # ax.plot(xplot2, spline2D1[flav]((np.log10(xplot2), np.repeat([jeteta_bins.centres[etaidx]],len(xplot2)))),
        #       '-.', markersize=0, **color_scheme[flav], linewidth=1.0)
        ax.plot(xplot2, spline0[flav][etaidx](np.log10(xplot2)),
                      '-.', markersize=0, **color_scheme[flav], linewidth=1.0)        
        # ax.plot(xplot2, spline2D0[flav]((np.log10(xplot2), np.repeat([jeteta_bins.centres[etaidx]],len(xplot2)))),
        #               '-.', markersize=0, **color_scheme[flav], linewidth=1.0)
        ax.plot(xplot2, spline1[flav][etaidx](np.log10(xplot2)),
              '-.', markersize=0, **color_scheme[flav], linewidth=1.0)
        # ax.plot(xplot2, spline2D1[flav]((np.log10(xplot2), np.repeat([jeteta_bins.centres[etaidx]],len(xplot2)))),
        #       '-.', markersize=0, **color_scheme[flav], linewidth=1.0)
# interp((np.log(np.arange(20,60,2)),[1]*20))
        if list(qfracs0.keys())[0] == flav:
            points_ls.append(points[0])
            points_ls.append(points2[0])
        
    
    ax.set_xscale('log')
    ax.set_xlabel('$p_{T,ptcl}$ (GeV)')
    ax.set_ylabel("Flavor fraction")

    xlims = ax.get_xlim()

    ax.set_xticks([])
    ax.hlines(0,8,10000, linestyles='--',color="black", linewidth=1,)
    ax.set_xticks([10, 20, 50, 100, 200, 500, 1000, 2000, 5000])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    # ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    legend_labs = [legenddict[samples[0]], legenddict[samples[1]] ] if legenddict is not None else [samples[0], samples[1]]
    legend1 = ax.legend(points_ls, legend_labs, loc="upper left", bbox_to_anchor=(0.56, 1))
    leg2 = ax.legend(ncol=1, loc='upper left', bbox_to_anchor=(0.42, 1))
    ax.add_artist(legend1)

    ylims = ax.get_ylim()
    ax.set_xlim(xlims)
    ax.set_ylim(ylims[0], ylims[1]*1.25)

    # ax.yaxis.get_ticklocs(minor=True)
    ax.minorticks_on()
    hep.cms.label(hep_label, loc=0, data=False, ax=ax, rlabel='')
    hep.label.exp_text(text=jeteta_bins.idx2plot_str(etaidx), loc=2)

    if saveplot:
        if not os.path.exists("fig/fractions"):
            os.mkdir("fig/fractions")

        etastr = jeteta_bins.idx2str(etaidx)
        fig_name = 'fig/fractions/fraction_'+etastr+"_"+"_".join(samples)
        print("Saving plot with the name = ", fig_name)
        plt.savefig(fig_name+'.pdf');
        plt.savefig(fig_name+'.png');
    # fig.close()

from helpers import hist_div, hist_add, hist_mult
from helpers import hist_div, hist_add, hist_mult
def plot_Efractions_ratio(sampledict, samples, etaidx, jeteta_bins, ptbins, legenddict=None, legenddict2 = None, saveplot=False,
                          legend1_loc=(0.42, 1), legend2_loc=(0.56, 1), ratio_title="Her7/Py8", ratio_lim=(0.5,1.5), fig_name=None, ylab_name=None, ylim=(-0.1, 1.25*1.3)):
    # samples = list(sampledict.keys())
    ptbins_c = ptbins.centres
    ptbins_e = ptbins.edges

    # ### Check that Herwig is the first sample and Pythia the second
    # if not ('Her' in samples[1] and 'Py' in samples[0]):
    #     raise ValueError('key in the dictionary happened to get reversed')
    
    
    qfracs0, qfrac_var0, spline0, spline2D0 = sampledict[samples[0]]
    qfracs1, qfrac_var1, spline1, spline2D1 = sampledict[samples[1]]
    
    plot_range = range(0, np.searchsorted(ptbins_c,1250)) if 'DY' in "".join(samples) else range(0, np.searchsorted(ptbins_c,2750))
    ptbins_c_plot = ptbins_c[plot_range]
    
    fig, (ax_main, ax_ratio) = plt.subplots(nrows=2, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0})
    xplot = np.geomspace(ptbins_c_plot.min() - (1), ptbins_c_plot.max(),1000)
    xplot2 = np.geomspace(ptbins_c_plot.min(), ptbins_c_plot.max(),1000)
    points_ls = []
    for flav in qfracs0.keys():
        if not legenddict2 is None:
            lab = legenddict2[flav]
        else:
            lab = legend_dict_short[flav]

        points = ax_main.errorbar(ptbins_c_plot, qfracs0[flav][plot_range, etaidx],
                             yerr=np.sqrt(qfrac_var0[flav][plot_range, etaidx]),
                             linestyle='none', label=lab,  **color_scheme[flav], capsize=1.6, capthick=0.7, linewidth=1.0)
        points2 = ax_main.errorbar(ptbins_c_plot, qfracs1[flav][plot_range, etaidx],
                              yerr=np.sqrt(qfrac_var1[flav][plot_range, etaidx]),
                              linestyle='none', mfc='none', markeredgewidth=1.2, **color_scheme[flav], capsize=1.6, capthick=0.7, linewidth=1.0)

        valid_fit_val = ~(np.isnan(qfracs1[flav]) | np.isinf(qfracs1[flav]) | (qfracs1[flav]==0))
        ax_main.plot(xplot2, spline0[flav][etaidx](np.log10(xplot2)),
                      '-.', markersize=0, **color_scheme[flav], linewidth=1.0)        
        # ax_main.plot(xplot2, spline2D0[flav]((np.log10(xplot2), np.repeat([jeteta_bins.centres[etaidx]],len(xplot2)))),
        #               '-.', markersize=0, **color_scheme[flav], linewidth=1.0)
        ax_main.plot(xplot2, spline1[flav][etaidx](np.log10(xplot2)),
              '-.', markersize=0, **color_scheme[flav], linewidth=1.0)
        # ax_main.plot(xplot2, spline2D1[flav]((np.log10(xplot2), np.repeat([jeteta_bins.centres[etaidx]],len(xplot2)))),
        #       '-.', markersize=0, **color_scheme[flav], linewidth=1.0)
# interp((np.log(np.arange(20,60,2)),[1]*20))
        if list(qfracs0.keys())[0] == flav:
            points_ls.append(points[0])
            points_ls.append(points2[0])
        
    
    ax_main.set_xscale('log')
    # ax_main.set_xlabel('$p_{T,ptcl}$ (GeV)')
    if ylab_name == None:
        ax_main.set_ylabel("Flavor fraction")
    else:
        ax_main.set_ylabel(ylab_name)

    xlims = ax_main.get_xlim()

    ax_main.hlines(0,8,10000, linestyles='--',color="black", linewidth=1,)
    ax_main.set_xticks([])
    ax_main.set_xticks([10, 20, 50, 100, 200, 500, 1000, 2000, 5000])
    ax_main.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    # ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    legend_labs = [legenddict[samples[0]], legenddict[samples[1]] ] if legenddict is not None else [samples[0], samples[1]]
    legend1 = ax_main.legend(points_ls, legend_labs, loc="upper left", bbox_to_anchor=legend2_loc)
    leg2 = ax_main.legend(ncol=1, loc='upper left', bbox_to_anchor=legend1_loc)
    ax_main.add_artist(legend1)
    # ax.add_artist(leg2)

    ylims = ax_main.get_ylim()
    ax_main.set_xlim(xlims)
    # ax_main.set_ylim(ylims[0], ylims[1]*1.3)
    ax_main.set_ylim(ylim)

    # ax.yaxis.get_ticklocs(minor=True)
    ax_main.minorticks_on()
    hep.cms.label(hep_label, loc=0, data=False, ax=ax_main, rlabel='')
    hep.label.exp_text(text=jeteta_bins.idx2plot_str(etaidx), loc=2, ax=ax_main)

    #### Ratio plot
    wd = np.diff(ptbins_e[range(plot_range[0], plot_range[-1]+2)])
    ax_ratio.hlines(1,-10, 10000, linestyles='--',color="black", 
               linewidth=1,)
    
    ratio = hist_div(qfracs1, qfracs0)
    ratio_unc_central = hist_div(qfrac_var0, hist_mult(qfracs0, qfracs0))
    ratio_unc_points = hist_div(qfrac_var1, hist_mult(qfracs0, qfracs0))
    ratio_unc = hist_add(ratio_unc_central, ratio_unc_points)
    # ratio_unc_points = hist_div(qfrac_var1, qfracs0)
    # data_model_ratio = yvals/yvals[0]
    # data_model_ratio_unc = stds / yvals[0]

    # for flav in qfracs0.keys():
    #     non_nan_ratio = ~np.isnan(ratio_unc_central[flav][plot_range, etaidx])
    #     ax_ratio.bar(
    #         ptbins_c_plot[non_nan_ratio],
    #         2 * np.sqrt(ratio_unc_central[flav][plot_range, etaidx][non_nan_ratio]),
    #         width=wd[non_nan_ratio],
    #         bottom=1.0 - np.sqrt(ratio_unc_central[flav][plot_range, etaidx][non_nan_ratio]),
    #         fill=False,
    #         linewidth=0,
    #         edgecolor=color_scheme[flav]['color'],
    #         hatch=10 * "/",
    #         # **color_scheme[flav],
    #     )

    for flav in qfracs0.keys():
        ax_ratio.errorbar(
            ptbins_c_plot,
            ratio[flav][plot_range, etaidx], #[nonzero_model_yield],
            yerr=np.sqrt(ratio_unc[flav][plot_range, etaidx]), #[nonzero_model_yield],
            linestyle="none",
            capsize=1.6, capthick=0.7, linewidth=1.0,
            mfc='none', markeredgewidth=1.2,
            **color_scheme[flav],
            #fmt=marker,
        )

    ax_ratio.set_ylim(ratio_lim)
    ### make the y-axis ticks in the ratio plot look nice: add a decent amount of major and minor ticks
    ax_ratio.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))
    ax_ratio.yaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins=25, steps=[1, 2, 5, 10])) #mpl.ticker.LinearLocator(numticks=25)
    ### remove the highest tick lavel from the ratio plot as it overlaps with the lowest label from the main plot 
    tick_labels = ax_ratio.get_yticks() 
    tick_labels = [f'{tick:.10g}' for tick in tick_labels]  ### remove floating point digits
    tick_labels = tick_labels[:-1]
    ax_ratio.set_yticks(ax_ratio.get_yticks()[:-1])
    ax_ratio.set_yticklabels(tick_labels)
    ax_ratio.set_xlabel('$p_{T,ptcl}$ (GeV)')
    ax_ratio.set_ylabel(ratio_title)
#     ax_ratio.set_ylabel()
    if saveplot:
        if not os.path.exists("fig/fractions"):
            os.mkdir("fig/fractions")

        etastr = jeteta_bins.idx2str(etaidx)
        if fig_name==None:
            fig_name = 'fig/fractions/fraction_'+etastr+"_"+"_".join(samples)
        if ('.pdf' in fig_name) or ('.png' in fig_name):
            fig_name = fig_name[:-4]
        print("Saving plot with the name = ", fig_name, ".pdf /.png")
        plt.savefig(fig_name+'.pdf')
        plt.savefig(fig_name+'.png')
    # fig.close()

# def plot_Efractions_ratio_antiflav(sampledict, samples, etaidx, jeteta_bins, ptbins, legenddict=None, saveplot=False,
#                           legend1_loc=(0.42, 1), legend2_loc=(0.56, 1), ratio_title="Her7/Py8", ratio_lim=(0.5,1.5)):
#     # samples = list(sampledict.keys())
#     ptbins_c = ptbins.centres
#     ptbins_e = ptbins.edges
    
#     qfracs0, qfrac_var0, spline0, spline2D0 = sampledict[samples[0]]
#     qfracs1, qfrac_var1, spline1, spline2D1 = sampledict[samples[1]]
    
#     plot_range = range(0, np.searchsorted(ptbins_c,1250)) if 'DY' in "".join(samples) else range(0, np.searchsorted(ptbins_c,2750))
#     ptbins_c_plot = ptbins_c[plot_range]
    
#     fig, (ax_main, ax_ratio) = plt.subplots(nrows=2, sharex=True, gridspec_kw={'height_ratios': [2, 2], 'hspace': 0})
#     xplot = np.geomspace(ptbins_c_plot.min() - (1), ptbins_c_plot.max(),1000)
#     xplot2 = np.geomspace(ptbins_c_plot.min(), ptbins_c_plot.max(),1000)
#     points_ls = []
#     for flav in qfracs0.keys():
#         lab = legend_dict_q_vs_qbar[flav]
#         # print('flav:', flav, ', lab:', lab)
#         mc1 = next(ax_main._get_lines.prop_cycler)

#         points = ax_main.errorbar(ptbins_c_plot, qfracs0[flav][plot_range, etaidx],
#                              yerr=np.sqrt(qfrac_var0[flav][plot_range, etaidx]),
#                              linestyle='none', label=lab,  **mc1, capsize=1.6, capthick=0.7, linewidth=1.0)
#         mc2 = next(ax_main._get_lines.prop_cycler)
#         points2 = ax_main.errorbar(ptbins_c_plot, qfracs1[flav][plot_range, etaidx],
#                               yerr=np.sqrt(qfrac_var1[flav][plot_range, etaidx]),
#                               linestyle='none', mfc='none', markeredgewidth=1.2  **mc2, capsize=1.6, capthick=0.7, linewidth=1.0)

#         valid_fit_val = ~(np.isnan(qfracs1[flav]) | np.isinf(qfracs1[flav]) | (qfracs1[flav]==0))
        
#         ax_main.plot(xplot2, spline2D0[flav]((np.log10(xplot2), np.repeat([jeteta_bins.centres[etaidx]],len(xplot2)))),
#                       '-.', markersize=0, **mc1, linewidth=1.0)
#         ax_main.plot(xplot2, spline2D1[flav]((np.log10(xplot2), np.repeat([jeteta_bins.centres[etaidx]],len(xplot2)))),
#               '-.', markersize=0, **mc2, linewidth=1.0)
# # interp((np.log(np.arange(20,60,2)),[1]*20))
#         if list(qfracs0.keys())[0] == flav:
#             points_ls.append(points[0])
#             points_ls.append(points2[0])
        
    
#     ax_main.set_xscale('log')
#     # ax_main.set_xlabel('$p_{T,ptcl}$ (GeV)')
#     ax_main.set_ylabel("Flavor fraction")

#     xlims = ax_main.get_xlim()

#     ax_main.hlines(0,8,10000, linestyles='--',color="black", linewidth=1,)
#     ax_main.set_xticks([])
#     ax_main.set_xticks([10, 20, 50, 100, 200, 500, 1000, 2000, 5000])
#     ax_main.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
#     # ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
#     legend_labs = [legenddict[samples[0]], legenddict[samples[1]] ] if legenddict is not None else [samples[0], samples[1]]
#     legend1 = ax_main.legend(points_ls, legend_labs, loc="upper left", bbox_to_anchor=legend2_loc)
#     leg2 = ax_main.legend(ncol=1, loc='upper left', bbox_to_anchor=legend1_loc)
#     ax_main.add_artist(legend1)
#     # ax.add_artist(leg2)

#     ylims = ax_main.get_ylim()
#     ax_main.set_xlim(xlims)
#     # ax_main.set_ylim(ylims[0], ylims[1]*1.3)
#     ax_main.set_ylim(-0.1, 1.25*1.3)

#     # ax.yaxis.get_ticklocs(minor=True)
#     ax_main.minorticks_on()
#     hep.cms.label(hep_label, loc=0, data=False, ax=ax_main, rlabel='')
#     hep.label.exp_text(text=jeteta_bins.idx2plot_str(etaidx), loc=2, ax=ax_main)

#     #### Ratio plot
#     wd = np.diff(ptbins_e[range(plot_range[0], plot_range[-1]+2)])
#     ax_ratio.hlines(1,-10, 10000, linestyles='--',color="black", 
#                linewidth=1,)
    
#     ratio = hist_div(qfracs1, qfracs0)
#     ratio_unc_central = hist_div(qfrac_var0, hist_mult(qfracs0, qfracs0))
#     ratio_unc_points = hist_div(qfrac_var1, hist_mult(qfracs0, qfracs0))
#     ratio_unc = hist_add(ratio_unc_central, ratio_unc_points)

#     for flav in qfracs0.keys():
#         ax_ratio.errorbar(
#             ptbins_c_plot,
#             ratio[flav][plot_range, etaidx], #[nonzero_model_yield],
#             yerr=np.sqrt(ratio_unc[flav][plot_range, etaidx]), #[nonzero_model_yield],
#             linestyle="none",
#             capsize=1.6, capthick=0.7, linewidth=1.0,
#             mfc='none', markeredgewidth=1.2,
#             **color_scheme[flav],
#             #fmt=marker,
#         )

#     ax_ratio.set_ylim((0.7,1.3))
#     ### make the y-axis ticks in the ratio plot look nice: add a decent amount of major and minor ticks
#     ax_ratio.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))
#     ax_ratio.yaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins=25, steps=[1, 2, 5, 10])) #mpl.ticker.LinearLocator(numticks=25)
#     ### remove the highest tick lavel from the ratio plot as it overlaps with the lowest label from the main plot 
#     tick_labels = ax_ratio.get_yticks() 
#     tick_labels = [f'{tick:.10g}' for tick in tick_labels]  ### remove floating point digits
#     tick_labels = tick_labels[:-1]
#     ax_ratio.set_yticks(ax_ratio.get_yticks()[:-1])
#     ax_ratio.set_yticklabels(tick_labels)
#     ax_ratio.set_xlabel('$p_{T,ptcl}$ (GeV)')
#     ax_ratio.set_ylabel(ratio_title)
# #     ax_ratio.set_ylabel()
#     if saveplot:
#         if not os.path.exists("fig/fractions"):
#             os.mkdir("fig/fractions")

#         etastr = jeteta_bins.idx2str(etaidx)
#         fig_name = 'fig/fractions/fraction_'+etastr+"_"+"_".join(samples)
#         print("Saving plot with the name = ", fig_name)
#         plt.savefig(fig_name+'.pdf');
#         plt.savefig(fig_name+'.png');
#     # fig.close()

from uncertainty_helpers import get_ratio, read_data2

def poly4(x, *p):
    c0, c1, c2, c3, c4 = p
    xs = np.log10(x)
    res = c0+c1*xs+c2*xs**2+c3*xs**3+c4*xs**4
    return res

def poly4lims(x, xmin, xmax, *p):
    xcp = x.copy()
    lo_pos = xcp<xmin
    hi_pos = xcp>xmax
    xcp[lo_pos] = xmin
    xcp[hi_pos] = xmax
    return poly4(xcp, *p)

color_scheme2 = color_scheme.copy()
color_scheme2['QCD, 3 jets'] = {'color': 'brown', 'marker': 'o'}
color_scheme2['DY, 2 jets'] = {'color': 'cyan', 'marker': 'o'}

from correction_fitter_helpers import find_stationary_pnt_poly

def plot_ratio_comparisons_samples(flav, etaidx, jeteta_bins, ptbins_c,
                                   eta_binning_str, 
                                   evaluator,
                                   evaluator_names:dict,
                                   divide:bool=False,
                                   inverse:bool=False,
                                   use_recopt:bool=False,
                                   maxlimit_static_pnt:bool=True,
                                   max_point_fit_idx=3,
                                   plotnewfit:bool=True,
                                   plotcorrectionratios:bool=False,
                                   inflate_smallest_std_bool:bool=True,
                                   show_original_uncertainties:bool=True,
                                   denom_samples:list=['_QCD-MG-Py', '_Pythia-TTBAR', '_DY-MG-Py'],
                                    samples:list=['_QCD-MG-Her', '_Herwig-TTBAR', '_DY-MG-Her'],
                                    sample_lab:list=['QCD', 'TTBAR', 'DY'],
                                    plot_corr_fits:dict={'J':'QCD', 'T': 'TTBAR', 'S': 'Sim.'},
                                    draw_DY200line:bool = False,
                                    fit_simultaneously:bool = True,
                                   ):
    ''' Put ratio plots of many all flavors at the same place. Reproduce Fig. 31 in arXiv:1607.03663
    Output, polynomial coeficients of the data ratio fit
    divide: True if divide Herwig by Pythia, False if subtract Pythia from Herwig
    inverse: True if plot corrections, False if plot responses
    use_recopt:  True if use reco pt, False if use gen pt
    plotcorrectionratios: True if plot the curves obtained from fitting the individual corrections and then taking ratios of them
    '''        
    mean_name = "Median"
    mean_name_std = mean_name+'Std'

    ### Set plotting range (can be different from fitting range)
    start = np.searchsorted(ptbins_c, 16, side='left')
    end = 27

    #### Read median response/correction data
    yvals = np.array([read_data2(mean_name, samp, flav, eta_binning_str)[start:end,etaidx] for samp in samples])
    stds  = np.array([read_data2(mean_name_std, samp, flav, eta_binning_str)[start:end,etaidx] for samp in samples])
    xvals = np.array([read_data2("MeanRecoPt", samp, flav, eta_binning_str)[start:end,etaidx] for samp in samples])
    
    yvals_d = np.array([read_data2(mean_name, samp, flav, eta_binning_str)[start:end,etaidx] for samp in denom_samples])
    stds_d  = np.array([read_data2(mean_name_std, samp, flav, eta_binning_str)[start:end,etaidx] for samp in denom_samples])


        
    #### Clean and set up the data for plotting
    yvals[(yvals==0) | (np.abs(yvals)==np.inf)] = np.nan
    yvals_d[(yvals_d==0) | (np.abs(yvals_d)==np.inf)] = np.nan
    
    ratios = get_ratio(yvals, yvals_d, divide)
    if divide==True:
        ratio_unc = ((stds / yvals_d)**2 + (yvals/yvals_d**2 * stds_d)**2)**(1/2)
    else:
        ratio_unc = (stds**2+stds_d**2)**(1/2)

    ratio_unc_plot = ratio_unc.copy()
    if inflate_smallest_std_bool:
        ratio_unc = inflate_smallest_std(ratio_unc)
    if not show_original_uncertainties:
        ratio_unc_plot = ratio_unc.copy()

    if not use_recopt:
        xvals = ptbins_c[start:end]    
        

    fig, ax = plt.subplots()
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_minor_locator(mpl.ticker.AutoMinorLocator())
        
    #### Plot the points
    for yval, std, samp in zip(ratios, ratio_unc_plot, sample_lab):
        ax.errorbar(xvals, yval, yerr=std,
                    linestyle="none", label=legend_dict_short[samp], **color_scheme2[samp],
                    capsize=1.6, capthick=0.7, linewidth=1.0)
       
    #### Plot pre-fitted curves
    xvals_cont = np.geomspace(np.min(xvals), np.max(xvals), 200)
    etaval = jeteta_bins.centres[etaidx]
    legend_dict_short['Sim.'] = 'Sim.'
    color_scheme['Sim.'] = {'color': 'purple', 'marker': 'o'}

    flav_unc = flav if not flav == 'all' else 'a'
    if plotcorrectionratios:
        for fit_samp in plot_corr_fits: # , lab in zip(tags, labs):
            lab = plot_corr_fits[fit_samp]
            eva = evaluator[f'{evaluator_names["Sum20Her"]}_{flav_unc}{fit_samp}']
            eva_d = evaluator[f'{evaluator_names["Sum20Py"]}_{flav_unc}{fit_samp}']
            corr_etabins = eva._bins['JetEta'] 
            corr_bin_idx = np.searchsorted(corr_etabins, etaval, side='right')-1
            ptmax = list(eva._eval_clamp_maxs.values())[0][corr_bin_idx]
            ptmax_d = list(eva_d._eval_clamp_maxs.values())[0][corr_bin_idx]
            ptmax = min([ptmax, ptmax_d])

            yvals_cont = eva(np.array([etaval]),xvals_cont)
            yvals_cont_d = eva_d(np.array([etaval]),xvals_cont)

            if inverse==True:
                yvals = 1/yvals
                yvals_d = 1/yvals_d
                ### Error propagation
                stds = yvals**2*stds
                stds_d = yvals_d**2*stds_d

            if inverse==False:
                yvals_cont = 1/yvals_cont
                yvals_cont_d = 1/yvals_cont_d
            
            ratios_cont = get_ratio(yvals_cont, yvals_cont_d, divide)
            # xspacing = xvals_cont[1]-xvals_cont[1]
            # breakpoint()
            ratios_cont[xvals_cont>(ptmax)] = ratios_cont[np.searchsorted(xvals_cont, (ptmax))-1]
            ax.plot(xvals_cont, ratios_cont, markersize=0, **color_scheme[lab], label=legend_dict_short[lab]+' fit')

    ax.set_xscale('log')
    xlims = ax.get_xlim()
    
    ax.hlines(1,1, 10000, linestyles='--',color="black", linewidth=1,)
    
    ####################### Fit ####################
    # if fit_simultaneously:
    #     fit_minx = np.searchsorted(ptbins_c, ptmin_global, side='left') - 1
    #     fit_maxx = np.searchsorted(ptbins_c, ptmax_global, side='left')
        
    #     xval4fit = np.tile(xvals[fit_minx:fit_maxx], len(sample_lab))
    #     yval4fit = np.concatenate(ratios[:,fit_minx:fit_maxx])
    #     ratio_unc4fit = np.concatenate(ratio_unc[:,fit_minx:fit_maxx])
    #     validpt_mask = ~(np.isnan(yval4fit) | np.isinf(yval4fit) | (yval4fit==0))
    #     xval4fit = np.array([xval4fit[validpt_mask]])
    #     yval4fit = np.array([yval4fit[validpt_mask]])
    #     ratio_unc4fit = np.array([ratio_unc4fit[validpt_mask]])
    # else:
    #     # xval4fit = xvals
    #     # yval4fit = ratios
    #     # ratio_unc4fit = ratio_unc
    #     validpt_mask = ~(np.isnan(ratios) | np.isinf(ratios) | (ratios==0))
    #     xval4fit = np.array([ xvals[validpt_mask[ii]] for ii in range(len(validpt_mask))])
    #     yval4fit = np.array([ ratios[ii][validpt_mask[ii]] for ii in range(len(validpt_mask))]) #yval4fit[validpt_mask]
    #     ratio_unc4fit = np.array([ ratio_unc[ii][validpt_mask[ii]] for ii in range(len(validpt_mask))]) # ratio_unc4fit[validpt_mask]
    # ### Put the minimum limit on the relative uncertainty to min_rel_uncert
    min_rel_uncert = 0.001
    
    polys = []
    xfitmins = []
    xfitmaxs = []
    if fit_simultaneously:
        len_fits = 1
    else:
        len_fits = len(ratios)
    for ii in range(len_fits):
        fit_label = r'Poly, n=4'
        if fit_simultaneously:
            fit_minx = np.searchsorted(ptbins_c, ptmin_global, side='left') - 1
            fit_maxx = np.searchsorted(ptbins_c, ptmax_global, side='left')
            
            xval4fit = np.tile(xvals[fit_minx:fit_maxx], len(sample_lab))
            yval4fit = np.concatenate(ratios[:,fit_minx:fit_maxx])
            ratio_unc4fit = np.concatenate(ratio_unc[:,fit_minx:fit_maxx])
        else:
            xval4fit = xvals
            yval4fit = ratios[ii]
            ratio_unc4fit = ratio_unc[ii]
            fit_label+=f', {sample_lab[ii]}'

        validpt_mask = ~(np.isnan(yval4fit) | np.isinf(yval4fit) | (yval4fit==0))
        xval4fit = xval4fit[validpt_mask]
        yval4fit = yval4fit[validpt_mask]
        ratio_unc4fit = ratio_unc4fit[validpt_mask]

        if divide == True:
            where_limit_std = (ratio_unc4fit/yval4fit)<min_rel_uncert
            ratio_unc4fit[where_limit_std] = min_rel_uncert*yval4fit[where_limit_std]
        else:
            where_limit_std = ratio_unc4fit<min_rel_uncert
            ratio_unc4fit[where_limit_std] = min_rel_uncert
        p_poly4_1, arr = curve_fit(poly4, xval4fit, yval4fit, p0=[ 1, 1, 1, 1, 1])
        p_poly4, arr = curve_fit(poly4, xval4fit, yval4fit, p0=p_poly4_1, sigma=ratio_unc4fit)
    #     p_poly4_1, arr = curve_fit(np.tile(xvals,len(sample_lab)), np.concatenate(ratios), means2fit, p0=[ 1, 1, 1, 1, 1])
        xfitmin = xval4fit.min()
        xfitmax = xval4fit.max()
        polys.append(p_poly4)
        xfitmins.append(xfitmin)
        xfitmaxs.append(xfitmax)
        poly4fun = lambda x, p: poly4lims(x, xfitmin, xfitmax, *p)
        y_poly4 = poly4fun(xvals_cont, p_poly4)

        if maxlimit_static_pnt:
            fit_max_lim_new = find_stationary_pnt_poly(xfitmin, xfitmax, *p_poly4, degree=4)
        else:
            fit_max_lim_new = xfitmax

        fit_max_lim_idx = np.searchsorted(np.sort(xval4fit), fit_max_lim_new, side="right")
        if maxlimit_static_pnt & (fit_max_lim_idx==len(xval4fit)) | (fit_max_lim_idx<=len(xval4fit)-max_point_fit_idx):
            # static point is too low or the last point that usually fluctuates out
            fit_max_lim_idx = len(xval4fit)-max_point_fit_idx
            fit_max_lim_new = np.sort(xval4fit)[fit_max_lim_idx]
        xplot_max_new = np.searchsorted(xvals_cont, fit_max_lim_new)
        y_poly4[xplot_max_new:] = y_poly4[xplot_max_new]
        # y_poly4_now = poly4fun(xvals_cont, p_poly4_1)
        if plotnewfit:
            mc_new_fit= next(ax._get_lines.prop_cycler)
            ax.plot(xvals_cont, y_poly4, label=fit_label ,linewidth=1.5, markersize=0, **mc_new_fit)

        if draw_DY200line:
            if not 'DY' in sample_lab:
                raise Exception("DY not in the sample_lab")
            else:
                if not 'DY' in sample_lab[ii]:
                    pt200idx = np.searchsorted(xvals_cont, 200, side='left')
                    # y_poly4[pt200idx]
                    ax.hlines(y_poly4[pt200idx], 10, 10000, linestyles='--', linewidth=1, color=mc_new_fit['color'], label='DY, 200 GeV mix')
    ####################### End fit ####################

    ####################### Calculate resonable limits excluding the few points with insane errors
    recalculate_limits=True
    if recalculate_limits:
        yerr_norm = np.concatenate(ratio_unc)
        y_norm = np.concatenate(ratios)
        norm_pos = (yerr_norm<0.01) &  (yerr_norm != np.inf) & (y_norm>-0.1)  
        if ~np.any(norm_pos):
            print("Cannot determine ylimits")
            norm_pos = np.ones(len(yerr_norm), dtype=int)
            raise Exception("Cannot determine ylimits")
        left_lim = np.min((y_norm-yerr_norm)[norm_pos])
        right_lim = np.max((yerr_norm+y_norm)[norm_pos])
        lim_pad = (right_lim - left_lim)/20
        ax.set_ylim(left_lim-lim_pad, right_lim+lim_pad*10)
    

    ####################### Formalities and save plot ####################`
    xlabel = r'$p_{T,reco}$ (GeV)' if use_recopt else r'$p_{T,ptcl}$ (GeV)'
    ax.set_xlabel(xlabel);
    ylab_pre = 'R(Her7)/R(Py8)' if divide else 'R(Her7)-R(Py8)'
    ylabel = r' (correction)' if inverse else r' (median response)'
    ax.set_ylabel(ylab_pre+ylabel);
    
    ax.set_xticks([10, 20, 50, 100, 500, 1000, 5000])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    hep.cms.label(hep_label, loc=0, data=False, ax=ax, rlabel='')
    hep.label.exp_text(text=jeteta_bins.idx2plot_str(etaidx)+f'\n{flav} jets', loc=2)
        

    leg1 = ax.legend(ncol=1, loc='upper right', bbox_to_anchor=(0.96, 1), handlelength=1.1)
    ax.set_xlim(xlims)
        
    figdir = "fig/uncertainty"
    if not os.path.exists(figdir):
        os.mkdir(figdir)
    etastr = jeteta_bins.idx2str(etaidx)
    fig_name = f'fig/uncertainty/Pythia_Herwig_all_samples_{flav}_jets_{etastr}'
    print("Saving plot with the name = ", fig_name)
    plt.savefig(fig_name+'.pdf');
    plt.savefig(fig_name+'.png');
    plt.show()
    return [[poly, xfitmin, xfitmax] for poly, xfitmin, xfitmax in zip(polys, xfitmins, xfitmaxs)]
    # return [p_poly4, xfitmin, xfitmax]

# def plot_ratio_all(flav, etaidx, jeteta_bins, ptbins_c,
#                                    eta_binning_str, 
#                                    evaluator,
#                                    evaluator_names:dict,
#                                    divide:bool=False,
#                                    inverse:bool=False,
#                                    use_recopt:bool=False,
#                                    maxlimit_static_pnt:bool=True,
#                                    max_point_fit_idx=3,
#                                    plotsimfit:bool=False,
#                                    plotnewfit:bool=True,
#                                    plotcorrectionratios:bool=False,
#                                    inflate_smallest_std_bool:bool=True,
#                                    show_original_uncertainties:bool=True,
#                                    denom_samples:list=['_QCD-MG-Py', '_Pythia-TTBAR', '_DY-MG-Py'],
#                                     samples:list=['_QCD-MG-Her', '_Herwig-TTBAR', '_DY-MG-Her'],
#                                     sample_lab:list=['QCD', 'TTBAR', 'DY'],
#                                     plot_corr_fits:dict={'J':'QCD'},
#                                    ):
#     ''' Fit ratio of responses for the inclusive flavor to obtain the normalization factors
#     Output, polynomial coeficients of the data ratio fit
#     '''        
#     mean_name = "Median"
#     mean_name_std = mean_name+'Std'

#     ### Set plotting range (can be different from fitting range)
#     start = np.searchsorted(ptbins_c, 16, side='left')
#     end = 27
    
#     yvals = np.array([read_data2(mean_name, samp, flav, eta_binning_str)[start:end,etaidx] for samp in samples])
#     stds  = np.array([read_data2(mean_name_std, samp, flav, eta_binning_str)[start:end,etaidx] for samp in samples])
#     xvals = np.array([read_data2("MeanRecoPt", samp, flav, eta_binning_str)[start:end,etaidx] for samp in samples])
    
#     yvals_d = np.array([read_data2(mean_name, samp, flav, eta_binning_str)[start:end,etaidx] for samp in denom_samples])
#     stds_d  = np.array([read_data2(mean_name_std, samp, flav, eta_binning_str)[start:end,etaidx] for samp in denom_samples])
        
#     #### Clean and set up the data for plotting
#     yvals[(yvals==0) | (np.abs(yvals)==np.inf)] = np.nan
#     yvals_d[(yvals_d==0) | (np.abs(yvals_d)==np.inf)] = np.nan
    
#     ratios = get_ratio(yvals, yvals_d, divide)
#     if divide==True:
#         ratio_unc = ((stds / yvals_d)**2 + (yvals/yvals_d**2 * stds_d)**2)**(1/2)
#     else:
#         ratio_unc = (stds**2+stds_d**2)**(1/2)

#     ratio_unc_plot = ratio_unc.copy()
#     if inflate_smallest_std_bool:
#         ratio_unc = inflate_smallest_std(ratio_unc)
#     if not show_original_uncertainties:
#         ratio_unc_plot = ratio_unc.copy()

#     if not use_recopt:
#         xvals = ptbins_c[start:end]    
        

#     fig, ax = plt.subplots()
#     for axis in [ax.xaxis, ax.yaxis]:
#         axis.set_minor_locator(mpl.ticker.AutoMinorLocator())
        
#     #### Plot the points
#     for yval, std, samp in zip(ratios, ratio_unc_plot, sample_lab):
#         ax.errorbar(xvals, yval, yerr=std,
#                     linestyle="none", label=legend_dict_short[samp], **color_scheme2[samp],
#                     capsize=1.6, capthick=0.7, linewidth=1.0)
       
#     #### Plot pre-fitted curves
#     xvals_cont = np.geomspace(np.min(xvals), np.max(xvals), 200)
#     etaval = jeteta_bins.centres[etaidx]
#     if plotsimfit:
#         tags, labs = ['J', 'T', 'S'], ['QCD', 'TTBAR', 'two fit difference']
#         legend_dict_short['two fit difference'] = 'Sim.'
#         color_scheme['two fit difference'] = {'color': 'purple', 'marker': 'o'}
#     else:
#         tags, labs = ['J', 'T'], ['QCD', 'TTBAR']

#     flav_unc = flav if not flav == 'all' else 'a'
#     if plotcorrectionratios:
#         for fit_samp, lab in zip(tags, labs):
#             eva = evaluator[f'{evaluator_names["Sum20Her"]}_{flav_unc}{fit_samp}']
#             eva_d = evaluator[f'{evaluator_names["Sum20Py"]}_{flav_unc}{fit_samp}']
#             corr_etabins = eva._bins['JetEta'] 
#             corr_bin_idx = np.searchsorted(corr_etabins, etaval, side='right')-1
#             ptmax = list(eva._eval_clamp_maxs.values())[0][corr_bin_idx]
#             ptmax_d = list(eva_d._eval_clamp_maxs.values())[0][corr_bin_idx]
#             ptmax = min([ptmax, ptmax_d])

#             yvals_cont = eva(np.array([etaval]),xvals_cont)
#             yvals_cont_d = eva_d(np.array([etaval]),xvals_cont)

#             if inverse==True:
#                 yvals = 1/yvals
#                 yvals_d = 1/yvals_d
#                 ### Error propagation
#                 stds = yvals**2*stds
#                 stds_d = yvals_d**2*stds_d

#             if inverse==False:
#                 yvals_cont = 1/yvals_cont
#                 yvals_cont_d = 1/yvals_cont_d
            
#             ratios_cont = get_ratio(yvals_cont, yvals_cont_d, divide)
#             # xspacing = xvals_cont[1]-xvals_cont[1]
#             # breakpoint()
#             ratios_cont[xvals_cont>(ptmax)] = ratios_cont[np.searchsorted(xvals_cont, (ptmax))-1]
#             ax.plot(xvals_cont, ratios_cont, markersize=0, **color_scheme[lab], label=legend_dict_short[lab]+' fit')

#     ax.set_xscale('log')
#     xlims = ax.get_xlim()
    
#     ax.hlines(1,1, 10000, linestyles='--',color="black", linewidth=1,)
    
#     ####################### Fit ####################
#     fit_minx = np.searchsorted(ptbins_c, ptmin_global, side='left') - 1
#     fit_maxx = np.searchsorted(ptbins_c, ptmax_global, side='left')
    
#     xval4fit = np.tile(xvals[fit_minx:fit_maxx], len(sample_lab))
#     yval4fit = np.concatenate(ratios[:,fit_minx:fit_maxx])
#     ratio_unc4fit = np.concatenate(ratio_unc[:,fit_minx:fit_maxx])
#     validpt_mask = ~(np.isnan(yval4fit) | np.isinf(yval4fit) | (yval4fit==0))
#     xval4fit = xval4fit[validpt_mask]
#     yval4fit = yval4fit[validpt_mask]
#     ratio_unc4fit = ratio_unc4fit[validpt_mask]
#     ### Put the minimum limit on the relative uncertainty to min_rel_uncert
#     min_rel_uncert = 0.001
#     if divide == True:
#         where_limit_std = (ratio_unc4fit/yval4fit)<min_rel_uncert
#         ratio_unc4fit[where_limit_std] = min_rel_uncert*yval4fit[where_limit_std]
#     else:
#         where_limit_std = ratio_unc4fit<min_rel_uncert
#         ratio_unc4fit[where_limit_std] = min_rel_uncert
    
#     p_poly4_1, arr = curve_fit(poly4, xval4fit, yval4fit, p0=[ 1, 1, 1, 1, 1])
#     p_poly4, arr = curve_fit(poly4, xval4fit, yval4fit, p0=p_poly4_1, sigma=ratio_unc4fit)
# #     p_poly4_1, arr = curve_fit(np.tile(xvals,len(sample_lab)), np.concatenate(ratios), means2fit, p0=[ 1, 1, 1, 1, 1])
#     xfitmin = xval4fit.min()
#     xfitmax = xval4fit.max()
#     poly4fun = lambda x, p: poly4lims(x, xfitmin, xfitmax, *p)
#     y_poly4 = poly4fun(xvals_cont, p_poly4)

#     if maxlimit_static_pnt:
#         fit_max_lim_new = find_stationary_pnt_poly(xfitmin, xfitmax, *p_poly4, degree=4)
#     else:
#         fit_max_lim_new = xfitmax

#     fit_max_lim_idx = np.searchsorted(np.sort(xval4fit), fit_max_lim_new, side="right")
#     if maxlimit_static_pnt & (fit_max_lim_idx==len(xval4fit)) | (fit_max_lim_idx<=len(xval4fit)-max_point_fit_idx):
#         # static point is too low or the last point that usually fluctuates out
#         fit_max_lim_idx = len(xval4fit)-max_point_fit_idx
#         fit_max_lim_new = np.sort(xval4fit)[fit_max_lim_idx]
#     xplot_max_new = np.searchsorted(xvals_cont, fit_max_lim_new)
#     y_poly4[xplot_max_new:] = y_poly4[xplot_max_new]
#     # y_poly4_now = poly4fun(xvals_cont, p_poly4_1)
#     if plotnewfit:
#         ax.plot(xvals_cont, y_poly4, label=r'Poly, n=4' ,linewidth=1.5, markersize=0);
#     ####################### End fit ####################

#     ####################### Calculate resonable limits excluding the few points with insane errors
#     recalculate_limits=True
#     if recalculate_limits:
#         yerr_norm = np.concatenate(ratio_unc)
#         y_norm = np.concatenate(ratios)
#         norm_pos = (yerr_norm<0.01) &  (yerr_norm != np.inf) & (y_norm>-0.1)  
#         if ~np.any(norm_pos):
#             print("Cannot determine ylimits")
#             norm_pos = np.ones(len(yerr_norm), dtype=int)
#             raise Exception("Cannot determine ylimits")
#         left_lim = np.min((y_norm-yerr_norm)[norm_pos])
#         right_lim = np.max((yerr_norm+y_norm)[norm_pos])
#         lim_pad = (right_lim - left_lim)/20
#         ax.set_ylim(left_lim-lim_pad, right_lim+lim_pad*10)
    

#     ####################### Formalities and save plot ####################`
#     xlabel = r'$p_{T,reco}$ (GeV)' if use_recopt else r'$p_{T,ptcl}$ (GeV)'
#     ax.set_xlabel(xlabel);
#     ylab_pre = 'R(Her7)/R(Py8)' if divide else 'R(Her7)-R(Py8)'
#     ylabel = r' (correction)' if inverse else r' (median response)'
#     ax.set_ylabel(ylab_pre+ylabel);
    
#     ax.set_xticks([10, 20, 50, 100, 500, 1000, 5000])
#     ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
#     hep.cms.label(hep_label, loc=0, data=False, ax=ax, rlabel='')
#     hep.label.exp_text(text=jeteta_bins.idx2plot_str(etaidx)+f'\n{flav} jets', loc=2)

#     leg1 = ax.legend(ncol=1, loc='upper right', bbox_to_anchor=(0.92, 1))
#     ax.set_xlim(xlims)
        
#     figdir = "fig/uncertainty"
#     if not os.path.exists(figdir):
#         os.mkdir(figdir)
#     etastr = jeteta_bins.idx2str(etaidx)
#     fig_name = f'fig/uncertainty/Pythia_Herwig_all_samples_{flav}_jets_{etastr}'
#     print("Saving plot with the name = ", fig_name)
#     plt.savefig(fig_name+'.pdf');
#     plt.savefig(fig_name+'.png');
#     plt.show()
#     return [p_poly4, xfitmin, xfitmax]

from matplotlib.legend_handler import HandlerLine2D

class HandlerLine2D_numpoints(HandlerLine2D):
    def __init__(self, marker_pad=0.3, numpoints=1, *args, **kwargs):
        self._marker_pad = marker_pad
        self._numpoints = numpoints
        HandlerLine2D.__init__(self, *args, **kwargs)

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        xdata, xdata_marker = self.get_xdata(legend, xdescent, ydescent,
                                             width, height, fontsize)
        ydata = ((height - ydescent) / 2.) * np.ones(xdata.shape, float)
        legline, = plt.plot(xdata, ydata, color=orig_handle.get_color(),
                            linestyle=orig_handle.get_linestyle(),
                            linewidth=orig_handle.get_linewidth())

        legline_marker, = plt.plot(xdata_marker, ydata[:len(xdata_marker)],
                                   marker=orig_handle.get_marker(),
                                   color=orig_handle.get_color(),
                                   linestyle=orig_handle.get_linestyle(),
                                   linewidth=orig_handle.get_linewidth(),
                                   markersize=orig_handle.get_markersize())

        return legline, legline_marker

def plot_uncertainty_antiflav(ptvals, etavals, additional_uncertainty_curves, uncertainties, ptoretastr, flavors, plotvspt=False):
    addc = additional_uncertainty_curves
    fig, ax = plt.subplots()

    xvals = ptvals if plotvspt else etavals
    flav_labs = []
    antiflav_labs = []

    for flav in flavors:
        lab = legend_dict_short[flav]
        norm_factor = 0               # no normalization for flal/antiflav uncertainty

        linestyle = '-.' if 'bar' in flav else '-'
        line = ax.plot(xvals, (addc[f'{flav}100']-norm_factor)*100, label=lab, markersize=0, linewidth=1.2, linestyle=linestyle,
                **color_scheme[flav])
        if 'bar' in flav:
            antiflav_labs.append(line[0])
        else:
            flav_labs.append(line[0])
    
    ax.hlines(0, ax.get_xlim()[0], ax.get_xlim()[1],color="gray",
        linewidth=1, alpha=0.4)

    smaller_spacing = plt.rcParams['legend.labelspacing']*0.45
    larger_spacing = plt.rcParams['legend.labelspacing']*2.0
    legend1 = ax.legend(handles=antiflav_labs, loc='upper right', bbox_to_anchor=(0.73, 0.978), handlelength=1.5, handleheight=1.55, labelspacing = smaller_spacing)
    leg2 = ax.legend(handles=flav_labs, ncol=1, loc='upper left', bbox_to_anchor=(0.69, 0.974), handlelength=0.9, handleheight=0.3, labelspacing = larger_spacing)#, title='assembled\nfrom QCD', title_fontsize=10)
    ax.add_artist(legend1)
#     ax.annotate("flavor", xy=(0.72,0.99), textcoords='axes fraction')
#     ax.annotate("antiflavor", xy=(0.6,0.99), textcoords='axes fraction')
    xlabel = r'$p_{T}$ (GeV)' if plotvspt else r'$\eta$'
    ax.set_xlabel(xlabel);
    ylabel = 'JEC uncertainty (%)'
    ax.set_ylabel(ylabel);
    ylim_old = ax.get_ylim()

    if plotvspt:
        ax.set_xscale('log')
        ax.set_xticks([10, 20, 50, 100, 500, 1000, 5000])
        ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
        ax.set_xlim(15,1000)
    else:        
        for HCal_border in JERC_Constants.etaBinsEdges_Win14():
            ax.vlines(HCal_border, ylim_old[0]*3, ylim_old[1]*3, linestyles='--',color="gray",
                linewidth=1,)
        ax.set_xlim(0,5)

#     ax.set_ylim(0.9885,1.0205)
    ylim_pad = (ylim_old[1]-ylim_old[0])*0.42 if plotvspt else (ylim_old[1]-ylim_old[0])*0.62
    ax.set_ylim(ylim_old[0],ylim_old[1]+ylim_pad)
    labtxt = f'{ptoretastr}' #if plotvspt else f'{ptoretastr}'
    labtxt = labtxt.replace("TTBAR", legend_labels["ttbar"]["lab"])
#     labtxt = f'$\eta$ = {etabins_abs[ptoretaidx]}' if plotvspt else f'$p_T$ = {ptbins_c[ptoretaidx]} GeV'
    hep.cms.label(hep_label, loc=0, data=False, ax=ax, rlabel='')
    hep.label.exp_text(text=labtxt, loc=2)

    ax.text(0.72, 0.95,'flavor',transform=ax.transAxes, fontsize=9.5)
    ax.text(0.48, 0.95,'antiflavor',transform=ax.transAxes, fontsize=9.5)
#     ax.annotate("flavor", xy=(0.05,0.05), textcoords='axes fraction')
#     ax.annotate("antiflavor", xy=(0,0), textcoords='axes fraction')
    figdir = "fig/uncertainty"
    if not os.path.exists(figdir):
        os.mkdir(figdir)

    if plotvspt:
        fig_name = figdir+f"/JECuncertainty_vs_pt_eta_{ptoretastr}".replace('.','')
    else:
        fig_name = figdir+f"/JECuncertainty_vs_pt_pt_{ptoretastr}".replace('.','_')
    fig_name = fig_name.replace(', ', '_').replace(' ', '_').replace('$', '').replace('=', '_').replace('\eta', 'eta').replace('|', '').replace('<', '').replace('\n', '_')
    print("Saving plot with the name = ", fig_name+".pdf / .png")
    plt.savefig(fig_name+'.pdf');
    plt.savefig(fig_name+'.png');
    plt.show()

def plot_uncertainty(ptvals, etavals, HerPy_differences, additional_uncertainty_curves, uncertainties, ptoretastr, pltstr, flavors, plotvspt=False, plot_qcd_DY=True):
    addc = additional_uncertainty_curves
    fig, ax = plt.subplots()

    xvals = ptvals if plotvspt else etavals
    old_uncs = []
    if plot_qcd_DY:
        for samp in ['QCD', 'DY']:    
            old_unc = ax.plot(xvals, (uncertainties[samp](etavals, ptvals)[:,0]-1)*100, '-.', markersize=0, linewidth=1.0,
                    **color_scheme[samp], alpha=0.6)
            ax.plot(xvals, (HerPy_differences[samp][0]-addc['Rref'])*100, linestyle=(2, (4, 2)), label=samp, markersize=0,
                    linewidth=1.2, **color_scheme[samp])
            old_uncs.append(old_unc[0])

    for flav in ['g', 'c', 'b', 'q']:
        color = color_scheme[flav] #if flav!='q' else color_scheme['ud']
        old_unc = ax.plot(xvals, (uncertainties[flav](etavals, ptvals)[:,0]-1)*100, '-.', markersize=0, linewidth=1.0,
                **color, alpha=0.6)
        old_uncs.append(old_unc[0])

    for flav in flavors:
        lab = legend_dict[flav]
        ax.plot(xvals, (addc[f'{flav}100']-addc['Rref'])*100, label=lab, markersize=0, linewidth=1.2,
                **color_scheme[flav])
    
    ax.hlines(0, ax.get_xlim()[0], ax.get_xlim()[1],color="gray",
        linewidth=1, alpha=0.4)

    legend1 = ax.legend(old_uncs, ['']*len(old_uncs), loc='upper right', bbox_to_anchor=(0.72, 1), handlelength=1.5, title='Run 1', title_fontsize=10)
    leg2 = ax.legend(ncol=1, loc='upper left', bbox_to_anchor=(0.67, 1), handlelength=0.9, title='Run 2' , title_fontsize=10)#, title='assembled\nfrom QCD', title_fontsize=10)
    ax.add_artist(legend1)
    xlabel = r'$p_{T}$ (GeV)' if plotvspt else r'$\eta$'
    ax.set_xlabel(xlabel)
    ylabel = 'JEC uncertainty (%)'
    ax.set_ylabel(ylabel)
    ylim_old = ax.get_ylim()
    if plotvspt:
        ax.set_xscale('log')
        ax.set_xticks([10, 20, 50, 100, 500, 1000, 5000])
        ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
        ax.set_xlim(15,1000)
    else:        
        for HCal_border in JERC_Constants.etaBinsEdges_Win14():
            ax.vlines(HCal_border, ylim_old[0]*3, ylim_old[1]*3, linestyles='--',color="gray",
                linewidth=1,)
        ax.set_xlim(0,5)

    # ax.set_ylim(0.9885,1.0205)
    ylim_pad = (ylim_old[1]-ylim_old[0])*0.4 if plotvspt else (ylim_old[1]-ylim_old[0])*0.62
    ax.set_ylim(ylim_old[0],ylim_old[1]+ylim_pad)
    labtxt = f'{ptoretastr}' #if plotvspt else f'{ptoretastr}'
    labtxt = labtxt.replace("TTBAR", legend_labels["ttbar"]["lab"])
#     labtxt = f'$\eta$ = {etabins_abs[ptoretaidx]}' if plotvspt else f'$p_T$ = {ptbins_c[ptoretaidx]} GeV'
    hep.label.exp_text(text=labtxt, loc=2)
    hep.cms.label(hep_label, loc=0, data=False, ax=ax, rlabel='')
    figdir = "fig/uncertainty"
    if not os.path.exists(figdir):
        os.mkdir(figdir)
    pltstr = pltstr.replace(', ', '_').replace(' ', '_').replace('$', '').replace('=', '_').replace('\eta', 'eta').replace('|', '').replace('<', '').replace('\n', '_')
    figdir = "fig/uncertainty"+"/uncertainty"+pltstr
    if not os.path.exists(figdir):
        os.mkdir(figdir)

    if plotvspt:
        fig_name = figdir+f"/JECuncertainty_vs_pt_eta_{ptoretastr}".replace('.','')
    else:
        fig_name = figdir+f"/JECuncertainty_vs_pt_pt_{ptoretastr}".replace('.','_')
    fig_name = fig_name.replace(', ', '_').replace(' ', '_').replace('$', '').replace('=', '_').replace('\eta', 'eta').replace('|', '').replace('<', '').replace('\n', '_')
    if not plot_qcd_DY:
        fig_name += '_noQCD-DY'

    print("Saving plot with the name = ", fig_name+".pdf / .png")
    plt.savefig(fig_name+'.pdf');
    plt.savefig(fig_name+'.png');
    plt.show()

def plot_HerPydiff(ptvals, HerPy_differences, additional_uncertainty_curves, divideHerPy, etaidx, jeteta_bins, pt_bins, pltstr, flavors, combine_antiflavour):
    addc = additional_uncertainty_curves
    fig, ax = plt.subplots()

    lines = []
    markers = []
    for samp in ['QCD', 'DY', 'TTBAR']:    
#         mc = next(ax._get_lines.prop_cycler)
        line = ax.plot(ptvals, HerPy_differences[samp][0], linestyle=(0, (3.3, 2)), markersize=0, **color_scheme[samp], linewidth=1.2)
        marker = ax.errorbar(pt_bins.centres, HerPy_differences[samp][1], yerr=HerPy_differences[samp][2],
                        linestyle='none', **color_scheme[samp], capsize=1.6, capthick=0.7, linewidth=1.0)
        markers.append(marker[0])

        lines.append(line[0])
        
    # pointsg20 = ax.plot(ptvals, addc['g20q80'], label='DY at 200 GeV', markersize=0, linewidth=1.2, **color_scheme["DY200"])
    for flav in flavors:
        linestyle = '-.' if 'bar' in flav else '-'
        if combine_antiflavour:
            lab = legend_dict[flav]
        else:
            lab = legend_dict_short[flav]
        ax.plot(ptvals, addc[f'{flav}100'], label=lab, markersize=0, linewidth=1.2, linestyle=linestyle, **color_scheme[flav])

    vlinecoord = 1 if divideHerPy else 0
    ax.hlines(vlinecoord ,1, 10000,color="gray",
        linewidth=1, alpha=0.4)

    ax.hlines(addc['g20q80_fixed'], 1, 10000, linestyles='--',color=color_scheme["DY200"]['color'],
        linewidth=1, alpha=0.9, label='DY at 200 GeV')

    leg1_handles = [(ai,bi) for ai, bi, in zip(lines,markers)]
    legend1 = ax.legend(leg1_handles, [legend_dict_short['QCD'], legend_dict_short['DY'], legend_dict_short['TTBAR']], loc="upper left", bbox_to_anchor=(0.02, 0.8), handlelength=1.5) # seg.len=5) #, title='correction', title_fontsize=10)
#     assert False
    leg2 = ax.legend(ncol=1, loc='upper left', bbox_to_anchor=(0.48, 1), handlelength=1.3)#, title='assembled\nfrom QCD', title_fontsize=10)
    ax.add_artist(legend1)
    xlabel = r'$p_{T}$ (GeV)'
    ax.set_xlabel(xlabel);
    ylab_pre = 'R(Her7)/R(Py8)' if divideHerPy else 'R(Her7)-R(Py8)'
    ylabel = r' (median response)'
    ax.set_ylabel(ylab_pre+ylabel);
    ax.set_xscale('log')
    ax.set_xticks([10, 20, 50, 100, 500, 1000, 5000])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    # ax.set_ylim(0.9885,1.0205)
    ax.set_xlim(15,1000)
    ylim_old = ax.get_ylim()
    ylim_pad = (ylim_old[1]-ylim_old[0])*0.3
    ax.set_ylim(ylim_old[0],ylim_old[1]+ylim_pad)
    labtxt = jeteta_bins.idx2plot_str(etaidx)+pltstr.replace("TTBAR", legend_labels["ttbar"]["lab"])
    hep.label.exp_text(text=labtxt, loc=2)
    hep.cms.label(hep_label, loc=0, data=False, ax=ax, rlabel='')

    figdir = "fig/uncertainty"
    if not os.path.exists(figdir):
        os.mkdir(figdir)
    pltstr = pltstr.replace(', ', '_').replace(' ', '_').replace('$', '').replace('=', '_').replace('\eta', 'eta').replace('|', '').replace('<', '').replace('\n', '_')
    figdir = "fig/uncertainty"+"/uncertainty"+pltstr
    if not os.path.exists(figdir):
        os.mkdir(figdir)
        
    add_name = '/Herwig_Pythia_ratio' if divideHerPy else '/Herwig_Pythia_difference'
    if not os.path.exists(figdir+add_name+'/'):
        os.mkdir(figdir+add_name+'/')
    fig_name = figdir+add_name+'/'+add_name+pltstr+'_'+jeteta_bins.idx2str(etaidx)
    print("Saving plot with the name = ", fig_name+".pdf / .png")
    plt.savefig(fig_name+'.pdf');
    plt.savefig(fig_name+'.png');
    plt.show()

# from uncertainty_plotters import color_scheme
# from fileNames.available_datasets import legend_labels
# hep_label = "Private work"
def convert_xpos(x, lims):
    ax, bx = lims
    return (np.log(x) - np.log(ax))/(np.log(bx) - np.log(ax))

def convert_ypos(y, lims):
    ay, by = lims
    return (y - ay)/(by - ay)

def plot_Rref(ptvals, Rdijet0, Rdijet, DY200, Rtot, Rtot_smooth, jeteta_bins, etaidx, pltstr, g_unc):
    fix, ax = plt.subplots()
    
    ax.plot(ptvals, Rdijet0, markersize=0, label=f"dijet, ${jeteta_bins.idx2plot_str(0)[5:]}",color=color_scheme["QCD"]['color'], linewidth=1)
    ax.plot(ptvals, Rdijet, markersize=0, label=f"dijet, {jeteta_bins.idx2plot_str(etaidx)}",color=color_scheme["QCD"]['color'], linestyle='--')
    ax.hlines(0, 1, 10000, linestyles='--', color='black',
        linewidth=0.8, alpha=0.9)
    ax.plot(ptvals, Rtot, markersize=0, label="$R_{ref}$")
    ax.plot(ptvals, Rtot_smooth, markersize=0, label="$R_{ref}$ smoothed")
    ax.plot(ptvals, g_unc, markersize=0, label="$R_{g, Her} - R_{g, Py}$", color=color_scheme["g"]['color'])
    ax.plot(ptvals, g_unc-Rtot_smooth, markersize=0, label="g uncertainty", color=color_scheme["g"]['color'], linestyle='--')
    ax.hlines(DY200, 1, 10000, linestyles='--',color=color_scheme["DY200"]['color'],
        linewidth=1, alpha=0.9, label='DY at 200 GeV')
    
    ax.set_xscale('log')
    ax.set_xticks([10, 20, 50, 100, 500, 1000, 5000])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_xlim(15,1000)
    ax.legend(handlelength=1.3)
    ax.set_xlabel(r'$p_{T}$ (GeV)')
    ax.set_ylabel('R(Her7) - R(Py8)')
    hep.cms.label(hep_label, loc=0, data=False, ax=ax, rlabel='')
    hep.label.exp_text(text=pltstr.replace("TTBAR", legend_labels["ttbar"]["lab"]), loc=2)
    pltstr = pltstr.replace(', ', '_').replace(' ', '_').replace('$', '').replace('=', '_').replace('\eta', 'eta').replace('|', '').replace('<', '').replace('\n', '_')
    figdir = f"fig/uncertainty/uncertainty_{pltstr}"
    if not os.path.exists(figdir):
        os.mkdir(figdir)
    figdir = figdir+"/Rref"
    if not os.path.exists(figdir):
        os.mkdir(figdir)

    ylim_old = ax.get_ylim()
    ylim_pad = (ylim_old[1]-ylim_old[0])*0.2
    ax.set_ylim(ylim_old[0],ylim_old[1]+ylim_pad)
    fig_name = figdir+f"/Rref_{jeteta_bins.idx2str(etaidx)}"
    print("Saving plot with the name = ", fig_name+".pdf / .png")
    plt.savefig(fig_name+'.pdf')
    plt.savefig(fig_name+'.png')
    
#     ax.set_xlim(10,40)
#     print("DY:", DY200)
#     hax = plt.subplot(1,2,2)
#     ax.set_position(good_pos)
#     hax.set_position(good_pos)
#     plt.annotate('', xy=(25, DY200[0]), xytext=(0, 5), 
#             arrowprops=dict(facecolor='black', shrink=0.),
#                 )
#     np.searchsorted(30)
    pt1_idx = np.searchsorted(ptvals, 35)
    pt1 = ptvals[pt1_idx]
#     Rtot_smooth
    xarrow = convert_xpos(pt1, ax.get_xlim())
    yarrow = convert_ypos(Rtot_smooth[pt1_idx], ax.get_ylim())
    dyarrow = convert_ypos(0, ax.get_ylim()) - yarrow
    mc = next(ax._get_lines.prop_cycler)
    ax.arrow(xarrow, yarrow, 0, dyarrow, width=0.01, transform=ax.transAxes, length_includes_head=True, color=mc['color'])
    
    yarrow2 = convert_ypos(g_unc[pt1_idx], ax.get_ylim())
    dyarrow2 = convert_ypos((g_unc-Rtot_smooth)[pt1_idx], ax.get_ylim()) - yarrow2
    ax.arrow(xarrow, yarrow2, 0, dyarrow2, width=0.01, transform=ax.transAxes, length_includes_head=True, color=mc['color'])

    pt3_idx = np.searchsorted(ptvals, 20)
    pt3 = ptvals[pt3_idx]
    xarrow3 = convert_xpos(pt3, ax.get_xlim())
    yarrow3 = convert_ypos(Rdijet0[pt3_idx], ax.get_ylim())
    dyarrow3 = convert_ypos(Rdijet[pt3_idx], ax.get_ylim()) - yarrow3
    mc = next(ax._get_lines.prop_cycler)
    ax.arrow(xarrow3, yarrow3, 0, dyarrow3, width=0.007, transform=ax.transAxes, length_includes_head=True, color=mc['color'])

    yarrow4 = convert_ypos(DY200[0], ax.get_ylim())
    ax.arrow(xarrow3, yarrow4, 0, dyarrow3, width=0.007, transform=ax.transAxes, length_includes_head=True, color=mc['color'])
    fig_name+='_arrows'
    print("Saving plot with the name = ", fig_name+".pdf / .png")
    plt.savefig(fig_name+'.pdf')
    plt.savefig(fig_name+'.png')
    #     hax.set_axis_off()
#     plt.plot(25, DY200[0], 'ko', marker=r'$\downarrow$', markersize=20)
#     arrow = mpatches.FancyArrow(0, 0, 2, 1,
#                                  width=2, length_includes_head=True)
#     ax.add_patch(arrow)
#     ax.arrow(25, DY200[0], 0, -DY200[0], width = 2, length_includes_head=True, head_width = 1.2)

    plt.show()
    
        # ax.plot(mov,markersize=0, label="After smoothing")


# def plot_all_flavor_comparison(num_sample_name,
#                          denom_sample_name, jeteta_bins, ptbins_c, eta_binning_str, fit_samp='J', etaidx=0):
#     ''' Put ratio plots of many all flavors at the same place. Reproduce Fig. 31 in arXiv:1607.03663
#     '''

#     inverse=False   #True if plot corrections, False if plot responses
#     use_recopt=False   #True if use reco pt, False if use gen pt
#     flavors = ['g', 'q' ,'c', 'b'] #, 'unmatched']
    
#     mean_name = "Median"
#     mean_name_std = mean_name+'Std'
#     start = np.searchsorted(ptbins, 15, side='left')
# #     etaidx = np.searchsorted(jeteta_bins_abs, 0, side='left')
    
#     yvals = np.array([read_data2(mean_name, num_sample_name, flav, eta_binning_str)[start:,etaidx] for flav in flavors])
#     stds  = np.array([read_data2(mean_name_std, num_sample_name, flav, eta_binning_str)[start:,etaidx] for flav in flavors])
#     xvals = np.array([read_data2("MeanRecoPt", num_sample_name, flav, eta_binning_str)[start:,etaidx] for flav in flavors])
    
#     yvals_d = np.array([read_data2(mean_name, denom_sample_name, flav, eta_binning_str)[start:,etaidx] for flav in flavors])
#     stds_d  = np.array([read_data2(mean_name_std, denom_sample_name, flav, eta_binning_str)[start:,etaidx] for flav in flavors])
#     xvals_d = np.array([read_data2("MeanRecoPt", denom_sample_name, flav, eta_binning_str)[start:,etaidx] for flav in flavors])
# #     print('etaidx = ', etaidx)

#     corr_loc_Sum20_Py = [f"* * Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs{eta_binning_str}.txt"]
#     corr_loc_Sum20_Her = [f"* * Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Her{eta_binning_str}.txt"]
#     ext = extractor()
#     ext.add_weight_sets(corr_loc_Sum20_Py+corr_loc_Sum20_Her)
#     ext.finalize()
#     evaluator = ext.make_evaluator()
        
#     yvals[(yvals==0) | (np.abs(yvals)==np.inf)] = np.nan
#     yvals_d[(yvals_d==0) | (np.abs(yvals_d)==np.inf)] = np.nan
    
#     ratios = yvals/yvals_d
#     ratio_unc = ((stds / yvals_d)**2 + (yvals/yvals_d**2 * stds_d)**2)**(1/2)
    
#     if not use_recopt:
#         xvals = ptbins_c[start:]    
        
#     etaval = jeteta_bins.centres[etaidx]
#     xvals_cont = np.geomspace(np.min(xvals), np.max(xvals), 100)
#     yvals_cont = np.array([evaluator[f'Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Her{eta_binning_str}_{flav}{fit_samp}'](np.array([etaval]),xvals_cont)
#                            for flav in flavors])
#     yvals_cont_d = np.array([evaluator[f'Summer20UL18_V2_MC_L5Flavor_AK4PFchs{eta_binning_str}_{flav}{fit_samp}'](np.array([etaval]),xvals_cont)
#                        for flav in flavors])
#     if inverse==True:
#         yvals = 1/yvals
#         yvals_d = 1/yvals_d
#         ### Error propagation
#         stds = yvals**2*stds
#         stds_d = yvals_d**2*stds_d
        
#     if inverse==False:
#         yvals_cont = 1/yvals_cont
#         yvals_cont_d = 1/yvals_cont_d
    
    
#     fig, ax = plt.subplots()
#     for axis in [ax.xaxis, ax.yaxis]:
#         axis.set_minor_locator(mpl.ticker.AutoMinorLocator())
        
# #     assert False
#     for yval, std, flav in zip(ratios, ratio_unc, flavors):
#         ax.errorbar(xvals, yval, yerr=std,
#                     linestyle="none", label=legend_dict[flav], **color_scheme[flav],
#                     capsize=1.6, capthick=0.7, linewidth=1.0)
# #         assert not lab == 'unmatched'
       
#     ratios_cont = yvals_cont/yvals_cont_d
# #     ax.set_prop_cycle(None)
#     for yval, flav in zip(ratios_cont, flavors):
#         ax.plot(xvals_cont, yval, markersize=0, **color_scheme[flav])
    
#     ax.set_xscale('log')
#     xlims = ax.get_xlim()
    
#     ax.hlines(1,1, 10000, linestyles='--',color="black", linewidth=1,)
#     ######################## Calculate resonable limits excluding the few points with insane errors
#     recalculate_limits=True
#     if recalculate_limits:
#         yerr_norm = np.concatenate(ratio_unc)
#         y_norm = np.concatenate(ratios)
#         norm_pos = (yerr_norm<0.01) &  (yerr_norm != np.inf) & (y_norm>-0.1)  
#         if ~np.any(norm_pos):
#             print("Cannot determine ylimits")
#             norm_pos = np.ones(len(yerr_norm), dtype=int)
#             raise Exception("Cannot determine ylimits")
#         left_lim = np.min((y_norm-yerr_norm)[norm_pos])
#         right_lim = np.max((yerr_norm+y_norm)[norm_pos])
#         lim_pad = (right_lim - left_lim)/20
#         ax.set_ylim(left_lim-lim_pad, right_lim+lim_pad*8)
    
#     xlabel = r'$p_{T,reco}$ (GeV)' if use_recopt else r'$p_{T,ptcl}$ (GeV)'
#     ax.set_xlabel(xlabel);
#     ylab_pre = 'R(Her7)/R(Py8)' 
#     ylabel = r' (correction)' if inverse else r' (median response)'
#     ax.set_ylabel(ylab_pre+ylabel);
    
#     ax.set_xticks([10, 20, 50, 100, 500, 1000, 5000])
#     ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
#     leg1 = ax.legend(ncol=1)
#     ax.set_xlim(xlims)
    
#     title_name = 'QCD' if fit_samp=='J' else 'ttbar'
#     hep.label.exp_text(text=jeteta_bins.idx2plot_str(eta_idx)+f', {title_name}', loc=0)
    
#     etastr = jeteta_bins.idx2str(eta_idx)
#     fig_name = f'fig/uncertainty/Pythia_Herwig_ratio_{etastr}_using_{fit_samp}_fits'
#     print("Saving plot with the name = ", fig_name)
#     plt.savefig(fig_name+'.pdf');
#     plt.savefig(fig_name+'.png');
#     plt.show();