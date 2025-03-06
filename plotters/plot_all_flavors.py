''' The script compares all the flavour responses vs pt in one plot for Herwig and Pythia.
Run using `python plotters/plot_all_flavors.py`
'''


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from coffea.lookup_tools import extractor
from JetEtaBins import JetEtaBins, PtBins
from collections.abc import Iterable
import mplhep as hep
from uncertainty_plotters import legend_dict_short, hep_label
# from correction_fitter import correction_fitter


# color_scheme_antiflav = {key: cycler_vals
#     for cycler_vals, key in zip(plt.rcParams['axes.prop_cycle'], ['g', 'udbar', 'cbar', 'bbar', 'QCD', 'DY', 'TTBAR', 'DY200', 'unmatched', 'sbar', 'qbar'])}
# color_scheme = color_scheme | color_scheme_antiflav



from data_tools import read_or_recreate_data, read_or_recreate_data_txt
out_txt_path = 'out_txt'

def read_data2(mean_name, samp, tag1):
    return read_or_recreate_data_txt(mean_name, samp, tag1, out_txt_path)

from fileNames.available_datasets import legend_labels
ttbarlab = legend_labels['ttbar']['lab']

from pltStyle import pltStyle
# from scipy.interpolate import CubicSpline
pltStyle(style='hep', font_frac=1.05)

# from uncertainty_plotters import color_scheme
color_scheme = {key: cycler_vals
    for cycler_vals, key in zip(plt.rcParams['axes.prop_cycle'], ['g', 'ud', 'c', 'b', 'QCD', 'DY', 'TTBAR', 'DY200', 'unmatched', 's', 'q', 'cs', 'u', 'd'])}
# color_scheme_antiflav = {key: cycler_vals
#     for cycler_vals, key in zip(plt.rcParams['axes.prop_cycle'], ['g', 'udbar', 'cbar', 'bbar', 'QCD', 'DY', 'TTBAR', 'DY200', 'unmatched', 'sbar', 'qbar'])}
# color_scheme = color_scheme | color_scheme_antiflav

# color_scheme = {key: cycler_vals
#     for cycler_vals, key in zip(plt.rcParams['axes.prop_cycle'], ['g', 'ud', 'c', 'b', 'd', 'u', 's', 'q', 'cs'])}

pltStyle('hep')

def read_data4plot(tag, closure=1, path=out_txt_path):
    '''Read the Mean, MeanStd, Median, MedianStd and RecoPt values of the data with tag `tag`.
    If closure==1, there is no clusure, otherwise it has to be of the same shape as the data read
    '''
    
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
        

def draw_all_flavors(data_dict,
                              function_dict,
                              samples = ['Py', 'Her'],
                              etabins=JetEtaBins("HCalPart", absolute=True), #np.array(JERC_Constants.etaBinsEdges_CaloTowers_full()),
                              ptbins=PtBins("MC_truth"), #np.array(JERC_Constants.ptBinsEdgesMCTruth()),
                              binidx=0, 
                              pt_min = 17,
                              inverse = True,
                              flavors = ['b', 'g', 'u', 'd', 'c', 's', 'all', 'q', 'unmatched']):
    ''' make a plot of the response_vs_pt for all the flavours for the bin given by the index `binidx`.
    data_dict: input data dictionary. The order of the keys for plotting is given by `samples`
    function_dict: splines to be drawn over the points
    samples: list of keys of the data_dict to be plotted
    etabins: the eta binning
    ptbins: the pt binning
    binidx: the index of the eta bin to be plotted
    pt_min: the minimum pt value to be plotted
    inverse: if True, the response is inverted to correction
    flavors: the list of flavours to be plotted
    '''

    plotvspt = True
    use_recopt = True

    ptbins_c = ptbins.centres
    plot_range = range(0, np.searchsorted(ptbins_c,1250)) if 'DY' in "".join(samples) else range(0, np.searchsorted(ptbins_c,2750))
    ptbins_c_plot = ptbins_c[plot_range]
    xvals_c = np.geomspace(ptbins_c_plot.min(), ptbins_c_plot.max(),1000)
    fig, ax = plt.subplots()

    points_ls = []
    for flav in flavors:
        lab = legend_dict_short[flav]

        if len(samples)>1:
            data2 = data_dict[samples[1]][flav]
            yvals, stds, xvals = prepare_points(data2, pt_min=pt_min, binidx=binidx, etabins=etabins, ptbins=ptbins, inverse=inverse, use_recopt=use_recopt)
            # points = ax.errorbar(xvals, yvals,
            #                 yerr=stds,
            #                 linestyle='none', label=lab,  **color_scheme[flav], capsize=1.6, capthick=0.7, linewidth=1.0)
            # breakpoint()
            points2 = ax.errorbar(xvals, yvals,
                            yerr=stds,
                            linestyle='none', mfc='white', markeredgewidth=1.2, **color_scheme[flav], capsize=1.6, capthick=0.7, linewidth=0.8)

        data1 = data_dict[samples[0]][flav]
        yvals, stds, xvals = prepare_points(data1, pt_min=pt_min, binidx=binidx, etabins=etabins, ptbins=ptbins, inverse=inverse, use_recopt=use_recopt)

        points = ax.errorbar(xvals, yvals,
                        yerr=stds,
                        linestyle='none', label=lab,  **color_scheme[flav], capsize=1.6, capthick=0.7, linewidth=0.8)

        if flavors[0] == flav:
            points_ls.append(points[0])
            if len(samples)>1:
                points_ls.append(points2[0])

        yvals_cont = prepare_splines(function_dict[samples[0]][flav], xvals_c, etabins, binidx, inverse)
        # print(yvals_cont)
        ax.plot(xvals_c, yvals_cont, markersize=0, **color_scheme[flav], linewidth=1.0)
        if len(samples)>1:
            yvals_cont = prepare_splines(function_dict[samples[1]][flav], xvals_c, etabins, binidx, inverse)
            ax.plot(xvals_c, yvals_cont, '-.', markersize=0, **color_scheme[flav], linewidth=1.0)
    bins = etabins
    if plotvspt:
        xlabel = r'$p_{T,reco}$ (GeV)' if use_recopt else r'$p_{T,ptcl}$ (GeV)'
    else:
        xlabel = r'$|\eta|$'
    ylabel = r'correction (1/median)' if inverse else r'median response'
    ax.set_ylabel(ylabel)
    if plotvspt:
        ax.set_xscale('log')

    xlims = ax.get_xlim()
    ax.hlines(1,-10, 10000, linestyles='--',color="black",
              linewidth=1,)

    ax.set_xlabel(xlabel)
    if plotvspt:
        ax.set_xticks([10, 20, 50, 100, 500, 1000, 5000]) 
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels([])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

    legend_labs = samples # [samples[0], samples[1]]
    legend1 = ax.legend(points_ls, legend_labs, loc="upper left", bbox_to_anchor=(0.54, 1))
    leg2 = ax.legend(ncol=3, loc='upper left', bbox_to_anchor=(0.54, 0.84))
    ax.add_artist(legend1)
    # leg1 = ax.legend(loc="upper right", ncol=1)

    if not plotvspt:
        xlims = (-0.2, 5.3)
    ax.set_xlim(xlims)
    left_lim, right_lim = ax.get_ylim()
    lim_pad = (right_lim - left_lim)/4.5
    ax.set_xlim(xlims)
    ax.set_ylim(left_lim, right_lim+lim_pad)


    ############################ Adding the CMS labels and saving the plots ######################################
    eta_string = bins.idx2str(binidx) #'_eta'+str(etabins_abs[etaidx])+'to'+str(etabins_abs[etaidx+1])
    fig_corr_name = 'corr' if inverse else 'med_resp'
    fig_x_str = 'pt' if plotvspt else 'eta'
    run_name =  f'{fig_corr_name}_vs_{fig_x_str}'
    dir_name1 = f'fig/{fig_corr_name}_vs_{fig_x_str}_comparisons_all_flav/'
    if not os.path.exists(dir_name1):
        os.mkdir(dir_name1)
        print("Creating directory ", dir_name1)
    if not os.path.exists(dir_name1):
        os.mkdir(dir_name1)
        print("Creating directory ", dir_name1)

    hep.cms.label(hep_label, loc=0, data=False, ax=ax, rlabel='')
    hep.label.exp_text(text=f'{bins.idx2plot_str(binidx)}', loc=2
                       , ax=ax)
    fig_name = dir_name1+'/'+run_name+'_'+eta_string
    print("Saving plot for eta = ", eta_string)
    print("Saving plot with the name = ", fig_name+".pdf / .png")
    plt.savefig(fig_name+'.pdf')
    plt.savefig(fig_name+'.png')
    plt.show()



def prepare_points(data_dict, pt_min=17, binidx=0,
                              etabins=JetEtaBins("HCalPart", absolute=True), #np.array(JERC_Constants.etaBinsEdges_CaloTowers_full()),
                              ptbins=PtBins("MC_truth"), inverse=True, plotvspt = True, use_recopt = False):
    
    start = ptbins.get_bin_idx(pt_min) if plotvspt else 0
    end = ptbins.nbins if plotvspt else etabins.nbins

    data_range = tuple([range(start,end),binidx]) if plotvspt else tuple([binidx, range(start,end)])
    yvals = data_dict["Median"][data_range]
    # end = np.min([len(yv) for yv in yvals])+start #cut, so that all yvals are the same
    stds = data_dict["MedianStd"][data_range]
    reco_pts = data_dict["MeanRecoPt"][data_range] 
    # reco_pts  = np.array([key[2][data_range] if len(key[2].shape)==2 else key[2][data_range] for key in data_dict.values()])
    ### Replacing response values to corrections
    use_recopt=inverse
    if not plotvspt:
        use_recopt=False
        
    yvals[(yvals==0) | (np.abs(yvals)==np.inf)] = np.nan

    if plotvspt:
        xvals = reco_pts if use_recopt else ptbins.centres[start:end]
    else:
        xvals = etabins.centres[start:end]
    validx = (xvals>0)*(yvals>0)

    if inverse==True:
        yvals = 1/yvals
        ### Error propagation
        stds = yvals**2*stds

    return yvals, stds, xvals #, validx    

def prepare_splines(function, xvals_c, bins, binidx, inverse):
    correction_fnc, closure = function
    xv_cont = xvals_c
    if closure is None or closure==1:
        def closure(a,b):
            return np.ones_like(a*b)
    #to ensure that the correction is applied from the right side of the bin border
    binval = bins.edges[binidx]+0.001
    vals_cont = (np.array([binval]), xv_cont)
    if inverse:
        yvals_cont = correction_fnc(*vals_cont)/closure(*vals_cont)
    else:
        yvals_cont = closure(*vals_cont)/correction_fnc(*vals_cont)
    return yvals_cont


# correction_fitter(saveplots = True, do_simfit = False, do_Mikkofit = False, correction_for = 'Py')
# # for reading older older corrections files
# Aut18_samples = ['all', 'b', 'c', 's', 'ud', 'g' ]
# Sum16_samples = ['b', 'c', 's', 'ud', 'g' ]


# list_corr_Sum16 = ["Summer16_07Aug2017_V15_Flavor_Pythia8_MC_"+samp+"_L5Flavor_AK4PFchs.txt" for samp in Sum16_samples]
# list_corr_Sum16.append("Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L2Relative_AK4PFchs.txt")
# list_corr_Sum16.append("Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L3Absolute_AK4PFchs.txt")
# corr_loc_Sum16 = ["* * ../Summer16_07Aug2017_V15_Flavor_Pythia8_MC/"+corr for corr in list_corr_Sum16]
# list_corr_Aut18 = ["Autumn18_V3_MC_Pythia8_"+samp+"_L2Relative_AK4PFchs.txt" for samp in Aut18_samples]
# corr_loc_Aut18 = ["* * ../Autumn18_V3_MC_Pythia8/"+corr for corr in list_corr_Aut18]
# corr_loc_Winter14 = ["* * ../Winter14_V8_MC_L5Flavor/Winter14_V8_MC_L5Flavor_AK5PFchs.txt"]
corr_loc_Sum20Her_Summer20Flavor = ["* * Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_MGQCD_Her.txt"]
corr_loc_Sum20_Summer20Flavor    = ["* * Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_MGQCD.txt"]
# corr_loc_Sum20Her_HCAL = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_simfit_Her_HCalPart.txt"]
# corr_loc_Sum20_HCAL    = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_standalonePy_simfit_HCalPart.txt"]
# corr_loc_Sum20_HCAL    = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_MGQCD_HCalPart.txt"]

# list_corr_Aut18 = ["Autumn18_V3_MC_Herwig7_"+samp+"_L2Relative_AK4PFchs.txt" for samp in Aut18_samples]
# corr_loc_Aut18_Her = ["* * ../Autumn18_V3_MC_Pythia8/"+corr for corr in list_corr_Aut18]
# corr_loc_Sum20_Her = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Her-etaAut18.txt"]

ext = extractor()
ext.add_weight_sets(corr_loc_Sum20Her_Summer20Flavor+corr_loc_Sum20_Summer20Flavor) #+corr_loc_Sum20_Her)
ext.finalize()

evaluator = ext.make_evaluator()

eta_binning  = "Summer20Flavor" #"CoarseCalo"  ### HCalPart, JERC, CoarseCalo, CaloTowers, Summer20Flavor, onebin;
eta_binning_str = '_'+eta_binning if eta_binning != "HCalPart" else ''
eta_binning_str_corr = '_'+eta_binning if eta_binning != "Summer20Flavor" else ''
# load_fit_res=True
# subsamples = ['all', 'b', 'c', 'd', 'u', 's', 'g']
flavors = ['b', 'c', 'd', 'u', 's', 'g']
# flavors = ['ud','b', 's', 'u', 'd', 'g']
# flavors = ['g']
# etabins = np.array(JERC_Constants.etaBinsEdges_Aut18_full())
# etabins = np.array(JERC_Constants.etaBinsEdges_CaloTowers_full())
# etabins_abs = etabins[(len(etabins)-1)//2:]

# flavors = ['ud', 'b']
# flavors = ['all']
tag1 = '_L5_QCD-Py'+eta_binning_str
# tag1_genwt = '_L5_QCD-Py_genwt'+eta_binning_str
tag2 = '_L5_QCD-MG-Py'+eta_binning_str
tag2Her = '_L5_QCD-MG-Her'+eta_binning_str
tag3 = '_L5_Pythia-TTBAR'+eta_binning_str
tag3Her = '_L5_Herwig-TTBAR'+eta_binning_str
tag4 = '_L5_DY-MG-Py'+eta_binning_str
tag4Her = '_L5_DY-MG-Her'+eta_binning_str

# tag1 = '_L5_QCD-MG-Her'+eta_binning_str
# tag2 = '_L5_QCD-Py'+eta_binning_str
# tag3 = '_L5_QCD-MG-Py'+eta_binning_str
# tag3 = '_L5_QCD-divided'


    
if __name__ == "__main__":
    closure_corr = read_data4plot(tag1)['all']['Median'] #divide by Pythia-standalone QCD

    mean_name = "Median"
    mean_name_std = mean_name+'Std'

    etabins = JetEtaBins(eta_binning, absolute=True)
    ptbins = PtBins("MC_truth")


    # data = {
    #         f"QCD MG+Py8": read_data4plot(tag2, closure_corr),
    #         f"QCD MG+Her7": read_data4plot(tag2Her, closure_corr),
    #         }

    # functions = {
    #         f"QCD MG+Py8":    {flav: [evaluator[f'Summer20UL18_V2_MC_L5Flavor_AK4PFchs_MGQCD{eta_binning_str_corr}_{flav}J'], None]  for flav in flavors},
    #         f"QCD MG+Her7":    {flav: [evaluator[f'Summer20UL18_V2_MC_L5Flavor_AK4PFchs_MGQCD{eta_binning_str_corr}_Her_{flav}J'], None]  for flav in flavors}
    #         }

    data = {
            f"{ttbarlab} Pow+Py8": read_data4plot(tag3, closure_corr),
            # f"{ttbarlab} Pow+Her7": read_data4plot(tag3Her, closure_corr),
            }

    functions = {
            f"{ttbarlab} Pow+Py8":    {flav: [evaluator[f'Summer20UL18_V2_MC_L5Flavor_AK4PFchs_MGQCD{eta_binning_str_corr}_{flav}T'], None]  for flav in flavors},
            # f"{ttbarlab} Pow+Her7":    {flav: [evaluator[f'Summer20UL18_V2_MC_L5Flavor_AK4PFchs_MGQCD{eta_binning_str_corr}_Her_{flav}T'], None]  for flav in flavors}
            }



    for k in range(etabins.nbins):
    # for k in etabins.get_bin_idx([0, 1.305, 2.5, 4]):
        # data = {tag:data_to_read[tag][samp] for tag in data_to_read}
        # data = {key:np.array([data[key][mean_name], data[key][mean_name_std], data[key]["MeanRecoPt"]]) for key in data}
    #     for k in range(1):
    #     for k in ptbins.get_bin_idx([20, 35, 150, 400]):
        print('Plotting eta: ', etabins.idx2str(k))

        draw_all_flavors(data,
                                functions,
                                samples = list(data.keys()),
                                etabins=etabins, #np.array(JERC_Constants.etaBinsEdges_CaloTowers_full()),
                                ptbins=ptbins, #np.array(JERC_Constants.ptBinsEdgesMCTruth()),
                                binidx=k, 
                                pt_min = 17,
                                inverse = False,
                                flavors = flavors)