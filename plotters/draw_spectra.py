from fileNames.available_datasets import legend_labels

include_unmatched = True
from uncertainty_helpers import get_output, sum_output, combine_flavors
from uncertainty_plotters import plot_spectra
from helpers import get_xsec_dict, sum_neg_pos_eta, rebin_hist, hist_mult
from fileNames.available_datasets import dataset_dictionary
from JetEtaBins import JetEtaBins, PtBins
from plotters.pltStyle import pltStyle
import matplotlib.pyplot as plt
# import mplhep as hep
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
