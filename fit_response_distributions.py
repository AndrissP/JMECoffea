#!/usr/bin/env python
    # coding: utf-8

"""
Analyses the output from  CoffeaJERCProcessor_L5 and fits the response distributions. 

Author(s): Andris Potrebko (RTU)
"""

# import sys
# coffea_path = '/afs/cern.ch/user/a/anpotreb/top/JERC/coffea/'
# if coffea_path not in sys.path:
#     sys.path.insert(0,coffea_path)

# ak_path = '/afs/cern.ch/user/a/anpotreb/top/JERC/local-packages/'

# if ak_path not in sys.path:
#     sys.path.insert(0,ak_path)

from coffea import util
import numpy as np
import argparse

# import inspect
import matplotlib.pyplot as plt
import hist
import warnings

from plotters.pltStyle import pltStyle
pltStyle(style='hep') #, font_frac=1.40
plt.rcParams['figure.subplot.left'] = plt.rcParams['figure.subplot.left']*1.4
# plt.rcParams['font.size'] = plt.rcParams['font.size']/0.98
plt.rcParams['figure.dpi'] = 150
import os

### import subpackages
import helpers as h
import plotters.plot_makers as plot_makers
# from common_binning import JERC_Constants
from fileNames.available_datasets import dataset_dictionary, dataset_labels
from plotters.plot_cutflow import plot_cutflow
from save_json import save_json

def fit_response_distributions(data_tag='Pythia-TTBAR', config=None):
    ''' The script fits the response histograms (or calculates the medians) and creates the `txt` files with the fit results
    (one file for each `Mean`, `MeanStd`, `Median`, `MedianStd`, `MeanRecoPt`)
    '''
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.realpath(__file__))

    if config is None:
        config = {}

    # Use the dictionary values, or the defaults if they are not in the dictionary
    test_run = config.get('test_run', False)
    saveplots = config.get('saveplots', False)
    save_fit_results = config.get('save_fit_results', True)
    fit_results_format = config.get('fit_results_format', 'json')
    if fit_results_format not in ('json', 'csv'):
        raise ValueError("fit_results_format must be either 'json' or 'csv'")
    combine_antiflavour = config.get('combine_antiflavour', False)
    eta_binning = config.get('eta_binning', 'HCalPart')
    pt_binning = config.get('pt_binning', 'MC_truth')
    sum_neg_pos_eta_bool = config.get('sum_neg_pos_eta_bool', True)
    tag_Lx = config.get('tag_Lx', '_L5')
    add_tag = config.get('add_tag', '')
    fit_tag = config.get('fit_tag', '')
    flavors = config.get('flavors', ['b', 'ud', 'all', 'g', 'c', 's', 'q', 'u', 'd', 'unmatched'])
    pt_to_fit = config.get('pt_to_fit', None)
    eta_to_fit = config.get('eta_to_fit', None)

    ################ End of the parameters of the run and switches  #########################
    # ### Do some logic with the input partameters and the rest of parameters of the run
    
    tag_full = tag_Lx+'_'+data_tag+add_tag
    if test_run:
        tag_full = tag_full+'_test'
    outname = os.path.join(script_dir, 'out', 'CoffeaJERCOutputs'+tag_full+'.coffea')
    
    tag_fit_res = tag_full
    
    if eta_binning != "HCalPart":
        tag_fit_res=tag_full+'_'+eta_binning
    if pt_binning != "MC_truth":
        tag_fit_res=tag_fit_res+'_pt-'+pt_binning
    combine_antiflavour_txt = '_split_antiflav' if not combine_antiflavour else ''
    tag_fit_res += combine_antiflavour_txt+fit_tag

    # if not os.path.exists("out"):
    #     os.mkdir("out")
            
    fig_path = os.path.join(script_dir, 'fig')
    if not os.path.exists(fig_path):
        os.mkdir(fig_path)
    cutflow_path = os.path.join(fig_path, 'cutflow')
    os.makedirs(cutflow_path, exist_ok=True)
        
    if test_run and not os.path.exists(script_dir+"/test"):
        os.mkdir(script_dir+"/test/")
        os.mkdir(script_dir+"/test/fig")
        
    out_txt_path = script_dir+"/out_txt" if not test_run else script_dir+"/test/out_txt"
    if not os.path.exists(out_txt_path):
        os.mkdir(out_txt_path)

    # ### End of do some logic with the input partameters and the rest of parameters of the run
    ################ Load the histograms and scale them according to their cross-sections  #########################
    output = util.load(outname)
    print("Loaded histograms from: ", outname)
    xsec_dict, legend_label = h.get_xsec_dict(data_tag, dataset_dictionary)
    
    keys = output.keys()
    try: ## in older files cutflow was not split into cutflow for jets and events. Should be removed in the future
        Nev = {key: output[key]['cutflow_events']['all_events'].value for key in keys}
    except KeyError:
        Nev = {key: output[key]['cutflow']['all_events'].value for key in keys}
    scale_factors = h.hist_div(xsec_dict, Nev)
    all_histo_keys = output[next(iter(output.keys()))].keys()
    hists_merged = {histo_key:h.sum_subhist(output, histo_key, scale_factors) for histo_key in all_histo_keys }
    
    # ### Fit responses
    
    # Define some global variables for the fit
    from JetEtaBins import JetEtaBins, PtBins
    
    jeteta_bins = JetEtaBins(eta_binning)
    pt_bins = PtBins(pt_binning)
    fiteta_bins = JetEtaBins(eta_binning, absolute=True) if sum_neg_pos_eta_bool else jeteta_bins

    if pt_to_fit is None:
        pt_bins_to_fit = range(pt_bins.nbins)
    elif len(pt_to_fit)==1:
        pt_bins_to_fit = [pt_bins.get_bin_idx(pt_to_fit[0])]
    elif len(pt_to_fit)==2:
        pt_bins_to_fit = range(pt_bins.get_bin_idx(pt_to_fit[0]), pt_bins.get_bin_idx(pt_to_fit[1])+1)

    if eta_to_fit is None:
        eta_bins_to_fit = range(fiteta_bins.nbins)
    elif len(eta_to_fit)==1:
        eta_bins_to_fit = [fiteta_bins.get_bin_idx(eta_to_fit[0])]
    elif len(eta_to_fit)==2:
        eta_bins_to_fit = range(fiteta_bins.get_bin_idx(eta_to_fit[0]), fiteta_bins.get_bin_idx(eta_to_fit[1])+1)
    
    def fit_responses(hists, flavor='b', saveplots = None, scaled_hist=None):
        ''' Extract the jet flavor `flavor` from the histogram dictionary `hists` and fit in all the eta and pt bins.
        Add `scaled_hist` if to produce the response distributions with all the samples stacked up.
        Return a dictionary of ["Mean", "MeanStd", "Median", "MedianStd", "MeanRecoPt"] values.
        
        '''
        warnings.filterwarnings('ignore')  ### filter out the many fit warnings
        if saveplots==None:
            saveplots = False if test_run or eta_binning != "HCalPart" else True
            
        response_hists = {}
        recopt_hists = {}
        if not scaled_hist==None:
            for sample in scaled_hist:
                response_hist, recopt_hist = h.add_flavors(scaled_hist[sample], flavor, combine_antiflavour) 
                response_hist = h.rebin_hist(response_hist, 'jeteta' , jeteta_bins.edges)
                recopt_hist = h.rebin_hist(recopt_hist, 'jeteta' , jeteta_bins.edges)

                if sum_neg_pos_eta_bool==True:
                    response_hist = h.sum_neg_pos_eta(response_hist)
                    recopt_hist = h.sum_neg_pos_eta(recopt_hist)
                response_hists[sample] = response_hist
                recopt_hists[sample] = recopt_hist
            
        response_hist, recopt_hist = h.add_flavors(hists, flavor, combine_antiflavour)
        
        # breakpoint()
        # print("response hist, values = ", response_hist.values()[5,45:55,jeteta_bins.get_bin_idx(4.8)])
        # response_hist.values()[:,:,1] = 0
        # recopt_hist.values()[:,1] = 0
        # response_hist.values()[:,:,-1] = 0
        # recopt_hist.values()[:,-1] = 0
        response_hist = h.rebin_hist(response_hist, 'jeteta' , jeteta_bins.edges)
        recopt_hist = h.rebin_hist(recopt_hist, 'jeteta' , jeteta_bins.edges)

        response_hist = h.rebin_hist(response_hist, 'pt_gen' , pt_bins.edges)
        recopt_hist = h.rebin_hist(recopt_hist, 'pt_gen' , pt_bins.edges)

        if sum_neg_pos_eta_bool==True:
            response_hist = h.sum_neg_pos_eta(response_hist)
            recopt_hist = h.sum_neg_pos_eta(recopt_hist)
        # print("response hist, values = ", response_hist[33j,45:55,1.305j:1.566j].values().flatten())
        # print("response hist, values = ", response_hist[33j,45:55,1.566j:1.74j].values().flatten())
        results = {key:np.zeros((pt_bins.nbins, fiteta_bins.nbins))
                      for key in ["Mean", "MeanStd", "Median", "MedianStd", "MeanRecoPt"]  }
                                
        N_converged = 0
        N_little_ev = 0
        N_failed = 0
    
        FitFigDir1 = fig_path+'/responses/responses'+tag_full
        if saveplots and not os.path.exists(fig_path+'/responses/'):
            os.mkdir(fig_path+'/responses/')
        if saveplots and not os.path.exists(FitFigDir1):
            os.mkdir(FitFigDir1)
        
        FitFigDir = FitFigDir1+'/response_pt_eta_'+flavor+tag_full
        if saveplots:
            if not os.path.exists(FitFigDir):
                os.mkdir(FitFigDir)
            print("Response fit plots will be saved under ", FitFigDir)
        elif not saveplots:
            print("Response fit plots won't be saved")
        
        for i in pt_bins_to_fit:
            for k in eta_bins_to_fit:
                if not scaled_hist==None:
                    histos = {sample: response_hists[sample][i, :, k] for sample in response_hists}
                    histos2plot = {key[10:]:histos[key] for key in histos.keys()}
                    h_stack = hist.Stack.from_dict(histos2plot)
                
                histo = response_hist[i, :, k]
                histopt = recopt_hist[i, k]            
                try:
                    Neff = histo.sum().value**2/(histo.sum().variance)
                except ZeroDivisionError:
                    Neff = histo.sum().value**2/(histo.sum().variance+1e-20)
    
                median, medianstd = h.get_median(histo, Neff) #, x_range=[0, 2]
                
                ##################### Mean of the pt_reco  ######################
                ### (The mean includes events that potentially had ptresponse in the second peak at low pt)
                ### No way to distinguish it if only x*weights are saved instead of the whole histogram.
                mean_reco_pt = histopt.value/np.sum(histo.values())
                
                ####################### Fitting ############################
                p2, cov, chi2, Ndof, status, fitlims = h.fit_response(histo, Neff, Nfit=3, sigma_fit_window=1.5)
                if status == 1:
                    N_converged += 1
                elif status == -1:
                    N_little_ev += 1
                else:
                    N_failed += 1
                    
                
                ####################### Store the results ############################
                results["Mean"][i,k] = p2[1]
                results["MeanStd"][i,k] = np.sqrt(np.abs(cov[1,1]))
                results["Median"][i,k] = median
                results["MedianStd"][i,k] = medianstd
                results["MeanRecoPt"][i,k] = mean_reco_pt
    
        ####################### Plotting ############################
                if  saveplots:
                    figName = FitFigDir+'/ptResponse'+pt_bins.idx2str(i)+fiteta_bins.idx2str(k)
                    hep_txt = pt_bins.idx2plot_str(i)+'\n'+fiteta_bins.idx2plot_str(k)+'\n'+f'{flav} jet' 
            
                    txt2print = ('\n'+r'Mean = {0:0.3f}$\pm${1:0.3f}'.format(p2[1], np.sqrt(cov[1,1]))
                                     + '\nWidth = {0:0.3f}$\pm${1:0.3f}'.format(np.abs(p2[2]), np.sqrt(cov[2,2]))
                                     + '\n'+r'Median = {0:0.3f}$\pm${1:0.3f}'.format(median, medianstd)
                                     + '\n'+r'$\chi^2/ndof$ = {0:0.2g}/{1:0.0f}'.format(chi2, Ndof)
                                     + '\n'+r'Neff = {0:0.3g}'.format(Neff))
                    
                    # txt2print = ('')
                    # hep_txt=''
                    plot_makers.plot_response_dist(histo, p2, fitlims,
                                       figName, dataset_name=legend_label, hep_txt=hep_txt, txt2print=txt2print, print_txt=True, simple_depiction=False)              
                    if not scaled_hist==None:
                        plot_makers.plot_response_dist_stack(h_stack, p2, fitlims,
                                                 figName+'stack', hep_txt=hep_txt, print_txt=False )
    
        print("fit summary: ")
        print(f"N bins converged = {N_converged}; N bins not fit because of too little data = {N_little_ev}; N bins not converged = {N_failed}")
        warnings.filterwarnings('default')
        
        return results  
    
    # ### Run fitting for each sample
    
    medians = []
    medianstds = []
    if not combine_antiflavour:
        flavors = np.concatenate([[flav, flav+'bar'] if flav in h.barable_flavors else [flav] for flav in flavors ])
    print('-'*25)
    print('-'*25)
    print(f'Starting to fit each flavor in: {flavors}')
    result_each_flav = {}
    for flav in flavors:
        print('-'*25)
        print('-'*25)
        print('Fitting flavor: ', flav)
        result = fit_responses(hists_merged, flav, saveplots=saveplots) #scaled_hist
        result_each_flav[flav] = result
        medians.append(result["Median"][0][0])
        medianstds.append(result["MedianStd"][0][0])
        if save_fit_results and fit_results_format == 'csv':
            for key in result:
                h.save_data(result[key], key, flav, tag_fit_res, pt_bins.centres, fiteta_bins.edges, out_txt_path)

    #     print("result = ", result)
        # median = result["Median"]
        # medianStd = result["MedianStd"] 
        # meanstd = np.sqrt(result["MeanStd"])
                
        if eta_binning=="onebin": #or fine_etabins:
            plot_makers.plot_corrections_eta(result["Median"], result["MedianStd"], pt_bins, fiteta_bins.centres, tag_fit_res, flav, plotptvals=[20, 35, 150, 400])
        else:
            plot_makers.plot_corrections(result, pt_bins.centres, fiteta_bins, tag_fit_res, flav, plotetavals=[0, 1.305, 2.5, 3.139], plotmean=True)
    #         plotters.plot_corrections_eta(result["Median"], result["MedianStd"], pt_bins, fiteta_bins.centres, tag_fit_res, flav, plotptvals=[20, 35, 150, 400])

    if save_fit_results and fit_results_format == 'json':
        save_json(result_each_flav, pt_bins, fiteta_bins, out_txt_path+'/response_fit_results'+tag_fit_res+'.json')


    print('-'*25)
    print('-'*25)
    print('Saving cutflow')
    rc_bottom_def = plt.rcParams['figure.subplot.bottom']
    plt.rcParams['figure.subplot.bottom'] = 0.39
    tag_cutflow = tag_full[4:]
    # hist1 = output[list(keys)[0]] ### plotting only the first

    scale_cutflow = {key: 1 for key in output}
    sum_cutflow = {histo_key:h.sum_subhist(output, histo_key, scale_cutflow) for histo_key in ['cutflow_events', 'cutflow_jets'] }

    ### normalize the cutflow histograms to the number of events
    def normalize_cutflow(hists_cutflow):
        for key in ['cutflow_events', 'cutflow_jets']:
            hists_cutflow[key].variances()[:] = hists_cutflow[key].values()
        return hists_cutflow
    
    sum_cutflow = normalize_cutflow(sum_cutflow)

    dataset_label = dataset_labels[tag_cutflow] if tag_cutflow in dataset_labels else tag_cutflow 
    print("dataset label: ", dataset_label)
    # sum_cutflow['cutflow_events'].variances()[:] = sum_cutflow['cutflow_events'].values()
    # sum_cutflow['cutflow_jets'].variances()[:] = sum_cutflow['cutflow_jets'].values()
    jets_pet_ev = sum_cutflow['cutflow_jets']/sum_cutflow['cutflow_events']['all_events'].value
    cutflow_tag_path = os.path.join(cutflow_path, tag_cutflow)
    plot_cutflow(sum_cutflow['cutflow_events'], tag_cutflow, ylab='N events', fig_name='cutflow_Nevents', title_name=dataset_label, figdir=cutflow_tag_path)
    plot_cutflow(sum_cutflow['cutflow_jets'], tag_cutflow, ylab='N jets', fig_name='cutflow_Njets', title_name=dataset_label, figdir=cutflow_tag_path)
    plot_cutflow(jets_pet_ev, tag_cutflow, ylab='N jets/N events', fig_name='cutflow_Njets_per_ev', title_name=dataset_label, figdir=cutflow_tag_path)
    # if 'cutflow_events' in hists_cutflow:
        # drawing only for the first sample as in the hist_merged, the total number of events are normalized
    if len(list(keys))>1:
        for key in keys:
            cutflow = normalize_cutflow(output[key])
            dataset_label = dataset_labels[key] if key in dataset_labels else key
            dataset_label.replace('_', ' ')

            jets_pet_ev = cutflow['cutflow_jets']/cutflow['cutflow_events']['all_events'].value
            plot_cutflow(cutflow['cutflow_events'], key, ylab='N events', fig_name='cutflow_Nevents', title_name=dataset_label, figdir=cutflow_tag_path)
            plot_cutflow(cutflow['cutflow_jets'], key, ylab='N jets', fig_name='cutflow_Njets', title_name=dataset_label, figdir=cutflow_tag_path)
            plot_cutflow(jets_pet_ev, key, ylab='N jets/N events', fig_name='cutflow_Njets_per_ev', title_name=dataset_label, figdir=cutflow_tag_path)

    # else:
    #     print("cutflow histograms cannot be drawn because the cutflow isn't split into events and jets. Potentially not all events are selected")
    #     cutflow = hist1['cutflow'][['all_events', 'selected_events',  'events passing the lepton selection',  'events, alpha cut' ]]
    #     plot_cutflow(cutflow, list(keys)[0], ylab='N events', fig_name='cutflow_Nevents')

    #     cutflow = hist1['cutflow'][[ 'all_jets', 'gen_matched', 'jets, tight lepton id',
    #                                     'jets, dR cut with leptons', 'jetpt cut',
    #                                     'alpha cut; leading jets','iso jets']]
    #     plot_cutflow(cutflow, list(keys)[0], ylab='N jets', fig_name='cutflow_Njets')
    plt.rcParams['figure.subplot.bottom'] = rc_bottom_def
    print('-----'*10)
    print("All done. Congrats!")
  
def main(argv=None):
    def parse_bool(value):
        normalized = value.lower()
        if normalized in ("true", "yes", "1", "on"):
            return True
        if normalized in ("false", "no", "0", "off"):
            return False
        raise argparse.ArgumentTypeError(
            f"expected a boolean value (true/false), got {value!r}"
        )

    def add_boolean_option(option, dest, default, help_text, add_negative_option=True):
        parser.add_argument(
            option, dest=dest, nargs="?", const=True, type=parse_bool,
            default=default, metavar="BOOL", help=help_text,
        )
        if add_negative_option:
            parser.add_argument(
                option.replace("--", "--no-", 1), dest=dest, action="store_false",
                default=default, help=argparse.SUPPRESS,
            )

    def split_list_values(values, option):
        parsed_values = []
        for value in values:
            parts = [part.strip() for part in value.split(",")]
            if any(not part for part in parts):
                parser.error(f"--{option} contains an empty list item")
            parsed_values.extend(parts)
        return parsed_values

    parser = argparse.ArgumentParser(
        description="Fit response distributions from CoffeaJERCProcessor output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-d", "--data", dest="data_tags", nargs="+", action="extend",
        default=argparse.SUPPRESS,
        help="dataset tag(s) from the available datasets; may provide several (default: QCD-Py-Run3-2024-Summer24)",
    )
    add_boolean_option(
        "--test-run", "test_run", False,
        "read processor output produced with test_run=True",
    )
    add_boolean_option(
        "--saveplots", "saveplots", True,
        "save response-distribution plots (many bins can take time and space)",
        add_negative_option=False,
    )
    add_boolean_option(
        "--save-fit-results", "save_fit_results", True,
        "save the fit results",
    )
    parser.add_argument(
        "--fit-results-format", choices=("json", "csv"),
        default="json",
        help="format used to save fit results",
    )
    add_boolean_option(
        "--combine-antiflavour", "combine_antiflavour", True,
        "combine particle and antiparticle flavors",
    )
    parser.add_argument(
        "--eta-binning", default="HCalPart",
        help=("eta binning for response fits: HCalPart, CaloTowers, CoarseCalo, "
              "JERC, Summer20Flavor, or onebin"),
    )
    parser.add_argument(
        "--pt-binning", default="MC_truth",
        help="pt binning for response fits (for example MC_truth, Uncert, Coarse, onebin)",
    )
    add_boolean_option(
        "--sum-neg-pos-eta", "sum_neg_pos_eta_bool", True,
        "combine negative- and positive-eta bins",
    )
    parser.add_argument(
        "--tag-lx", dest="tag_Lx", default="_L5",
        help="correction level tag (L23 is not currently supported)",
    )
    parser.add_argument(
        "--add-tag", default="",
        help="extra tag for output files when run settings differ",
    )
    parser.add_argument(
        "--fit-tag", default="",
        help="extra tag for fit strategy variations",
    )
    parser.add_argument(
        "--pt-to-fit", nargs="+", type=float, default=None, metavar="PT",
        help="fit one pt value or a pt range (low high); omit to fit all pt bins",
    )
    parser.add_argument(
        "--eta-to-fit", nargs="+", type=float, default=None, metavar="ETA",
        help="fit one eta value or an eta range (low high); omit to fit all eta bins",
    )
    parser.add_argument(
        "--flavors", nargs="+",
        default=['b', 'ud', 'all', 'g', 'c', 's', 'q', 'u', 'd', 'unmatched'],
        help="flavors to fit",
    )
    args = parser.parse_args(argv)
    config = vars(args).copy()
    data_tags = config.pop("data_tags", None)
    if data_tags is None:
        data_tags = ["QCD-Py-Run3-2024-Summer24"]
    else:
        data_tags = split_list_values(data_tags, "data")
    config["flavors"] = split_list_values(config["flavors"], "flavors")
    for option in ("pt_to_fit", "eta_to_fit"):
        values = getattr(args, option)
        if values is not None and len(values) not in (1, 2):
            parser.error(f"--{option.replace('_', '-')} expects one value or a low/high pair")
    for data_tag in data_tags:
        fit_response_distributions(data_tag=data_tag, config=config)


if __name__ == "__main__":
    main()
