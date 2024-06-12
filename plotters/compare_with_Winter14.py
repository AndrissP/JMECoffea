import numpy as np
import matplotlib.pyplot as plt
# import matplotlib as mpl
from make_comparison_plot import make_comparison_plot
from coffea.lookup_tools import extractor
from JetEtaBins import JetEtaBins, PtBins
from collections.abc import Iterable
from correction_fitter import correction_fitter


from data_tools import read_or_recreate_data, read_or_recreate_data_txt
out_txt_path = 'out_txt'

def read_data2(mean_name, samp, tag1):
    return read_or_recreate_data_txt(mean_name, samp, tag1, out_txt_path)

# from common_binning import JERC_Constants

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

def read_data4plot(tag, closure=1, path=out_txt_path):
    '''Read the Mean, MeanStd, Median, MedianStd and RecoPt values of the data with tag `tag`.
    If closure==1, there is no clusure, otherwise it has to be of the same shape as the data read
    '''
#     file_path = f'../out_txt/fit_results_L5_{tag}.json'
#     with open(file_path, 'r') as json_file:
#         json_data = json.load(json_file)
    
    data = read_or_recreate_data(tag, path)['data']

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
    

# for reading older older corrections files
Aut18_samples = ['all', 'b', 'c', 's', 'ud', 'g' ]
Sum16_samples = ['b', 'c', 's', 'ud', 'g' ]


list_corr_Sum16 = ["Summer16_07Aug2017_V15_Flavor_Pythia8_MC_"+samp+"_L5Flavor_AK4PFchs.txt" for samp in Sum16_samples]
list_corr_Sum16.append("Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L2Relative_AK4PFchs.txt")
list_corr_Sum16.append("Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L3Absolute_AK4PFchs.txt")
corr_loc_Sum16 = ["* * Summer16_07Aug2017_V15_Flavor_Pythia8_MC/"+corr for corr in list_corr_Sum16]
list_corr_Aut18 = ["Autumn18_V3_MC_Pythia8_"+samp+"_L2Relative_AK4PFchs.txt" for samp in Aut18_samples]
corr_loc_Aut18 = ["* * Autumn18_V3_MC_Pythia8/"+corr for corr in list_corr_Aut18]
corr_loc_Winter14 = ["* * Winter14_V8_MC_L5Flavor/Winter14_V8_MC_L5Flavor_AK5PFchs.txt"]
# corr_loc_Sum20Her_Summer20Flavor = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_MGQCD_Her.txt"]
# corr_loc_Sum20_Summer20Flavor    = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_MGQCD.txt"]
corr_loc_Sum20Her_CoarseCalo = ["* * Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_CoarseCalo.txt"]
# corr_loc_Sum20_HCAL    = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_standalonePy_simfit_HCalPart.txt"]
# corr_loc_Sum20_HCAL    = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_MGQCD_HCalPart.txt"]

list_corr_Aut18 = ["Autumn18_V3_MC_Herwig7_"+samp+"_L2Relative_AK4PFchs.txt" for samp in Aut18_samples]
corr_loc_Aut18_Her = ["* * Autumn18_V3_MC_Pythia8/"+corr for corr in list_corr_Aut18]
# corr_loc_Sum20_Her = ["* * ../Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs_Her-etaAut18.txt"]

rerun_fits = False
import os

rerun_fit = rerun_fits
save_fit_plots = False
filename = corr_loc_Sum20Her_CoarseCalo[0][4:]
if (not rerun_fits) and (not os.path.exists(filename)):
    print(f"The evaluator file {filename} does not exist.")
    create_file = input("Do you want to rerun the fit? (yes/no): ")
    if create_file.lower() == 'yes' or create_file.lower() == '' or create_file.lower() == 'y':
        rerun_fit = True
        plots_input = input("Do you want to save the plots of the fits? (no/yes): ")
        save_fit_plots = not (plots_input.lower() == 'no' or plots_input.lower() == '' or plots_input.lower() == 'n')

if rerun_fit:
    print("Rerunning the fit")
    correction_fitter(saveplots = save_fit_plots, do_simfit = False, do_Mikkofit = False, correction_for = 'Py', eta_binning='CoarseCalo')
    # correction_fitter(saveplots = save_fit_plots, do_simfit = False, do_Mikkofit = False, correction_for = 'Py', eta_binning='CoarseCalo')

ext = extractor()
ext.add_weight_sets(corr_loc_Sum16+corr_loc_Aut18+corr_loc_Winter14
                    # +corr_loc_Sum20+corr_loc_Sum20Her
                    # +corr_loc_Sum20_fineta+corr_loc_Sum20Her_fineta
                    # +corr_loc_Sum20_coarseCalo+corr_loc_Sum20Her_coarseCalo
                    # +corr_loc_Sum20_JERC+corr_loc_Sum20Her_JERC
                    +corr_loc_Aut18_Her
                    # +corr_loc_Sum20standPy_JERC+corr_loc_Sum20standPy_CoarseCalo
                    +corr_loc_Sum20Her_CoarseCalo) #+corr_loc_Sum20_Her)
ext.finalize()

evaluator = ext.make_evaluator()

eta_binning  = "CoarseCalo" #"CoarseCalo"  ### HCalPart, JERC, CoarseCalo, CaloTowers, Summer20Flavor, onebin;
eta_binning_str = '_'+eta_binning if eta_binning != "HCalPart" else ''
eta_binning_str_corr = '_'+eta_binning if eta_binning != "Summer20Flavor" else ''
# load_fit_res=True
# subsamples = ['all', 'b', 'c', 'd', 'u', 's', 'g']
flavors = ['b', 'c', 'd', 'u', 's', 'g', 'ud', 'q', 'all']
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

closure_corr = read_data4plot(tag1)['all']['Median'] #divide by Pythia-standalone QCD

mean_name = "Median"
mean_name_std = mean_name+'Std'

closure_Aut18 = evaluator['Autumn18_V3_MC_Pythia8_all_L2Relative_AK4PFchs']
closure_Sum16 = evaluator['Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L2Relative_AK4PFchs']
etabins = JetEtaBins(eta_binning, absolute=True)
ptbins = PtBins("MC_truth")


for samp in flavors:
    samp_Aut18 = samp
#     samp_Sum20 = '_'+samp
    samp_Aut18 = '_ud' if samp_Aut18=='u' or samp_Aut18=='d' or samp_Aut18=='q' else '_'+samp_Aut18
    
    evo_Her = evaluator[f'Autumn18_V3_MC_Herwig7{samp_Aut18}_L2Relative_AK4PFchs']
    evo = evaluator[f'Autumn18_V3_MC_Pythia8{samp_Aut18}_L2Relative_AK4PFchs']
    if samp_Aut18=='_all':
        evo2 = evaluator['Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L2Relative_AK4PFchs']
    else:
        evo2 = evaluator[f'Summer16_07Aug2017_V15_Flavor_Pythia8_MC{samp_Aut18}_L5Flavor_AK4PFchs']
        
    if samp_Aut18!='_all':
        evo3 = evaluator[f'Winter14_V8_MC_L5Flavor_AK5PFchs{samp_Aut18}J']
        evo4_Py = evaluator[f'Summer20UL18_V2_MC_L5Flavor_AK4PFchs{eta_binning_str}{samp_Aut18}J']
    else:
        evo3 = evaluator[f'Winter14_V8_MC_L5Flavor_AK5PFchs_aJ']
        evo4_Py = evaluator[f'Summer20UL18_V2_MC_L5Flavor_AK4PFchs{eta_binning_str}_aJ']

    data_to_read = {
            "QCD Py8": read_data4plot(tag1, closure_corr),
           }
    
    functions = {
            "Summer20, Py":    [evo4_Py, None],

            "Autumn18":    [evo, closure_Aut18],
            "Run 1":       [evo3, None],
           }
    
    for k in range(etabins.nbins):
    # for k in etabins.get_bin_idx([0, 1.305, 2.5, 4]):
        data = {tag:data_to_read[tag][samp] for tag in data_to_read}
        data = {key:np.array([data[key][mean_name], data[key][mean_name_std], data[key]["MeanRecoPt"]]) for key in data}
#     for k in range(1):
#     for k in ptbins.get_bin_idx([20, 35, 150, 400]):
        print('Plotting subsample: ', samp)
        print('Eta: ', k)
        if not np.any(data[list(data.keys())[0]][2][:,k]>-0.1):
#             continue
#         if not np.any(median2[:,k]>-0.1):
            print("All median values are none")
            continue
        
        make_comparison_plot(data, 
                                  functions,
                                  etabins, ptbins,
                                  binidx=k, flav=samp, ratio_name='*/ QCD Py8', ratio_ylim=[0.98,1.02], plotvspt=True, reset_colors=True,
                            pt_min=15)
        
3;
