''' Reads json files with the results for the median responses and it's uncertainties and stores the results as a root file
'''

import json
from JetEtaBins import JetEtaBins, PtBins
from save_json import save_root

samps = ['DY-MG-Her', 'DY-MG-Py', 'QCD-MG-Her', 'QCD-MG-Py', 'Herwig-TTBAR', 'Pythia-TTBAR', 'QCD-Py']

eta_binning  = "Summer20Flavor"  ### HCalPart, CoarseCalo, CaloTowers, one_bin, Summer20Flavor;       
                         ### HCalPart: bin in HCal sectors, CaloTowers: the standard JERC binning,
                         ### CoarseCalo: like 'CaloTowers' but many bins united;    

jeteta_bins = JetEtaBins(eta_binning, absolute=True)
pt_bins = PtBins("MC_truth")
eta_binning_str = '_'+eta_binning if eta_binning != "HCalPart" else ''

all_data = {}
for samp in samps:
    with open(f'out_txt/response_fit_results_L5_{samp}{eta_binning_str}.json', 'r') as json_file:
        all_data[samp] = json.load(json_file)['data']


save_root(all_data, pt_bins, jeteta_bins, 'out_txt/response_fit_results_L5.root')