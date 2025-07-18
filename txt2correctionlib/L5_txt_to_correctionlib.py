''' Converting L5 text files to correctionlib format
'''


import json
import correctionlib
correctionlib.register_pyroot_binding()
import gzip

import numpy as np

import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

from coffea.jetmet_tools import JetCorrectionUncertainty #FactorizedJetCorrector
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory

from coffea.lookup_tools import extractor

def get_description(corr):
    '''Given the corr string (e.g. bT) obtain a description for the L5 corrections
    '''
    flav = corr[-2]
    corr_type = corr[-1]
    if flav=='a':
        flav_txt = 'all'
    else:
        flav_txt = flav
    
    if corr_type=='T':
        corr_type_txt = 'ttbar'
    elif corr_type=='J':
        corr_type_txt = 'QCD'
    description = f"L5Flavour for AK4PFchs {flav_txt} jets obtained from {corr_type_txt} events"
    return description

def main():
    ######### Open the text file using coffea and then convert and use it's parser
    ######### then format it according to correctionlib
    txtpath = 'Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs.txt'
    extL5 = extractor()
    extL5.add_weight_sets([
                f"* * {txtpath}"
    ])
    extL5.finalize()
    evaluatorL5 = extL5.make_evaluator()
    print(f"Loaded coffea L5 evaluator from {txtpath}")

    # Append the corrections on top of an existing correctionlib file
    correctionsjsonUL18 = '/afs/cern.ch/user/a/anpotreb/top/JERC/JMECoffea/jet_jerc.json.gz' # for UL18
    new_json_name = '/afs/cern.ch/user/a/anpotreb/top/JERC/JMECoffea/jet_jerc_withL5.json.gz'

    with gzip.open(correctionsjsonUL18,mode='rb') as f4:
        data_bytes4 = f4.read()
    json_str4 = data_bytes4.decode('utf-8')
    data_UL18 = json.loads(json_str4)
    f4.close()
    print(f"Loaded existing correctionlib json from {correctionsjsonUL18}")
        
    # evalu = evaluatorL5["Summer20UL18_V2_MC_L5Flavor_AK4PFchs_bT"]
    for corr in dir(evaluatorL5):
        ### copy an already existing entry, replace the relevent info and append to the list
        L5_entry = data_UL18['corrections'][17].copy() #select Summer19UL18_V5_MC_L2Relative_AK4PFchs
        L5_entry['name'] = corr
        L5_entry['description'] = get_description(corr)
        
        
        ### use coffea evaluator to parse the txt file and convert the input into the relevent data
        evalu = evaluatorL5[corr]
        ### in coffea paramter 'x' in formula is replaced with JetPt, so we have to replace it back
        ### we also have to set by hand the upper and lower limits (the first two parameters)
        formula = evalu._formula_str
        formula = formula.replace('JetPt','max(min(x,[1]),[0])')
        for ii in range(10): ## large enough number to cover all the arguments;
            ## to do: implement a more logical upper limit
            ## the first two parameters are limits, shift the formula accordingly
            formula = formula.replace(f'p{ii}', f'[{ii+2}]')

        parms = (np.array(evalu._parms)[:,:,0])
        mins_maxs = np.array([evalu._eval_clamp_mins['JetPt'],evalu._eval_clamp_maxs['JetPt']])[:,:,0]
        parms_ordered = (np.concatenate([mins_maxs, parms]).T).tolist()
        
        ### replace the data
        data_new = L5_entry['data'].copy()
        new_edges = (evalu._bins['JetEta']).tolist()
        data_new['edges'] = new_edges
        
        new_content = []
        for parms_ii in parms_ordered:
            new_content_tmp = L5_entry['data']['content'][0].copy()
            new_content_tmp['expression'] = formula
            new_parms = list(parms_ii)
            new_content_tmp['parameters'] = new_parms
            new_content.append(new_content_tmp)
        
        data_new['content'] = new_content
        
        ### append the current json
        L5_entry['data'] = data_new
        data_UL18['corrections'].append(L5_entry)
        
    with gzip.open(new_json_name, mode='w') as fout4: # mode='wb'
        fout4.write(json.dumps(data_UL18, indent=4).encode('utf-8'))
    fout4.close()
    print(f"New json file with L5 corrections created: {new_json_name}")

if __name__ == "__main__":
    main()
    import test_correctionlib_L5 as test_correctionlib_L5
    test_correctionlib_L5.main()  # Run the test script to verify the corrections
    