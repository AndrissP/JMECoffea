### test L5 corrections on correctionlib by comparing them to the results from coffea and txt files
import correctionlib

import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

from coffea.jetmet_tools import JetCorrectionUncertainty #FactorizedJetCorrector
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory

from coffea.lookup_tools import extractor

import numpy as np

def main():
    print("Testing L5 corrections from correctionlib against coffea evaluator")
    ### test results on jets from a NanoAOD file
    filename = 'file:///eos/user/a/anpotreb/lxplus/JERC/PFNano_QCD_2000toInf_MGHer/job_13_Nano_herwig7_QCD_HT2000toInf_PF_NANO.root'

    events = NanoEventsFactory.from_root(
        filename,
        schemaclass=NanoAODSchema.v6,
        entry_start=0,
        entry_stop=10000,
    #      entrystart=195440, entrystop=202420
    ).events()

    ### load the txt file L5 corrections
    txtpath = 'Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs.txt'
    extL5 = extractor()
    extL5.add_weight_sets([
                f"* * {txtpath}"
    ])
    extL5.finalize()
    evaluatorL5 = extL5.make_evaluator()
    print(f"Loaded coffea L5 evaluator from {txtpath}")

    ### load the correctionlib file with L5 corrections
    # Append the corrections on top of an existing correctionlib file
    new_json_name = '/afs/cern.ch/user/a/anpotreb/top/JERC/JMECoffea/jet_jerc_withL5.json.gz'
    evaluatorCL = correctionlib.CorrectionSet.from_file(new_json_name)
    print(f"Loaded correctionlib evaluator from {new_json_name}")

    ### compare the correction factor from both sources on all the jets
    print("Testing the L5 corrections from correctionlib against coffea evaluator")
    jets = events.Jet
    for key in evaluatorCL.keys():
        if 'L5Flavor' in key:
            print("key:", key)
            evaluatorL5_corlib = evaluatorCL[key]
            corr_corlib = evaluatorL5_corlib.evaluate(ak.flatten(jets.eta),ak.flatten(jets.pt))
            
            corr_txt = ak.flatten(evaluatorL5[key](jets.eta,jets.pt))
            
    #         stats.describe(corr_txt)
    #         stats.describe(corr_corlib)
            # print("corr_txt:", corr_txt.to_numpy())
            # print("corr_corlib:", corr_corlib)
            agreeing_thres = 1e-5
            print("N of disagreeing corrections: ", np.sum((corr_corlib-corr_txt)/corr_txt > agreeing_thres))
            print("N of agreeing corrections: ", np.sum((corr_corlib-corr_txt)/corr_txt < agreeing_thres))

    # uncomment to print out the cocontents of the correctionlib file
    # for corr in evaluatorCL.values():
    #     print(f"Correction {corr.description}")
    #     print(f"  It's name {corr.name} has {len(corr.inputs)} inputs")
    #     for ix in corr.inputs:
    #         print(f"   Input {ix.name} ({ix.type}): {ix.description}")

    compareL2L3 = False
    if compareL2L3:
        ######### compare L1-L3 corrections from correctionlib and coffea #####################
        ### apply L1-L3 corrections from correctionlib
        ## apply L1
        evaluatorL1 = evaluatorCL["Summer19UL18_V5_MC_L1FastJet_AK4PFchs"]
        for ix in evaluatorL1.inputs:
            print(f"   Input {ix.name} ({ix.type}): {ix.description}")
            
        jets = events.Jet
        jets['pt_raw'] = (1 - jets['rawFactor']) * jets['pt']     #raw pt. pt before the corrections are applied to data
        jets['mass_raw'] = (1 - jets['rawFactor']) * jets['mass']
        # jets['pt_gen'] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
        jets['rho'] = ak.broadcast_arrays(events.fixedGridRhoFastjetAll, jets.pt)[0]

        corrL1 = evaluatorL1.evaluate(ak.flatten(jets.area), ak.flatten(jets.eta), ak.flatten(jets.pt_raw), ak.flatten(jets.rho))
        jets["pt_postL1"] = ak.unflatten(ak.flatten(jets.pt_raw)*corrL1, ak.num(jets))
        jets["mass_postL1"] = ak.unflatten(ak.flatten(jets.mass_raw)*corrL1, ak.num(jets))
        # jets["pt_postL1"] = jets.pt_raw

        ## apply L2
        evaluatorL2L3Res = evaluatorCL["Summer19UL18_V5_MC_L2L3Residual_AK4PFchs"]
        for ix in evaluatorL2L3Res.inputs:
            print(f"   Input {ix.name} ({ix.type}): {ix.description}")

        corr = evaluatorL2L3Res.evaluate(ak.flatten(jets.eta),ak.flatten(jets.pt_postL1)) # corrections are 1, so do not apply

        evaluatorL2Rel = evaluatorCL["Summer19UL18_V5_MC_L2Relative_AK4PFchs"]
        for ix in evaluatorL2Rel.inputs:
            print(f"   Input {ix.name} ({ix.type}): {ix.description}")

        corr = evaluatorL2Rel.evaluate(ak.flatten(jets.eta),ak.flatten(jets.pt_postL1))

        jets["pt_correctionlib"] = ak.unflatten(ak.flatten(jets.pt_postL1)*corr, ak.num(jets)) #corrected_pt

        evaluatorL3 = evaluatorCL["Summer19UL18_V5_MC_L3Absolute_AK4PFchs"]
        for ix in evaluatorL3.inputs:
            print(f"   Input {ix.name} ({ix.type}): {ix.description}")
            
        #for MC corrections are 1 so do not apply
        corr2 = evaluatorL3.evaluate(ak.flatten(jets.eta),ak.flatten(jets.pt)) 

        ### apply L1-L3 corrections from coffea text files

        ext = extractor()
        ext.add_weight_sets([
            "* * Summer19UL18_V5_MC/Summer19UL18_V5_MC_L1FastJet_AK4PFchs.txt",
            "* * Summer19UL18_V5_MC/Summer19UL18_V5_MC_L2Relative_AK4PFchs.txt",
            "* * Summer19UL18_V5_MC/Summer19UL18_V5_MC_L3Absolute_AK4PFchs.txt",
        #             "* * Summer20UL18_V2_MC/Summer20UL18_V2_MC_L5Flavor_AK4PFchs.txt",
        ])
        ext.finalize()

        jec_stack_names = [ "Summer19UL18_V5_MC_L1FastJet_AK4PFchs",
                        "Summer19UL18_V5_MC_L2Relative_AK4PFchs", 
                        "Summer19UL18_V5_MC_L3Absolute_AK4PFchs",
        #                    "Summer20UL18_V2_MC_L5Flavor_AK4PFchs",
                        ]

        evaluator = ext.make_evaluator()        
        jec_inputs = {name: evaluator[name] for name in jec_stack_names}
        jec_stack = JECStack(jec_inputs)

        name_map = jec_stack.blank_name_map
        name_map['JetPt'] = 'pt'
        name_map['JetMass'] = 'mass'
        name_map['JetEta'] = 'eta'
        name_map['JetA'] = 'area'
        name_map['ptGenJet'] = 'pt_gen'
        name_map['ptRaw'] = 'pt_raw'
        name_map['massRaw'] = 'mass_raw'
        name_map['Rho'] = 'rho'

        jet_factory = CorrectedJetsFactory(name_map, jec_stack)

        selected_jets = events.Jet
        selected_jets['pt_raw'] = (1 - selected_jets['rawFactor']) * selected_jets['pt']     #raw pt. pt before the corrects applied to data
        selected_jets['mass_raw'] = (1 - selected_jets['rawFactor']) * selected_jets['mass']
        # selected_jets['pt_gen'] = ak.values_astype(ak.fill_none(selected_jets.matched_gen.pt, 0), np.float32)
        selected_jets['rho'] = ak.broadcast_arrays(events.fixedGridRhoFastjetAll, selected_jets.pt)[0]
        events_cache = events.caches[0]

        reco_jets = jet_factory.build(selected_jets, lazy_cache=events_cache)

        print("Comparing the results from correctionlib and coffea")
        print("jet pt from correctionlib:", ak.flatten(reco_jets.pt).to_numpy())
        print("jet pt from coffea:", jets.pt_correctionlib)

        '''
        Checked that:
            - [x] Correction is the same for L2L3
            - [ ] Jets are not the same after both of the corrections after L1
            - [x] The correction factor after L1 is the same
            - [x] Checked that in CorrectedJetsFactory.py l.201, the correction factor multiplies ptRaw, the same as I do
            - [ ] My corrections are applied on 20UL18, but the they are not available yet. People should use 19UL18
        '''

### function to apply the L5 corrections from the txt files based on the partonFlavour
### todo: to be finished
# dictionary showing what correction to apply on what partonFlavour. In the real life L5 should be applied on jets tagged experimentally
L5_dict = { 5: 'b',
        21: 'g',
            'rest':'q'
        }

def apply_L5(evaluatorL5, jets, L5_dict:dict, corr_name:str='Summer20UL18_V2_MC_L5Flavor_AK4PFchs', corr_type:str='T', debug:bool=False):
    '''Apply the L5 corrections on the jets based on the partonFlavour of the jet. Return: corrected jets.
    '''
    if debug:
        print("Check non-corrected pt: ", jets.pt)
    ### it is neccessary to flatten the jets to substitute only a subarray of jets.pt
    num_jets = ak.num(jets)
    len_jets = sum(num_jets)
    jets_flat = ak.flatten(jets)
    ## placeholder for correction factor
    corr_factor = np.ones(len_jets)

    
    for flavNr in L5_dict:
        flav = L5_dict[flavNr]
        if flavNr == 'rest' :
            ### needs to be finished to add all the other flavour
            flavNr = 1
            
        ## get jets of the specific flavour
        flavour_sel = np.abs(jets_flat.partonFlavour)==flavNr
        jets_sel = jets_flat[flavour_sel]

        corr_factor[flavour_sel] = evaluatorL5[corr_name+f'_{flav}{corr_type}'](jets_sel.eta,jets_sel.pt)
        if debug:
            print(f"Adding flavour: {flav}, len selected jets: {len(jets_sel)}")
            print("Check the correction factor:", corr_factor[flavour_sel])
    corr_factor = ak.unflatten(corr_factor, num_jets)

    new_pt = jets.pt*corr_factor
    new_mass = jets.mass*corr_factor
    jets["pt"] = new_pt
    jets["mass"] = new_mass
    if debug:
        print("Check corrected pt:", jets.pt)
    return jets
if __name__ == "__main__":
    main()





