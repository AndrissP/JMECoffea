#CoffeaJERCProcessor_L5.py
''' coffea processor for calculating the jet energy response in bins of pt_ptcl and jet_eta.
    The processor makes a separate histogram for each jet flavor.
output: a dictionary over datasets of dictionaries over histograms.
output histograms: ptresponse histogram, pt_reco histogram for each flavor and the cutflow
''' 

from memory_profiler import profile
from common_binning import JERC_Constants
import JERCProcessorcuts as cuts
from JERCProcessorcuts import apply_jetNevent_cuts

# workaround to get a locally installed coffea and awkwrd version using lch on lxplus
# comment out or replace the path if I happened to forget to remove these lines before pushing:
import sys
import os
coffea_path = '/afs/cern.ch/user/a/anpotreb/top/JERC/coffea/'
if not os.path.exists(coffea_path):
    raise ValueError(f"The path to the coffea installation does not exist. Please supply the correct path or comment out this line if using the environment path. The provided path is: {coffea_path}.")
if coffea_path not in sys.path:
    sys.path.insert(0,coffea_path)
# 
# ak_path = '/afs/cern.ch/user/a/anpotreb/top/JERC/local-packages/'
# if ak_path not in sys.path:
#         sys.path.insert(0,ak_path)
# sys.path.insert(0,'/afs/cern.ch/user/a/anpotreb/top/JERC/JMECoffea')
# print("sys path = ", sys.path)
# from os import listdir
# listdir('.')
# listdir('./coffea')

from coffea import processor
import numpy as np
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory
from coffea.lookup_tools import extractor
import correctionlib

from LHE_flavour import get_LHE_flavour, get_LHE_flavour2
import hist
import awkward as ak

# from coffea import some_test_func
# some_test_func.test_func()

# manual_bins = [400, 500, 600, 800, 1000, 1500, 2000, 3000, 7000, 10000]
ptbins = np.array(JERC_Constants.ptBinsEdgesMCTruth())
etabins = np.array(JERC_Constants.etaBinsEdges_CaloTowers_full())

class CutMaker():
    def __init__(self, processor_config):   
        self.cfg = processor_config
        self.jetflavour = processor_config["jetflavour"]

        ext = extractor()
        ext.add_weight_sets([
            "* * Summer20UL18_V2_MC/Summer20UL18_V2_MC_L1FastJet_AK4PFchs.txt",
            "* * Summer20UL18_V2_MC/Summer20UL18_V2_MC_L2Relative_AK4PFchs.txt",
            "* * Summer20UL18_V2_MC/Summer20UL18_V2_MC_L3Absolute_AK4PFchs.txt",
#             "* * Summer20UL18_V2_MC/Summer19UL18_V5_MC_L2L3Residual_AK4PFchs.txt", #Doesn't do anything but for transparancy I add it
        ])
        ext.finalize()

        jec_stack_names = ["Summer20UL18_V2_MC_L1FastJet_AK4PFchs",
                           "Summer20UL18_V2_MC_L2Relative_AK4PFchs", 
                           "Summer20UL18_V2_MC_L3Absolute_AK4PFchs",
#                            "Summer19UL18_V5_MC_L2L3Residual_AK4PFchs",
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
        
        
        self.jet_factory = CorrectedJetsFactory(name_map, jec_stack)
        self.flavor2partonNr = {'b':5,
                           'c':4,
                           's':3,
                           'u':2,
                           'd':1,
                           'bbar':-5,
                           'cbar':-4,
                           'sbar':-3,
                           'ubar':-2,
                           'dbar':-1,
                           'g':21,
                        #    'FSR_gluon':21,
                        #    'ISR_gluon':21, ## will be split into FSR/ME gluons and ISR later
                           'unmatched':0,
                           }
        
        self.flavors = self.flavor2partonNr.keys() #['b', 'c', 'u', 'd', 's', 'g', 'bbar', 'cbar', 'ubar', 'dbar', 'sbar', 'untagged']

        path_to_PU_weights = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/2018_UL/puWeights.json.gz"
        self.pucorr = correctionlib.CorrectionSet.from_file(path_to_PU_weights)
            
    # @property
    # def accumulator(self):
    #     return self._accumulator

    # @profile    
    # def for_memory_testing(self):
    #     a=1
        
    # @profile    
    def process(self, events, dataset=''):

        ############ Define the histograms ############
        flavors = self.flavors        
        # flavour_axis = hist.axis.StrCategory(flavors, growth=False, name="jet_flav", label=r"jet_flavour")  ###not completelly sure if defining an axis is better than doing through a dictionary of flavors. See, https://github.com/CoffeaTeam/coffea/discussions/705
        pt_gen_axis = hist.axis.Variable(ptbins, name="pt_gen", overflow=True, underflow=True, label=r"$p_{T,gen}$")
        ptresponse_axis = hist.axis.Regular( 100, 0, 2.5, overflow=True, underflow=True, name="ptresponse", label="RECO / GEN response")
        jeteta_axis = hist.axis.Variable(etabins, name="jeteta", label=r"Jet $\eta$")

        #self.for_memory_testing()
        output = {'ptresponse_'+samp:hist.Hist(pt_gen_axis, ptresponse_axis, jeteta_axis, storage="weight", name="Counts")
                  for samp in flavors}
        #self.for_memory_testing()
        # To calculate the mean recopt, store only the sums of values for each bin.
        # Thus it takes much less space than storing the whole reco_pt distribution.  
        for samp in flavors:
            output['reco_pt_sumwx_'+samp] = hist.Hist(pt_gen_axis, jeteta_axis, storage="weight", name="Counts")
        #self.for_memory_testing()
        
        cutflow_axis = hist.axis.StrCategory([], growth=True, name="cutflow", label="Cutflow Scenarios")
        output['cutflow_events'] = hist.Hist(cutflow_axis, storage="weight", label="N events")
        output['cutflow_jets'] = hist.Hist(cutflow_axis, storage="weight", label="N jets")
        output['sum_weights'] = hist.Hist(cutflow_axis, storage="weight", label="sum of weights")

    
        selectedEvents, reco_jets, cutflow_evts, cutflow_jets = apply_jetNevent_cuts(events, self.cfg, output['cutflow_events'], output['cutflow_jets'], self)
        output['cutflow_events'] = cutflow_evts
        output['cutflow_jets'] = cutflow_jets
        gen_jets = reco_jets.matched_gen

        jet_flavour = reco_jets[self.jetflavour] 

        ########### Split the samples into jet flavours ###############
        shapes_jets = ak.num(gen_jets.pt) #for event weights
        gen_jetpt  = ak.flatten(gen_jets.pt).to_numpy( allow_missing=True)
        gen_jeteta = ak.flatten(gen_jets.eta).to_numpy( allow_missing=True)
        jetpt      = ak.flatten(reco_jets.pt).to_numpy( allow_missing=True)

        ptresponse_np = jetpt / gen_jetpt #/ self.closure_corr[correction_pos_pt, correction_pos_eta]
        
        if 'LHEWeight' not in selectedEvents.fields: ### no LHEWeight.originalXWGTUP stored in standalone Pythia8 but Generator.weight instead
            gen_weights = selectedEvents.Generator.weight
        else:
            gen_weights = selectedEvents.LHEWeight.originalXWGTUP
        
        if self.cfg["use_gen_weights"]:
            weights = gen_weights
        else:
            weights = np.ones(len(selectedEvents))

        if self.cfg["use_pu_weights"]:
            weights = weights*self.pucorr['Collisions18_UltraLegacy_goldenJSON'].evaluate(selectedEvents.Pileup.nTrueInt, "nominal")
        weights_jet = np.repeat(weights, shapes_jets)


        masks = {flav: ak.flatten((jet_flavour == self.flavor2partonNr[flav] )).to_numpy( allow_missing=True)
                 for flav in flavors if 'unmatched' not in flav}
        from functools import reduce
        ## find the jets that are not taggeed as any of the flavours
        # masks['FSR_gluon'] = ak.flatten((reco_jets["partonFlavour"] == 21) & (reco_jets["LHE_flavour2"] != 21)).to_numpy( allow_missing=True)
        # masks['ISR_gluon'] = ak.flatten((reco_jets["partonFlavour"] == 21) & (reco_jets["LHE_flavour2"] == 21)).to_numpy( allow_missing=True)
        masks['unmatched'] = reduce(lambda x, y: x+y, masks.values()) == 0


        ptresponses     = { flav: ptresponse_np[masks[flav]]        for flav in flavors }
        gen_jetpts      = { flav: gen_jetpt[masks[flav]]            for flav in flavors }
        gen_jetetas     = { flav: gen_jeteta[masks[flav]]           for flav in flavors }
        jetpts          = { flav: jetpt[masks[flav]]                for flav in flavors }
        # if self.cfg["use_weights"]==True:
        weights_jet     = { flav: weights_jet[masks[flav]]             for flav in flavors }


        ########### Filling of the histograms ###############
        for flav in flavors:
            output['ptresponse_'+flav].fill(pt_gen=gen_jetpts[flav],
                                              jeteta=gen_jetetas[flav],
                                              ptresponse=ptresponses[flav],
                                              weight=weights_jet[flav]
                                             )
            
            output['reco_pt_sumwx_'+flav].fill(pt_gen=gen_jetpts[flav],
                                                 jeteta=gen_jetetas[flav],
                                                 weight=jetpts[flav]*weights_jet[flav]
                                                )
        #self.for_memory_testing()
        output['sum_weights'].fill(cutflow='sum_weights', weight=ak.sum(gen_weights))

        return {dataset: output}, reco_jets

from CoffeaJERCProcessor_L5_config import processor_config, processor_dependencies
cut_maker = CutMaker(processor_config)

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

events_scale1 = NanoEventsFactory.from_root(
    '/eos/cms/store/user/anpotreb/pi_vs_pibar_all/job_391_scale_1_JME-RunIISummer20UL18NanoAODv9-00002.root', #fileslist[0],
    schemaclass=NanoAODSchema.v6,
    entry_start=0,
    entry_stop=1000,
).events()


events = NanoEventsFactory.from_root(
    '/eos/cms/store/user/anpotreb/pi_vs_pibar_all/job_0_scale_0_JME-RunIISummer20UL18NanoAODv9-00002.root', #fileslist[0],
    schemaclass=NanoAODSchema.v6,
    entry_start=0,
    entry_stop=1000,
).events()

events_scale5 = NanoEventsFactory.from_root(
    '/eos/cms/store/user/anpotreb/pi_vs_pibar_all/job_1173_scale_5_JME-RunIISummer20UL18NanoAODv9-00002.root', #fileslist[0],
    schemaclass=NanoAODSchema.v6,
    entry_start=0,
    entry_stop=1000,
).events()
events.metadata['dataset'] = 'not_scaled_pion'
# breakpoint()
# events_scale1.metadata['dataset'] = 'not_scaled_pion'
# events_scale5.metadata['dataset'] = 'not_scaled_pion'
_, jets = cut_maker.process(events)
_, jets_scaled1 = cut_maker.process(events_scale1)
_, jets_scaled5 = cut_maker.process(events_scale5)

print('jetspt = ', jets.pt)

matchedJets = ak.cartesian([jets, jets_scaled5])
deltaR = matchedJets.slot0.delta_r(matchedJets.slot1)
matchedJets = matchedJets[deltaR < 0.2]
jets = matchedJets.slot0
jets_scaled = matchedJets.slot1

import hist
from hist import Hist

from common_binning import JERC_Constants
ptbins = np.array(JERC_Constants.ptBinsEdgesMCTruth())
pt_gen_axis = hist.axis.Variable(ptbins, name="pt_gen", overflow=True, underflow=True, label=r"$p_{T,gen}$")
pt_gen_scaled_axis = hist.axis.Variable(ptbins, name="pt_gen_scaled", overflow=True, underflow=True, label=r"$p_{T,gen, sc}$")
h_jetpt = Hist(pt_gen_axis, pt_gen_scaled_axis)
h_genjetpt = Hist(pt_gen_axis, pt_gen_scaled_axis)
h_jetpt.fill(ak.flatten(jets.pt), ak.flatten(jets_scaled.pt))
h_genjetpt.fill(ak.flatten(jets.matched_gen.pt), ak.flatten(jets_scaled.matched_gen.pt))

h = h_jetpt
proj = h.project("pt_gen_scaled").values()
scale = np.array([proj]*len(proj))

from plotters.pltStyle import pltStyle
import mplhep
import matplotlib.pyplot as plt
pltStyle(style='hep') #, font_frac=1.40
plt.rcParams['figure.subplot.right'] = plt.rcParams['figure.subplot.right']-0.04
plt.rcParams['figure.subplot.left'] = plt.rcParams['figure.subplot.left']*0.7

fig, ax = plt.subplots()
mplhep.hist2dplot(h/scale, ax=ax)
ax.set_xscale('log')
ax.set_yscale('log')
print(f'Saving figure to test.pdf')
fig.savefig('test.pdf')


resp_axis = hist.axis.Regular( 100, 0, 2.5, overflow=True, underflow=True, name="ptresponse", label="R")
resp_axis_scaled = hist.axis.Regular( 100, 0, 2.5, overflow=True, underflow=True, name="ptresponse_scaled", label="$R_{scaled}$")
# resp_axis = hist.axis.Variable(ptresponse_axis, name="response", overflow=True, underflow=True, label=r"$R$")
# resp_scaled_axis = hist.axis.Variable(ptresponse_axis, name="response_scaled", overflow=True, underflow=True, label=r"$R_{scaled}}$")
h_resp_2d = Hist(resp_axis, resp_axis_scaled)
# h_genjetpt = Hist(resp_axis, resp_axis_scaled)
h_resp_2d.fill(ak.flatten(jets.pt/jets.matched_gen.pt), ak.flatten(jets_scaled.pt/jets_scaled.matched_gen.pt))
# h_genjetpt.fill(ak.flatten(jets.matched_gen.pt), ak.flatten(jets_scaled.matched_gen.pt))


print("corr coef = ", np.corrcoef(ak.flatten(jets.pt/jets.matched_gen.pt), ak.flatten(jets_scaled.pt/jets_scaled.matched_gen.pt)))

h = h_resp_2d
proj = h.project("ptresponse_scaled").values()
scale = np.array([proj]*len(proj))

fig, ax = plt.subplots()
mplhep.hist2dplot(h/scale, ax=ax)
# ax.set_xscale('log')
# ax.set_yscale('log')
print(f'Saving figure to response_2d.pdf')
fig.savefig('response_2d.pdf')