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
from helpers import legend_str_to_filename

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

        #self.for_memory_testing()
        output = {}
        
        cutflow_axis = hist.axis.StrCategory([], growth=True, name="cutflow", label="Cutflow Scenarios")
        output['cutflow_events'] = hist.Hist(cutflow_axis, storage="weight", label="N events")
        output['cutflow_jets'] = hist.Hist(cutflow_axis, storage="weight", label="N jets")
        output['sum_weights'] = hist.Hist(cutflow_axis, storage="weight", label="sum of weights")

        selectedEvents, reco_jets, cutflow_evts, cutflow_jets, leptons = apply_jetNevent_cuts(events, self.cfg, output['cutflow_events'], output['cutflow_jets'], self, dataset)
        output['cutflow_events'] = cutflow_evts
        output['cutflow_jets'] = cutflow_jets
        gen_jets = reco_jets.matched_gen

        return reco_jets, cutflow_evts, leptons

from CoffeaJERCProcessor_L5_config import processor_config, processor_dependencies

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from plotters.pltStyle import pltStyle
import mplhep as hep
import matplotlib.pyplot as plt
pltStyle(style='hep') #, font_frac=1.40
plt.rcParams['figure.subplot.right'] = plt.rcParams['figure.subplot.right']-0.04
plt.rcParams['figure.subplot.left'] = plt.rcParams['figure.subplot.left']*0.7

# partonNr2legend = { 5:'b',
#                     4:'c',
#                     3:'s',
#                     2:'u',
#                     1:'d',
#                     -5:'\overline{b)',
#                     -4:'\overline{c)',
#                     -3:'\overline{s)',
#                     -2:'\overline{u)',
#                     -1:'\overline{d)',
#                     21:'g',
#                 #    'FSR_gluon':21,
#                 #    'ISR_gluon':21, ## will be split into FSR/ME gluons and ISR later
#                     0:'unmatched',
#                     }

flavor2partonNr = {'b':5,
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

events = NanoEventsFactory.from_root(
    'root://cmsxrootd.fnal.gov///store/mc/RunIISummer20UL18NanoAODv9/DYJetsToLL_M-50_TuneCH3_13TeV-madgraphMLM-herwig7/NANOAODSIM/20UL18JMENano_HerwigJetPartonBugFix_106X_upgrade2018_realistic_v16_L1v1-v1/50000/A5653C07-C418-2442-B1F2-55688FC0CBED.root',
    # 'root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18NanoAODv9/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/20UL18JMENano_106X_upgrade2018_realistic_v16_L1v1-v1/40000/BCB3E2FC-D575-0341-A211-5C9A8D8798B9.root', #fileslist[0],
    schemaclass=NanoAODSchema.v6,
    entry_start=0,
    entry_stop=10000,
).events()

cut_maker = CutMaker(processor_config)
reco_jets, cutflow_evts, leptons = cut_maker.process(events, dataset='Pythia-TTBAR')
# processor_config["jetflavour"] = 'LHE_flavour2'
# cut_maker2 = CutMaker(processor_config)
# reco_jets_LHE, cutflow_evts = cut_maker2.process(events, dataset='Pythia-TTBAR')
# breakpoint()

all_ev = cutflow_evts['all_events'].value
dir_name = 'plotters/fig/misc'
if not os.path.exists(dir_name):
    os.mkdir(dir_name)

dphi = hist.new.Reg(20,-4,4).Double()
Z_px = leptons[:, 0].px + leptons[:, 1].px
Z_py = leptons[:, 0].py + leptons[:, 1].py
lepton_phi = np.arctan2(Z_py, Z_px)

delta_phi = reco_jets[:, 0].phi - lepton_phi
delta_phi = (delta_phi + np.pi) % (2 * np.pi) - np.pi  # ensure the result is between -pi and pi
dphi.fill(delta_phi.flatten())
fig, ax = plt.subplots()
dphi.plot1d(histtype='fill', alpha=0.75) #, label='all jets')
hep.cms.label("Private work", loc=0, data=False, ax=ax, rlabel='')
hep.label.exp_text(text=f'$\Delta\phi$ between leading jet and Z', loc=2, ax=ax)
ax.set_xlabel(f'$\Delta\phi$')
ax.set_ylabel('Events')
# ax.legend()
figname = dir_name+f'/dphi_leading_jet_Z'
plt.savefig(figname+'.pdf')
plt.savefig(figname+'.png')
print(f'Figure saved: {figname}.pdf /.png')




# ran = list(range(0,6))
# ran = [0]
# from helpers import composite_flavor_dict
# flavors = ['b', 'c', 's', 'ud', 'g', 'unmatched']
# ran.append(21)
# for flav in flavors:
#     fig, ax = plt.subplots()
#     n_beforecut = hist.new.Var(range(8)).Double()
#     n_aftercut = hist.new.Var(range(8)).Double()
#     n_LHE_jet = hist.new.Var(range(8)).Double()

#     if flav in composite_flavor_dict:
#         flavs = composite_flavor_dict[flav]
#     else:
#         flavs = [flav]
    
#     for flavii in flavs:
#         flav_nr = flavor2partonNr[flavii]
#         n_beforecut.fill(ak.sum(np.abs(events.Jet.partonFlavour)==flav_nr,axis=1))
#         n_aftercut.fill(ak.sum(np.abs(reco_jets.partonFlavour)==flav_nr,axis=1))
#         n_LHE_jet.fill(ak.sum(np.abs(reco_jets_LHE["LHE_flavour2"])==flav_nr,axis=1))

#     n_beforecut = n_beforecut/all_ev
#     n_aftercut = n_aftercut/all_ev
        
#     n_beforecut.plot1d(histtype='fill', alpha=0.75, label='all jets')
#     n_aftercut.plot1d(histtype='fill', alpha=0.75, label='selected jets')

#     # if 21 == flav:
#         # n_LHE_glu = hist.new.Var(range(8)).Double()
#         # LHEPart = events.LHEPart
#         # absLHEid = np.abs(LHEPart.pdgId)
#         # LHE_outgoing = LHEPart[(LHEPart.status==1) & ((absLHEid < 6) | (absLHEid == 21))]
#         # n_LHE_glu.fill(ak.sum(LHE_outgoing.pdgId == 21, axis=1))
#         # n_LHE_glu = n_LHE_glu/all_ev
#         # n_LHE_glu.plot1d(histtype='fill', alpha=0.8, label='LHE partons')

#     n_LHE_jet = n_LHE_jet/all_ev
#     n_LHE_jet.plot1d(histtype='fill', alpha=0.75, label='LHE jets')

#     hep.cms.label("Private work", loc=0, data=False, ax=ax, rlabel='')
#     hep.label.exp_text(text=f'{flav} jets', loc=2, ax=ax)
#     ax.set_xticks(range(8))
#     ax.set_xlabel(f'number of {flav} jets')
#     ax.set_ylabel('Events')
#     ax.legend()

#     figname = dir_name+f'/n_{legend_str_to_filename(flav)}_jets_after_cut_ttbar'
#     plt.savefig(figname+'.pdf')
#     plt.savefig(figname+'.png')
#     print(f'Figure saved: {figname}.pdf /.png')