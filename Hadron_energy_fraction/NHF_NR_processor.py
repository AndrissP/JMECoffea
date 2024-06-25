'''
A little messy processo to get the 
'''

from coffea import processor
import hist
import awkward as ak
import numpy as np
from hist import Hist

from common_binning import JERC_Constants
ptbins = np.array(JERC_Constants.ptBinsEdgesMCTruth())

nums = [-3334, -3322, -3312, -3222, -3122, -3112, -2212,
              -2112, -321, -211, -13, -11, 11, 13, 22, 130, 211,
              310, 321, 2112, 2212, 3112,
              3122, 3222, 3312, 3322, 3334]


from coffea import processor
class Processor(processor.ProcessorABC):
    def __init__(self, processor_config):
        pass
    
    @property
    def accumulator(self):
        return self._accumulator
    
    def process(self, events):
        particle_ax = hist.axis.IntCategory(nums, name="had_type")
        neutHF_ax  = hist.axis.Regular(20, 0, 1, name="neutHF")
        ptresponse_axis = hist.axis.Regular( 100, 0, 2.5, overflow=True, underflow=True, name="ptresponse", label="RECO / GEN response")

        pt_gen_axis = hist.axis.Variable(ptbins, name="pt_gen", overflow=True, underflow=True, label=r"$p_{T,gen}$")

        # h_bad = Hist(particle_ax, flavour_ax, neutHF_ax, storage="weight")
        # h_good = Hist(particle_ax, flavour_ax, neutHF_ax, storage="weight")
        
        h_bad_NHF = Hist(particle_ax, pt_gen_axis, storage="weight")
        h_good_NHF = Hist(particle_ax, pt_gen_axis, storage="weight")

        h_bad_NHF_reco = Hist(neutHF_ax, pt_gen_axis, storage="weight")
        h_good_NHF_reco = Hist(neutHF_ax, pt_gen_axis, storage="weight")
        
        gen_jets = events.GenJet
        jets = events.Jet

        jet_gen_match_mask = ~ak.is_none(jets.matched_gen,axis=1)
        jets = jets[jet_gen_match_mask]

        ### genjets don't have the parton flavor stored, so we copy jets.partonFlavor to gen_gets.partonFlavour2.
        ### The index from jets to gen_jets is available but here also the index from gen_jets to jets is needed.
        genjet_shape = ak.num(gen_jets.partonFlavour)
        partonFlavour_flat = ak.flatten(ak.zeros_like(gen_jets.partonFlavour)).to_numpy()
        toJetIdx_flat = ak.flatten(ak.ones_like(gen_jets.partonFlavour)*-1).to_numpy()

        genjet_shapecum = np.concatenate([[0], np.cumsum(genjet_shape)[:-1]])
        genjetidx = ak.flatten(genjet_shapecum+jets.genJetIdx).to_numpy()
        partonFlavour_flat[genjetidx] = ak.flatten(jets.partonFlavour).to_numpy()
        toJetIdx_flat[genjetidx] = ak.flatten(np.argsort(jets.genJetIdx) ).to_numpy()
        gen_jets["partonFlavour2"] = ak.unflatten(partonFlavour_flat, genjet_shape)
        gen_jets["jetIdx"] = ak.unflatten(toJetIdx_flat, genjet_shape)

        ### select the gen particles
        if 'genCandsIdx' in events.GenJetCands.fields:
            genCand = events.GenCands[events.GenJetCands.genCandsIdx]
        else:
            genCand = events.GenCands[events.GenJetCands.pFCandsIdx]
        
        genCand["jetIdx"] = events.GenJetCands.jetIdx

        ### select the gen particles
        gen_jets_vCand = gen_jets[events.GenJetCands.jetIdx]
        has_mathed_reco = gen_jets_vCand["jetIdx"]>-1
        gen_jets_vCand = gen_jets_vCand[has_mathed_reco]
        gen_jets_vCand = gen_jets_vCand[np.abs(gen_jets_vCand.partonFlavour2)==5]

        jets_cand = jets[gen_jets_vCand.jetIdx]
        jetpt = jets_cand.pt

        # ptmask = (gen_jets_vCand.pt>400) & (gen_jets_vCand.pt<1000)
        ### if masking jets with a good response and bad response
        jet_response = jetpt/gen_jets_vCand.pt
        badjet_mask = jet_response>1.1
        goodjet_mask = (jet_response<1.1) & (jet_response>0.9)

        total_good_mask = goodjet_mask
        total_bad_mask = badjet_mask
        gen_jets_vCand_good = gen_jets_vCand[total_good_mask]
        genCand_good = genCand[total_good_mask]
        jets_good = jets_cand[total_good_mask]
        
        gen_jets_vCand_bad = gen_jets_vCand[total_bad_mask]
        genCand_bad = genCand[total_bad_mask]
        jets_bad = jets_cand[total_bad_mask]
#         jetpt_bad = jetpt[total_mask]
#         resp_hist.fill(ak.flatten(jetpt/gen_jets_vCand.pt).to_numpy())

        h_good_NHF.fill(ak.flatten(genCand_good.pdgId), ak.flatten(gen_jets_vCand_good.pt), weight=ak.flatten(genCand_good.pt))
        h_bad_NHF.fill(ak.flatten(genCand_bad.pdgId), ak.flatten(gen_jets_vCand_bad.pt), weight=ak.flatten(genCand_bad.pt))

        jets = jets[np.abs(jets.partonFlavour)==5]
        # ptmask = (jets.matched_gen.pt>400) & (jets.matched_gen.pt<1000)
        ### if masking jets with a good response and bad response
        jet_response = jets.pt/jets.matched_gen.pt
        # print("jet_response", jet_response)
        badjet_mask = jet_response>1.1
        # print("badjet_mask", badjet_mask)
        goodjet_mask = (jet_response<1.1) & (jet_response>0.9)
        # print("goodjet_mask", goodjet_mask)
        total_good_mask = goodjet_mask
        total_bad_mask = badjet_mask

        jets_good = jets[total_good_mask]
        jets_bad = jets[total_bad_mask]

        h_good_NHF_reco.fill(ak.flatten(jets_good["neHEF"]+jets_good["neEmEF"]), ak.flatten(jets_good.matched_gen.pt))
        h_bad_NHF_reco.fill(ak.flatten(jets_bad["neHEF"]+jets_bad["neEmEF"]), ak.flatten(jets_bad.matched_gen.pt))

        # h_good_NHF_reco.fill(ak.flatten(jets_good["neHEF"]+jets_good["neEmEF"]), gen_jets_vCand_good.pt)
        # h_bad_NHF_reco.fill(ak.flatten(jets_bad["neHEF"]+jets_bad["neEmEF"]), gen_jets_vCand_bad.pt)

        return {"h_good":h_good_NHF, "h_bad":h_bad_NHF, "h_good_NHF_reco":h_good_NHF_reco, "h_bad_NHF_reco":h_bad_NHF_reco}
    
    def postprocess(self, accumulator):
        return accumulator