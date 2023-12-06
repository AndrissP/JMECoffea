#!/usr/bin/env python
# JERCProcessorcuts.py
"""
Defines the cuts that can be used in CoffeaJERCProcessor_L5 

Author(s): Andris Potrebko (RTU)
"""

import awkward as ak
import numpy as np
from memory_profiler import profile
from LHE_flavour import get_LHE_flavour, get_LHE_flavour2

### input numbers
mZpdg = 91.1876

def inv_mass_plus(lepton_pairs):
    ''' Compute the invariant mass of the leptons in `lepton pairs`
    '''
    return np.sqrt(np.abs(np.sum(lepton_pairs.E, axis=1)**2 - np.sum(lepton_pairs.px, axis=1)**2
                    - np.sum(lepton_pairs.py, axis=1)**2 - np.sum(lepton_pairs.pz, axis=1)**2))

def jet_iso_cut(reco_jets, dr_cut=0.8):
    ''' Jet isolation cut'''
    drs, _ = reco_jets.metric_table(reco_jets, return_combinations=True, axis=1)
    jet_iso_mask = ~ ak.any((1e-10<drs) & (drs<dr_cut), axis=2 )
    return reco_jets[jet_iso_mask]
# @profile
def leading_jet_and_alpha_cut(reco_jets, leptons, events, dataset, alphaQCD, alphaDY, NjetsQCD, NjetsDY):
    '''
    Selects the leading generator jets and performs the alpha cut
    Alpha cut = cut on the additional jet activity: to avoid effects due to a non-physical jet spectrum in the MC
    Output reco_jets are sorted according to the matched_gen pt.
    '''
    if "QCD" in dataset:
        if NjetsQCD!=-1:
            leading_gen_jet_indices = ak.argsort(reco_jets.matched_gen.pt, axis=1, ascending=False)[:,:NjetsQCD]
        else:
            leading_gen_jet_indices = ak.argsort(reco_jets.matched_gen.pt, axis=1, ascending=False)[:,:]
        reco_jets = reco_jets[leading_gen_jet_indices]

        if NjetsQCD>2 or NjetsQCD==-1:
            # To correctly/safely treat the cases where there are less then 3(2) jets in QCD (DY) left after the cuts
            # pad nons to the correct size of the jets and then accept all the events where there are less than 3(2) jets assuming that all the bad jets where already cut out
            gen_jetpt = ak.pad_none(reco_jets.matched_gen.pt, 3, axis=1, clip=False)
            alpha = gen_jetpt[:,2]*2/(gen_jetpt[:,0]+gen_jetpt[:,1])
            alpha = ak.fill_none(alpha,0)
            reco_jets = reco_jets[alpha<alphaQCD] #[:,:NjetsQCD]
            events = events[alpha<alphaQCD]
            leptons = leptons[alpha<alphaQCD]
        
    elif 'DY' in dataset:
        if NjetsDY!=-1:
            leading_gen_jet_indices = ak.argsort(reco_jets.matched_gen.pt, axis=1, ascending=False)[:,:NjetsDY]
        else:
            leading_gen_jet_indices = ak.argsort(reco_jets.matched_gen.pt, axis=1, ascending=False)[:,:]
        reco_jets = reco_jets[leading_gen_jet_indices]

        if NjetsDY>1 or NjetsDY==-1:
            gen_jetpt = ak.pad_none(reco_jets.matched_gen.pt, 2, axis=1, clip=False)
            alpha = gen_jetpt[:,1]/ak.sum(leptons.pt,axis=1)
            alpha = ak.fill_none(alpha,0)
            reco_jets = reco_jets[alpha<alphaDY] #[:,:NjetsDY]
            events = events[alpha<alphaDY]
            leptons = leptons[alpha<alphaQCD]

    return reco_jets, events, leptons

def dphi_separation(reco_jets, events, leptons, dataset, dphi_DY, dphi_QCD, debug=False):
    ''' Select back-to-back jets
    '''
    non_triv_jets = ak.num(reco_jets)>1
    reco_jets = reco_jets[non_triv_jets]
    events = events[non_triv_jets]
    leptons = leptons[non_triv_jets]

    # print("num(jets) = ", ak.num(reco_jets))
    if "QCD" in dataset:
        delta_phi = reco_jets[:,0].delta_phi(reco_jets[:,1])
        if debug:
            from plotters.plot_dphi import plot_dphi
            plot_dphi(delta_phi, dphi_QCD, dataset)
        dphi_pass = np.abs(delta_phi)>dphi_QCD

    elif "DY" in dataset:        
        # Calculate the phi of the di-lepton system
        Z_px = leptons[:, 0].px + leptons[:, 1].px
        Z_py = leptons[:, 0].py + leptons[:, 1].py
        lepton_phi = np.arctan2(Z_py, Z_px)

        delta_phi = reco_jets[:, 0].phi - lepton_phi
        delta_phi = (delta_phi + np.pi) % (2 * np.pi) - np.pi  # ensure the result is between -pi and pi
        if debug:
            from plotters.plot_dphi import plot_dphi
            plot_dphi(delta_phi, dphi_DY, dataset)
        dphi_pass = np.abs(delta_phi)>dphi_DY
    else:
        dphi_pass = np.array([True]*len(events))

    reco_jets = reco_jets[dphi_pass]
    events = events[dphi_pass]
    leptons = leptons[dphi_pass]
    return reco_jets, events, leptons


def jet_pt_cut(reco_jets, mingenjetpt, apply=False):
    if not apply:
        return reco_jets
    jet_pt_mask = reco_jets.matched_gen.pt>mingenjetpt
    ## funny workaround to change the ak.type of jet_pt_mask from '10 * var * ?bool' to '10 * var * bool'
    ## otherwise after the correction .matched_gen field is not found.
    jet_pt_mask_shape = ak.num(jet_pt_mask)
    jet_pt_mask_np = ak.flatten(jet_pt_mask).to_numpy()
    jet_pt_mask = ak.unflatten(jet_pt_mask_np.data, jet_pt_mask_shape)
    reco_jets = reco_jets[jet_pt_mask]
    return reco_jets

def good_lepton_cut(reco_jets, events, dataset, leptons, tightelectrons, tightmuons):
    '''Comparing to number of generated prompt leptons can deal with all the sample cases (2 for DY, 0,1,2 for TTBAR, 0 for QCD)
    Cuts on DY and ttbar based on L3Res selections https://twiki.cern.ch/twiki/bin/view/CMS/L3ResZJet
    '''
    events_with_good_lep = ((ak.num(tightmuons) == ak.num(leptons))
                    | (ak.num(tightelectrons) == ak.num(leptons) )
                    )        

    DYcond = np.array([True]*len(events))

    if 'DY' in dataset:
        DYcond = DYcond * (
            (np.sum(tightelectrons.pt, axis=1)>15) | (np.sum(tightmuons.pt, axis=1)>15)
        )
        DYcond = DYcond * (
            (np.abs(inv_mass_plus(tightelectrons) - mZpdg) < 20)
            | (np.abs(inv_mass_plus(tightmuons) - mZpdg) < 20)
        )

    events = events[events_with_good_lep*DYcond]
    reco_jets = reco_jets[events_with_good_lep*DYcond]
    leptons = leptons[events_with_good_lep*DYcond]
    tightelectrons = tightelectrons[events_with_good_lep*DYcond]
    tightmuons = tightmuons[events_with_good_lep*DYcond]
    return(events, reco_jets, leptons, tightelectrons, tightmuons)

def select_leptons(selectedEvents):
    ''' Select reconstructed electrons and muons according to the L3Res selections https://twiki.cern.ch/twiki/bin/view/CMS/L3ResZJet
    select generated leptons
    '''
    muon = selectedEvents.Muon
    tight_mu_cut = (muon.tightId) & (muon.pfIsoId>=4) & (np.abs(muon.eta)<2.3) & (muon.pt>20)
    tightmuons = muon[tight_mu_cut]

    electron = selectedEvents.Electron
    tight_ele_cut = (electron.cutBased==4) &(np.abs(electron.eta)<2.4) & (electron.pt>25)
    tightelectrons = electron[tight_ele_cut]
    
    genpart = selectedEvents.GenPart
    lepton_mask = (
            ((np.abs(genpart.pdgId) == 11) | (np.abs(genpart.pdgId) == 13) | (np.abs(genpart.pdgId) == 15 ))
            & (genpart.statusFlags>>13&1 == 1) 
            & (genpart.statusFlags&1 == 1)
    )
    leptons = genpart[lepton_mask]
    return leptons, tightelectrons, tightmuons

def recolep_drcut(reco_jets, tightelectrons, tightmuons, dR=0.2):
    ''' Additional dR cut on the jets to not overlap with leptons
    (tight lepton veto id does not seem to cut all the leptons)
    '''
    drs = reco_jets.metric_table(tightelectrons, return_combinations=False, axis=1 )
    overlappng_reco_lep_mask = np.all((drs>dR),axis=2)

    drs = reco_jets.metric_table(tightmuons, return_combinations=False, axis=1 )
    overlappng_reco_lep_mask = overlappng_reco_lep_mask*np.all((drs>dR),axis=2)
    reco_jets = reco_jets[overlappng_reco_lep_mask]
    return reco_jets

def wrong_recolep_drcut(reco_jets, cut_prompt_lep):
    ''' Additional dR cut on not overlapping with leptons
    cut_prompt_lep = True: cut on only prompt leptons, else cut on all leptons (including from semileptonic jet decays -> wrong)
    '''
    if cut_prompt_lep:
        ele_partFlav = reco_jets.matched_electrons.genPartFlav
        mu_partFlav = reco_jets.matched_muons.genPartFlav
        dressed_electron_mask = np.logical_not(np.sum((ele_partFlav == 1) | (ele_partFlav == 15),axis=2))
        dressed_muon_mask = np.logical_not(np.sum((mu_partFlav == 1) | (mu_partFlav == 15),axis=2))
    else:
        dressed_electron_mask = ak.sum(ak.is_none(reco_jets.matched_electrons,axis=2), axis=2)==2
        dressed_muon_mask     = ak.sum(ak.is_none(reco_jets.matched_muons,axis=2), axis=2)==2

    jet_mask = dressed_electron_mask & dressed_muon_mask
    return reco_jets[jet_mask]

def select_Nth_jet(reco_jets, selectedEvents, N):
    ''' Select the Nth jet.
    For example, when N=3: take the jet that comes from ME in MG+Pythia8 and from the shower in Pythia8
    '''
    N_jets_exist = ak.num(reco_jets)>=N
    reco_jets = reco_jets[N_jets_exist]
    selectedEvents = selectedEvents[N_jets_exist]
    reco_jets = reco_jets[:,N-1:N]
    return reco_jets, selectedEvents

def jetMCmatching(reco_jets, dR=0.2, apply=False):
    # MC jet matching
    if not apply:
        return reco_jets
    matchedJets = reco_jets[reco_jets.delta_r(reco_jets.matched_gen) < dR]
    # matchedJets = ak.cartesian([reco_jets.matched_gen, reco_jets])
    # deltaR = matchedJets.slot0.delta_r(matchedJets.slot1)
    # matchedJets = matchedJets[deltaR < dR]
    # matchedJets = matchedJets[ak.num(matchedJets) > 0]
    return matchedJets

def remove_apply(cut_tmp):
    return {key: cut_tmp[key] for key in cut_tmp if key!='apply'}

def correct_jets(selected_jets, selectedEvents, processor):
    ############ Apply Jet energy corrections on the jets ###########
    # define variables needed for corrected jets
    # https://coffeateam.github.io/coffea/notebooks/applying_corrections.html#Applying-energy-scale-transformations-to-Jets
    ## raw - subtracting back the corrections applying when generating the NanoAOD
    selected_jets['pt_raw'] = (1 - selected_jets['rawFactor']) * selected_jets['pt']     #raw pt. pt before the corrections are applied to data
    selected_jets['mass_raw'] = (1 - selected_jets['rawFactor']) * selected_jets['mass']
    selected_jets['pt_gen'] = ak.values_astype(ak.fill_none(selected_jets.matched_gen.pt, 0), np.float32)
    selected_jets['rho'] = ak.broadcast_arrays(selectedEvents.fixedGridRhoFastjetAll, selected_jets.pt)[0]
    events_cache = selectedEvents.caches[0]

    reco_jets = processor.jet_factory.build(selected_jets, lazy_cache=events_cache)
    return reco_jets

def apply_jetNevent_cuts(events, cfg, cutflow_evts, cutflow_jets, processor, dataset):
    
    ############ Event Cuts ############
    # apply npv cuts
    cutflow_evts.fill(cutflow='all_events', weight=len(events))
    
    npvCut = (events.PV.npvsGood > 0)
    pvzCut = (np.abs(events.PV.z) < 24)
    rxyCut = (np.sqrt(events.PV.x*events.PV.x + events.PV.y*events.PV.y) < 2)
    eventsCut = npvCut & pvzCut & rxyCut
    if cfg["gen_vtx_dz_cut"]["apply"]==True:
        gen_vtx_dz_cut = (np.abs(events.GenVtx.z-events.PV.z)<0.2)
        eventsCut = eventsCut & gen_vtx_dz_cut
    
    selectedEvents = events[eventsCut] 
    cutflow_evts.fill(cutflow='gen vertex cut', weight=len(selectedEvents))
    # get GenJets and Jets
    jets = selectedEvents.Jet
    cutflow_jets.fill(cutflow='all_jets', weight=ak.sum(ak.num(jets)))

    ########### Redo the flavour tagging if neccesarry. LHE Flavour2 derivation has to be done before the jet cuts  ###########
    #### Some samples have a missing LHE flavour infomration ####
    if (not 'LHEPart' in selectedEvents.fields) and ('LHE_flavour' in processor.jetflavour):
        raise ValueError(f"jet flavour is chosen as {processor.jetflavour}, but the sample does not contain 'LHEPart' "+
                                ", so the jet flavour cannot be recalculated.")
            
    if processor.jetflavour=='LHE_flavour2' or cfg["split_ISR_FSR_gluons"]:
    # if True: ### for testing the ISR/ FSR splitting
        jets = get_LHE_flavour2(jets, selectedEvents)

    ############ Jet selection ###########
    # Require that at least one gen jet is matched
    # jet_gen_match_mask = ~ak.is_none(jets.matched_gen,axis=1)
    selected_jets = jets[~ak.is_none(jets.matched_gen,axis=1)]
    # del jet_gen_match_mask
    cutflow_jets.fill(cutflow='gen matched', weight=ak.sum(ak.num(selected_jets)))

    reco_jets = correct_jets(selected_jets, selectedEvents, processor)
    # ############ Apply Jet energy corrections on the jets ###########
    # # define variables needed for corrected jets
    # # https://coffeateam.github.io/coffea/notebooks/applying_corrections.html#Applying-energy-scale-transformations-to-Jets
    # ## raw - subtracting back the corrections applying when generating the NanoAOD
    # selected_jets['pt_raw'] = (1 - selected_jets['rawFactor']) * selected_jets['pt']     #raw pt. pt before the corrections are applied to data
    # selected_jets['mass_raw'] = (1 - selected_jets['rawFactor']) * selected_jets['mass']
    # selected_jets['pt_gen'] = ak.values_astype(ak.fill_none(selected_jets.matched_gen.pt, 0), np.float32)
    # selected_jets['rho'] = ak.broadcast_arrays(selectedEvents.fixedGridRhoFastjetAll, selected_jets.pt)[0]
    # events_cache = selectedEvents.caches[0]

    # reco_jets = processor.jet_factory.build(selected_jets, lazy_cache=events_cache)

    leptons, tightelectrons, tightmuons = select_leptons(selectedEvents)
    # ###### Event selection based on leptons: 2 (0/1/2) reco leptons for DY (TTBAR semilep) ######
    # # cuts on DY and ttbar based on L3Res selections https://twiki.cern.ch/twiki/bin/view/CMS/L3ResZJet
    if cfg["good_lepton_cut"]["apply"]==True:
        selectedEvents, reco_jets, leptons, tightelectrons, tightmuons = good_lepton_cut(reco_jets, selectedEvents, dataset, leptons, tightelectrons, tightmuons)
    cutflow_evts.fill(cutflow='lepton selection', weight=len(selectedEvents))

    # Require tight lepton veto id on jets = no matched (dressed) leptons in the jet;
    # Leptons are also reconstructed as jets with just one (or more) particle, so it is important to remove them
    if cfg["tight_lepton_veto_id"]["apply"]==True:
        reco_jets[(reco_jets.jetId >> 2 & 1)==1] ### tight lepton veto id
    cutflow_jets.fill(cutflow='tight lep. id', weight=ak.sum(ak.num(reco_jets)))

    if cfg["recolep_drcut"]["apply"]==True:
        reco_jets = recolep_drcut(reco_jets, leptons, leptons)
    cutflow_jets.fill(cutflow=r'$\Delta R$'+' cut with leptons', weight=ak.sum(ak.num(reco_jets)))

    cut_tmp = cfg["jet_pt_cut"]
    if cut_tmp["apply"]==True:
        reco_jets = jet_pt_cut(reco_jets, cut_tmp["mingenjetpt"])
    cutflow_jets.fill(cutflow='jetpt cut', weight=ak.sum(ak.num(reco_jets)))

    # redo the jet matching with potentially lower dr cut than matched automatically
    cut_tmp = cfg["reco_jetMCmatching"]
    if cut_tmp["apply"]==True:
        reco_jets = jetMCmatching(reco_jets, **cut_tmp)
    cutflow_jets.fill(cutflow='matched gen cut', weight=ak.sum(ak.num(reco_jets)))
    ######### Alpha cut = cut on the additional jet activity  ############    
    # Not used since run 2 because the large pileup causes a bias    
    cut_tmp = cfg["leading_jet_and_alpha_cut"]
    if cut_tmp["apply"]==True:
        reco_jets, selectedEvents, leptons = leading_jet_and_alpha_cut(reco_jets, leptons, selectedEvents, dataset, **remove_apply(cut_tmp))

    cut_tmp = cfg["select_Nth_jet"]
    if cut_tmp["apply"]==True:
        reco_jets, selectedEvents = select_Nth_jet(reco_jets, selectedEvents, cut_tmp["N"])
    cutflow_jets.fill(cutflow=r'$\alpha$ cut'+'\nleading jets', weight=ak.sum(ak.num(reco_jets)))
    cutflow_evts.fill(cutflow=r'$\alpha$ cut',       weight=len(selectedEvents))

    cut_tmp = cfg["dphi_separation"]
    if cut_tmp["apply"]==True:
        reco_jets, selectedEvents, leptons = dphi_separation(reco_jets, selectedEvents, leptons, dataset, dphi_DY=cut_tmp["dphi_DY"], dphi_QCD=cut_tmp["dphi_QCD"], debug=False)
    cutflow_jets.fill(cutflow=r'$\Delta \phi$ cut', weight=ak.sum(ak.num(reco_jets)))
    cutflow_evts.fill(cutflow=r'$\Delta \phi$ cut',       weight=len(selectedEvents))

    # # Cut on overlapping jets
    cut_tmp = cfg["jet_iso_cut"]
    if cut_tmp["apply"]==True:
        reco_jets = jet_iso_cut(reco_jets, **remove_apply(cut_tmp))
    cutflow_jets.fill(cutflow='iso jets', weight=ak.sum(ak.num(reco_jets)))
    return selectedEvents, reco_jets, cutflow_evts, cutflow_jets