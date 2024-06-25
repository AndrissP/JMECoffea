#!/usr/bin/env python
    # coding: utf-8
### Hadron_energy_fractions.py
### File automatically converted using ConvertJupyterToPy.ipynb from Hadron_energy_fractions.ipynb
### No comments or formatting is preserved by the transfer.
def main():
    
    import awkward as ak
    import numpy as np
    import matplotlib.pyplot as plt
    
    import sys
    top_path = '../'
    if top_path not in sys.path:
        sys.path.append(top_path)
    # coffea_path = '/afs/cern.ch/user/a/anpotreb/top/JERC/coffea/'
    # if coffea_path not in sys.path:
    #     sys.path.insert(0,coffea_path)
    
    # ak_path = '/afs/cern.ch/user/a/anpotreb/top/JERC/local-packages/'
    
    # if ak_path not in sys.path:
    #     sys.path.insert(0,ak_path)
    
    import hist
    from hist import Hist
    
    import pickle
    
    from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
    
    from plotters.pltStyle import pltStyle
    pltStyle(style='hep')
    plt.rcParams['figure.dpi'] = 150
    
    # events = NanoEventsFactory.from_root(
    # #     'file:///afs/cern.ch/work/m/mseidel/public/forAndris/testNano_PFNANO.root',
    # #     'file:///afs/cern.ch/work/m/mseidel/public/forAndris/testNanoCustomise_herwig7_NANO.root',
    #     'file:///eos/cms/store/group/cmst3/group/top/TOP21009/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/mc_2016_ULPreVFP/230330_162057/0000/nano_mc2016pre_27.root',
    #     schemaclass=NanoAODSchema.v6,
    # #     entry_stop=5000,
    #     metadata={
    # #             "dataset": "TTBAR",
    # #               "xsec":semilepxsec,
    #              },
    # ).events()
    
    import os
    figdir = 'fig/hadrons/'
    if not os.path.exists(figdir):
        os.mkdir(figdir)
    
    names = ["Omega+", "a-Theta0", "Theta+", "a-Sigma+=-", "a-Lambda=0", "a-Sigma-=+", "p-",  "a-n",
    "KS-","pi-", "mu-" ,"el-", "el", "mu", "gamma", "K0L", "pi+", "KS0", "K+",
    "n", "p", "Sigma-", "Lambda", "Sigma+", "Theta-", "Theta0", "Omega-" ]
    
    nums = [-3334, -3322, -3312, -3222, -3122, -3112, -2212,
                  -2112, -321, -211, -13, -11, 11, 13, 22, 130, 211,
                  310, 321, 2112, 2212, 3112,
                  3122, 3222, 3312, 3322, 3334]
    # for a, b in zip(names, nums):
    #     print(a+': ', b)
    
    particle_ax = hist.axis.IntCategory(nums, name="had_type")
    
    '''
    Omega+, a-Theta0, Theta+, a-Sigma+=-, a-Lambda=0, a-Sigma-=+, p-,  a-n,
    KS-,pi-, mu- ,el-, el, mu, gamma, K0L, pi+, KS0, K+,
    n, p, Sigma-, Lambda, Sigma+, Theta-, Theta0, Omega- 
    '''
    
    jetpt_ax = hist.axis.Regular(30, 1, 500, overflow=False, underflow=False, name="jet_pt")
    had_pt_ax = hist.axis.Regular(10, 1, 20, overflow=False, underflow=False, name="had_pt")
    # Hist.new.Reg(100, 5000, 100, name="B").Double()
    h = Hist(particle_ax, jetpt_ax) #, had_pt_ax)
    
    # ### split in jet flavour, inclusive in pt
    
    ### define different hadron categories
    pos_mesons = [211, 321 ]
    neg_mesons = [-211, -321 ]
    neut_mesons = [130, 310]
    leptons = [-13, -11, 11, 13]
    protons = [2212]
    aprotons = [-2212]
    neutrons = [2112]
    aneutrons = [-2112]
    sPbaryons = [3222, -3112, -3312, -3334]
    sNbaryons = [-3222, 3112, 3312, 3334]
    sbaryons = [3122, 3322]
    asbaryons = [-3122, -3322]
    gamma = [22]
    
    ### define hadron categories
    categories = {"pos_mesons": [211, 321 ],
                    "neg_mesons": [-211, -321 ],
                    "neut_mesons": [130, 310],
                    "leptons": [-13, -11, 11, 13],
                    "protons": [2212],
                  "a-protons": [-2212],
                  "neutrons": [2112],
                  "a-neutrons": [-2112],
                  "sPbaryons": [3222, -3112, -3312, -3334],
                  "sNbaryons": [-3222, 3112, 3312, 3334],
                  "sbaryons": [3122, 3322],
                  "a-sbaryons": [-3122, -3322],
                  "gamma": [22]}
    
    categories = {"pos_mesons": [211, ],
                  "pos_s_mesons": [321, ],
                    "neg_mesons": [-211, ],
                  "neg_s_mesons": [-321, ],
                    "neut_mesons": [130, 310],
                    "leptons": [-13, -11, 11, 13],
                    "protons": [2212],
                  "a-protons": [-2212],
                  "neutrons": [2112],
                  "a-neutrons": [-2112],
                  "pos_s_baryons": [3222, -3112, -3312, -3334],
                  "neg_s_baryons": [-3222, 3112, 3312, 3334],
                  "neut_s_baryons": [3122, 3322],
                  "neut_s_anti_baryons": [-3122, -3322],
                  "gamma": [22]}
    
    ### check for consistency
    set(np.concatenate([ii for ii in categories.values()])) == set(nums)
    
    particle_ax = hist.axis.IntCategory(nums, name="had_type")
    flavour_ax = hist.axis.IntCategory([5, -5, 4, -4, 3, -3, 2, -2, 1, -1, 21, 0, -999], name="parton_flavour")
    flavour_ax = hist.axis.IntCategory([5, 4, 3, 2, 1, 21, 0, -999], name="parton_flavour")
    
    # jetpt_ax = hist.axis.Regular(30, 1, 500, overflow=False, underflow=False, name="jet_pt")
    # jeta_eta_ax = hist.axis.Variable(np.array([0, 1.3, 2.5, 3, 5]), overflow=False,
    #                                  underflow=False, name="jet_eta")
    had_pt_ax = hist.axis.Regular(10, 1, 20, overflow=False, underflow=False, name="had_pt")
    h = Hist(particle_ax, flavour_ax, storage="weight")
    
    # ### Reco gen matching for samples (e.g Herwig) where gen jets do not contain partonFlavor (all 0) but reco jets do
    # ### Then fill the histogram
    # ### Question: Why do I reco the gen matching instead of using events.jet.matchedGen?
    
    ### To check the response
    resp_hist = Hist.new.Reg(100, 0, 5, overflow=True, underflow=True, name="response").Double()
    
    from coffea import processor
    class Processor(processor.ProcessorABC):
        def __init__(self):
            pass
        
        @property
        def accumulator(self):
            return self._accumulator
        
        def process(self, events):
            particle_ax = hist.axis.IntCategory(nums, name="had_type")
            flavour_ax = hist.axis.IntCategory([5, -5, 4, -4, 3, -3, 2, -2, 1, -1, 21, 0, -999], name="parton_flavour")
            flavour_ax = hist.axis.IntCategory([5, 4, 3, 2, 1, 21, 0, -999], name="parton_flavour")
    
            # jetpt_ax = hist.axis.Regular(30, 1, 500, overflow=False, underflow=False, name="jet_pt")
            # jeta_eta_ax = hist.axis.Variable(np.array([0, 1.3, 2.5, 3, 5]), overflow=False,
            #                                  underflow=False, name="jet_eta")
            had_pt_ax = hist.axis.Regular(10, 1, 20, overflow=False, underflow=False, name="had_pt")
            h = Hist(particle_ax, flavour_ax, storage="weight")
            
            gen_jets = events.GenJet
            jets = events.Jet
    
            jet_gen_match_mask = ~ak.is_none(jets.matched_gen,axis=1)
            jets = jets[jet_gen_match_mask]
    
            # in the case of Herwig, the gen jets have empty partonFlavour, so it has to be rematched
            toJetIdx = ak.zeros_like(gen_jets.partonFlavour)
            genjet_shape = ak.num(gen_jets.partonFlavour)
            jet_shape = ak.num(jets.partonFlavour)
            partonFlavour_flat = ak.flatten(ak.zeros_like(gen_jets.partonFlavour)).to_numpy()
            toJetIdx_flat = ak.flatten(ak.ones_like(gen_jets.partonFlavour)*-1).to_numpy()
    
            genjet_shapecum = np.concatenate([[0], np.cumsum(genjet_shape)[:-1]])
            jet_shapecum = np.concatenate([[0], np.cumsum(jet_shape)[:-1]])
            genjetidx = ak.flatten(genjet_shapecum+jets.genJetIdx).to_numpy()
            partonFlavour_flat[genjetidx] = ak.flatten(jets.partonFlavour).to_numpy()
            toJetIdx_flat[genjetidx] = ak.flatten(np.argsort(jets.genJetIdx) ).to_numpy()
            gen_jets["partonFlavour2"] = ak.unflatten(partonFlavour_flat, genjet_shape)
            gen_jets["jetIdx"] = ak.unflatten(toJetIdx_flat, genjet_shape) # index of the matched reco jet
    
            if 'genCandsIdx' in events.GenJetCands.fields:
                genCand = events.GenCands[events.GenJetCands.genCandsIdx]
            else:
                genCand = events.GenCands[events.GenJetCands.pFCandsIdx]
            gen_jets_vCand = gen_jets[events.GenJetCands.jetIdx]
            jetpt = jets[gen_jets_vCand.jetIdx].pt
    
            has_mathed_reco = gen_jets_vCand["jetIdx"]>-1
            ptmask = (gen_jets_vCand.pt>400) & (gen_jets_vCand.pt<1000)
            ### if masking jets with a good response and bad response
            funnyjet_mask = jetpt>gen_jets_vCand.pt*1.1
            goodjet_mask = (jetpt<gen_jets_vCand.pt*1.1) & (jetpt>gen_jets_vCand.pt*0.9)
    
            total_mask = has_mathed_reco & ptmask & goodjet_mask # funnyjet_mask #
            gen_jets_vCand = gen_jets_vCand[total_mask]
            genCand = genCand[total_mask]
            jetpt = jetpt[total_mask]
    #         resp_hist.fill(ak.flatten(jetpt/gen_jets_vCand.pt).to_numpy())
    
            h.fill(ak.flatten(genCand.pdgId), ak.flatten(np.abs(gen_jets_vCand.partonFlavour2)), weight=ak.flatten(genCand.pt));
        
            return h
        
        def postprocess(self, accumulator):
            return accumulator
    
    xrootdstr = 'file:///eos/cms/store/group/cmst3/group/top/TOP21009/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/mc_2016_ULPreVFP/230330_162057/0000/'
    def txt2filesls(dataset_name):
        with open(dataset_name) as f:
            rootfiles = f.read().split()
            fileslist = [xrootdstr + file for file in rootfiles]
        return fileslist
    
    Nfiles = 100
    dataset = 'fileNames/PFNanoTTToSemiLepPowPy.txt'
    fileslist = txt2filesls(dataset)[:Nfiles]
    filesets = {'dataset1': {"files": fileslist, "metadata": {"xsec": 1}}}
    
    UsingDaskExecutor = True
    CoffeaCasaEnv = False
    CERNCondorCluster = True
    
    if(UsingDaskExecutor and not CoffeaCasaEnv):
        from dask.distributed import Client 
     # Dask set up for LPC only 
        if not CERNCondorCluster:
            client = Client()
            client.get_versions(check=True)
    #         client.nanny = False
    
        else:
            from dask_lxplus import CernCluster
            import socket
    
            cluster = CernCluster(
    # #             memory=config.run_options['mem_per_worker'],
    # #             disk=config.run_options.get('disk_per_worker', "20GB"),
    #             env_extra=env_extra,
                cores = 1,
                memory = '4000MB',
                disk = '2000MB',
                death_timeout = '60',
                lcg = True,
                nanny = False,
                container_runtime = 'none',
                log_directory = '/eos/user/a/anpotreb/condor/log',
                scheduler_options = {
                    'port': 8786,
                    'host': socket.gethostname(),
                },
                job_extra = {
                    'MY.JobFlavour': '"longlunch"',
    #                 'transfer_input_files': '/afs/cern.ch/user/a/anpotreb/top/JERC/JMECoffea/count_2d.py',
                },
            )
            cluster.adapt(minimum=2, maximum=200)
            cluster.scale(8)
            client = Client(cluster)
        
    #     client.upload_file('CoffeaJERCProcessor'+tag_Lx+'.py')
    #     client.upload_file('count_2d.py')
    
        client
    
    # socket.gethostname()
    
    import time
    from numpy.random import RandomState
    load_preexisting = False
    from coffea import processor, util
    
    outname = 'HadronEfractions.coffea'
    
    if UsingDaskExecutor:
        client.close()
        time.sleep(5)
        if CERNCondorCluster or CoffeaCasaEnv:
            cluster.close()
    
    # filesets
    
    tstart = time.time()
    
    outputs_unweighted = {}
    
    seed = 1234577890
    prng = RandomState(seed)
    chunksize = 10000
    maxchunks = None
    
    if not load_preexisting:
        if not UsingDaskExecutor:
            chosen_exec = 'futures'
            output = processor.run_uproot_job(filesets,
                                              treename='Events',
                                              processor_instance=Processor(),
                                              executor=processor.iterative_executor,
        #                                        executor=processor.futures_executor,
                                              executor_args={
                                                  'skipbadfiles':True,
                                                  'schema': NanoAODSchema, #BaseSchema
                                                  'workers': 2},
                                              chunksize=chunksize,
                                              maxchunks=maxchunks)
        else:
            chosen_exec = 'dask'
            output = processor.run_uproot_job(filesets,
                                              treename='Events',
                                              processor_instance=Processor(),
                                              executor=processor.dask_executor,
                                              executor_args={
                                                  'client': client,
                                                  'skipbadfiles':True,
                                                  'schema': NanoAODSchema, #BaseSchema
                                                  'xrootdtimeout': 60,
                                                  'retries': 2,
    #                                               'workers': 2
                                              },
                                              chunksize=chunksize,
                                              maxchunks=maxchunks)
    
        elapsed = time.time() - tstart
        print("Processor finished. Time elapsed: ", elapsed)
    #     outputs_unweighted[name] = output
        print("Saving the output histograms under: ", outname)
        util.save(output, outname)
    #     outputs_unweighted[name] = output
    else:
        output = util.load(outname)
        print("Loaded histograms from: ", outname)
    
    #### Attempt to prevent the error when the cluster closes. Doesn't always work.
    if UsingDaskExecutor:
        client.close()
        time.sleep(5)
        if CERNCondorCluster or CoffeaCasaEnv:
            cluster.close()
    
if __name__ == "__main__":
    main()