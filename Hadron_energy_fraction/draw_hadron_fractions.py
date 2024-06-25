''' Draw the hadron energy fractions in jets (jet hadron content plots) from the coffea output files processed in `Hadron_energy_fractions.ipynb` or using `Hadron_energy_fracions.py`.

'''

import numpy as np
import matplotlib.pyplot as plt
import hist
from collections import OrderedDict
import mplhep as hep

import sys
top_path = '../'
if top_path not in sys.path:
    sys.path.append(top_path)
    
from plotters.pltStyle import pltStyle
pltStyle(style='hep')
plt.rcParams['figure.dpi'] = 150

plt.figure(num=None, figsize=(2, 2), dpi=80)
plt.plot([1,2,3],[1,3,3])
import matplotlib.pyplot as plt
pltStyle('hep')

nums = [-3334, -3322, -3312, -3222, -3122, -3112, -2212,
              -2112, -321, -211, -13, -11, 11, 13, 22, 130, 211,
              310, 321, 2112, 2212, 3112,
              3122, 3222, 3312, 3322, 3334]

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

categories = {"$\pi^+$": [211, ],
              "$\pi^-$": [-211, ],
              "$K^+$": [321, ],
              "$K^-$": [-321, ],
              "$[\pi^0; K^0_S]$": [130, 310],
              "$p$": [2212],
              "$\overline{p}$": [-2212],
              "$n$": [2112],
              "$\overline{n}$": [-2112],
              "$l^\pm$": [-13, -11, 11, 13],
              "$\gamma$": [22],
              "$[\Sigma^+; \Sigma^+; \Theta^+; \Omega^+]$": [3222, -3112, -3312, -3334],
              "$[\Sigma^-; \Sigma^-; \Theta^-; \Omega^-]$": [-3222, 3112, 3312, 3334],
              "$[\Lambda; \Theta^0]$": [3122, 3322],
              "$[\overline{\Lambda}; \overline{\Theta}^0]$": [-3122, -3322],
             }

### check for consistency
set(np.concatenate([ii for ii in categories.values()])) == set(nums)

def plot_stack_all_hads(h, hist_name='tot', figdir='fig/'):

    # rc_old = plt.rcParams.copy()
    # rc_top_old = plt.rcParams['figure.subplot.top']
    # plt.rcParams['figure.subplot.top'] = 0.63
    # plt.rcParams['figure.figsize'] = [4.2, 4.3] 

    s = h[:,:].stack("had_type")
    sumVals = h[sum,:].values()
    sumVals[sumVals==0] = 1
    # s*norms
    for ii in range(len(s)):
    #     hist_scale = s[ii].sum()
        hist_scale= 1/sumVals
        s[ii] = s[ii]*hist_scale # if hist_scale!=0 else s[ii];

    fig, ax = plt.subplots()
    s.plot(stack=True, histtype="fill")
    plt.ylabel("Energy fraction")
    plt.legend(bbox_to_anchor=(0, 1, 1, 0), loc="lower left", mode="expand",ncol=3); #, labels = had_labels
    fig.savefig(figdir+hist_name+'.png')
    fig.savefig(figdir+hist_name+'.pdf')
    print("Figure saved = ", figdir+hist_name+'.png')
    # plt.rcParams = rc_old

def plot_stack_merged(h, hist_name='tot', figdir='fig/', sample='MG+Py8', antiflav=True):
    print(plt.rcParams['figure.subplot.top'])
    rc_old = plt.rcParams.copy()
    plt.rcParams['figure.subplot.top'] = 0.70
    plt.rcParams['figure.subplot.left'] = 0.12
    plt.rcParams['figure.figsize'] = [4.2, 4.3] 
    print(plt.rcParams['figure.subplot.top'])
    my_stack = {}
    cat_keys = categories.keys()
    for key in cat_keys:
    #     complex_cat = map(lambda a: complex(a)*1j, [3112, 3122])
        complex_cat = [complex(valii)*1j for valii in categories[key]]
        hist_cat = h[complex_cat,:].project('parton_flavour')
        my_stack[key] = hist_cat
        
    my_stack = OrderedDict(reversed(list(my_stack.items())))
    hist_stack = hist.Stack.from_dict(my_stack)

    # s = h[:,:].stack("had_type")
    s = hist_stack
    sumVals = h[sum,:].values()
    sumVals[sumVals==0] = 1
    # s*norms
    for ii in range(len(s)):
    #     hist_scale = s[ii].sum()
        hist_scale= 1/sumVals
        s[ii] = s[ii]*hist_scale # if hist_scale!=0 else s[ii];

    fig, ax = plt.subplots()
    s.plot(stack=True, histtype="fill")
    ax.set_xlabel("Parton flavour")
    plt.ylabel("Energy fraction")
    ax.set_ylim(0,1)
    plt.legend(bbox_to_anchor=(0, 1.043, 1, 0), loc="lower left", mode="expand",ncol=3); #, labels = had_labels
    if antiflav:
        ax.set_xticklabels(['$b$', '$\overline{b}$', '$c$', '$\overline{c}$', '$s$', '$\overline{s}$', '$u$', '$\overline{u}$', '$d$', '$\overline{d}$', '$g$', "un.\njet", "un.\n had."])
        ## move the bbar and dbar labels lower to make them centered 
        xticks = ax.get_xticks()
        xticklabels = ax.get_xticklabels()
        # Adjust the position of specific tick labels
        for i, label in enumerate(xticklabels):
        #     xticks[i]
            if not xticks[i] in [1.5, 9.5]:  # Specify the x-tick values you want to move
                label.set_y(-0.01)  # Adjust this value as needed        
        # Redraw the plot with adjusted tick labels
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_xlim(-0,13)
    else:
        ax.set_xticklabels(['$b$', '$c$', '$s$', '$u$', '$d$', '$g$', 0, -100])

    hep.cms.label("Private work", loc=0, data=False, ax=ax, rlabel='')
    hep.label.exp_text(text=sample, loc=2, ax=ax)
    fig.savefig(figdir+hist_name+'.png')
    fig.savefig(figdir+hist_name+'.pdf')
    print("Figure saved = ", figdir+hist_name+'.png')
    # breakpoint()
    plt.rcParams = rc_old
#     print("scaled sum = ", s[ii].sum())
    
from coffea import util
output = util.load('HadronEfractions_QCD_Herwig.coffea')

import os
figdir = 'fig/hadrons/'
if not os.path.exists(figdir):
    os.makedirs(figdir)

for key in ['h_bad', 'h_good']:
    print("drawing for key = ", key)
    h = output[key]
    plot_stack_all_hads(h, hist_name="all_Her_"+key, figdir='fig/hadrons/')
    print("Outside funcs: ", plt.rcParams['figure.subplot.top'])
    plot_stack_merged(h, hist_name="merged_Her_"+key, figdir='fig/hadrons/', antiflav=False)

output = util.load('HadronEfractions_QCD_Pythia.coffea')

for key in ['h_bad', 'h_good']:
    print("drawing for key = ", key)
    h = output[key]
    plot_stack_all_hads(h, hist_name="all_Py_"+key, figdir='fig/hadrons/')
    plot_stack_merged(h, hist_name="merged_Py_"+key, figdir='fig/hadrons/', antiflav=False)

import pickle 

with open("hadrons/Pythia.pkl", "rb") as f:
    h = pickle.load(f)
# Py = True
plot_stack_all_hads(h, hist_name="all_Py", figdir='fig/hadrons/')
plot_stack_merged(h, hist_name="merged_Py", figdir='fig/hadrons/', sample='MG+Py8')

with open("hadrons/Herwig.pkl", "rb") as f:
    h = pickle.load(f)
# Py = False
plot_stack_all_hads(h, hist_name="all_Her", figdir='fig/hadrons/')
plot_stack_merged(h, hist_name="merged_Her", figdir='fig/hadrons/', sample='MG+Her7')
