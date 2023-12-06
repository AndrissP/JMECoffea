import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import hist
import os
from plotters.pltStyle import pltStyle
pltStyle('hep')

dir_name = 'plotters/fig/misc'
if not os.path.exists(dir_name):
    os.mkdir(dir_name)

def plot_dphi(dphi, dphi_cut, dataset):
    h_dphi = hist.new.Reg(16, -np.pi, np.pi).Double()
    # breakpoint()
    h_dphi.fill(dphi)
    fig, ax = plt.subplots()
    h_dphi.plot1d(histtype='fill', alpha=0.75) #, label='all jets')
    hep.cms.label("Private work", loc=0, data=False, ax=ax, rlabel='')
    # hep.label.exp_text(text=f'$\Delta\phi$ between leading jet and Z', loc=2, ax=ax)
    ax.set_xlabel(f'$\Delta\phi$')
    ax.set_ylabel('Events')
    ax.vlines(dphi_cut, 0, 1.1*h_dphi.values().max(), color='b', linestyle='--')
    ax.vlines(-dphi_cut, 0, 1.1*h_dphi.values().max(), color='b', linestyle='--')
    # ax.legend()
    if 'DY' in dataset:
        figname = dir_name+f'/dphi_leading_jet_DY'
    elif 'QCD' in dataset:
        figname = dir_name+f'/dphi_leading_jet_QCD'
    plt.savefig(figname+'.pdf')
    plt.savefig(figname+'.png')
    print(f'Figure saved: {figname}.pdf /.png')