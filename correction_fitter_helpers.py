header_txt = '''# L5 flavor corrections for IC5 algorithm
 [*J] (flavor * from diJet mixture)
 [*T] (flavor * from ttbar sample)
 [*S] (flavor * from simultaneous fit of diJet, ttbar and Drell-Yan)
 [*M] (flavor * from a combination of Drell-Yan (low pt) and diJet )
 [g*] (gluonss )
 [b*] (b quarks )
 [c*] (c quarks )
 [q*] (uds quarks )
 [ud*] (ud quarks )
 [a*] (all quarks )
 parametrization: Pt<p0 or Pt>p1, p2,... : fitting constants
 etamin  etamax  #ofparameters  ptmin  ptmax    p0         p1        p2        p3        p4
'''

def save_correction_txt_file(txtfile_outname, fit_res_all_tags):
    '''
    Saves the corrections in the txt file with the name `txtfile_outname`.
    fit_res_all_tags: for each correction tag (e.g., T, J, S, M), a dictionary of corrections for each flavor
    '''
    with open(txtfile_outname, 'w') as file:
        file.write(header_txt+'\n')
        for tag in fit_res_all_tags:
            fit_res_all = fit_res_all_tags[tag]
            for key in fit_res_all.keys():
                fit_res = fit_res_all[key]
                file.write(f'[{key}]\n')
                file.write('{1 JetEta 1 JetPt ([0]+[1]*log10(x)+[2]*pow(log10(x),2)+[3]*pow(log10(x),3)+[4]*pow(log10(x),4)) Correction L5Flavor}\n')
                ### copy from the positive eta region into the negative
                fit_res = np.vstack([np.hstack([np.flip(fit_res[:,0:2]*-1), np.flip(fit_res[:,2:], 0)]), fit_res])
                for row in fit_res:
                    row[2] = row[2]+2  #+2 because of the pt lower/higher limits that are not accounted into the # parameters before
                    line2write = ('{:>11} '*5+'{:>13} '*(int(row[2])-2)).format(*row[:2], int(row[2]), *np.round(row[3:], 7))+'\n'
                    file.write(line2write);
    print("Saving the corrections with the name = ", txtfile_outname)

def save_correction_txt_file_Mikko(txtfile_outname, fit_res_all_tags):
    '''
    Saves the corrections in the txt file with the name `txtfile_outname`.
    fit_res_all_tags: for each correction tag (e.g., T and J), a dictionary of corrections for each flavor
    '''
    with open(txtfile_outname, 'w') as file:
        # str_poly='([0]+[1]*log10(x)+[2]*pow(log10(x),2)+[3]*pow(log10(x),3)+[4]*1/pow(log10(x),[5]))'
        str_poly='([0]+[1]*log10(x)+[2]*pow(log10(x),2)+[3]*pow(log10(x),3)+[4]*pow(log10(x),4))'
        file.write(header_txt+'\n')
        for tag in fit_res_all_tags:
            fit_res_all = fit_res_all_tags[tag]
            for key in fit_res_all.keys():
                fit_res = fit_res_all[key]
                file.write(f'[{key}]\n')
#                 n_add_parm = 2 #+2 because of the pt lower/higher limits that are not accounted into the # parameters before
                if 'sM' in key:
                    fun_tmp = Mikkofun_ud
                    for ii in range(6,-1,-1):
                        fun_tmp = fun_tmp.replace(f'[{ii}]', f'[{ii+5}]')
                    fun_str = f'1/({fun_tmp})/({Mikkofun_s})'
#                     n_add_parm+=6
                elif 'udM' in key:
                    fun_str = f'1/({Mikkofun_ud})'
                elif 'M' in key:
                    fun_tmp = Mikkofun_ud
                    for ii in range(6,-1,-1):
                        fun_tmp = fun_tmp.replace(f'[{ii}]', f'[{ii+4}]')
                    fun_str = f'1/({fun_tmp})/({Mikkofun})'
#                     n_add_parm+=6
                else:
                    fun_str = str_poly
                file.write('{1 JetEta 1 JetPt '+fun_str+' Correction L5Flavor}\n')
                ### copy from the positive eta region into the negative
                fit_res = np.vstack([np.hstack([np.flip(fit_res[:,0:2]*-1), np.flip(fit_res[:,2:], 0)]), fit_res])
                for row in fit_res:
                    row[2] = row[2]+2  #+2 because of the pt lower/higher limits that are not accounted into the # parameters before
                    line2write = ('{:>11} '*5+'{:>13} '*(int(row[2])-2)).format(*row[:2], int(row[2]), *np.round(row[3:], 7))+'\n'
                    file.write(line2write);
    print("Saving the corrections with the name = ", txtfile_outname)

# from scipy.optimize import fsolve
from scipy.optimize import brentq
# scipy.optimize.brentq

def find_stationary_pnt_poly(xmin, xmax, *p, degree=4):
    '''Finds the last x within the limits [xmin, xmax], where the derivative of the n-th order polynomial with
    coefficients p, changes sign. n is allowed to be 4 or 3 at the moment. If there is no such point in outputs xmax.  '''
    if degree not in [3,4]:
        raise ValueError(f"Degree can be either 3 or 4. The value given is {degree}")
    if degree==3:
        c0, c1, c2, c3 = p
    elif degree==4:
        c0, c1, c2, c3, c4 = p
    xmin_l = np.log10(xmin)
    xmax_l = np.log10(xmax)
    xs = np.linspace(xmin_l, xmax_l, 1000)
    if degree==3:
        deriv = lambda xs: c1+2*c2*xs+c3*3*xs**2
    elif degree==4:
        deriv = lambda xs: c1+2*c2*xs+c3*3*xs**2+4*c4*xs**3
    signx = np.sign(deriv(xs))
    changes_sign = signx[1:]*signx[:-1]
    last_change_idx = np.where(changes_sign==-1)[0]
    if len(last_change_idx)==0:
        last_change_idx = 0
        return xmax
    else:
        last_change_idx = last_change_idx[-1]
        root = brentq(deriv, xs[last_change_idx], xs[last_change_idx+1] )
        return 10**root

#### initial values borrowed from Winter14 data
#### https://github.com/cms-jet/JECDatabase/blob/master/textFiles/Winter14_V8_MC/Winter14_V8_MC_L5Flavor_AK5Calo.txt/
#### used for fit `response_fnc`
init_vals_2014 = {
    'b':
    [[0.540002, 13.8495, 17.8549, -0.215711, 0.576285, 1.42258],
    [0.73119, 7.52866, 17.3681, -0.078402, 1.21665, 1.69878],
    [0.999952, 0.0322738, -1.05606, -19.6994, 0.720321, -1.58314],
    [0.135913, 7.92441, 3.85698, -0.804604, 1.11911, 0.732041]],
    'c' :
    [[ 0.940259, 0.705481, 0.23917, -0.826926, 0.311473, -0.514041],
    [0.982083, 0.238007, 4.35924, -0.0314618, 5.91028, 1.67749],
    [0.733505, 7.26794, 12.2028, -0.756302, 0.0895257, -1.96324],
    [0.932305, 1.15954, 17.1731, -0.471313, 2.58424, 0.254917]],
    'g' :
    [[0.877892, 3.10194, 1.16568, -677.876, 0.0325026, -12.9485],
    [0.983775, 0.247943, 1.55373, -0.0254802, 3.35748, 1.71263],
    [-0.972548, 38.8683, 2.47151, -44.0233, 0.0901665, -3.15495],
    [1.0655, -0.0680325, -0.509038, -8.59434e+06, 42.6162, 0.357177]],
    'd':
    [[1.28488, -46.3648, 151.749, -0.0108461, 15.4256, 1.63377],
    [ 1.50931, -118.71, 224.19, -0.0196468, 4.62655, 1.51581],
    [0.692016, 8.26488, 11.1655, -0.802769, 0.116182, -1.16094],
    [1.01244, -0.0926519, -0.12138, -3.69494e+07, 7.15634, -0.625288]],  
    'u':
    [[1.28488, -46.3648, 151.749, -0.0108461, 15.4256, 1.63377],
    [ 1.50931, -118.71, 224.19, -0.0196468, 4.62655, 1.51581],
    [0.692016, 8.26488, 11.1655, -0.802769, 0.116182, -1.16094],
    [1.01244, -0.0926519, -0.12138, -3.69494e+07, 7.15634, -0.625288]],  
    's':
    [[1.28488, -46.3648, 151.749, -0.0108461, 15.4256, 1.63377],
    [ 1.50931, -118.71, 224.19, -0.0196468, 4.62655, 1.51581],
    [0.692016, 8.26488, 11.1655, -0.802769, 0.116182, -1.16094],
    [1.01244, -0.0926519, -0.12138, -3.69494e+07, 7.15634, -0.625288]],
    'all':
    [[0.540002, 13.8495, 17.8549, -0.215711, 0.576285, 1.42258],
    [0.73119, 7.52866, 17.3681, -0.078402, 1.21665, 1.69878],
    [0.999952, 0.0322738, -1.05606, -19.6994, 0.720321, -1.58314],
    [0.135913, 7.92441, 3.85698, -0.804604, 1.11911, 0.732041]],  
    'ud':
    [[1.28488, -46.3648, 151.749, -0.0108461, 15.4256, 1.63377],
    [ 1.50931, -118.71, 224.19, -0.0196468, 4.62655, 1.51581],
    [0.692016, 8.26488, 11.1655, -0.802769, 0.116182, -1.16094],
    [1.01244, -0.0926519, -0.12138, -3.69494e+07, 7.15634, -0.625288]],  
    'q':
    [[1.28488, -46.3648, 151.749, -0.0108461, 15.4256, 1.63377],
    [ 1.50931, -118.71, 224.19, -0.0196468, 4.62655, 1.51581],
    [0.692016, 8.26488, 11.1655, -0.802769, 0.116182, -1.16094],
    [1.01244, -0.0926519, -0.12138, -3.69494e+07, 7.15634, -0.625288]],
    
}

## Initial values for the two gaussian fit
init_two_gaus = [3,0,1,2,0,1,2,3,4]
# Better starting fit values I found
init_vals_2014['b'][0] = [ 9.81014871e-01, -6.46744813e-03, -1.05658840e+00,  5.35445486e+03, 2.99200015e+01,  1.21399356e+02]
init_vals_2014['b'][3] = [ 9.81014871e-01, -6.46744813e-03, -1.05658840e+00,  5.35445486e+03, 2.99200015e+01,  1.21399356e+02]


from JetEtaBins import JetEtaBins, PtBins
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

def poly4(x, *p):
    c0, c1, c2, c3, c4 = p
    xs = np.log10(x)
    res = c0+c1*xs+c2*xs**2+c3*xs**3+c4*xs**4
    return res

def poly4lims(x, xmin, xmax, *p):
    xcp = x.copy()
    lo_pos = xcp<xmin
    hi_pos = xcp>xmax
    xcp[lo_pos] = xmin
    xcp[hi_pos] = xmax
    return poly4(xcp, *p)

def poly3(x, *p):
    c0, c1, c2, c3 = p
    xs = np.log10(x)
    res = c0+c1*xs+c2*xs**2+c3*xs**3
    return res

def poly3_1_x(x, *p):
    c0, c1, c2, c3, c4 = p
    xs = np.log10(x)
    res = c0+c1*xs+c2*xs**2+c3*xs**3+c4/xs**10
    return res

def poly3_1_xn(x, *p):
    c0, c1, c2, c3, c4, c5 = p
    xs = np.log10(x)
    res = c0+c1*xs+c2*xs**2+c3*xs**3+c4/xs**c5
    return res

def poly3_1_x3(x, *p):
    c0, c1, c2, c3, c4, c5 = p
    xs = np.log10(x)
    res = c0+c1*xs+c2*xs**2+c3*xs**3+c4/(xs+c5)
    return res

from helpers import gauss
Mikkofun_ud = "[0]+[1]*pow(x,[2])+[3]*pow(x/1000.,[4])+[5]/x"
Mikkofun = ("1+[0]*pow(x,[1])"
		    +"+[2]/x"
		    +"+[3]*log(x)/x"
		    # +"+[4]*pow(x,fabs([5]))"
            )
Mikkofun_s = ("1+[0]*pow(x,[1])"
		     +"+[2]*1/(2.5066*[4])*exp(-0.5*((log(x)-[3])/[4])*((log(x)-[3])/[4]))")

def Mikko_fun1(x, *p):
    c0, c1, c2, c3 = p
    xs = np.log10(x)
    res = (1+c0*pow(xs,c1)
	+ c2/xs
	+ c3*np.log10(xs)/xs
	# + c4*pow(xs,np.abs(c5))
    )
    return res

def Mikko_fun2(x, *p):
    c0, c1, c2, c3, c4 = p
    xs = np.log10(x)
    res = (1+c0*pow(xs,c1)
           + gauss(xs, *[c2, c3, c4]) )
    return res

def Mikko_fun_ud(x, *p):
    c0, c1, c2, c3, c4, c5 = p
    res = c0+c1*x**c2 +c3*(x/1000)**c4+c5/x
    return res

def poly3lims(x, xmin, xmax, *p):
    xcp = x.copy()
    lo_pos = xcp<xmin
    hi_pos = xcp>xmax
    xcp[lo_pos] = xmin
    xcp[hi_pos] = xmax
    return poly3(xcp, *p)

def response_fnc_simp(x, *p):
    p0, p3, p4, p5 = p
    logx = np.log10(x)
    return p0 + (p3*np.exp(-p4*((logx-p5)*(logx-p5))))

def response_fnc(x, *p):
    p0, p1, p2, p3, p4, p5 = p
    logx = np.log10(x)
    return p0+(p1/((logx**2)+p2)) + (p3*np.exp(-p4*((logx-p5)*(logx-p5))))

def two_gaus_fnc(x, *p):
    p0, p1, p2, p3, p4, p5, p6, p7, p8 = p
    return (  p0
            + (p1/((np.log10(x)**2)+p2))
            + (p3*np.exp(-p4*((np.log10(x)-p5)*(np.log10(x)-p5))))
            + (p6*np.exp(-p7*((np.log10(x)-p8)*(np.log10(x)-p8))))
           )

def response_fnc_raw(x, p0, p1, p2, p3, p4, p5):
    response_fnc(x, *[p0, p1, p2, p3, p4, p5])

import ROOT

def draw_result_root(fun, xvals, graph, num_parameters, fit_func, p0):
    # Create canvas for plotting
    canvas = ROOT.TCanvas("fit_canvas", "Fit Canvas", 800, 600)

    # Plot the data points
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.0)
    graph.SetMarkerColor(ROOT.kBlack)
    graph.Draw("AP")

    # Plot the initial fit with dashed line
    initial_fit_func = ROOT.TF1("initial_fit", fun, min(xvals), max(xvals), num_parameters)
    initial_fit_func.SetParameters(*p0)
    initial_fit_func.SetLineColor(ROOT.kBlue)
    initial_fit_func.SetLineStyle(ROOT.kDashed)
    initial_fit_func.Draw("same")

    # Plot the final fit with solid line
    fit_func.SetLineColor(ROOT.kRed)
    fit_func.Draw("same")

    # Add legend
    legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
    legend.AddEntry(graph, "Data", "p")
    legend.AddEntry(initial_fit_func, "Initial Fit", "l")
    legend.AddEntry(fit_func, "Final Fit", "l")
    legend.Draw()

    # Show the canvas
    canvas.Draw()
    canvas.SaveAs('test.pdf')

def fit_with_root(fun, xvals, yvals, p0, sigma, bounds=None):
    is_good = yvals != 0
    xvals = xvals[is_good]
    yvals = yvals[is_good]
    sigma = sigma[is_good]
    num_parameters=len(p0)
    graph = ROOT.TGraphErrors(len(xvals), np.array(xvals, dtype=np.double), np.array(yvals, dtype=np.double),
                             np.array([0.0]*len(xvals), dtype=np.double), np.array(sigma, dtype=np.double))
    fit_func = ROOT.TF1("fit_function", fun, min(xvals), max(xvals), 6)
#     print(p0)
    fit_func.SetParameters(*p0)
    if bounds is not None:
        for ii, bound in enumerate(bounds):
            if None==bound or np.inf in bound or -np.inf in bound:
                continue
            fit_func.SetParLimits(ii, *bound)
    graph.Fit(fit_func, "R")
    fit_params = [fit_func.GetParameter(i) for i in range(num_parameters)]
    return fit_params

from coffea.lookup_tools.jme_standard_function import wrap_formula
def root_to_fun(string, n_parms):
    ps = []
    string = string.replace('pi', 'numpy.pi')
    for ii in range(n_parms):
        string = string.replace(f'[{ii}]',f'p{ii}')
        ps.append(f'p{ii}')
    return wrap_formula(string, ['x']+ps)

from scipy.stats import chi2
import mplhep as hep
from fileNames.available_datasets import legend_labels
figdir = "fig/median_correction_fits"


fits2plot = { ##name: [function, initial values, # parameters]
             "MC truth": [response_fnc, init_vals_2014, 6],
             "MC truth extended": [two_gaus_fnc, init_two_gaus, 9],
             "Poly, n=4": [poly4, [ 1, 1, 1, 1, 1], 5],
             "Poly, n=3": [poly3, [ 1, 1, 1, 1], 4],
             }

main_fit = "Poly, n=4"

def inflate_smallest_std(stds):
    stds[np.isnan(stds)] = np.inf
    min_positions = np.argmin(stds, axis=0)
    second_min_values = np.sort(stds, axis=0)[1, :]
    inf_pos = (second_min_values == np.inf) & (np.min(stds, axis=0)!=np.inf)
    second_min_values[inf_pos] = np.min(stds, axis=0)[inf_pos] 
    ax2 = np.arange(len(min_positions))
    stds[min_positions, ax2] = second_min_values
    stds[stds==np.inf] = np.nan
    return stds

def fit_corrections(etaidx, data_dict, flav, data_tags,
                           fits2plot, main_fit,
                    figdir2=figdir+'correction_fit',
                    jeteta_bins=JetEtaBins(), pt_bins=PtBins(),
                    plot_initial_val:list=[],
                    use_recopt=True,
                    maxlimit_static_pnt=True,
                    inverse=True,
                    max_point_fit_idx=3,
                    min_pt_val = None,
                    max_ptval=None,
                    min_rel_uncert=0.0,
                    min_rel_uncert_relative=0,
                    inflate_smallest_std_bool=True,
                    show_original_errorbars=False,
                    ncoefs_out=5,
                    saveplots=True,
                    colors = None,
                    fit_sample=None,
                    custom_jet_legend=None ):
    """ fit the data and plot
    etaidx: index of the eta bin from the data in the `data_dict` to fit
    data_dict: dictionary of the data to fit
    flav: flavor of the jets to fit 
    data_tags: tags of the data to fit. The same order as in data_dict.keys()
    fits2plot: dictionary of the fit functions to use for fitting and their initial values and number of input parameters. The keys are the names of the fits.
    main_fit: the name of the fit from `fits2plot` to use for the output. The output will be the coefficients of this fit.
    figdir2: directory to save the plots
    jeteta_bins: JetEtaBins object, eta bins used for the data
    pt_bins: PtBins object, pt bins used for the data
    plot_initial_val: list of functions for which to plot the initial values. Default: empty list
    use_recopt: if True, use the reco pt as the x values for the fit. Otherwise, use the pt bin centres (thus, same values for all the data samples).
    maxlimit_static_pnt: if True, find the last point where the derivative of the fit function changes sign and use it as the upper limit of the fit range.
    inverse: if True, fit the inverse of the data (to get the correction).
    max_point_fit_idx: if `maxlimit_static_pnt` is True but the static point is too far to the left, use the `max_point_fit_idx`-th point from the end as the upper limit of the fit range.
    max_ptval: if `maxlimit_static_pnt` is False, use this value as the upper limit of the fit range.
    min_rel_uncert: minimum relative uncertainty to use for the fit. If 0, use the `min_rel_uncert_relative` to define the minimum uncertainty relative to the range.
    min_rel_uncert_relative: minimum relative uncertainty to use for the fit relative to the range. If 0, use the `min_rel_uncert` to define the minimum uncertainty.
    show_original_errorbars: if True, show the original errorbars of the data points.
    ncoefs_out: number of coefficients of the output function. If less than the number of coefficients of the main fit, pad with zeros.
    saveplots: if True, save the plots
    colors: dictionary of the colors to use for the data samples. The keys are the data tags.
    fit_sample: if not None, use all samples to fit, otherwise index of the sample to fit
    """
    ###################### Logistics with the input ######################
    keys = [key for key in data_dict.keys()]
    if maxlimit_static_pnt:
        if max_ptval==None:
            max_ptval=5000
        else:
            raise ValueError(f"Remove `max_ptval` if using `maxlimit_static_pnt`. `max_ptval` given as {max_ptval} but `maxlimit_static_pnt` is set to {maxlimit_static_pnt}.")
    elif max_ptval==None:
        max_ptval=5000

    if min_pt_val is None:
        min_pt_val = np.min(pt_bins.centres)

    if len(data_dict)<2:
        inflate_smallest_std_bool=False

    ### pt limits for the fit
    ptmin_idx = np.searchsorted(pt_bins.centres, min_pt_val, side='left') #-1: so that the first bin includes the value
    ptmax_idx = np.searchsorted(pt_bins.centres, max_ptval, side='right')
    data_range = tuple([range(ptmin_idx,ptmax_idx), etaidx])

    yvals = np.array([data_dict[key]["Median"][data_range] for key in keys])
    stds  = np.array([data_dict[key]["MedianStd"][data_range] for key in keys])
    reco_pts  = np.array([data_dict[key]["MeanRecoPt"][data_range] if len(data_dict[key]["MeanRecoPt"].shape)==2 else data_dict[key]["MeanRecoPt"][data_range] for key in keys])

    yvals[yvals==0] = np.nan
    if inverse==True:
        yvals = 1/yvals
        ### Error propagation
        stds = yvals**2*stds

    validpt_mask = ~(np.isnan(yvals) | np.isinf(yvals) | (yvals==0))
    if not fit_sample is None:
        ## don't use the samples of indices that don't match fit_sample for fitting 
        validpt_mask[np.logical_not(np.arange(len(validpt_mask)) == fit_sample)] = False
    xvals = reco_pts if use_recopt else np.array([pt_bins.centres[data_range]]*len(yvals))

    ### Put the minimum limit on the relative uncertainty to min_rel_uncert
    # only the first case makes sence here imo, because it defines the minimum uncertainty relative to the range
    # but the second case kept for a while to be consistent with the not
    if min_rel_uncert<=0:
        min_rel_uncert_tmp = min_rel_uncert_relative*(np.nanmax(yvals)-np.nanmin(yvals))
        print('men_rel_uncert_tmp:', min_rel_uncert_tmp)
    else:
        min_rel_uncert_tmp = min_rel_uncert
    stds_orig = stds.copy()
    if inflate_smallest_std_bool:
        stds = inflate_smallest_std(stds)
    where_limit_std = (stds/yvals)<min_rel_uncert_tmp
    stds[where_limit_std] = min_rel_uncert_tmp*yvals[where_limit_std]


    if np.sum(validpt_mask)==0:
        fit_res_new = np.concatenate([[jeteta_bins.edges[etaidx], jeteta_bins.edges[etaidx+1],
                               ncoefs_out, 
                               pt_bins.centres[ptmin_idx], pt_bins.centres[ptmax_idx-1]],
                              [1,0,0,0,0] ])
        print("No points to fit. Returning a unity function.")
        return fit_res_new

    means2fit = yvals[validpt_mask]
    means_unc2fit = np.abs(stds[validpt_mask])
    ptbins2fit = xvals[validpt_mask]
    fit_min_lim = np.min(ptbins2fit)
    fit_max_lim = np.max(ptbins2fit)

    xplot = np.linspace(fit_min_lim, fit_max_lim, 1000)

    ###################### Fits ######################
    
    fitres = {}
    for fit in fits2plot:
        if len(means2fit)>(fits2plot[fit][2]+2):
            init_vals = fits2plot[fit][1]
            fit_kwargs = {}
            if type(init_vals) is dict:
                fit_kwargs["p0"]= init_vals[flav][etaidx]
            else:
                fit_kwargs["p0"]= init_vals

            if len(fits2plot[fit])>3 and fits2plot[fit][3] is not None:
                fit_kwargs["bounds"] = fits2plot[fit][3]

            if len(fits2plot[fit]) == 5: #if the fit is defined with a string, use root to fit, otherwise use scipy
                fit_kwargs["sigma"] = means_unc2fit
                p_err = fit_with_root(fits2plot[fit][4], ptbins2fit, means2fit, **fit_kwargs)
            else:
                try:
                    p, arr = curve_fit(fits2plot[fit][0], ptbins2fit, means2fit, **fit_kwargs)
                    fit_kwargs["sigma"] = means_unc2fit
                    fit_kwargs["p0"] = p
                    p_err, arr = curve_fit(fits2plot[fit][0], ptbins2fit, means2fit, **fit_kwargs)
                except(RuntimeError):
                    print(f"{fit} fit failed")
                    p, p_err = [[np.nan]*fits2plot[fit][2]]*2
            # except(TypeError):
        else:
            print(f"Too little points for {fit} fit")
            p, p_err = [[np.nan]*fits2plot[fit][2]]*2
        fitres[fit] = p_err

    # convert strings for the root fits to python lamda functions
    for fit in fits2plot:
        if len(fits2plot[fit]) == 5:
            fits2plot[fit][0] = root_to_fun(fits2plot[fit][4], fits2plot[fit][2])
    
    chi2s = {fit: np.sum((fits2plot[fit][0](ptbins2fit, *fitres[fit]) - means2fit)**2/means_unc2fit**2)
                 for fit in fits2plot}
    Ndofs = {fit: len(ptbins2fit) - fits2plot[fit][2] for fit in fits2plot}

    ## if chi2 of n=3 polynomial outside the 2 one-sided std of the chi2 distribution, use n=3 polynomial.
    if main_fit == "Poly, n=4" and "Poly, n=3" in fits2plot and chi2.ppf(1-0.158, Ndofs["Poly, n=3"])>chi2s["Poly, n=3"]:
        main_fit = "Poly, n=3"
    
    print(f"Using the {main_fit} fit results ")
    if maxlimit_static_pnt:
        if main_fit=="Poly, n=3":
            fit_degree = 3 
        elif main_fit=="Poly, n=4":
            fit_degree = 4
        else:
            raise ValueError(f"Main fit is {main_fit} but the derivative for the static point is not defined for this fit.")
        fit_max_lim_new = find_stationary_pnt_poly(fit_min_lim, fit_max_lim, *fitres[main_fit], degree=fit_degree)
    else:
        fit_max_lim_new = fit_max_lim
    fit_max_lim_idx = np.searchsorted(np.sort(ptbins2fit), fit_max_lim_new, side="right")
    if maxlimit_static_pnt & (fit_max_lim_idx==len(ptbins2fit)) | (fit_max_lim_idx<=len(ptbins2fit)-max_point_fit_idx):
        # static point is too low or the last point that usually fluctuates out
        fit_max_lim_idx = len(ptbins2fit)-max_point_fit_idx
        fit_max_lim_new = np.sort(ptbins2fit)[fit_max_lim_idx]
    xplot_max_new = np.searchsorted(xplot, fit_max_lim_new)

    main_fit_res = fitres[main_fit]
    if ncoefs_out<len(main_fit_res):
        raise ValueError(f"ncoefs is smaller than the number of coefficients of the main fit."
                        +f"ncoefs_out={ncoefs_out}, len(main_fit_res)={len(main_fit_res)}."
                        +f"Either raise the number of coefficients of the output or choose a different output function.")
    fit_res_new = np.concatenate([[jeteta_bins.edges[etaidx], jeteta_bins.edges[etaidx+1],
                                ncoefs_out, 
                                fit_min_lim, fit_max_lim_new],
                                np.pad(main_fit_res, (0, ncoefs_out-len(main_fit_res))) ])

    curve_yvals = {fit: fits2plot[fit][0](xplot, *fitres[fit]) for fit in fits2plot}


    ###################### Plotting ######################
    fig, ax = plt.subplots()
    legend_original_initiated = False
    for yval, std, reco_pt, std_orig, data_tag in zip(yvals, stds, xvals, stds_orig, data_tags):
        if colors is not None:
            color = colors[data_tag]
        else:
            color = next(ax._get_lines.prop_cycler)['color']
        if "merged" in data_tag:
            fun_kwargs = {'mfc':'none', 'markeredgewidth':1.2}
        else:
            fun_kwargs = {}
        plt.errorbar(reco_pt, yval, yerr=std, marker='o',
                    linestyle="none", label=data_tag, #{jeteta_bins.idx2plot_str(etaidx)}',
                    capsize=1.7, capthick=0.9, linewidth=1.0, color=color, **fun_kwargs)
        if show_original_errorbars:
            if legend_original_initiated:
                plt.errorbar(reco_pt, yval, yerr=std_orig, marker='o',
                    linestyle="none",
                    capsize=1.5, capthick=0.7, linewidth=0.9, markersize=0, color='black', **fun_kwargs) #, color='blue')
            else:
                legend_original_initiated = True
                plt.errorbar(reco_pt, yval, yerr=std_orig, marker='o',
                    linestyle="none", label=f'Original errorbars',
                    capsize=1.5, capthick=0.7, linewidth=0.9, markersize=0, color='black', **fun_kwargs)
    
    for fit in fits2plot:
        if np.isnan(chi2s[fit]): 
            lab = f'{fit} failed'
        else:
            polytxt3 = ', selected' if fit==main_fit and len(fits2plot)>1 else ''
            lab= fit+r'; $\chi^2/N_{dof} = $'+r' {0:0.3g}/{1:0.0f}'.format(chi2s[fit], Ndofs[fit])+polytxt3
        ax.plot(xplot, curve_yvals[fit], label=lab, markersize=0);
        if maxlimit_static_pnt and fit==main_fit:
            ax.plot(xplot[:xplot_max_new], curve_yvals[fit][:xplot_max_new], label=f'{fit}; range chosen', markersize=0, linewidth=0.8); #
        # if plot_initial_val and fit=="MC truth":
        if fit in plot_initial_val:
            yvals_init = fits2plot[fit][0](xplot, *fits2plot[fit][1])
            ax.plot(xplot, yvals_init, label=f"Initial values for {fit}", markersize=0);
        
    ###################### Plot formalities ######################
    ylim_tmp = ax.get_ylim()
    ylim_pad = (ylim_tmp[1] - ylim_tmp[0])/1.6
    ax.set_ylim(ylim_tmp[0], ylim_tmp[1]+ylim_pad)

    ax.set_xlabel(r'$p_{T,reco}$ (GeV)')
    ylabel = r'correction (1/median)' if inverse else r'median response'
    ax.set_ylabel(ylabel)
    ax.set_xscale('log')

    ax.set_xticks([])
    ax.set_xticks([20, 50, 100, 200, 500, 1000])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    # hep.cms.label("Preliminary", loc=0, data=False, ax=ax)
    hep.cms.label("Private work", loc=0, data=False, ax=ax, rlabel='')
    jet_legend = custom_jet_legend if custom_jet_legend is not None else f'{flav} jets'
    hep.label.exp_text(text=f'{jeteta_bins.idx2plot_str(etaidx)}; '+jet_legend, loc=2, fontsize=mpl.rcParams["font.size"]/1.15)

    ### hack put errorbars before the curves in the legend
    #get handles and labels
    handles, labels = plt.gca().get_legend_handles_labels()

    #specify order of items in legend
    nfits = len(fits2plot)+maxlimit_static_pnt
    order = np.concatenate([np.arange(nfits,nfits+len(keys)+show_original_errorbars), np.arange(nfits)]) 

    #add legend to plot
    odered_handles = [handles[idx] for idx in order]
    ordered_labels = [labels[idx] for idx in order]
 
    ax.legend(odered_handles,
              ordered_labels,
              loc="upper left", bbox_to_anchor=(0.01, 0.90)) #, prop={'size': 10}
    
    if saveplots:
        figdir2 = (figdir+'/'+data_tag.replace(legend_labels["ttbar"]["lab"], 'ttbar').replace(', ', '-')
                    .replace(" ", "_").replace("+", "_").replace('(', '').replace(')', '').replace('/', '').replace('\n', '')
                )
        if not os.path.exists(figdir2):
            os.mkdir(figdir2)
            
        add_name = f'correction_fit_{flav}_'+jeteta_bins.idx2str(etaidx)
        fig_name = figdir2+'/'+add_name
            
        print("Saving plot with the name = ", fig_name+".pdf / .png")
        plt.savefig(fig_name+'.pdf');
        plt.savefig(fig_name+'.png');

    plt.show();
    plt.close();
    
    return fit_res_new