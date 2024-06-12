'''Extract the overall pion uncertainty factor from the difference between the scaled 10 and not scaled samples.
The scaled x10 sample is used since the not scaled sample only shows some statistical fluctuations
x10 on the other hand shows the uncertainties more prominently. This maybe overestimates One can divide them /10 afterwards.
'''

from data_tools import read_or_recreate_data
import numpy as np

eta_binning  = "HCalPart"  ### HCalPart, JERC, CoarseCalo, CaloTowers, Summer20Flavor, onebin;
pt_binning  = "onebin"  ### MC_truth, Uncert, Coarse, onebin;
eta_binning_str = '_'+eta_binning if eta_binning != "HCalPart" else ''
pt_binning_str = '_pt-'+pt_binning if pt_binning != "MC_truth" else ''
# flav = 'b'
# ratio_x10 = np.array(read_or_recreate_data('_L5_scaled_times10_pion'+eta_binning_str+pt_binning_str    +'_split_antiflav', '../out_txt')['data'][flav]['Median'])/np.array(data['data'][flav]['Median'])
# ratio = np.array(read_or_recreate_data('_L5_scaled_times10_pion'+eta_binning_str+pt_binning_str    +'_split_antiflav', '../out_txt')['data'][flav]['Median'])/np.array(data_x10['data'][flav]['Median'])
# uncertainty = ratio_x10-ratio

flavor_unc_name_dict = {'c': 'PureCharm', 'b':'PureBottom', 'ud':'PureUpDown',
                        'q':'PureQuark', 's':'PureStrange',
                        'cbar': 'PureAntiCharm', 'bbar':'PureAntiBottom', 'udbar':'PureAntiUpDown',
                        'qbar':'PureAntiQuark', 'sbar':'PureAntiStrange'}
flavors = ['c', 'b', 'ud', 'q', 's']

lines = {}
data_x10 = read_or_recreate_data('_L5_scaled_pion'+eta_binning_str+pt_binning_str    +'_split_antiflav', '../out_txt')
data = read_or_recreate_data('_L5_not_scaled_pion'+eta_binning_str+pt_binning_str    +'_split_antiflav', '../out_txt')
for flav in flavors:
    aflav = flav+'bar'
    flav_name = 'Pion'+flavor_unc_name_dict[flav]
    aflav_name = 'Pion'+flavor_unc_name_dict[aflav]
    ratio_x10 = np.array(data_x10['data'][flav]['Median'])/np.array(data['data'][flav]['Median'])
    ratio = np.array(data_x10['data'][flav]['Median'])/np.array(data_x10['data'][flav]['Median'])
    shift = (ratio_x10-ratio)[0]

    ratio_x10 = np.array(data_x10['data'][aflav]['Median'])/np.array(data['data'][aflav]['Median'])
    ratio = np.array(data_x10['data'][aflav]['Median'])/np.array(data_x10['data'][aflav]['Median'])
    shift_bar = (ratio_x10-ratio)[0]
    uncertainty = (shift-shift_bar)/2

    lines[flav_name] = []
    lines[aflav_name] = []
    lines[flav_name].append('{1 JetEta 1 JetPt "" Correction JECSource} \n')
    bins_txt_rev = ["-5.4 -3.2", "-3.2 -2.5", "-2.5 -1.3", "-1.3 0"]
    bins_txt = ["0 1.3", "1.3 2.5", "2.5 3.2", "3.2 5.4"]
    for binii in range(len(bins_txt_rev)):
        unc = np.round(uncertainty[3-binii],6)
        lines[flav_name].append(f'{bins_txt_rev[binii]} 6 9.0 {unc} {unc} 6538.0 {unc} {unc}\n')
    for binii in range(len(bins_txt)):
        unc = np.round(uncertainty[binii],6)
        lines[flav_name].append(f'{bins_txt[binii]} 6 9.0 {unc} {unc} 6538.0 {unc} {unc}\n')

    uncertainty2 = -uncertainty
    lines[aflav_name].append('{1 JetEta 1 JetPt "" Correction JECSource} \n')
    bins_txt_rev = ["-5.4 -3.2", "-3.2 -2.5", "-2.5 -1.3", "-1.3 0"]
    bins_txt = ["0 1.3", "1.3 2.5", "2.5 3.2", "3.2 5.4"]
    for binii in range(len(bins_txt_rev)):
        unc = np.round(uncertainty2[3-binii],6)
        lines[aflav_name].append(f'{bins_txt_rev[binii]} 6 9.0 {unc} {unc} 6538.0 {unc} {unc}\n')
    for binii in range(len(bins_txt)):
        unc = np.round(uncertainty2[binii],6)
        lines[aflav_name].append(f'{bins_txt[binii]} 6 9.0 {unc} {unc} 6538.0 {unc} {unc}\n')
        # lines[flav_name].append(f'-5.4 5.4 6 {uncertainty} {uncertainty} 6538.0 {uncertainty} {uncertainty}')




def append_uncertainties(file_path, newlines, tag="_pion"):
    ''' Append the flavor-antiflavor uncertainty file with the pion uncertainty
    '''
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # first real all the lines and copy to the new file. Check if there is no 'PionPureAntiQuark' etc.. The write new lines with the new uncertainties 
    output_lines = []
#     matched_section = False
    for line in lines:
        if line.startswith("["):
            section_title = line.strip()[1:-1]
            if 'PionPure' in section_title:
                raise ValueError(f"Section '{section_title}' already exists in the file '{file_path}'")
            else:
                output_lines.append(line)
        else:
            output_lines.append(line)

    sections_titles = ['PionPureQuark', 'PionPureAntiQuark', 'PionPureStrange', 'PionPureAntiStrange', 'PionPureCharm', 'PionPureAntiCharm', 'PionPureBottom', 'PionPureAntiBottom', 'PionPureUpDown', 'PionPureAntiUpDown']
    for section_title in sections_titles:
        output_lines.append(f"[{section_title}]\n")
        for line in newlines[section_title]:
            output_lines.append(line)

    outfilename = file_path[:-4]+tag+'.txt'
    print(f"Writing to a file '{outfilename}'")
    with open(outfilename, 'w') as file:
        file.writelines(output_lines)
    
append_uncertainties("../Summer19UL18_V5_MC/Summer19UL18_V5_MC_UncertaintySources_AK4PFchs_run2flavor_antiflavor.txt", lines)
    # for section in 
#             if line.startswith("{"):
#                 output_lines.append(line)
#             else:
#                 elements = line.split()
#                 updated_vars = [[float(elements[i]), mult*float(elements[i+1]), mult*float(elements[i+2])] for i in range(3, len(elements), 3)]
#                 updated_line = ' '.join(elements[:3]) + ' ' + ' '.join([f'{x[0]} {x[1]} {x[2]}' for x in updated_vars]) + ' \n'
#                 output_lines.append(updated_line)

        
