legend_labels = {"ttbar": {"lab":"$t\overline{\, t\!}$ ",
                            "short": "ttbar"},
                "QCD": {"lab":"QCD",
                        "short": "QCD"},
                "DY": {"lab":"DY",
                        "short": "DY"}
                            }

'''
Lift of available datasets.

Datasets can either given by a path to a file storing file names or by a path to a file storing datasets
and their cross-sections

data_tag: [file_name_path, xsec_path, legend_label]
'''

dataset_labels = {
    "Pythia-TTBAR": legend_labels["ttbar"]["lab"]+'Pow+Py8',
    "Pythia-TTBAR_TTbarSemiLep": legend_labels["ttbar"]["lab"]+', l+jets Pow+Py8',
    "Pythia-TTBAR_TTbarDilep": legend_labels["ttbar"]["lab"]+r', 2l2$\nu$ Pow+Py8',
    "Pythia-TTBAR_TTbarFullHad": legend_labels["ttbar"]["lab"]+', all had Pow+Py8',
    "Herwig-TTBAR": legend_labels["ttbar"]["lab"]+'Pow+Her7',
    "DY-MG-Py": 'ZJets MG+Py8',
    "DY-MG-Her": 'ZJets MG+Her7',
    "QCD-MG-Py": 'QCD MG+Py8',
    "QCD-MG-Her": 'QCD MG+Her7',
    "QCD-Py": 'QCD Py8',
    "DY-FxFx": 'ZJets FxFx',
    "scaled_pion_kaon": 'pion up/ kaon up',
    "scaled_pion": 'pion up',
    "scaled_times2_pion": 'pion up x2',
    "scaled_times5_pion": 'pion up x5',
    "scaled_times10_pion": 'pion up x10',
    "scaled_times100_pion": 'pion up x100',
    "not_scaled_pion": 'pion central'
}

# putting labels in the dataset_dictionary was an older way of specifying the labels but it does not work if e.g. Pythia-TTBAR consists of three channels and one needs labels for each
dataset_dictionary = {
    "Pythia-TTBAR": [None, 'fileNames/TTBAR_Pythia_20UL18/xsecs_TTBAR_Pow-Py8.txt', legend_labels["ttbar"]["lab"]+'Pow+Py8'],
    "Pythia-semilep-TTBAR": ['fileNames/TTBAR_Pythia_20UL18/TTToSemi20UL18_JMENano.txt', 1, legend_labels["ttbar"]["lab"]+', l+jets Pow+Py8'],
    "Pythia-dilep-TTBAR": ['fileNames/TTBAR_Pythia_20UL18/TTToDilep20UL18_JMENano.txt', 1, legend_labels["ttbar"]["lab"]+r', 2l2$\nu$ Pow+Py8'],
    "Pythia-fullhad-TTBAR": ['fileNames/TTBAR_Pythia_20UL18/TTToHad20UL18_JMENano.txt', 1, legend_labels["ttbar"]["lab"]+', all had Pow+Py8'],
    "Pythia-semilep-TTBAR_0-500": ['fileNames/TTBAR_Pythia_20UL18/TTToSemi20UL18_JMENano_redi_0-500.txt', 1, legend_labels["ttbar"]["lab"]+'Pow+Py8'],
    "Pythia-semilep-TTBAR_501-end": ['fileNames/TTBAR_Pythia_20UL18/TTToSemi20UL18_JMENano_redi_501-end.txt', 1, legend_labels["ttbar"]["lab"]+'Pow+Py8'],
    "Pythia-non-semilep-TTBAR": [None, 'fileNames/TTBAR_Pythia_20UL18/xsecs_TTBAR_Pow-Py8-non-semilep.txt', legend_labels["ttbar"]["lab"]+'Pow+Py8'],
    "Herwig-TTBAR": ['fileNames/TT20UL18_JMENano_Herwig.txt', 1, legend_labels["ttbar"]["lab"]+'Pow+Her7'],
    "DY-MG-Py":     ['fileNames/DYJets_MG-Py.txt', 1, 'ZJets MG+Py8'],
    "DY-MG-Her":    ['fileNames/DYJets_MG-Her.txt', 1, 'ZJets MG+Her7'],
    "QCD-MG-Py":    [None, 'fileNames/QCD_MG_Py8_20UL18/xsecs_QCD_MG_py8.txt', 'QCD MG+Py8'],
    # "QCD-MG-Her":   [None, 'fileNames/QCD_MG_Py8_20UL18/xsecs_QCD_MG_py8.txt', 'QCD MG+Her7'],
    "QCD-MG-Her":   [None, 'fileNames/QCD_Herwig_20UL18/xsecs_QCD_Herwig_corrected.txt', 'QCD MG+Her7'],
    "QCD-Py":       ['fileNames/QCD20UL18_JMENano.txt', 1373000000, 'QCD Py8'],
    "DY-FxFx":      ['fileNames/DYJets.txt', 1, 'ZJets FxFx'],
    "noJME-QCD-Py_pu": ['fileNames/fileNames_QCD20UL18.txt', 1, 'QCD Py8'],
    "scaled_pion_kaon": ['fileNames_pion_response/fileNames_scaled_pion_kaon.txt', 1, 'pion up/ kaon up'],
#     "scaled_pion": ['fileNames_pion_response/fileNames_scaled_pion.txt', 1, 'pion up'],
    "scaled_pion": ['fileNames_pion_response/fileNames_scaled_times1_pion.txt', 1, 'pion up'],
    "scaled_times2_pion": ['fileNames_pion_response/fileNames_scaled_times2_pion.txt', 1, 'pion up'],
    "scaled_times5_pion": ['fileNames_pion_response/fileNames_scaled_times5_pion.txt', 1, 'pion up'],
    "scaled_times10_pion": ['fileNames_pion_response/fileNames_scaled_times10_pion.txt', 1, 'pion up'],
    "scaled_times100_pion": ['fileNames_pion_response/fileNames_scaled_times100_pion.txt', 1, 'pion up'],
    "not_scaled_pion": ['fileNames_pion_response/fileNames_not_scaled_pion.txt', 1, 'pion central'],}