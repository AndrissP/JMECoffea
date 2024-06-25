import json
import numpy as np

class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                            np.int16, np.int32, np.int64, np.uint8,
                            np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32,
                              np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def save_json(df, pt_bins, jeteta_bins, filename):
    binning = {"pt_edges": pt_bins.edges,
           "eta_edges": jeteta_bins.edges}
    df_to_json2 = {"binning": binning,
              "data": df}
    json_string = json.dumps(df_to_json2, cls=NumpyEncoder)
    # Write the JSON string to the file
    with open(filename, 'w') as json_file:
        json_file.write(json_string)
    print(f"Saved {filename}")

def save_json_fractions(df, pt_bins, jeteta_bins, filename):
    df2 = {samp:{"mean": df[samp][0], "std": df[samp][1]} for samp in df}
    save_json(df2, pt_bins, jeteta_bins, filename)

import ROOT

def save_fractions_root(df, pt_bins, jeteta_bins, filename):
    x_values = pt_bins.centres

    # Create a ROOT file
    root_file = ROOT.TFile(filename, "RECREATE")

    # Loop over the keys in the df dictionary
    for key in df.keys():
        # Create a directory for each key in the ROOT file
        root_dir = root_file.mkdir(key)
        root_dir.cd()

        # Loop over the flavors
        for flav in df[key][0].keys():

            # Create a TGraphErrors for each column (bin)
            for i in range(jeteta_bins.nbins):
                # Get the data from the df dictionary
                mean_values = df[key][0][flav][:, i]
                std_values = np.sqrt(df[key][1][flav][:, i])

                # Filter out NaN values
                valid_indices = ~np.isnan(mean_values)
                x_vals_loc = np.array(x_values)[valid_indices]
                mean_values = mean_values[valid_indices]
                std_values = std_values[valid_indices]

                # Create a TGraphErrors for each bin
                graph_name = f"fractions_{flav}_{jeteta_bins.idx2str(i)}"
                graph_title = f"jet energy fractions, flavor: {flav}; eta bin: {jeteta_bins.idx2plot_str(i)[1:-1]}"
                graph = ROOT.TGraphErrors(len(x_vals_loc), x_vals_loc, np.array(mean_values),
                                          np.zeros(len(x_vals_loc)), np.array(std_values))
                # Set graph properties (optional)
                graph.SetName(graph_name)
                graph.SetTitle(graph_title)

    #             print('name = ', graph_name)
    #             print('tit = ', graph_title)
                # Write the TGraphErrors to the ROOT file
                graph.Write()

    # Close the ROOT file
    root_file.Close()
    print(f"Saved {filename}")

def save_root(df, pt_bins, jeteta_bins, filename):
    x_values = pt_bins.centres

    # Create a ROOT file
    root_file = ROOT.TFile(filename, "RECREATE")

    # Loop over the keys in the df dictionary
    for key in df.keys():
        print("Storing the histogram:", key)
        # Create a directory for each key in the ROOT file
        root_dir = root_file.mkdir(key)
        root_dir.cd()

        # Loop over the flavors
        for flav in df[key].keys():

            # Create a TGraphErrors for each column (bin)
            for i in range(jeteta_bins.nbins):
                # Get the data from the df dictionary
                mean_values = np.array(df[key][flav]["Median"])[:, i]
                std_values = np.array(df[key][flav]["MedianStd"])[:, i]

                # Filter out NaN values
                valid_indices = ~np.isnan(mean_values) & np.logical_not(mean_values==0)
                x_vals_loc = np.array(x_values)[valid_indices]
                mean_values = mean_values[valid_indices]
                std_values = std_values[valid_indices]

                # Create a TGraphErrors for each bin
                graph_name = f"response_{flav}_{jeteta_bins.idx2str(i)}"
                graph_title = f"response, flavor: {flav}; eta bin: {jeteta_bins.idx2plot_str(i)[1:-1]}"
                # print("x_vals_loc: ", x_vals_loc)
                if len(x_vals_loc) > 0:
                    graph = ROOT.TGraphErrors(len(x_vals_loc), x_vals_loc, np.array(mean_values),
                                            np.zeros(len(x_vals_loc)), np.array(std_values))
                else:
                    print(f"The histogram for flavor {flav} and eta bin {jeteta_bins.idx2str(i)} is empty. Saving as a zero graph.")
                    graph = ROOT.TGraphErrors(len(x_values), x_values, np.zeros(len(x_values)),
                            np.zeros(len(x_values)), np.zeros(len(x_values)))
                    # graph = ROOT.TGraphErrors(2, np.array([0,100]), np.array([0,0]), np.array([0,0]), np.array([0,0]))
                # Set graph properties (optional)
                graph.SetName(graph_name)
                graph.SetTitle(graph_title)

    #             print('name = ', graph_name)
    #             print('tit = ', graph_title)
                # Write the TGraphErrors to the ROOT file
                graph.Write()

    # Close the ROOT file
    root_file.Close()
    print(f"Saved {filename}")