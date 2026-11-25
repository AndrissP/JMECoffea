# JMECoffea

This code contains a columnar-based calculation of MC jet energy corrections (JEC) and flavor uncertainties, following the perscriptions of [arxiv:1607.03663](https://arxiv.org/pdf/1607.03663). At the moment, the code (L5_flavour_dep_corr branch) is adapted for the calculation of L5 flavor JEC and L2/L3 correction part is not tested.

The main structure of the code:
- `run_processor.py` runs the suplied coffea processor (e.g., `CoffeaJERCProcessor_L5.py`) and saves the results in .coffea files.
- `fit_response_distributions.py` fits the response distributions and stores the results in .txt/.json files.
- `correction_fitter.ipynb` or `correction_fitter.py` reads the response fit results, fits them vs the reco pt and stores JEC as `.txt` files.
- `flavour_fractions.py` reads the response histograms and produces flavor fractions and their plots.
- `flavour_uncertainties.py` combines the flavor fractions and response corrections to calculate flavor uncertainties. The original `flavor_fractions_and_uncertainties.ipynb` is retained as a reference.
- `plotters/Plotting_comparison.ipynb` allows to plot different responses and corrections vs pt or jet_eta.


## Set-up
Several set-up options exist depending if you are running on [lxplus](#set-up-coffea-on-the-lxplus-with-the-lcg-nightlies-environment-and-dask-on-htcondor), [lpc](#instructions-for-the-lpc-and-singularity) or [coffea casa](#instructions-for-coffea-casa). Running on any lxpus or lpc using singularity is possible but not stated here.
### Set up coffea on the **lxplus** with the lcg environment and running using dask on HTCondor:

Load the `lcg` environment with the relevant packages including `coffea`
```
source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh
```

Follow the instruction of [installing](#install-the-repo) the repo.

### Instructions for the **LPC** (and singularity)

1. Inside the JMECoffea directory:
```
wget https://raw.githubusercontent.com/CoffeaTeam/lpcjobqueue/main/bootstrap.sh
bash bootstrap.sh
./shell
```
2. Install the repo

### Instructions for **Coffea Casa**

1. Log into https://coffea.casa/hub/login?next=%2Fhub%2F with your CERN account
2. Go to the "Git" menu and select "Clone a repository" 

## Install the repo

Clone the repo
```
git clone git@github.com:AndrissP/JMECoffea.git
```

Install the package in the editable state in a virtual environment. Explanation of this is in the [stackoverflow page](https://stackoverflow.com/questions/714063/importing-modules-from-parent-folder/50194143#50194143):
```
python -m venv venv
pip install -e .
```

Every time logging in, activate the virtual environment and after that load the lcg
```
. venv/bin/activate
source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh
```

Use coffea update with a fix for the grid crashes (unless coffea with the fix is used)
```
cd ..
git clone git@github.com:AndrissP/coffea.git
pip install -e .
realpath
### load the path before running `import coffea` in python
```

## ssh port forwarding in to use `jupyter notebooks`:
To use jupyter notebooks you have to log into lxplus using port forwarding
```
ssh -L <local_port>:localhost:<remote_port> <username>@lxplus.cern.ch
jupyter notebook --no-browser --port=<remote_port>
```
e.g
```
ssh -L 8099:localhost:8095 <username>@lxplus.cern.ch
source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh
jupyter notebook --no-browser --port=8095
```
Then copy the link in the output and connect in to your broswer, replacing the remote_port valur with local_port. Here you are mapping the port local_port (e.g 8099) of the localhost (your machine) to the port remote_port (e.g. 8095) of the remote server (lxplus.cern.ch). We assume that the two ports are free and therefore available.

## Running the code:

### Running the L5 (or L2L3Rel) processors
Run the histogram creation
```
python run_processor.py -d QCD-Py-Run3-2024-Summer24 -e iterative
```
By default it runs the `CoffeaJERCProcessor_L5.py` processor.
The results are stored in .coffea files in `out/`.
Run `python run_processor.py --help` for help, run `python run_processor.py --test-run` to test.

### Fitting and plotting the response distributions
Run the fits with the default dataset:
```
python fit_response_distributions.py --data QCD-Her-Run3-2024-Summer24
python fit_response_distributions.py --data QCD-Py-Run3-2024-Summer24,QCD-Her-Run3-2024-Summer24
```

Run `python run_processor.py --help` for help.

Fit results are saved in `out_txt/` by default as a combined JSON file ; select CSV with `--fit-results-format csv` to save the per-flavor, per-result CSV files. Response plots are saved in `fig/` by default.

### Fitting respponses and creating output txt files
```
python correction_fitter.py
```
Or follow `correction_fitter.ipynb`.

### Fitting the flavor uncertainties
Calculate the flavor fractions (defaults match the Run 2 notebook setup):
```
python flavour_fractions.py
```
The script saves a JSON cache under `out_txt/`. It can optionally save plots and ROOT output; use `--help` to see the available settings.

Calculate the flavor uncertainties:
```
python flavour_uncertainties.py
```
The uncertainty script reuses a fraction cache only when its input files and configuration match; otherwise it recomputes the fractions. Pass `--recompute-fractions` to force a fresh calculation. By default it saves uncertainty results as JSON and plots, without changing any correction or uncertainty text payload. To create updated payload copies, provide source files with `--uncertainty-files`; output files are written alongside those sources with a `_run2flavor` suffix and are not overwritten unless `--overwrite-output` is specified. Use `--help` for the full command-line options.
