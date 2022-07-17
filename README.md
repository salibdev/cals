# Cals

An open-source Python package used to calculate stellar chromospheric activity parameters and generate spectrum diagram of Ca II H&K lines with LAMOST (Large sky Area Multi-Object fiber Spectroscopic Telescope, also named Guoshoujing Telescope) LRS (Low-Resolution Spectroscopic Survey) spectra.

### Require the following Python libraries:
- astropy
- scipy
- numpy
- matplotlib

## Installation:
`python setup.py install`

## Brief example:

```
from cals import CaIISindex
data_path = "example.fits"
cs = CaIISindex(data_path)
cs.calcSindex()
cs.calcError()
cs.stellar_parameters = ['5534.22','4.423','-0.025']   # [teff, logg, feh]
cs.plotSpectrum()
cs.saveRecord()
```

## User's guide:
### Import package

`from cals import CaIISindex`

### Initialization

`data_path = "example.fits"`

`cs = CaIISindex(data_path)`

### Basic settings

#### Set csv file save path

```
cs.Sindex_savepath = "S_index_example.csv"
  # filled in with default path if not set explicitly
```

The csv file is for saving the calculated stellar chromospheric activity parameters and their error values. Default csv file save folder is current directory and default csv file name is `S_index_<primary name of the data file>.csv`. For instance, if the input data file name is `example.fits`, the default csv file name is `S_index_example.csv`.

#### Set figure file save path

```
cs.figure_savepath = "example.png"
  # filled in with default path if not set explicitly
```

The figure file is for saving the generated spectrum diagram. Default figure file save folder is current directory and default figure file name is `<primary name of the data file>.png`. For instance, if the input data file name is `example.fits`, the default figure file name is `example.png`.

#### Set stellar parameters

```
cs.stellar_parameters = ['5534.22', '4.423', '-0.025']
  # content of stellar parameter list: [teff (K), logg (dex), feh (dex)]
  # filled in with 'unknown' if not set explicitly
```

Because the stellar parameters 'teff', 'logg', and 'feh' are not included in LAMOST LRS fits file, if the three parameters are not provided via `cs.stellar_parameters`, they will be displayed as 'unknown' in spectrum diagram. Note that this setting is only used in spectrum diagram and does not affect the result of stellar chromospheric activity parameters calculated by this package.

### Functions

#### Calculate stellar chromospheric activity parameters

`cs.calcSindex()`

Return a dictionary of stellar activity parameters: `{'R_mean':R_mean, 'V_mean':V_mean, 'H_mean_tri':H_mean_tri, 'K_mean_tri':K_mean_tri, 'S_tri':S_tri, 'S_MWL':S_MWL, 'H_mean_rec':H_mean_rec,'K_mean_rec':K_mean_rec, 'S_rec':S_rec}`.

#### Calculate errors of stellar chromospheric activity parameters

`cs.calcError()`

Return a dictionary of activity parameter errors: `{'R_mean_err':R_mean_err, 'V_mean_err':V_mean_err, 'H_mean_tri_err':H_mean_tri_err, 'K_mean_tri_err':K_mean_tri_err, 'S_tri_err':S_tri_err, 'S_MWL_err':S_MWL_err, 'H_mean_rec_err':H_mean_rec_err, 'K_mean_rec_err':K_mean_rec_err, 'S_rec_err':S_rec_err}`.

#### Plot spectrum diagram

```
cs.plotSpectrum(figure_savepath=None,
                  stellar_parameters=None,
                  figshow=True)
```

#### Save results

```
cs.saveRecord(Sindex_savepath=None,
                figure_savepath=None,
                stellar_parameters=None,
                csv_header=True,
                figshow=False)
```

Save stellar chromospheric activity parameters (and their error values) to a csv file and spectrum diagram to a figure file.

## Terminal model

```
Start cals
usage: cals [-h] --cs CS [--savepath SAVEPATH] [--figsav {0,1}] [--figshow {0,1}]


optional arguments:
  -h, --help            show this help message and exit
  --cs CS, -c CS        Calculate the S by setting file path. | e.g. cals -c example.fits
  --savepath SAVEPATH, -sp SAVEPATH
                        Setting save path
  --figsav {0,1}, -fs {0,1}
                        Save/not save figure by setting figsav as 1/0, default=1.
  --figshow {0,1}, -fw {0,1}
                        Display/not display figure by setting figshow as 1/0, default=0.

```

#### terminal example

```
cd example/
cals -c ./example/example.fits -sp ./example/
```



