# EWAVE

EWAVE is a one-dimensional, time-averaged energy balance model implemented in MATLAB. It predicts the cross-shore transformation of short waves and infragravity (IG) waves, including incoming and outgoing IG wave heights and wave setup.

The model represents wave breaking, bottom friction, surface rollers, and wave reflection. It uses wave-induced water-mass exchange to parameterize low-frequency energy transfer from incoming short waves and their rollers. The transferred energy is partitioned between incoming IG waves and wave setup.

## Repository layout

| Directory | Contents |
| --- | --- |
| [`src/`](src/) | The main solver, `EWAVE_functions_2020d.m`, and its supporting MATLAB functions. |
| [`cases/`](cases/) | Validation scripts and their bathymetry, boundary conditions, and observations. |

## Run a validation case

1. Clone or download the repository. Keep `src/` and `cases/` at the same level.
2. In MATLAB, set the Current Folder to the desired case directory under `cases/`. Keep the entry script and its three `.inp` files in that directory.
3. Run the command in the table below. The entry scripts locate `src/` from their own location. For the Xisha Bay case, keep `plot_EWAVE_timeseries.m` alongside the entry script.

| Validation case | Site or experiment | Default case ID | MATLAB command | Reference |
| --- | --- | --- | --- | --- |
| GLOBEX | Gently sloping beach in a laboratory flume | `A3` | `main_function` | [Ruessink et al. (2013)](https://hydralab.eu/uploads/TAdocuments/170_Ruessink_Gerben.pdf) |
| Boers | Barred beach in a laboratory flume | `1c` | `main_function` | [Boers (1996)](https://resolver.tudelft.nl/uuid:22dfbc0a-adf8-44cb-9659-9d0a366c88d6) |
| DELILAH | Barred beach at Duck, North Carolina, USA | `delilah2` | `main_function` | [Van Dongeren et al. (2003)](https://doi.org/10.1029/2002JC001332) |
| Agate Beach | Low-slope beach in Oregon, USA | `case9` | `main_function` | [Fiedler et al. (2015)](https://doi.org/10.1002/2015GL066124) |
| Xisha Bay | Macrotidal low-tide terrace in Fujian, China | `xsw` | `main_function_YL_timeseries` | [Li et al. (2026)](https://doi.org/10.1016/j.oceaneng.2025.123998) |

## Case inputs and outputs

The `casename` variable in the entry script sets the case ID. To run another input set in the same profile case directory, set `casename` to the ID shared by its three input files.

| Input file | Contents |
| --- | --- |
| `be_<id>.inp` | Cross-shore coordinates and bed elevations in metres, starting at the offshore boundary. |
| `bc_<id>.inp` | Boundary forcing records, including still water level and incident wave conditions. |
| `obs_<id>.inp` | Observations at the validation gauges. |

The entry script defines the boundary and observation columns in `input.index`. Wave heights are expressed as root-mean-square (RMS) values internally. The `spectype` flag is `1` for RMS input and `2` for significant input.

The profile validation scripts use one boundary record and save `validation_<id>_paired.mat` in the case directory. The file contains paired observations and predictions at the gauges for total short-wave height, total IG wave height, and incoming and outgoing IG wave heights. The scripts report RMSE and signed relative bias.

The Xisha Bay script runs EWAVE for each boundary record and pairs the predictions with observations by time and gauge position. Use the same time units in the boundary and observation files. The time-series validation, selected cross-shore IG wave profile, and comparison at the offshore boundary gauge (G0) are saved under `validation_xsw_timeseries/`.

## Citation

When using EWAVE, cite the repository and report the release tag or commit used. The companion manuscript is *Parametric Modelling of Nearshore Infragravity Waves* by Yuan Li, Chi Zhang, and Jinhai Zheng. Cite the corresponding references in the case table when using the validation data.

## License

EWAVE source code and original case scripts are provided under the MIT License. See [LICENSE](LICENSE). The validation datasets and third-party materials are subject to the terms of their original providers.
