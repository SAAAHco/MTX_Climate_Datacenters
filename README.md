# Water and energy footprints of data centers in warm climates

Model and analyses of the article *Water and energy footprints of data centers in warm climates: thermodynamic modeling of cooling modes and climatic variability*, by Raúl Sánchez, Zainab Ashkanani, Reem Abdulrahman, and Rabi Mohtar (manuscript NEXUS-D-26-00499, *Energy Nexus*).

The code computes the cooling-related power usage effectiveness (PUE) and the consumptive evaporative water usage effectiveness (WUE) of a data center for every weather record. The cooling system has three modes:

- free cooling with a dry cooler,
- evaporative cooling with a cooling tower,
- chiller plus evaporative cooling.

The code applies the model to 14 years of half-hourly records at seven stations in Madrid (Spain) and 27 years of hourly records at four locations in Texas (USA). It then produces every number, table, and figure of the results.

## The model (`revmodel.py`)

- **Moist air.** Moist air is an ideal-gas mixture. The entropies of dry air and water vapor are evaluated at their partial pressures. The properties come from polynomial fits in `X = T[°C] + 273.16` (Table 2 of the article).
- **Heat load of the tower.** In the chilling mode, the tower rejects `P_ITE (1 + 1/COP)`, which includes the compressor work.
- **Tower state.** The operating state of the tower (and of the dry cooler) is closed with the counterflow limit of the Merkel theory plus a design approach `a`. The state is the minimum air flow with saturated outlet air, not warmer than `T1 - a`. The entropy generation is computed as a check.
- **Mode feasibility.** The conditions are analytical:
  - Free cooling is feasible when `T2 <= 18 - a_dry`.
  - Evaporative cooling is feasible when the wet-bulb temperature satisfies `Twb2 <= 18 - a_ev`.
  - Otherwise the chiller operates with the condenser-water set point `T4 = max(26, Twb2 + a_ev)` and `T1 = T4 + 24`.
- **Mode selection.** The mode minimizes `J = PUE + lambda * WUE`. With `lam = inf` (the default), the rule is free cooling whenever feasible, then evaporative cooling, then chilling. This rule is the minimum of `J` for every lambda between about 0.02 and 0.2 kWh/L; `sens.py` computes the exact interval (Section 3.2 of the article).
- **Pressure.** The atmospheric pressure of each location comes from its elevation (standard atmosphere, `stations.py`).
- **Fans and pumps.** Their power is `dp * V`, with 250 Pa for the air and 0.21 MPa for the water.

The design case uses `a_ev = 4 K` and `a_dry = 6 K`. The thermodynamic limit uses `a_ev = a_dry = 0`.

Optional switches of `revmodel.run` (keyword arguments) include:

- `cop_model='lift'`: the lift-based COP curve, Eq. (17).
- `dT_chiller`: the condenser range.
- `tower_lam`: a fan control that raises the air flow above the minimum when this lowers the fan energy plus `tower_lam` times the evaporation.
- `RH5`: the outlet relative humidity.
- `comp_heat=False`: the compressor work is not rejected in the tower.
- `lam`: the water value of the mode selection.

## Requirements

Python 3.11 or later. The article used Python 3.13.3 and the package versions of `requirements.txt`. Get the code and install the packages with

```
git clone https://github.com/SAAAHco/MTX_Climate_Datacenters.git
cd MTX_Climate_Datacenters
pip install -r requirements.txt
```

## Input data

The weather data are not included in this repository. They are publicly available from the sources below. Place the two files in `data/` under the names given here, or set the environment variables `MADRID_DATA` and `TEXAS_DATA` to their paths (see `config.py`).

**`data/todosLosDatos.txt`: Madrid.** Half-hourly records, January 1, 2006, to December 31, 2019, from seven stations of the SiAR network (Sistema de Información Agroclimática para el Regadío, Spanish Ministry of Agriculture, Fisheries and Food) in the province of Madrid:

| Station | Name |
|---|---|
| 1 | CENTER |
| 2 | Arganda |
| 3 | Aranjuez |
| 4 | Fuentidueña |
| 5 | San Martín |
| 6 | Chinchón |
| 102 | Villaprado |

The file is tab-separated and Latin-1 encoded. Its columns are `IdProvincia`, `IdEstacion`, `Fecha` (dd/mm/yyyy), `Año`, `Dia`, `HoraMin` (hhmm), `Temp Media (ºC)`, and `Hum Media (%)`, followed by wind, radiation, and precipitation, which are not used.

**`data/Texas_abreviatedData.txt`: Texas.** Hourly records, January 1, 1998, to December 31, 2024, of the National Solar Radiation Database (NSRDB, National Renewable Energy Laboratory) at the four locations of `stations.py` (El Paso, Lubbock, Dallas, and Houston). The file is tab-separated with the columns `Year`, `Month`, `Day`, `Hour`, `Minute`, `Temperature` (°C), `Relative Humidity` (%), and `Station` (`ElPaso`, `Lubbock`, `Dallas`, or `Houston`).

## Running the analyses

Run the scripts in this order from the repository folder, or run `python run_all.py`. The times were measured on a desktop computer; the complete run takes about 1 h 45 min.

| Step | Script | What it does | Main outputs | Time |
|---|---|---|---|---|
| 1 | `loaddata.py` | Reads and caches the two data files; reports the missing values | `madrid.pkl`, `texas.pkl` | a few seconds |
| 2 | `run_base.py` | Runs the design case and the thermodynamic limit for every record | `res_design.pkl`, `res_limit.pkl` | about 3 min |
| 3 | `analysis.py` | Descriptive statistics, nested variance decomposition, averaging convergence, mode occupancy, climate transfer (kernel densities, Hellinger distance), annual and monthly demands of a 100 MW facility, peaks, ET0 | `analysis_out.json`, `transfer_maps.pkl` | under 1 min |
| 4 | `sens.py` | One-at-a-time cases on the complete records and the Pareto front over the water value lambda, with the lambda interval of the selection rule | `sens_out.json` | about 51 min |
| 5 | `gsa.py` | Sobol indices (Saltelli sampling, N = 1024) and 95 % intervals from 4,096 Latin hypercube samples, on climate histograms of each location | `gsa_out.json`, `gsa_base.json`, `gsa_*.npy` | about 50 min |
| 6 | `extras.py` | Data completeness, drift share, Merkel numbers and water-to-air ratios, air-flow ratios, marginal water value of extra tower air, Carnot fractions, transfer errors of all donor-target pairs, ET0 totals | `extras_out.json` | under 1 min |
| 7 | `figs.py` | Figures 3 to 12 of the article | `figs/*.png`, `figs/eps/*.eps` | under 1 min |
| 8 | `fig2.py` | Figure 2, the schematic of the cooling circuits | `figs/fig02_schematic.png`, `figs/eps/fig02_schematic.eps` | seconds |
| 9 | `propcheck.py` | Checks the property fits and psychrometric relations against IAPWS-95, the reference equation of state of dry air (CoolProp), and ASHRAE (psychrolib) | printed report | under 1 min |

`resq.py` reads the result files for `figs.py`. `numbering.json` holds the equation numbers of the article, which the flowchart of Figure 3 cites.

## Reproducibility

The folder `results_article/` holds the result files used in the article. A complete run of steps 1 to 9 from a clean copy of this repository, with the package versions of `requirements.txt`, reproduced all five files exactly (`analysis_out.json`, `sens_out.json`, `gsa_base.json`, `gsa_out.json`, and `extras_out.json`, with a largest relative difference of zero). All random elements use fixed seeds: the Latin hypercube samples, the bootstrap confidence intervals of the Sobol indices, and the Saltelli sequence, which is deterministic. Other versions of NumPy, SciPy, or SALib can change the last digits.

## Contents

| File | Purpose |
|---|---|
| `revmodel.py` | The cooling-system model |
| `stations.py` | Coordinates, elevations, and standard-atmosphere pressure of the locations |
| `config.py` | Paths of the input data |
| `loaddata.py`, `run_base.py`, `analysis.py`, `sens.py`, `gsa.py`, `extras.py` | Analyses (steps 1 to 6) |
| `figs.py`, `fig2.py` | Figures |
| `propcheck.py` | Property checks |
| `resq.py` | Helper that reads the results |
| `run_all.py` | Runs all the steps in order |
| `results_article/` | Result files of the article |
| `numbering.json` | Equation, figure, and table numbers of the article |

## Contact Information
For questions email: Ashkanani@tamu.edu, Ashkananai@saaah.co
