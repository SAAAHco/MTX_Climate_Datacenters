# Input data

Place the two weather files here, or set the environment variables `MADRID_DATA` and `TEXAS_DATA` to their paths. They are not included in the repository, because of their size and because they belong to their providers. Both are publicly available.

| File | Region | Source | Period | Time step | Size |
|---|---|---|---|---|---|
| `todosLosDatos.txt` | Madrid, seven stations | SiAR network, Spanish Ministry of Agriculture, Fisheries and Food | 2006 to 2019 | 30 min | about 106 MB |
| `Texas_abreviatedData.txt` | El Paso, Lubbock, Dallas, Houston | National Solar Radiation Database (NSRDB), National Renewable Energy Laboratory | 1998 to 2024 | 1 h | about 33 MB |

## `todosLosDatos.txt`

This is a tab-separated file in Latin-1 encoding with one header row. Its first eight columns are listed below; the scripts use the second to the eighth, and the remaining columns (wind, radiation, and precipitation) are not used.

| Column | Content |
|---|---|
| `IdProvincia` | Province code (28, Madrid) |
| `IdEstacion` | Station code (1, 2, 3, 4, 5, 6, or 102) |
| `Fecha` | Date, dd/mm/yyyy |
| `Año` | Year |
| `Dia` | Day of the year |
| `HoraMin` | Time, hhmm |
| `Temp Media (ºC)` | Mean air temperature, °C |
| `Hum Media (%)` | Mean relative humidity, % |

## `Texas_abreviatedData.txt`

This is a tab-separated file with one header row and the columns `Year`, `Month`, `Day`, `Hour`, `Minute`, `Temperature` (°C), `Relative Humidity` (%), and `Station` (`ElPaso`, `Lubbock`, `Dallas`, or `Houston`). Each location has 236,687 hourly records, time-stamped at the half hour. The coordinates of the NSRDB data points are listed in `stations.py`.
