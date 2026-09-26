# Station metadata. Elevations (m a.s.l.) are used for the standard-atmosphere pressure (Reviewer 7.10).
import numpy as np
import revmodel as R

MADRID = [1, 2, 3, 4, 5, 6, 102]
TEXAS = ['ElPaso', 'Lubbock', 'Dallas', 'Houston']
LABEL = {1: 'M01', 2: 'M02', 3: 'M03', 4: 'M04', 5: 'M05', 6: 'M06', 102: 'M102',
         'ElPaso': 'El Paso', 'Lubbock': 'Lubbock', 'Dallas': 'Dallas', 'Houston': 'Houston'}
NAME = {1: 'CENTER', 2: 'Arganda', 3: 'Aranjuez', 4: 'Fuentidueña', 5: 'San Martín', 6: 'Chinchón', 102: 'Villaprado'}
LAT = {1: 40.4120, 2: 40.3111, 3: 40.0416, 4: 40.1071, 5: 40.2333, 6: 40.1922, 102: 40.2510,
       'ElPaso': 31.7791, 'Lubbock': 33.5281, 'Dallas': 32.8384, 'Houston': 29.4718}
LON = {1: -3.4966, 2: -3.4979, 3: -3.6304, 4: -3.1741, 5: -3.5599, 6: -3.4687, 102: -4.2730,
       'ElPaso': -106.3112, 'Lubbock': -101.8761, 'Dallas': -96.8358, 'Houston': -95.0832}
ELEV = {1: 553.0, 2: 532.0, 3: 486.0, 4: 552.0, 5: 516.0, 6: 532.0, 102: 466.0,
        'ElPaso': 1214.0, 'Lubbock': 984.0, 'Dallas': 147.0, 'Houston': 7.0}
# Madrid: official SiAR (MAPA) altitudes and positions; Texas: USGS 3DEP ground elevation at the NSRDB data points.


def pressure(st):
    return float(R.p_std(ELEV[st]))
