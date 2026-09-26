# Locations of the input weather data. Set the environment variables MADRID_DATA and TEXAS_DATA, or place the files in
# the data/ folder under the names below (see README.md for the sources and the column layout of each file).
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
MADRID_DATA = os.environ.get('MADRID_DATA', os.path.join(_HERE, 'data', 'todosLosDatos.txt'))
TEXAS_DATA = os.environ.get('TEXAS_DATA', os.path.join(_HERE, 'data', 'Texas_abreviatedData.txt'))
