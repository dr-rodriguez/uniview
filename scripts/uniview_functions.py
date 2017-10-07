"""
File to store some of the generic functions for use elsewhere
"""

#from astropy.table import Table
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from math import cos, sin, isnan

# ===================================================================================
# Prepare data for output
def parsedata(data):
    datatxt = ''
    labeltxt = ''

    # Re-index to avoid errors
    data = data.reindex(range(len(data)))

    for i in range(len(data)):
        if isnan(data.X[i]):
            continue

        temp = str(data.X[i]) + ' ' + str(data.Y[i]) + ' ' + str(data.Z[i]) + ' '
        datatxt += temp
        labeltxt += temp

        for column in data:
            if column in ('X','Y','Z'): continue
            if (column == 'Name') | (column == '_2MASS'):
                datatxt += ' # ' + data[column][i] + ' '
                labeltxt += 'text ' + data[column][i]
            else:
                datatxt += str(data[column][i]) + ' '

        datatxt += '\n'
        labeltxt += '\n'

    return datatxt, labeltxt

# ===================================================================================
# Transform coordinates
def transform(ra, dec, dist):
    # Precess
    c = SkyCoord(ra=ra, dec=dec, frame='icrs', equinox='J2000', unit='deg')
    #c = c.transform_to(FK5(equinox='B1950'))

    # Convert to XYZ
    l,b = c.galactic.l.radian,c.galactic.b.radian
    xgc = dist * map(cos,b) * map(cos,l)
    ygc = dist * map(cos,b) * map(sin,l)
    zgc = dist * map(sin,b)

    return xgc, ygc, zgc

# ====
# Get cmap colors
def rgb(x1,x2,x3):
    return list(np.array([x1,x2,x3])/255.)