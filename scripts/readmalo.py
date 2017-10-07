"""
Script to read the Malo 2013 table of bonafide members *and* candidates
"""

from astropy.table import Table
import pandas as pd
from astropy.coordinates import SkyCoord, FK5
from math import cos, sin

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

# ===================================================================================
# Read in the members table
def read_members(file):
    t = Table.read(file)

    # Grabbing RA, Dec, Distance
    ra = t['_RAJ2000']
    dec = t['_DEJ2000']
    dist = 1000./t['Plx']

    xgc, ygc, zgc = transform(ra, dec, dist)

    # Manually set the X Y Z columns
    data = {'X':pd.Series(xgc), 'Y':pd.Series(ygc), 'Z':pd.Series(zgc)}
    data = pd.DataFrame(data)

    # Import additional columns of interest
    cnames = ['YKG','Name']
    for col in cnames:
        data.loc[:,col] = pd.Series(t[col])

    # Converting groups to numbers (for coloring purposes)
    data = group_numbers(data, 'YKG')

    return data

# ===================================================================================
# Read in the candidates table
def read_candidates(file):
    t = Table.read(file)

    # Grabbing RA, Dec
    ra = t['_RAJ2000']
    dec = t['_DEJ2000']

    # Grabbing distance
    d_plx = t['Dp']
    plx_mask = t['Dp'].mask
    d_kin = t['Ds']

    # Remove ambiguous and previous members (lot of missing info on them)
    group = t['YKG']
    mask = ~((group=='Amb') | (group=='Prev'))

    ra, dec = ra[mask], dec[mask]
    d_kin, d_plx = d_kin[mask], d_plx[mask]
    plx_mask = plx_mask[mask]

    # Set the distance to be parallax if present, otherwise use kinematic (statistical)
    dist = d_plx
    dist[plx_mask] = d_kin[plx_mask]

    xgc, ygc, zgc = transform(ra, dec, dist)

    # Manually set the X Y Z columns
    data = {'X':pd.Series(xgc), 'Y':pd.Series(ygc), 'Z':pd.Series(zgc)}
    data = pd.DataFrame(data)

    # Import additional columns of interest
    cnames = ['YKG','_2MASS']
    for col in cnames:
        data.loc[:,col] = pd.Series(t[col][mask])

    # Converting groups to numbers (for coloring purposes)
    data = group_numbers(data, 'YKG')

    return data

# ===================================================================================
# Set groups defined in *column* as numbers
def group_numbers(data, column):
    gname = ['TWA','bPMG','THA','ABDMG','COL','ARG','CAR']
    ind = 0
    for elem in gname:
        data.loc[data[column]==elem, column] = ind
        ind += 1

    return data

# ===================================================================================
# Prepare data for output
def parsedata(data):
    datatxt = ''
    labeltxt = ''
    for i in range(len(data)):
        temp = str(data.X[i]) + ' ' + str(data.Y[i]) + ' ' + str(data.Z[i]) + ' '
        datatxt += temp
        labeltxt += temp

        for column in data:
            if column in ('X','Y','Z'): continue
            if (column == 'Name') | (column == '_2MASS'):
                datatxt += ' # ' + data[column][i]
                labeltxt += 'text ' + data[column][i]
            else:
                datatxt += str(data[column][i]) + ' '

        datatxt += '\n'
        labeltxt += '\n'

    return datatxt, labeltxt

# ===================================================================================
# Generate the colormap
def gen_colormap(filename):
    # Color map
    txt = """7
0.0 1.0 0.0 # green
1.0 0.0 0.0 # red
0.0 0.0 1.0 # blue
1.0 1.0 1.0 # white
1.0 1.0 0.0 # yellow
0.0 1.0 1.0 # cyan
1.0 0.0 1.0 # magenta"""

    with open(filename,'w') as f:
        f.write(txt)

# ===================================================================================
# Generate files for each group separately
def files_by_group(data):
    gname = ['TWA','bPMG','THA','ABDMG','COL','ARG','CAR']
    filenames = map(lambda x: str.lower(x), gname)

    for i in range(len(gname)):
        outname = 'nymg/'+filenames[i]

        newdata = data.loc[data['YKG']==i]
        newdata.index = range(len(newdata))
        maintext,labeltext = parsedata(newdata)

        with open(outname+'.speck','w') as f:
            f.write('datavar 0 groupid\n')
            f.write(maintext)

        with open(outname+'.label','w') as f:
            f.write('textcolor 1\n')
            f.write(labeltext)

# ===================================================================================

file = 'input_tables/malo13_members.vot'
outname = 'nymg/nymg_all'

data = read_members(file)
maintext,labeltext = parsedata(data)
print maintext

file = 'input_tables/malo13_candidates.vot'
data2 = read_candidates(file)
maintext2,labeltext2 = parsedata(data2)

with open(outname+'.speck','w') as f:
    f.write('datavar 0 groupid\n')
    f.write(maintext)
    f.write(maintext2)

with open(outname+'.label','w') as f:
    f.write('textcolor 1\n')
    f.write(labeltext)
    f.write(labeltext2)

# Create colormap
gen_colormap('nymg/nymg.cmap')

# Merge two data sets and produce individual files
data2.rename(columns = {'_2MASS':'Name'}, inplace = True)
data3 = data.append(data2)
files_by_group(data3)
