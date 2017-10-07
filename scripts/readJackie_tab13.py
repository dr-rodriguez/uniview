"""
Script to read in Tables from J. Faherty's paper and process them for input to Uniview/Partiview
"""

import pandas as pd
from uniview_functions import parsedata
from druvw import xyz, uvw

filename = 'input_tables/jackie_tab13.csv'
#colnames = ['Name','U','V','W','X','Y','Z']

t = pd.read_csv(filename, header=0) # header = 0 to grab column names

# Remove the UVW columns
#t.drop(t.columns[[1,2,3]], axis=1, inplace=True)

# Get distances
t['distance'] = 1000./t['parallax']

# Re-order to have name and group last
cols = t.columns.tolist()
cols.remove('Name')
cols.remove('Group')
cols = cols + ['Name','Group']
t = t[cols]

# Process to calculate XYZ
ind = pd.isnull(t['X']) & pd.notnull(t['distance'])
x,y,z = xyz(t['RA'][ind], t['Dec'][ind], t['distance'][ind])
loc = t[ind].index.tolist()

t['X'][loc] = pd.Series(x)
t['Y'][loc] = pd.Series(y)
t['Z'][loc] = pd.Series(z)

# Now for UVW
ind = pd.isnull(t['U']) & pd.notnull(t['distance']) & pd.notnull(t['RV'])
u,v,w = uvw(t['RA'][ind],t['Dec'][ind],t['distance'][ind],t['mura'][ind]*1000,t['mudec'][ind]*1000,t['RV'][ind])
loc = t[ind].index.tolist()

t['U'][loc] = pd.Series(u)
t['V'][loc] = pd.Series(v)
t['W'][loc] = pd.Series(w)

# Clean up table with relevant columns
cols = t.columns.tolist()
map(lambda x: cols.remove(x), ['RA','Dec','parallax'])
t = t[cols]

# Generate text versions
datatxt, labeltxt = parsedata(t)

outname = 'bdkp/tab13'

# Generate datavar IDs for the extra columns
cols = t.columns.tolist()
map(lambda x: cols.remove(x), ['X','Y','Z','Name','Group'])

with open(outname+'.speck','w') as f:
    for i in range(len(cols)):
        txt = 'datavar '+str(i)+' '+cols[i]+'\n'
        f.write(txt)
    f.write(datatxt)

with open(outname+'.label','w') as f:
    f.write('textcolor 1\n')
    f.write(labeltxt)
