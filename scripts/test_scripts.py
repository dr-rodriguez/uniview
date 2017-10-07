"""
HIP7:
 -12.2689    40.5845   -37.1473     0.540      0.34363      5.88      9.64     1     183.77
 1     17.740      7.33     60.567     30.129    -47.705     82.776       7 # HIP7 Gli
"""

import pandas as pd
from druvw import xyz, uvw

# Solar velocity, to be added to the UVW to convert to LSR frame
lsr_vel=(-8.5,13.38,6.49) # Coskunoglu et al. 2011 MNRAS
lsr_vwl[0] = lsr_vel[0]*-1

ra = 0.02254892
dec = 20.03660198
plx = 17.28
pmra = -208.26
pmdec = -200.65
d = 57.56
rv = 8.3

# Data from Extended Hipparcos Compilation (XHIP) (Anderson+, 2012)
#XYZ -12.5	41.4	-37.9
#UVW 71.7	2.1	-34.0

# These agree with the catalog, as they should. But they don't agree with stars.speck
print xyz(ra=ra, dec=dec, d=d)
print uvw(ra=ra, dec=dec, d=d, pmra=pmra, pmde=pmdec, rv=rv)

# Using the stars.speck parallax
d = 1000./17.74
print xyz(ra=ra, dec=dec, d=d) # this now agrees with the stars.speck file

#VX, VY, VZ, Speed: 60.567     30.129    -47.705     82.776
u,v,w = uvw(ra=ra, dec=dec, d=d, pmra=pmra, pmde=pmdec, rv=rv)
print u,v,w
print u+lsr_vel[0], v+lsr_vel[1], w+lsr_vel[2]


# Extra TWA members
df = pd.read_table('input_tables/extratwa.txt', sep=' ', header=0)
x,y,z = xyz(df['ra'], df['dec'], df['dist'])
new_df = pd.DataFrame({'X':x, 'Y':y, 'Z':z, 'name':df['name']})

# Test stars
df = pd.read_table('input_tables/testobjects.txt', sep=' ', header=0)
d = 1000./df['pi']
x,y,z = xyz(df['ra'], df['dec'], d)
new_df = pd.DataFrame({'X':x, 'Y':y, 'Z':z, 'Name':df['name']})
from uniview_functions import parsedata
maintext,labeltext = parsedata(new_df)
outname = 'teststars/test'
with open(outname+'.speck','w') as f:
    f.write(maintext)

with open(outname+'.label','w') as f:
    f.write('textcolor 1\n')
    f.write(labeltext)