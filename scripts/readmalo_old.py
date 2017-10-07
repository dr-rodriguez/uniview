"""
Script to read the Malo 2013 table of bonafide members and output in the proper format
First version, no coordinate transformations
"""

from astropy.table import Table
import pandas as pd
from astropy.coordinates import SkyCoord, FK5
#from astropy.io.votable import parse
#import numpy as np
#from druvw import xyz, uvw

def readmalo_members(file):
    t = Table.read(file)
    #votable = parse(file)
    #t = votable.get_first_table()

    # Creating data frame
    #t.colnames # to get all column names
    cnames = ['X','Y','Z','YKG','Uvel','Vvel','Wvel','Name']
    data = pd.DataFrame
    for col in cnames:
        try:
            data.loc[:,col] = pd.Series(t[col]) #old parse format uses t.array[col]
        except TypeError: # This is for the first entry
            data = pd.DataFrame(t[col])
            data.columns = [col]

    # Rename columns
    cnames[3] = 'group'
    data.columns = cnames

    # Converting groups to numbers (for coloring purposes)
    gname = ['TWA','bPMG','THA','ABDMG','COL','ARG','CAR']
    ind = 0
    for elem in gname:
        data.loc[data.group==elem,'group'] = ind
        ind += 1

    # data.loc[data.group=='TWA','group'] = 0
    # data.loc[data.group=='bPMG','group'] = 1
    # data.loc[data.group=='THA','group'] = 2
    # data.loc[data.group=='ABDMG','group'] = 3
    # data.loc[data.group=='COL','group'] = 4
    # data.loc[data.group=='ARG','group'] = 5
    # data.loc[data.group=='CAR','group'] = 6

    return data


# Convert RA/Dec to XYZ
# This is accurate when compared to the member table
def convertxyz(file):
    t = Table.read(file)

    ra = t['_RAJ2000']
    dec = t['_DEJ2000']
    d = 1000./t['Plx']

    # This part is a bit slow
    x,y,z = [],[],[]
    for i in range(len(ra)):
        x0,y0,z0 = xyz(ra[i], dec[i], d[i])
        x.append(x0)
        y.append(y0)
        z.append(z0)

    return x,y,z

# Pandas version
def parsedata(data):
    nrow = len(data)

    datatxt = ''
    labeltxt = ''
    for i in range(nrow):
        temp = str(data.X[i]) + ' ' + str(data.Y[i]) + ' ' + str(data.Z[i]) + ' '
        datatxt += temp
        labeltxt += temp

        for column in data:
            if column in ('X','Y','Z'): continue
            if column == 'Name':
                datatxt += ' # ' + data[column][i]
                labeltxt += 'text ' + data[column][i]
            else:
                datatxt += str(data[column][i]) + ' '

        datatxt += '\n'
        labeltxt += '\n'

    return datatxt, labeltxt

def createheader():
    return 0

file = 'malo13_members.vot'
outname = 'nymg'

data = readmalo_members(file)
maintext,labeltext = parsedata(data)
print maintext
#print labeltext

with open(outname+'.speck','w') as f:
    f.write('datavar 0 groupid\n')
    f.write('datavar 1 u\n')
    f.write('datavar 2 v\n')
    f.write('datavar 3 w\n')
    f.write(maintext)

with open(outname+'.label','w') as f:
    f.write('textcolor 1\n')
    f.write(labeltext)

# Color map
txt = """7
0.0 1.0 0.0 # green
1.0 0.0 0.0 # red
0.0 0.0 1.0 # blue
1.0 1.0 1.0 # white
1.0 1.0 0.0 # yellow
0.0 1.0 1.0 # cyan
1.0 0.0 1.0 # magenta
"""
with open('nymg.cmap','w') as f:
    f.write(txt)