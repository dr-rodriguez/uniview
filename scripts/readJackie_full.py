"""
Script to read in the full table from J. Faherty's paper and process it for input to Uniview/Partiview
"""

import pandas as pd
from uniview_functions import parsedata, transform

# Function to change group numbers. Does so in place
def group_numbers(data, column):
    gname = ['TWA','betaPic','TucHor','ABDor','Columba','Argus','CAR']
    ind = 0
    for elem in gname:
        data.loc[data[column]==elem, column] = ind
        ind += 1
    return

filename = 'input_tables/jackie_full_v2.csv'

t = pd.read_csv(filename, header=0)  # header = 0 to grab column names

# Re-order columns
cols = t.columns.tolist()
cols.remove('Name')
cols.remove('Young')
cols.remove('SpT-Either')
cols.remove('Group')
cols = cols + ['Group', 'Name', 'Young']
t = t[cols]

# Process to calculate XYZ
x, y, z = transform(t['RA'], t['Dec'], t['Distance'])

t['X'] = pd.Series(x)
t['Y'] = pd.Series(y)
t['Z'] = pd.Series(z)

# Replace Group with an ID number (same as NYMG)
group_numbers(t, 'Group')
# Replace NaN with 99
t.loc[t['Group'].isnull(), 'Group'] = 7

# Clean up table with relevant columns
t = t.drop(['RA', 'Dec'], axis=1)

# Generate text versions
t2 = t.copy().drop('Young', axis=1)
datatxt, labeltxt = parsedata(t2)

outname = 'bdkp/tab_full'

# Generate datavar IDs for the extra columns
cols = t.columns.tolist()
map(lambda x: cols.remove(x), ['X','Y','Z','Name','Young'])

with open(outname+'.speck','w') as f:
    for i in range(len(cols)):
        txt = 'datavar '+str(i)+' '+cols[i]+'\n'
        f.write(txt)
    f.write(datatxt)

with open(outname+'.label','w') as f:
    f.write('textcolor 1\n')
    f.write(labeltxt)

# ===============================================================
# Repeat, but eliminate the Young stars
ind = t['Young'] == 0
t_old = t[ind].copy()

# Remove Young column
t_old.drop(['Young', 'Group'], axis=1, inplace=True)

# Generate text versions
datatxt, labeltxt = parsedata(t_old)

outname = 'bdkp/tab_noyoung'

# Generate datavar IDs for the extra columns
cols = t_old.columns.tolist()
map(lambda x: cols.remove(x), ['X','Y','Z','Name'])

with open(outname+'.speck','w') as f:
    for i in range(len(cols)):
        txt = 'datavar '+str(i)+' '+cols[i]+'\n'
        f.write(txt)
    f.write(datatxt)

with open(outname+'.label','w') as f:
    f.write('textcolor 1\n')
    f.write(labeltxt)