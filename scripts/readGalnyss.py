"""
Script to read in the GALNYSS table and parse to Uniview
"""

from uniview_functions import transform, parsedata
import pandas as pd
from astropy.table import Table


filename = 'input_tables/galnyss_basic_v2.vot'
t = Table.read(filename, format='votable')
t = t.to_pandas()

df = t.loc[:, ['sname','ra','dec','d10']]

xgc, ygc, zgc = transform(df['ra'], df['dec'], df['d10'])

data = {'Name':df['sname'], 'X':xgc, 'Y':ygc, 'Z':zgc}
galnyss = pd.DataFrame(data)

datatxt, labeltxt = parsedata(galnyss)

outname = 'galnyss/galnyss'

with open(outname+'.speck','w') as f:
    #f.write('datavar 0 groupid\n')
    f.write(datatxt)

with open(outname+'.label','w') as f:
    f.write('textcolor 1\n')
    f.write(labeltxt)
