# Script to read the Tycho-2 data, match it against GAIA, and output the results to a speck file

import os
import pandas as pd
from uniview import Uniview
from astroML.crossmatch import crossmatch_angular
import numpy as np

# Load the data
path = '/Users/strakul/data/gaia/tycho2'
files = os.listdir(path)

cols = ['tyc', 'pflag', 'ra', 'dec', 'pmRA', 'pmDE', 'e_RAmdeg', 'e_DEmdeg', 'e_pmRA', 'e_pmDE',
        'EpRAm', 'EpDEm', 'Num', 'q_RAmdeg', 'q_DEmdeg', 'q_pmRA', 'q_pmDE', 'BTmag', 'e_BTmag', 'VTmag',
        'e_VTmag', 'prox', 'TYC', 'HIP-CCDM', 'RAdeg', 'DEdeg', 'EpRA-1990', 'EpDE-1990', 'e_RAdeg',
        'e_DEdeg', 'posflg', 'corr']

# file = 'tyc2.dat.00'
# df = pd.read_table(os.path.join(path, file), sep='|', header=None, names=cols)

# All data
count = 0
for file in files:
    if not file.startswith('tyc2.dat'):
        continue
    filename = os.path.join(path, file)
    temp = pd.read_table(filename, sep='|', header=None, names=cols)

    if count == 0:
        df = temp.copy()
    else:
        df = pd.concat([df, temp], axis=0, ignore_index=True)

    count += 1
    # if count > 3: break  # do only part

# Read GAIA speck file
uni = Uniview()
df_gaia = uni.load_speck('gaia_dr1/gaia_dr1.speck')

# Make certain columns numeric
df['ra'] = pd.to_numeric(df['ra'], errors='coerce')  # mean ra/dec
df['dec'] = pd.to_numeric(df['dec'], errors='coerce')
df['RAdeg'] = pd.to_numeric(df['RAdeg'], errors='coerce')  # actual ra/dec
df['DEdeg'] = pd.to_numeric(df['DEdeg'], errors='coerce')
df['BTmag'] = pd.to_numeric(df['BTmag'], errors='coerce')
df['VTmag'] = pd.to_numeric(df['VTmag'], errors='coerce')

# Fix missing ra/dec. Otherwise, 21614 are missing
bad = df['ra'].isnull()
df['ra'][bad] = df['RAdeg'][bad]
df['dec'][bad] = df['DEdeg'][bad]

for col in ['ra', 'dec', 'pmra', 'pmdec']:
    if col in ['comment', 'texture']:
        continue
    df_gaia[col] = pd.to_numeric(df_gaia[col], errors='coerce')

# Matches with astroML
T = np.empty((len(df), 2), dtype=np.float64)
T[:, 0] = df['ra']
T[:, 1] = df['dec']

ra_2000 = df_gaia['ra'] - ((df_gaia['pmra']*np.cos(df_gaia['dec']*np.pi/180))*1/1000.*1/3600.)*(2015.-2000.)
dec_2000 = df_gaia['dec'] - (df_gaia['pmdec']*1/1000.*1/3600.)*(2015.-2000.)
G = np.empty((len(df_gaia), 2), dtype=np.float64)
G[:, 0] = ra_2000
G[:, 1] = dec_2000

# Match takes about 16 minutes to complete
max_radius = 5. / 3600  # 5 arcsec
import datetime
print datetime.datetime.now().time()
dist, ind = crossmatch_angular(G, T, max_radius)  # ~2 mins for matching 1 mil vs 250,000, ~3 for vs 500,000
print datetime.datetime.now().time()
match = ~np.isinf(dist)

df_gaia['separation'] = dist * 3600
df_gaia.loc[match, 'match'] = ind[match]

col_subset = ['tyc', 'ra', 'dec', 'BTmag', 'VTmag']
gaia = pd.merge(df_gaia, df[col_subset], how='left', left_on='match', left_index=False, right_index=True, sort=False,
                    suffixes=('', '_t'), copy=True, indicator=False)

# Compute color
gaia['Tcolor'] = gaia['BTmag'] - gaia['VTmag']
# 744 with either no Tycho-2 match or not enough photometry info

# Drop/rename columns
gaia.drop(['parallax_error', 'absmag', 'distance', 'match', 'ra_t', 'dec_t',
           'separation', 'texture'], axis=1, inplace=True)
gaia.columns = ['Gmag' if x == 'phot_g_mean_mag' else x for x in gaia.columns]

# Replace missing info with 0
gaia.fillna(value={'Tcolor': 0, 'BTmag': 999, 'VTmag': 999}, inplace=True)
gaia.drop(['tyc'], axis=1, inplace=True)  # dropping tycho identifier

# Color correct
def correct(x):
    if x <= 0.5:
        BV = x - 0.006 - 1.069E-01*x + 1.459E-01*x**2
    else:
        BV = x - 7.813E-03*x - 1.489E-01*x**2 + 3.384E-02*x**3
    return BV

gaia['BVcolor'] = gaia['Tcolor'].map(correct)


# Write to uniview (~8 mins)
header = """# GAIA Data Release 1
# Including only 5-sigma parallaxes
# Matched against Tycho-2 for color information
# Colors corrected from Mamajek, Meyer, Liebert, 2002 AJ, 124, 1670 and Mamajek, Meyer, Liebert, 2006, AJ, 131, 2360
# Prepared by: David R. Rodriguez"""
uni.write_speck('gaia_dr1/gaia_dr1_tycho2.speck', gaia, header=header, texture=True, comment='comment')
