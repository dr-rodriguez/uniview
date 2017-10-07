# Loading up and matching the stars.speck (Hip) data against GAIA+Tycho2
# Goal is to examine the color mismatch and save a stars+GAIA subset for comparison

import pandas as pd
from uniview import Uniview
from astroML.crossmatch import crossmatch_angular
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

uni = Uniview()

stars = uni.load_speck('DU/stars.speck')
gaia = uni.load_speck('gaia_dr1/gaia_dr1_tycho2.speck')

# Drop the sun
stars.drop(0, axis=0, inplace=True)

# Make certain columns numeric
cols = ['ra', 'dec', 'BVcolor', 'Tcolor', 'pmra', 'pmdec', 'Gmag', 'parallax']
for col in cols:
    gaia[col] = pd.to_numeric(gaia[col], errors='coerce')

cols = ['X', 'Y', 'Z', 'colorb_v', 'absmag', 'lum']
for col in cols:
    stars[col] = pd.to_numeric(stars[col], errors='coerce')

# Convert XYZ to ra/dec
# in astropy coordinates, x=w, y=u, z=v when in the galactic frame
t = SkyCoord(w=stars['X'].tolist() * u.pc, u=stars['Y'].tolist() * u.pc, v=stars['Z'].tolist() * u.pc,
                   representation='cartesian', frame='galactic')
t.representation = 'spherical'
c_stars = t.transform_to('icrs')

# Matches with astroML
H = np.empty((len(stars), 2), dtype=np.float64)
H[:, 0] = c_stars.ra.value
H[:, 1] = c_stars.dec.value

ra_2000 = gaia['ra'] - ((gaia['pmra']*np.cos(gaia['dec']*np.pi/180))*1/1000.*1/3600.)*(2015.-2000.)
dec_2000 = gaia['dec'] - (gaia['pmdec']*1/1000.*1/3600.)*(2015.-2000.)
G = np.empty((len(gaia), 2), dtype=np.float64)
G[:, 0] = ra_2000
G[:, 1] = dec_2000

max_radius = 2. / 3600  # 2 arcsec
import datetime
print datetime.datetime.now().time()
dist, ind = crossmatch_angular(H, G, max_radius)  # flipped so stars are first (less entries to match)
print datetime.datetime.now().time()
match = ~np.isinf(dist)

stars['separation'] = dist * 3600
stars.loc[match, 'match'] = ind[match]

col_subset = ['ra', 'dec', 'Gmag', 'BTmag', 'VTmag', 'Tcolor', 'BVcolor', 'parallax', 'lum']
hip_gaia = pd.merge(stars, gaia[col_subset], how='left', left_on='match', left_index=False, right_index=True, sort=False,
                    suffixes=('', '_g'), copy=True, indicator=False)

# Drop those that didn't match
hip_gaia.dropna(inplace=True)

# Drop columns and output to file
output = hip_gaia.drop(['vx', 'vy', 'vz', 'speed', 'match', 'texnum', 'dcalc', 'plxerr'], axis=1, inplace=False)
output['hipnum'] = pd.to_numeric(output['hipnum'], errors='coerce')
header = """# TGAS HIP crossmatch
# Including only 5-sigma parallaxes
# Matched against Tycho-2 for color information and HIP for comparison with stars.speck
# Colors corrected from Mamajek, Meyer, Liebert, 2002 AJ, 124, 1670 and Mamajek, Meyer, Liebert, 2006, AJ, 131, 2360
# Prepared by: David R. Rodriguez"""
uni.write_speck('gaia_dr1/gaia_dr1_hip.speck', output, header=header, texture=True, comment='comment')

output[(output['hipnum'] < 49680) & (output['hipnum'] > 49660)][['colorb_v', 'BVcolor', 'hipnum', 'separation', 'comment']]
output['appmag'] = pd.to_numeric(output['appmag'], errors='coerce')
output.sort_values(by='appmag', ascending=True).head()

# Drop some columns
hip_gaia.drop(['texnum', 'distly', 'dcalc', 'match', 'separation', 'plxerr'], axis=1, inplace=True)

# Drop those with no color info in Tycho-2 (may want to look into why these were matched in TGAS)
# Also, drop those with minimum/maximum colorb_v values (looks like a floor/ceiling was set here)
bad = (hip_gaia['BTmag'] == '999.0') | (hip_gaia['colorb_v'] == -0.63) | (hip_gaia['colorb_v'] == 2.057)
hip_gaia = hip_gaia[~bad]





# Plot color differences
import seaborn as sns
import matplotlib.pyplot as plt

# Compare histograms of B-V
n, bins = np.histogram(gaia['BVcolor'], 30)
g = sns.distplot(gaia['BVcolor'], kde=False, norm_hist=True, label='GAIA', bins=bins)
g = sns.distplot(stars['colorb_v'], kde=False, norm_hist=True, ax=g, label='HIP', bins=bins, axlabel='B-V')
g.legend()
plt.savefig('colordist.png')

g = sns.jointplot('colorb_v', 'BVcolor', data=hip_gaia, kind='reg')
g.set_axis_labels('B-V (HIP)', 'B-V (GAIA)')
plt.savefig('colorcompare.png')

# Fit the data
from scipy.optimize import curve_fit

# Write a function for scipy
def myfunction(x, a, b):
    return a + b * x

res, cov = curve_fit(myfunction, hip_gaia['colorb_v'], hip_gaia['BVcolor'])
print(res)

plt.scatter(hip_gaia['colorb_v'], hip_gaia['BVcolor'], marker='.', c='blue')
xdata = np.arange(-1, 2.1, 0.1)
plt.plot(xdata, myfunction(xdata, *res), color='red', marker='', linestyle='-', lw=2)

# How many GAIA stars redder than 2.05 and bluer than -0.6
(gaia['BVcolor'] > 2.05).sum()  # 1183
(gaia['BVcolor'] < -0.6).sum()  # 294

# Compare luminosities
dist = 1000./hip_gaia['parallax']
gaia_abs = hip_gaia['Gmag'] - (5*np.log10(dist) - 5)
plt.scatter(hip_gaia['absmag'], gaia_abs, marker='.', c='blue')

res, cov = curve_fit(myfunction, hip_gaia['absmag'], gaia_abs)
print(res)

hip_gaia['G_abs'] = gaia_abs

# Calculating lum
z0 = 4.75  # may need to update this. GAIA G0=25.525 to convert to e/s (?)
# Update: may not matter as Partiview renormalizes from 0 to 1
lum = np.power(10, (gaia_abs - z0) / (-2.5))
lum.head()
hip_gaia[['lum', 'absmag', 'G_abs']].head()
(np.power(10, (hip_gaia['absmag'] - 4.73) / (-2.5))).head()


hip_gaia['lum'] = pd.to_numeric(hip_gaia['lum'], errors='coerce')
res, cov = curve_fit(myfunction, np.log10(hip_gaia['lum']), hip_gaia['absmag'])
print(res)

plt.scatter(hip_gaia['lum'], np.power(10, (hip_gaia['absmag'] - 4.72) / (-2.5)), marker='.', c='blue')

np.log10(hip_gaia['lum'] / np.power((hip_gaia['absmag'] - 4.72) / (-2.5), 10))

# Quick cleanup and output for luminosities
gaia.drop('texture', axis=1, inplace=True)
dist = 1000./gaia['parallax']
gaia_abs = gaia['Gmag'] - (5*np.log10(dist) - 5)
lum = np.power(10, (gaia_abs - z0) / (-2.5))
gaia['lum'] = lum

# Write to uniview (~8 mins)
header = """# GAIA Data Release 1
# Including only 5-sigma parallaxes
# Matched against Tycho-2 for color information
# Colors corrected from Mamajek, Meyer, Liebert, 2002 AJ, 124, 1670 and Mamajek, Meyer, Liebert, 2006, AJ, 131, 2360
# Prepared by: David R. Rodriguez"""
uni.write_speck('gaia_dr1/gaia_dr1_tycho2.speck', gaia, header=header, texture=True, comment='comment')


# Extra checks vs Tycho2
gaia[gaia['BTmag'] == '999.0'].head()
# TYC48-748-1 did not match
df[df['tyc'] == '0048 00748 1']  # has no ra, dec, but does have RAdeg, DEdeg
df[['tyc','ra','RAdeg','dec','DEdeg']][bad]
bad = df['ra'].isnull()

