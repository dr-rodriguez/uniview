# Produce a speck file of the GAIA DR1 data with parallax information

from uniview import Uniview
import numpy as np
from tqdm import tqdm
import os
import math
import pandas as pd
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord

# Load the data
gaia_path = '/Users/strakul/data/gaia/dr1-tgas/fits'
files = os.listdir(gaia_path)
columns = ['ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'phot_g_mean_mag', 'hip', 'tycho2_id']

# All data
count = 0
with tqdm(total=len(files)) as pbar:
    for file in files:
        if not file.endswith('.fits'):
            continue
        filename = os.path.join(gaia_path, file)
        t = Table.read(filename, format='fits')

        if count == 0:
            df = t[columns].to_pandas()
        else:
            temp = t[columns].to_pandas()
            df = pd.concat([df, temp], axis=0, ignore_index=True)

        count += 1
        pbar.update(1)

# Filter out parallaxes. Only 5 sigma or better
ind = df['parallax'] / df['parallax_error'] >= 5
df = df[ind]
df['distance'] = 1000./df['parallax']

# Set name columns, drop hip and tycho2_id
df['name'] = ['HIP{:.0f}'.format(x) if not math.isnan(x) else 'TYC{}'.format(y.strip())
              for x, y in zip(df['hip'], df['tycho2_id'])]
df.drop(['hip', 'tycho2_id'], axis=1, inplace=True)

# Convert magnitudes to absolute then to 'lum'
dm = 5 * np.log10(df['distance']) - 5
absmag = df['phot_g_mean_mag'] - dm
z0 = 4.75  # may need to update this. GAIA G0=25.525 to convert to e/s (?)
# Update: may not matter as Partiview renormalizes from 0 to 1
lum = np.power(10, (absmag - z0) / (-2.5))
df['lum'] = lum
df['absmag'] = absmag
df[['lum','absmag', 'phot_g_mean_mag','distance']].describe()

# Convert the units
# This fixes for the 2015 epoch
ra_2000 = df['ra'] - ((df['pmra']*np.cos(df['dec']*np.pi/180))*1/1000.*1/3600.)*(2015.-2000.)
dec_2000 = df['dec'] - (df['pmdec']*1/1000.*1/3600.)*(2015.-2000.)
c_gaia = SkyCoord(ra=ra_2000.tolist()*u.degree, dec=dec_2000.tolist()*u.degree)

factor = np.pi / 180.
df['X'] = df['distance'] * np.cos(c_gaia.galactic.b.value*factor) * np.cos(c_gaia.galactic.l.value*factor)
df['Y'] = df['distance'] * np.cos(c_gaia.galactic.b.value*factor) * np.sin(c_gaia.galactic.l.value*factor)
df['Z'] = df['distance'] * np.sin(c_gaia.galactic.b.value*factor)

# Drop some columns
# df.drop(['parallax_error', 'pmra', 'pmdec', 'absmag', 'distance', 'lum'], axis=1, inplace=True)
# df.reset_index(inplace=True)  # very important
#
# # Convert some columns to strings to fix the number of decimal points
# cols = ['parallax', 'phot_g_mean_mag']
# # Merge together
# df = pd.concat([df[['X', 'Y', 'Z', 'ra', 'dec']].round(5), df[cols].round(2), df['name']], axis=1, ignore_index=False)
# df.reset_index(inplace=True)
# df.drop('index', axis=1, inplace=True)

# Write to speck and label
uni = Uniview()
header = "# GAIA DR1\n# Only 5-sigma parallaxes"
uni.write_speck('gaia_dr1/gaia_dr1.speck', df, header=header, texture=True, comment='name')  # speck file
# uni.write_speck('gaia_dr1/gaia_dr1_notexture.speck', df, header=header, texture=False, comment='name')  # no texture
# uni.write_label('gaia_dr1/gaia_dr1.label', df, name='name', header=header)  # label file: CAUSES ERRORS IN UNIVIEW
