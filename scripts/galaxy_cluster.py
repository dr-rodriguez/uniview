
import matplotlib.pyplot as plt
from uniview import Uniview, transform
from sklearn.cluster import DBSCAN
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.misc import comb
from math import cos, sin, acos
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import astropy.units as u
from mpl_toolkits.mplot3d import Axes3D  # For 3D plot

uni = Uniview()
# df = uni.load_speck('2mass_extragal/2MASS_XSCz_0.00_0.01.speck')
df = uni.load_speck('2df/2dfgals.speck')
df = df.apply(lambda x: pd.to_numeric(x, errors='ignore'))  # convert all to numeric
if 'rahrs' in df.columns.tolist(): df['radeg'] = df.apply(lambda row: row['rahrs']*15, axis=1)

# Check if we need angular distance or just euclidean for XYZ  -> euclidian in XYZ is fine

# Compute pairwise angular distance
# http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist
def angle_dist(u, v):
    ra1 = u[0]
    ra2 = v[0]
    dec1 = u[1]
    dec2 = v[1]
    ra1 = np.deg2rad(ra1)
    dec1 = np.deg2rad(dec1)
    ra2 = np.deg2rad(ra2)
    dec2 = np.deg2rad(dec2)
    d = sin(dec1) * sin(dec2) + cos(dec1) * cos(dec2) * cos(ra1-ra2)
    # global pbar
    pbar.update()
    return np.rad2deg(acos(d))
with tqdm(total=int(comb(len(df), 2))) as pbar:
    dm = pdist(df[['radeg', 'decdeg']], angle_dist)  # Very time-consuming
dm_square = squareform(dm)


# c1 = SkyCoord(ra=53.602305*u.degree, dec=-21.235798*u.degree)
# c2 = SkyCoord(ra=170.780340*u.degree, dec=30.478548*u.degree)
# c1.separation(c2)

# cluster = DBSCAN(eps=5, metric='precomputed', min_samples=10)
# cluster.fit(dm_square)

cluster = DBSCAN(eps=5, metric='euclidean', min_samples=10)
cluster.fit(df[['X', 'Y', 'Z']])
core_samples_mask = np.zeros_like(cluster.labels_, dtype=bool)
core_samples_mask[cluster.core_sample_indices_] = True
df['cluster'] = cluster.labels_

# Converting coordinates for sky plots (this takes a while)
c = SkyCoord(ra=df['radeg']*u.deg, dec=df['decdeg']*u.deg, frame='icrs')
df['ra_rad'] = c.ra.wrap_at(180 * u.deg).radian
df['dec_rad'] = c.dec.radian


plt.subplot(111, projection="aitoff")
plt.scatter(df['ra_rad'], df['dec_rad'], c=df['cluster'], alpha=0.2)

# Creating a subset
df_c = df[df['cluster'] == 2]
# df_c = df[core_samples_mask & (df['cluster'] >= 1)]

fig = plt.figure()
plt.subplot(111, projection="aitoff")
# plt.scatter(df['ra_rad'], df['dec_rad'], c='grey', alpha=0.1)
plt.scatter(df_c['ra_rad'], df_c['dec_rad'], c=df_c['cluster'], alpha=0.8)

# 3D Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# ax.scatter(df['X'], df['Y'], df['Z'], c='grey', alpha=0.1)
ax.scatter(df_c['X'], df_c['Y'], df_c['Z'], c=df_c['cluster'], alpha=0.8)
ax.set_xlabel('X (Mpc)')
ax.set_ylabel('Y (Mpc)')
ax.set_zlabel('Z (Mpc)')

# Average distance
dist_avg = df_c['distMpc'].mean()
dist_std = df_c['distMpc'].std()
sqz_factor = 0.5  # How much to shrink standard deviation. Final one will be dist_std * sqz_factor
new_dist = (df_c['distMpc'] - dist_avg) * sqz_factor + dist_avg

new_dist.mean()
new_dist.std()
new_dist.std()/dist_std

df_n = df_c.copy()
df_n['new_dist'] = new_dist
df_n['X'], df_n['Y'], df_n['Z'] = transform(df_n['radeg'], df_n['decdeg'], df_n['new_dist'])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(df_c['X'], df_c['Y'], df_c['Z'], c='red', alpha=0.8)
ax.scatter(df_n['X'], df_n['Y'], df_n['Z'], c='blue', alpha=0.8)
ax.set_xlabel('X (Mpc)')
ax.set_ylabel('Y (Mpc)')
ax.set_zlabel('Z (Mpc)')
# ax.set_xlim(-100, 100)
# ax.set_ylim(-100, 100)
# ax.set_zlim(-100, 100)


sqz_factor = 0.5  # How much to shrink standard deviation. Final one will be dist_std * sqz_factor
df_new = df.copy()
for cl in set(cluster.labels_):
    df_c = df[df['cluster'] == cl]

    dist_avg = df_c['distMpc'].mean()
    dist_std = df_c['distMpc'].std()

    if cl < 0 or dist_std > 10:
        continue

    print("Cluster #{}: {} objects".format(cl, len(df_c)))
    print('Old Distances. Mean: {} StDev: {}'.format(dist_avg, dist_std))

    new_dist = (df_c['distMpc'] - dist_avg) * sqz_factor + dist_avg
    print('New Distances. Mean: {} StDev: {}'.format(new_dist.mean(), new_dist.std()))

    # Set new distances
    df_new.loc[df['cluster'] == cl, 'distMpc'] = new_dist

    # Same for lightyears
    dist_avg = df_c['distMly'].mean()
    dist_std = df_c['distMly'].std()
    new_dist = (df_c['distMly'] - dist_avg) * sqz_factor + dist_avg
    df_new.loc[df['cluster'] == cl, 'distMly'] = new_dist

cl=2
df_new[df['cluster']==cl]
df_new[df['cluster']==cl][['cluster','distMpc','distMly']]
df[df['cluster']==cl][['cluster','distMpc','distMly']]
df_new[df['cluster']==cl]['distMpc'] - df[df['cluster']==cl]['distMpc']

# Transform and output to speck
df_new['X'], df_new['Y'], df_new['Z'] = transform(df_new['radeg'], df_new['decdeg'], df_new['distMpc'])
uni.write_speck('2df/2dfgals_new.speck', df_new, '2df/2dfgals.speck')
