# Cluster compression script for 2df galaxies

from uniview import Uniview, transform
from sklearn.cluster import DBSCAN
import pandas as pd
import numpy as np
from tqdm import tqdm
from astropy import units as u
from astropy.coordinates import SkyCoord

uni = Uniview()
df = uni.load_speck('2df/2dFgals.speck')  # now using the new data Brian provided
df = df.apply(lambda x: pd.to_numeric(x, errors='ignore'))  # convert all to numeric
if 'rahrs' in df.columns.tolist(): df['radeg'] = df.apply(lambda row: row['rahrs']*15, axis=1)

cluster = DBSCAN(eps=5, metric='euclidean', min_samples=10)
cluster.fit(df[['X', 'Y', 'Z']])
core_samples_mask = np.zeros_like(cluster.labels_, dtype=bool)
core_samples_mask[cluster.core_sample_indices_] = True
df['cluster'] = cluster.labels_

sqz_factor = 0.2  # How much to shrink standard deviation. Final one will be dist_std * sqz_factor
stdcut = 10  # maximum distance standard deviation to consider; clusters more dispersed will be skipped
ang_stdcut = 1  # maximum standard deviation, in degrees on sky, to consider; clusters more dispersed will be skipped
df_new = df.copy()
verbose = True
usepbar = False
if usepbar: pbar = tqdm(total=len(set(cluster.labels_)))
for cl in set(cluster.labels_):
    df_c = df[df['cluster'] == cl]

    # Skip those not in a cluster
    if cl < 0:
        if usepbar: pbar.update(1)
        continue

    dist_avg = df_c['distMpc'].mean()
    dist_std = df_c['distMpc'].std()

    # Convert XYZ to ra/dec.
    # in astropy coordinates, x=w, y=u, z=v when in the galactic frame
    t = SkyCoord(w=df_c['X'].tolist() * u.Mpc, u=df_c['Y'].tolist() * u.Mpc, v=df_c['Z'].tolist() * u.Mpc,
                 representation='cartesian', frame='galactic')
    t.representation = 'spherical'
    coords = t.transform_to('icrs')
    # Find center and deviation from center
    center_coords = SkyCoord(ra=np.mean(coords.ra.value) * u.degree, dec=np.mean(coords.dec.value) * u.degree,
                             frame='icrs')
    ang_dist = center_coords.separation(coords)
    ang_std = np.std(ang_dist).value

    if dist_std > stdcut or ang_std > ang_stdcut:
        if usepbar: pbar.update(1)
        continue

    new_dist = (df_c['distMpc'] - dist_avg) * sqz_factor + dist_avg

    if verbose:
        print("Cluster #{}: {} objects".format(cl, len(df_c)))
        print('Old Distances. Mean: {} StDev: {}'.format(dist_avg, dist_std))
        print('New Distances. Mean: {} StDev: {}'.format(new_dist.mean(), new_dist.std()))
        print('Angular standard deviation: {}'.format(ang_std))

    # Set new distances
    df_new.loc[df['cluster'] == cl, 'distMpc'] = new_dist

    # Same for lightyears
    dist_avg = df_c['distGly'].mean()
    dist_std = df_c['distGly'].std()
    new_dist = (df_c['distGly'] - dist_avg) * sqz_factor + dist_avg
    df_new.loc[df['cluster'] == cl, 'distGly'] = new_dist

    # Transform to XYZ
    df_new.loc[df['cluster'] == cl, 'X'], df_new.loc[df['cluster'] == cl, 'Y'], df_new.loc[df['cluster'] == cl, 'Z'] = \
        transform(coords.ra.value, coords.dec.value, df_new.loc[df['cluster'] == cl, 'distMpc'])

    if usepbar: pbar.update(1)
if usepbar: pbar.close()

print('Writing file...')
uni.write_speck('2df/2dfgals_new_{}_std{}_ang{}.speck'.format(sqz_factor, stdcut, ang_stdcut), df_new, oldfile='2df/2dFgals.speck')
