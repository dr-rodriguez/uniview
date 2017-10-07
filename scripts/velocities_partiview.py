# Script to add partiview specific velocities

from astropy.table import Table
import numpy as np
from uniview import Uniview, transform
from astropy import units as u
from astropy.coordinates import SkyCoord


def space_velocity(ra, dec, dist, pmRA, pmDE, rvel):
    """Adapted from Brian's script emailed to Jackie"""
    ra = ra * np.pi/180
    dec = dec * np.pi/180

    # 1 km/sec = 1.0226903 pc/Myr.  Neat huh?  --STUART LEVY
    rv = rvel * 1.0226903

    # 1 mas/yr = .004848 pc/Myr / pc
    tvs = dist * 0.004848

    Rx = np.cos(ra) * np.cos(dec)
    Ry = np.sin(ra) * np.cos(dec)
    Rz = np.sin(dec)

    v_dec_x = -1*Rx
    v_dec_y = -1*Ry
    v_dec_z = 1 - Rz

    v_ra_x = -1*Ry
    v_ra_y = Rx
    v_ra_z = 0.0

    v_dec_mag = np.sqrt(v_dec_x * v_dec_x + v_dec_y * v_dec_y + v_dec_z * v_dec_z)
    v_ra_mag = np.sqrt(v_ra_x * v_ra_x + v_ra_y * v_ra_y + v_ra_z * v_ra_z)

    sv_dec = (pmDE * tvs) / v_dec_mag
    sv_ra = (pmRA * tvs) / v_ra_mag

    vx = rv * Rx + sv_dec * v_dec_x + sv_ra * v_ra_x
    vy = rv * Ry + sv_dec * v_dec_y + sv_ra * v_ra_y
    vz = rv * Rz + sv_dec * v_dec_z + sv_ra * v_ra_z

    # TRANSFORM TO GALACTIC COORDS
    vx_galactic = -0.05487554 * vx - 0.8734371 * vy - 0.483835 * vz
    vy_galactic = 0.4941095 * vx - 0.4448296 * vy + 0.7469823 * vz
    vz_galactic = -0.8676661 * vx - 0.1980764 * vy + 0.4559838 * vz

    speed = np.sqrt(vx_galactic * vx_galactic + vy_galactic * vy_galactic + vz_galactic * vz_galactic)

    return vx_galactic, vy_galactic, vz_galactic, speed


filename = 'input_tables/oh_table_matched_to_rave.fits'

t = Table.read(filename)

df = t.to_pandas()

columns = ['row_id', 'tgas_source_id', 'name', 'ra', 'dec', 'parallax_1', 'distance_1', 'G', 'rave_obs_id_1', 'rv',
           'group_id', 'group_size', 'pmRA_TGAS', 'pmDE_TGAS']

df[columns]

vx_galactic, vy_galactic, vz_galactic, speed = space_velocity(df['ra'], df['dec'], df['distance_1'], df['pmRA_TGAS'], df['pmDE_TGAS'], df['rv'])

df['vx'] = vx_galactic
df['vy'] = vy_galactic
df['vz'] = vz_galactic
df['speed'] = speed

x, y, z = transform(df['ra'], df['dec'], df['distance_1'])

df['X'] = x
df['Y'] = y
df['Z'] = z

columns = ['X','Y', 'Z', 'row_id', 'tgas_source_id', 'name', 'ra', 'dec', 'parallax_1', 'distance_1', 'G', 'rave_obs_id_1', 'rv',
           'group_id', 'group_size', 'pmRA_TGAS', 'pmDE_TGAS', 'vx', 'vy', 'vz', 'speed']
df[columns]

# Write out speck file
uni = Uniview()
header = """# GAIA Semyeong Oh Co-Moving stars with RAVE
# Matched against RAVE for radial velocity information
# Oh, S., Price-Whelan, A.~M., Hogg, D.W., Morton, T.D., & Spergel, D.N. 2017, AJ, 153, 257
# Prepared by: David R. Rodriguez"""
uni.write_speck('DU/oh_gaiadr1_rave.speck', df[columns], header=header, texture=False)

stars = uni.load_speck('DU/stars.speck')



# TODO: Cross match to check vx, vy, vz, speed
from astroML.crossmatch import crossmatch_angular
max_radius = 5. / 3600  # 5 arcsec
import datetime

G = np.empty((len(df), 2), dtype=np.float64)
G[:, 0] = df['ra']
G[:, 1] = df['dec']

t = SkyCoord(w=stars['X'].tolist() * u.pc, u=stars['Y'].tolist() * u.pc, v=stars['Z'].tolist() * u.pc,
                   representation='cartesian', frame='galactic')
t.representation = 'spherical'
c_stars = t.transform_to('icrs')

# Matches with astroML
H = np.empty((len(stars), 2), dtype=np.float64)
H[:, 0] = c_stars.ra.value
H[:, 1] = c_stars.dec.value

print datetime.datetime.now().time()
dist, ind = crossmatch_angular(G, H, max_radius)  # ~2 mins for matching 1 mil vs 250,000, ~3 for vs 500,000
print datetime.datetime.now().time()
match = ~np.isinf(dist)

df_gaia['separation'] = dist * 3600
df_gaia.loc[match, 'match'] = ind[match]