# Uniview object class
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from math import cos, sin, isnan
from tqdm import tqdm  # for progress bars


class Uniview:
    def __init__(self):
        return

    def load_speck(self, filename):
        """
        Load the speck file

        :param filename:
        :return: pandas DataFrame
        """

        columns = ['X', 'Y', 'Z']
        data = list()
        count = 0
        print('Loading file...')
        with open(filename, 'r') as f:
            for line in f:
                # Skip comments
                if line.startswith('#') or line in ['\n', '\r\n']:
                    continue

                # Prepare columns
                if line.startswith('datavar'):
                    elem = line.split()
                    columns.append(elem[-1])
                    continue

                if line.startswith('texture'):
                    continue

                # Comments column if present
                if '#' in line and count == 0:
                    columns.append('comment')
                    count += 1

                # Add entry
                elem = line.strip().split()
                datum = dict()
                for i, col in enumerate(columns):
                    x = elem[i]
                    if x == '#':
                        x = ' '.join(elem[i+1:])
                    datum[col] = x
                data.append(datum)

        df = pd.DataFrame(data, columns=columns)
        print('Data loaded.')
        return df

    def write_speck(self, filename, df, oldfile='', header='', texture=False, comment='name'):
        """
        Write the speck file

        :param filename:
        :param df:
        :param oldfile:
        :return:
        """

        # Add header from oldfile
        if oldfile and not header:
            header = ''
            with open(oldfile, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        header += line

        text = header + '\n'

        # Add datavar keys
        ignore = ['X', 'Y', 'Z', comment, 'index', 'level_0', 'texture']
        count = 0
        for i, col in enumerate(df.columns):
            if col in ignore: continue
            text += 'datavar {} {}\n'.format(count, col)
            count += 1

        if texture:
            text += 'datavar {0} texture\ntexturevar {0}\ntexture -M 1 halo.sgi\n'.format(count)

        # Add main text
        with tqdm(total=len(df)) as pbar:  # progress bar
            for i, row in df.iterrows():
                line = ' '.join(['{}'.format(x) for x in row[['X', 'Y', 'Z']]]) + ' '
                line += ' '.join(['{}'.format(x) for x,y in zip(row, row.keys()) if y not in ignore])
                if comment in row.keys() and not texture:
                    line += ' # {}'.format(row[comment])
                elif comment in row.keys() and texture:
                    line += ' 1 # {}'.format(row[comment])
                elif texture and comment not in row.keys():
                    line += ' 1 '
                line += '\n'

                text += line
                pbar.update(1)

        # Actually write the file
        with open(filename, 'w') as f:
            f.write(text)

        print('{} written.'.format(filename))

    def write_label(self, filename, df, name='name', oldfile='', header=''):
        """
        Write a label file

        :param filename:
        :param df:
        :param name:
        :return:
        """

        # Add header from oldfile
        if oldfile and not header:
            header = ''
            with open(oldfile, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        header += line

        text = header + '\n'

        text += 'textcolor 1\n'

        # Add main text
        with tqdm(total=len(df)) as pbar:  # progress bar
            for i, row in df.iterrows():
                line = ' '.join(['{}'.format(x) for x in row[['X', 'Y', 'Z']]])
                line += ' text {}'.format(row[name])
                line += '\n'

                text += line
                pbar.update(1)

        # Actually write the file
        with open(filename, 'w') as f:
            f.write(text)

        print('{} written.'.format(filename))


def transform(ra, dec, dist):
    c = SkyCoord(ra=ra, dec=dec, frame='icrs', equinox='J2000', unit='deg')
    #c = c.transform_to(FK5(equinox='B1950'))  # Precess

    # Convert to XYZ
    l, b = c.galactic.l.radian, c.galactic.b.radian
    xgc = dist * map(cos, b) * map(cos, l)
    ygc = dist * map(cos, b) * map(sin, l)
    zgc = dist * map(sin, b)

    return xgc, ygc, zgc


def uniview_rgb(x1, x2, x3):
    return list(np.array([x1,x2,x3])/255.)


def parsedata(data):
    datatxt = ''
    labeltxt = ''

    # Re-index to avoid errors
    data = data.reindex(range(len(data)))

    for i in range(len(data)):
        if isnan(data.X[i]):
            continue

        temp = str(data.X[i]) + ' ' + str(data.Y[i]) + ' ' + str(data.Z[i]) + ' '
        datatxt += temp
        labeltxt += temp

        for column in data:
            if column in ('X','Y','Z'): continue
            if (column == 'Name') | (column == '_2MASS'):
                datatxt += ' # ' + data[column][i] + ' '
                labeltxt += 'text ' + data[column][i]
            else:
                datatxt += str(data[column][i]) + ' '

        datatxt += '\n'
        labeltxt += '\n'

    return datatxt, labeltxt


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
