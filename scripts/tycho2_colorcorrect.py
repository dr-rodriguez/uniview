# Correct colors for the GAIA-Tycho2 data

import numpy as np
import pandas as pd
from uniview import Uniview

# Read GAIA-Tycho2 speck file
uni = Uniview()
df = uni.load_speck('gaia_dr1/gaia_dr1_tycho2.speck')

# Process columns
# df.columns = ['Tcolor' if x == 'color' else x for x in df.columns]
df['Tcolor'] = pd.to_numeric(df['Tcolor'], errors='coerce')
df.drop(['texture'], axis=1, inplace=True)
# Tcolor ranges from -1.88 to 4.8, which seems very red

def correct(x):
    if x <= 0.5:
        BV = x - 0.006 - 1.069E-01*x + 1.459E-01*x**2
    else:
        BV = x - 7.813E-03*x - 1.489E-01*x**2 + 3.384E-02*x**3
    return BV

df['BVcolor'] = df['Tcolor'].map(correct)

"""
for 0.5 < (BT-VT) < 2.0:
B-V = (BT-VT) - 7.813E-03(BT-VT) - 1.489E-01(BT-VT)^2 + 3.384E-02(BT-VT)^3


for -0.25 < (BT-VT) < 0.5:
B-V = (BT-VT) - 0.006 -1.069E-01(BT-VT) + 1.459E-01(BT-VT)^2

From:
Mamajek, Eric E.; Meyer, Michael R.; Liebert, James
Astronomical Journal, Volume 124, Issue 3, pp. 1670-1694
2002AJ....124.1670M
and
Mamajek, Meyer, Liebert, 2006, AJ, 131, 2360 (erratum)
"""

# Write to speck
header = """# GAIA Data Release 1
# Including only 5-sigma parallaxes
# Matched against Tycho-2 for color information
# Colors corrected from Mamajek, Meyer, Liebert, 2002 AJ, 124, 1670 and Mamajek, Meyer, Liebert, 2006, AJ, 131, 2360
# Prepared by: David R. Rodriguez"""
uni.write_speck('gaia_dr1/gaia_dr1_tycho2.speck', df, header=header, texture=True, comment='comment')


# CMD, for fun
# import matplotlib.pyplot as plt
# cols = ['parallax', 'BTmag', 'VTmag', 'Gmag']
# for col in cols:
#     df[col] = pd.to_numeric(df[col], errors='coerce')
# dist = 1000./df['parallax']
# DM = 5*np.log10(dist/10)
# plt.scatter(df['BVcolor'], df['Gmag']-DM, s=20, marker='.', c='blue')
