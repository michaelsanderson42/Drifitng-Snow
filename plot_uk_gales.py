import os
import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib.pyplot as plt


def read_gale_obs(datadir, filename):

    df = pd.read_csv(os.path.join(datadir, filename), header=0)
    dates = pd.to_datetime(dict(year=df['Year'], month=df['Month'], day=df['Day']))
    df.insert(0, 'date', dates)
    df.set_index('date', inplace=True)

    return df


def plot_gale_days(df):

    fig = plt.figure()
    nrows = 3
    ncols = 1

    years = list(set(df['Year'].tolist()))
    years.sort()

    gale_types = ['gales', 'severe gales', 'very severe gales']
    gales = []
    severe = []
    vsevere = []

    for year in years[:-1]:
# Select data recorded between October and May, the period used by the
# Snow Survey of Great Britain.
        d1 = "{:4d}-10-01".format(year)
        d2 = "{:4d}-05-31".format(year+1)
    
        df_tmp = df.loc[d1:d2]

        s = df_tmp['Gale Index exceedances'].value_counts()
        data_names = list(s.index)

        if gale_types[0] in data_names:
            gales.append(int(s['gales']))
        else:
            gales.append(0)
        if gale_types[1] in data_names:
            severe.append(int(s['severe gales']))
        else:
            severe.append(0)
        if gale_types[2] in data_names:
            vsevere.append(int(s['very severe gales']))
        else:
            vsevere.append(0)

# Plot the counts
    i = years.index(1900)
    j = years.index(2020)
    xmin = 1899.4
    xmax = 2020.5
    width = 0.7

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=True)

    plotnum = 0
#   ax = fig.add_subplot(nrows, ncols, plotnum)
    ax = axes[0]
    ax.bar(years[i:j+1], gales[i:j+1], width=width, color='grey')
    ax.set_xlim(xmin, xmax)
    ax.text(0.02, 0.88, gale_types[0], transform=ax.transAxes)
    print(linregress(years[i:j+1], gales[i:j+1]))
    plotnum = 2
#   bx = fig.add_subplot(nrows, ncols, plotnum)
    bx = axes[1]
    bx.bar(years[i:j+1], severe[i:j+1], width=width, color='grey')
    bx.set_xlim(xmin, xmax)
    bx.text(0.02, 0.88, gale_types[1], transform=bx.transAxes)
    print(linregress(years[i:j+1], severe[i:j+1]))
    plotnum = 3
#   cx = fig.add_subplot(nrows, ncols, plotnum)
    cx = axes[2]
    cx.bar(years[i:j+1], vsevere[i:j+1], width=width, color='grey')
    cx.set_xlim(xmin, xmax)
    cx.text(0.02, 0.88, gale_types[2], transform=cx.transAxes)
    cx.set_xlabel('Year')
    print(linregress(years[i:j+1], vsevere[i:j+1]))

    zx = fig.add_subplot(111, frame_on=False)
    zx.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    zx.set_ylabel('Frequency')

#   plt.subplots_adjust(hspace=0.3)
    plt.savefig('UK_Jenkinson_gales.png', dpi=300)
    plt.close()
#   plt.show()


if __name__ == '__main__':

    datadir = '/home/h03/hadmi/Python/Snow/'
    filename = "CRU_UK_jenkinson_gale_index.csv"
    df = read_gale_obs(datadir, filename)
    plot_gale_days(df)

