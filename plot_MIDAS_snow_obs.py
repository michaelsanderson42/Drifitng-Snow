import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def make_station_name(s):

    return ' '.join([c.capitalize() for c in s.split('_')])


def read_snow_obs(datadir, filename):

    nheaders = 9
    names = ['date', 'tmax', 'tmin', 'tmean', 'max gust', 'snow lying', 'snow falling']

    df = pd.read_csv(os.path.join(datadir, filename), skiprows=nheaders, header=0,
        names=names)
    df.set_index('date', inplace=True)

    return df


def plot_days_drifts(datadir, station_names):

    fig = plt.figure()
    nrows = 2
    ncols = 1
    plotnum = 0

# Wind speed minima for drifts and blizzards as used by Met Office
# are 12 knots and 17 knots respectively.
# 10 knots is equivalent to 5.14444 m/s.
    wind_drift_thresh = 12 * 0.514444

    for station_name in station_names:
        f = '_'.join(['0008594', station_name, 'results.csv'])
        df = read_snow_obs(datadir, f)
        if len(df) == 0:
            continue

# Count the number of non-missing values
        sr_tot = df.count()
        n_lying = int(sr_tot['snow lying'])
        n_falling = int(sr_tot['snow falling'])
        n_gusts = int(sr_tot['max gust'])
        print(f, n_lying, n_falling, n_gusts)

# If no or very few snow obs available, skip this station
        if n_lying < 5 and n_falling < 5:
            continue
        if n_gusts < 5:
            print(f'No wind gust data for {station_name}, skipping')
            continue

        years = list(set([int(s[:4]) for s in df.index.values]))
        years.sort()

        d_out = {}

        for year in years[:-1]:
# Select data recorded between October and May, the period used by the
# Snow Survey of Great Britain.
            d1 = "{:4d}-10-01".format(year)
            d2 = "{:4d}-05-31".format(year+1)
    
            df_tmp = df.loc[d1:d2]

            season_drift = 0

# Get the met data for this day
            for idx in df_tmp.index.values:
                tx = df_tmp.at[idx, 'tmax']
                tn = df_tmp.at[idx, 'tmin']
                tm = df_tmp.at[idx, 'tmean']
                gs = df_tmp.at[idx, 'max gust']
                sl = df_tmp.at[idx, 'snow lying']
                sf = df_tmp.at[idx, 'snow falling']

                if np.isnan(tm) and not (np.isnan(tx) or np.isnan(tn)):
                    tm = np.mean([tx, tn])
                if np.isnan(tn):
                    tn = 99.9
                if np.isnan(gs):
                    gs = 99.9
# If tmin > 0, assume too warm for snow to form drifts, as it would melt
                if tn > 0.0:
                    continue

                if gs > wind_drift_thresh and (sl > 0 or sf > 0):
                    season_drift += 1

            d_out[year] = season_drift
    
        sname = make_station_name(station_name)
        print(sname)
#       print(sname, years[0], years[-1])
        xdata = list(d_out.keys())
        season_drift = list(d_out.values())
        print(xdata)
        print(season_drift)

        plotnum += 1
        ax = fig.add_subplot(nrows, ncols, plotnum)
        ax.plot(xdata, season_drift, marker='s', color='grey', linestyle='None', label='Drifting snow')
#       ax.set_ylim(0, 100)
        if 'TULLOCH' in f:
            ax.set_xlim(xmin=2009.5, xmax=2021.5)
            ax.xaxis.get_major_locator().set_params(integer=True)
            ax.set_xlabel('Year')
            ax.set_ylabel('Number of days')
        lines, labels = ax.get_legend_handles_labels()
        ax.set_title(sname)

    plt.subplots_adjust(hspace=0.3)
    plt.savefig('Cruach_Scotland_snowdrift_days.png', dpi=300)
    plt.close()
#   plt.show()


def plot_days_snow(datadir, station_names):

    fig = plt.figure()
    nrows = 3
    ncols = 2
    plotnum = 0
    nsta = len(station_names)

    for station_name in station_names:
        f = '_'.join(['0008594', station_name, 'results.csv'])
        df = read_snow_obs(datadir, f)
        if len(df) < 2:
            print(f'No data for {station_name}, skipping')
            continue

# Count the number of non-missing values
        sr_tot = df.count()
        n_lying = int(sr_tot['snow lying'])
        n_falling = int(sr_tot['snow falling'])
        print(f, n_lying, n_falling)

# If no or very few snow obs available, skip this station
        if n_lying < 5 and n_falling < 5:
            print(f'No snow data for {station_name}, skipping')
            continue

        xdata = []
        ndays_lying = []
        ndays_falling = []
        plotnum += 1
        years = list(set([int(s[:4]) for s in df.index.values]))
        years.sort()

        for year in years[:-1]:
# Select data recorded between October and May, the period used by the
# Snow Survey of Great Britain.
            d1 = "{:4d}-10-01".format(year)
            d2 = "{:4d}-05-31".format(year+1)
    
            df_tmp = df.loc[d1:d2]
            nl = df_tmp['snow lying'].dropna().astype(int).sum()
            nf = df_tmp['snow falling'].dropna().astype(int).sum()
            if nl > 0 and nf > 0:
                xdata.append(year)
                nmiss_nl = df_tmp['snow lying'].isna().sum()
                nmiss_nf = df_tmp['snow falling'].isna().sum()
#               print(year, nl, nmiss_nl, nf, nmiss_nf)
                ndays_lying.append(nl)
                ndays_falling.append(nf)
    
        sname = make_station_name(f[8:-12])
#       print(sname, years[0], years[-1])
        ax = fig.add_subplot(nrows, ncols, plotnum)
        ax.plot(xdata, ndays_lying, marker='o', color='grey', linestyle='None', label='Lying snow')
        ax.plot(xdata, ndays_falling, 'k+', label='Falling snow')
        ax.set_ylim(0, 100)
        if 'TULLOCH' in f:
            ax.set_xlim(xmin=2011)
            ax.xaxis.get_major_locator().set_params(integer=True)
        if 'TYNDRUM' in f:
            ax.set_xlim(xmin=1989)
        if 'ARDTALNAIG' in f:
            ax.set_xlim(xmin=1969)
        if 'FASKALLY' in f:
            ax.set_xlim(xmin=1970)
        if plotnum == 6:
            ax.legend(loc=(0.4, 0.55), frameon=False, labelspacing=0.3)
        ax.set_title(sname)

    zx = fig.add_subplot(111, frame_on=False)
    zx.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    zx.set_ylabel('Number of days')

    plt.subplots_adjust(hspace=0.5)
    plt.savefig('Cruach_Scotland_snow_days.png', dpi=300)
    plt.close()
#   plt.show()


if __name__ == '__main__':

    datadir = '/home/h03/hadmi/Python/Snow/data/'

# Stations within a 50 km radius of the Cruach snow shed
    station_names = ['TULLOCH_BRIDGE', 'AONACH_MOR', 
                     'DALWHINNIE_NO_2', 'TYNDRUM_NO_3', 'GLEN_OGLE',
                     'ARDTALNAIG', 'ABERFELDY_DULL',
                     'FASKALLY'
    ]

#   plot_days_snow(datadir, station_names)
    plot_days_drifts(datadir, station_names)

