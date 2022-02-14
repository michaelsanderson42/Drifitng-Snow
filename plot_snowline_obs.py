import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def read_ssgb_data(station_name):

    ssgb_datadir = '/project/earthobs/SNOW/SSGB/data/'
    filename = f"SSGB_st_{station_name}.txt"

    df = pd.read_csv(os.path.join(ssgb_datadir, filename), delimiter="|", header=0,
        usecols=[0,2,3])

# Replace days without observations ('n') with a zero
    df.replace(to_replace="n", value="0", inplace=True)

# Convert snowline elevations to ints
    df['SnowlineElev'] = df['SnowlineElev'].astype(int)

    df['date'] = pd.to_datetime(df['Date']) 
    df.set_index(['date'], inplace=True)

    return df


def count_snowline_elevations(station_name):

    df = read_ssgb_data(station_name)

    snow_elevs = df.loc[df['SnowlineElev'] > 99, 'SnowlineElev'].tolist()
    unique_elevs = sorted(list(set(snow_elevs)))
    bin_edges = unique_elevs[:]
    bin_edges.insert(0, 99)
    bin_edges = np.array(bin_edges) + 1
    print(bin_edges)

    l_out = []

# Get the unique years - should be sorted
    years = sorted(list(set([int(d.split('-')[0]) for d in df['Date'].tolist()])))

    for year in years:
        d0 = '{:4d}-10-01'.format(year)
        d1 = '{:4d}-05-31'.format(year+1)

        s_elevs = df.loc[d0:d1, 'SnowlineElev'].to_numpy()
        snow_elevs = s_elevs[(s_elevs > 99)]
        n_per_elev_range, _ = np.histogram(snow_elevs, bins=bin_edges)
        n_per_elev_range = list(n_per_elev_range)
        n_per_elev_range.insert(0, year)
        l_out.append(n_per_elev_range)

    columns = ['year']
    columns.extend(unique_elevs)
    df_out = pd.DataFrame(l_out)
#   df_out = df_out.transpose()
    df_out.columns = columns
    print(df_out)

    return df_out


def plot_num_per_elev_range(df, station_name):

    columns = list(df.columns)
    elevs = [c for c in columns[1:]]
    years = df['year'].to_numpy()
    colours = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628']

    fig = plt.figure()
#   ax = fig.add_subplot(2,1,1)
    bx = fig.add_subplot(1,1,1)

    cumulative_counts = np.zeros(len(years), dtype=np.float)
    for i, elev in enumerate(elevs):
        num_days = df[elev].to_numpy()
#       ax.plot(years, num_days, linestyle='-', color=colours[i], marker='None', label=str(elev))
        cumulative_counts += num_days
        bx.plot(years, cumulative_counts, linestyle='-', color=colours[i], marker='None', label=str(elev))

#   ax.set_xlim(np.min(years)-1, np.max(years)+1)
#   ax.set_xlabel('Year')
#   ax.set_ylim(ymin=-1, ymax=80)
#   ax.set_ylabel('Number of days')
    handles, labels = bx.get_legend_handles_labels()
    dummy_patch = mpatches.Patch(color='w', label=' ', alpha=0)
    handles.insert(3, dummy_patch)
#   ax.legend(handles=handles, loc='upper right', ncol=2, frameon=False)
#   ax.legend(loc='upper right', ncol=2, frameon=False)

    bx.set_xlim(np.min(years)-1, np.max(years)+1)
    bx.set_xlabel('Year')
#   bx.set_ylim(ymin=-1, ymax=80)
    bx.set_ylabel('Number of days below elevation')
    bx.legend(handles=handles, loc='upper right', ncol=2, frameon=False)

    plt.savefig('SSGB_{}_snow_counts.png'.format(station_name), dpi=300)
    plt.close()
#   plt.show()


if __name__ == '__main__':

    station_name = 'RannochPS'
    df = count_snowline_elevations(station_name)

    plot_num_per_elev_range(df, station_name)
