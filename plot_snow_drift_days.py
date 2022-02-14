import os
import iris
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


def get_ukcp18_label(var_name):

    plot_labels = {'tasmin': 'Days below freezing',
                   'wsgmax2m': 'Gusts above 12 kn',
                   'snw': 'Lying snow',
                   'snow_frac': 'Cover > 50%',
                   'drift': 'Drifting snow'}

    if var_name in plot_labels:
        plot_label = plot_labels[var_name]
    else:
        plot_label = ''

    return plot_label


def read_snow_drift_data(datadir, ensemble_member, year_range, var_name):

    clist = iris.cube.CubeList()
    for the_year in list(range(year_range[0]+1, year_range[1])):
        filename = "drift_rcp85_land-cpm_uk_2.2km_{:02d}_{:4d}.nc".format(ensemble_member, the_year)
        cube = iris.load_cube(os.path.join(datadir, filename), var_name)
        clist.append(cube)

    return clist.merge_cube()


def read_and_plot_data(datadir, ensemble_member, year_ranges, var_names, i, j):

    colours = ['red', 'green', 'blue', 'orange', 'grey']
    nrows = len(var_names)
    ncols = len(year_ranges)
    fig = plt.figure(figsize=(15,15))

    for row, var_name in enumerate(var_names):
        ymin = 999.0
        ymax = -999.0
        axes = []
        for col, year_range in enumerate(year_ranges, start=1):
            plotnum = row*ncols + col
            print(plotnum, var_name)
            year_centre = (year_range[0] + year_range[1]) // 2
            cube = read_snow_drift_data(datadir, ensemble_member, year_range, var_name)

            ax = fig.add_subplot(nrows, ncols, plotnum)
            years = cube.coord('year').points
            data = cube.data[:,j,i]
            ax.plot(years, data, marker='None', linestyle='-', color=colours[row])
            plot_label = get_ukcp18_label(var_name)
            the_title = '{}: {:4d}-{:4d}'.format(plot_label, year_range[0], year_range[1])

            ax.set_title(the_title)
            ax.set_xlim(year_range[0], year_range[1])
            ax.set_xticks([year_range[0], year_centre, year_range[1]])
            if var_name == 'tasmin':
                ax.set_yticks([20, 60, 100, 140])
            ax_ymin = 10 * np.floor(np.min(data) / 10)
            ax_ymax = 10 * np.ceil(np.max(data) / 10)
            if ymin > ax_ymin:
                ymin = ax_ymin
            if ymax < ax_ymax:
                ymax = ax_ymax
            axes.append(ax)
            print(the_title, np.min(data), np.max(data))
# Set a common y-coordinate range for each variable
        for zx in axes:
            zx.set_ylim(ymin, ymax)

# Add a big subplot for common x- and y-axis labels
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel("Year")
    plt.ylabel("Number of days")
    plt.subplots_adjust(hspace=0.3)

    plt.savefig('projected_numbers_snow_days.png', dpi=300)
    plt.close()
#   plt.show()


def main():

    datadir = '/scratch/hadmi/UKCP18/processed/'
    ensemble_member = 1
    year_ranges = [[1980, 2000], [2020, 2040], [2060, 2080]]
    var_names = ['tasmin', 'wsgmax2m', 'snw', 'snow_frac', 'drift']
    i = 2
    j = 2

    read_and_plot_data(datadir, ensemble_member, year_ranges, var_names, i, j)


if __name__ == '__main__':
    main()
