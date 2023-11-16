from pyopenms import *
import pandas
import matplotlib as mpl
import matplotlib.pyplot as pyplot
import matplotlib.transforms as transforms
import matplotlib.ticker as mticker
import math
import seaborn
import os
import re

######## Edit these variables ########

# The hits table that you exported from XCMS
# Recommended: click the "hide isotope peaks" button on the XCMS webpage before exporting!
xcms_hits_table = "220308 Correlation Human Mouse RTs.csv"

# data_path is the folder the data files are in
data_path = os.getcwd()

# figure_path is where to save the figures
# The folder will be created if it doesn't exist already
figure_path = "Figures"

# DEFINE YOUR DATA SETS HERE
# Data groups are lists of filenames for the mzdata files you exported from Agilent's Qualitative Analysis
# They should be the same thing you fed to XCMS
# Each group can have a different color
# data_group_1 will be at the bottom of the EIC plot
data_group_1 = ["Con_24hr_1.mzdata", "Con_24hr_4.mzdata", "Con_24hr_7.mzdata", "Con_24hr_10.mzdata", "Con_24hr_13.mzdata", "Con_24hr_16.mzdata", "Con_24hr_19.mzdata", "Con_24hr_22.mzdata"]
data_group_2 = ["HK_24hr_2.mzdata", "HK_24hr_5.mzdata", "HK_24hr_8.mzdata", "HK_24hr_11.mzdata", "HK_24hr_14.mzdata", "HK_24hr_17.mzdata", "HK_24hr_20.mzdata", "HK_24hr_23.mzdata"]
data_group_3 = ["Live_24hr_3.mzdata", "Live_24hr_6.mzdata", "Live_24hr_12.mzdata","Live_24hr_24.mzdata"]
# data_group_4 etc...

# RGB values (0-255) that define the colors for each dataset as defined above
# Add or remove the groups in parentheses as needed -- make sure you have one per data_group_#
data_group_colors = [(147, 161, 173), (165, 28, 48), (77, 184, 72)]

# Name of each data group to display in the legend
data_group_names = ["Con", "HK", "Live"]

# ppm tolerance for the extracted ion chromatograms
# this is the maximum distance from the mass reported by XCMS, so 10 ppm makes for a 20 ppm window
eic_ppm_tolerance = 10

# width (in minutes) of the x-axis to plot
rt_window = 8

# distance the EICs are offset from each other (in points, 1 point = 1/72 inches)
eic_x_offset = 4
eic_y_offset = 4

# size of the figure in (width, height)
figure_size = (6, 3)

# resolution of the figure in dots per inch
figure_dpi = 300

# file format to save the figures as -- can be anything matplotlib supports, you probably want "png" or "pdf"
figure_file_format = "png"

######################################

seaborn.set_style("white")
seaborn.set_style("ticks")
seaborn.set_context("paper")

mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.titlesize'] = 12
mpl.rcParams['figure.titlesize'] = 12

mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['ytick.major.size'] = 4
mpl.rcParams['xtick.major.pad'] = 1
mpl.rcParams['ytick.major.pad'] = 1

mpl.rcParams['axes.edgecolor'] = 'k'
mpl.rcParams['axes.labelcolor'] = 'k'
mpl.rcParams['xtick.color'] = 'k'
mpl.rcParams['ytick.color'] = 'k'

mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['ytick.major.width'] = 0.5

mpl.rcParams['lines.markersize'] = 2
mpl.rcParams['axes.axisbelow'] = False
mpl.rcParams['savefig.pad_inches'] = 0

mpl.rcParams['figure.subplot.bottom'] = 0
mpl.rcParams['figure.subplot.hspace'] = 0
mpl.rcParams['figure.subplot.left'] = 0
mpl.rcParams['figure.subplot.right'] = 1
mpl.rcParams['figure.subplot.top'] = 1
mpl.rcParams['figure.subplot.wspace'] = 0

mpl.rcParams['font.sans-serif'] = ['Arial']
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['path.snap'] = False
mpl.rcParams['path.simplify'] = False
mpl.rcParams['mathtext.default'] = 'rm'

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['text.color'] = 'k'

# Returns the sum of intensities in a given spectrum within the specified ppm limit
def getIonSum(mz, ppm, spectrum):
    tolerance = mz * ppm / 1.e6
    return sum([peak.getIntensity() for peak in spectrum if peak.getMZ() >= mz - tolerance and peak.getMZ() <= mz + tolerance])

# Returns an extracted ion chromatogram as a tuple of (rt, intensity)
def getEIC(mz, ppm, experiment):
    return [(spectrum.getRT() / 60, getIonSum(mz, ppm, spectrum)) for spectrum in experiment]

# Returns the pyopenms spectrum that most closely matches the provided retention time
def getSpectrum(rt, experiment):
    spectrum1 = None
    spectrum2 = None

    for spectrum in experiment:
        if spectrum.getRT() / 60 > rt:
            spectrum2 = spectrum
            break
        else:
            spectrum1 = spectrum

    if spectrum1 is not None and spectrum2 is not None:
        d1 = rt - spectrum1.getRT() / 60
        d2 = rt - spectrum2.getRT() / 60

        if abs(d1) < abs(d2):
            return spectrum1
        else:
            return spectrum2

    elif spectrum2 is not None:
        return spectrum2

    else:
        return spectrum1

# Plots a set of EICs and associated spectrum side-by-side
def graph(datasets, colors, rt, spectrum_data, spectrum_color, mz):
    # Snap the retention time on the x-axis to the nearest 0.5 min
    rt = round(rt * 2) / 2
    rtmin = rt - rt_window / 2
    rtmax = rt + rt_window / 2

    fig, axs = pyplot.subplots(nrows = 1, ncols = 2, figsize = figure_size, dpi = figure_dpi)
    eic_ax = axs[0]
    spectrum_ax = axs[1]

    #### EXTRACTED ION CHROMATOGRAMS ####

    eic_ax.set_xlim([rtmin, rtmax])
    eic_ax.set_xlabel('Retention Time (min)')
    eic_ax.set_xticks(np.linspace(rtmin, rtmax, 5))

    plotmaxima = list()

    # This is the height (in inches) of the y-axis so we can fudge the scale later to accomodate the diagonal stack of EICs
    axis_height = eic_ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted()).height

    lines = []
    legend_lines = []

    # Plot all of the lines
    # Figure out which one is the highest on the screen
    # Rescale the plot so that the ylim matches the figure-space maximum
    data_index = 0
    for i, datagroup in enumerate(datasets):
        for d in datagroup:
            # Crop the datarange to plot to just within the retention time window we care about
            coords = [(x,y) for x,y in d if x >= rtmin and x <= rtmax]

            dx = (data_index + 1) * eic_x_offset / 72.
            dy = (data_index + 1) * eic_y_offset / 72.
            offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
            finaltransform = eic_ax.transData + offset
            line, = eic_ax.plot([x for x,y in coords], [y for x,y in coords], zorder=-data_index, color = colors[i], transform = finaltransform)
            line.set_clip_on(False)
            
            # Store the maximum scaled for the offset from the y-axis
            xmax,ymax = max(d, key = lambda xy: xy[1])
            scaled_ymax = ymax * (1 + dy/axis_height) * 1.1
            plotmaxima.append(scaled_ymax)
            
            lines.append(line)
            data_index = data_index + 1

        legend_lines.append(lines[-1])

    eic_ax.set_ylim([0, max(plotmaxima)])

    # Adjust ytick labels to be 1-10 with an exponent in the y-axis label
    yticks = eic_ax.get_yticks()
    tickpower = math.floor(math.log10(yticks[-1]))
    tickdivisor = math.pow(10, tickpower)
    yticksnormalized = [y / tickdivisor for y in yticks]
    eic_ax.set_yticks(yticks)
    eic_ax.set_yticklabels(yticksnormalized)
    eic_ax.set_ylabel(r'Intensity / 10$^\mathregular{' + str(tickpower) + r'}$')

    # Legend
    eic_ax.legend(reversed(legend_lines), reversed(data_group_names), loc = 'upper left', frameon = False, handlelength = 1.5, fontsize = 8)
    
    #### SPECTRUM ####
    spectrum_ax.set_xlabel("m/z")
    spectrum_ax.set_yticks([])
    spectrum_ax.vlines([x for x,y in spectrum_data], [0 for x,y in spectrum_data], [y for x,y in spectrum_data], color = spectrum_color)
    spectrum_ax.set_ylim(bottom = 0)

    min_x = math.floor(mz - 2)
    max_x = math.ceil(mz + 4)
    
    spectrum_ax.set_xlim([min_x, max_x])
    spectrum_ax.set_xticks(range(min_x, max_x + 1, 1))

    # Figure title
    fig.suptitle("m/z = " + str("%0.4f" % mz))

    seaborn.despine()
    pyplot.tight_layout()
    
    return fig

experiment_groups = list()
all_globals = globals()

datasets = [v for v in all_globals if re.search('data_group_\d+', v)]
for dataset in datasets:
    experiment_group = list()

    for d in all_globals[dataset]:
        experiment = MSExperiment()
        MzDataFile().load(os.path.join(data_path, d), experiment)
        experiment_group.append(experiment)

    experiment_groups.append(experiment_group)

# Load and sort the results
if not os.path.exists(figure_path):
    os.makedirs(figure_path)

results = pandas.read_csv(xcms_hits_table)
results = results.sort_values("p", ascending = True)

print(results)

# Normalize colors for 0-1 range instead of 0-255
data_group_colors = [tuple([c/255. for c in color]) for color in data_group_colors]

# Empirical conversion for human data to equivalent rt in the Balskus lab
def rt_to_balskus_rt(rt):
    return  -0.0252 * rt * rt + 1.6266 * rt + 0.708

for i, row in results.iterrows():
    if os.path.exists(os.path.join(figure_path, str(row['Compound_ID']) + ".png")):
        continue

    rt = rt_to_balskus_rt(row['RT'])
    mz = float(row['MZ'])

    datasets = list()
    datasets_maxima = list()

    for group in experiment_groups:
        data = [getEIC(mz, 10, e) for e in group]
        datasets.append(data)

    # Data trimmed to narrow window around rt
    # (to find the dataset with maximum intensity)
    for i, dataset in enumerate(datasets):
        rt_data = [[(x,y) for x,y in d if x >= rt - 3 and x <= rt + 3] for d in dataset]

        # Find whichever dataset has the highest intensity EIC near the retention time    
        maxima = [max(xy, key = lambda xy: xy[1]) for xy in rt_data]
        m = max(maxima, key = lambda xy: xy[1])
        mi = maxima.index(m)
        datasets_maxima.append((m, i, experiment_groups[i][mi]))

    max_signal_xy, max_signal_group, max_signal_experiment = max(datasets_maxima, key = lambda d: d[0])
    # Get the spectrum at the data file that had the highest signal
    spectrum = getSpectrum(max_signal_xy[0], max_signal_experiment)
    spectrum_xy = [(peak.getMZ(), peak.getIntensity()) for peak in spectrum if peak.getMZ() >= mz - 4 and peak.getMZ() <= mz + 6]

    fig = graph(datasets, data_group_colors, rt, spectrum_xy, data_group_colors[max_signal_group], mz)
    fig.savefig(os.path.join(figure_path, str(row['Compound_ID']) + ".png"))