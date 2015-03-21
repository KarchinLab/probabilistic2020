"""
The plot.py does the actual calling of plot commands from matplotlib.
Essentially, plot.py encapsulates all the minor tweaks needed in matplotlib
to make a reasonable looking plot.
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats


def barplot(df,
            file_path,
            kind='bar',
            yerr=None,
            xerr=None,
            ecolor=None,
            title='',
            xlabel='',
            ylabel='',
            stacked=False):
    """barplot generates/saves a bar plot from a pandas data frame.

    Parameters
    ----------
    df : pd.DataFrame
        data frame for bar plot
    file_path : str
        path to save bar plot figure
    kind : str, ['bar' | 'barh']
        vertical or horizontal bar chart
    stderror : list
        stdev of each bar
    Matplotlib options for plotting
    """
    if yerr is not None:
        # error bars for y-axis
        df.plot(kind=kind, yerr=yerr, ecolor=ecolor, stacked=stacked)
    elif xerr is not None:
        # error bars for x-axis
        df.plot(kind=kind, xerr=xerr, ecolor=ecolor, stacked=stacked)
    else:
        # normal bar plot
        df.plot(kind=kind, stacked=stacked)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(file_path)
    plt.clf()  # clear figure
    plt.close()


def histogram(df,
              save_path,
              bins=False,
              log=False,
              title='',
              xlabel='',
              ylabel=''):
    """Plots a histogram using matplotlib.

    Parameters
    ----------
    df : pd.DataFrame
        one dimensional data frame or series
    file_path : str
        path to save figure
    bins : list
        bin positions for histogram
    log : Bool
        boolean for log scaling y-axis
    title : str
        title of plot
    xlabel : str
        label on x-axis
    ylabel : str
        label on y-axis
    """
    if bins:
        df.hist(bins=bins, log=log)
    else:
        df.hist(log=log)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.clf()  # clear figure
    plt.close()


def line(data,
         file_path,
         style=[],
         title='',
         xlabel='',
         ylabel='',
         logx=False,
         logy=False,
         vlines=[]):
    """Plots a line plot using matplotlib.

    Parameters
    ----------
    data : pd.DataFrame
        two column df with x and y values
    file_path : str
        path to save figure
    title : str
        graph title
    xlabel : str
        x-axis label
    ylabel : str
        y-axis label
    vlines : list of int
        draw vertical lines at positions
    """
    # plot data
    data.plot(kind='line', style=style)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # log scale if neccessary
    if logx:
        plt.xscale('log')
    if logy:
        plt.yscale('log')

    # plot vertical lines
    ymin, ymax = plt.ylim()  # get plotting range of y-axis
    for l in vlines:
        plt.vlines(l, ymin=ymin,
                   ymax=ymax,
                   color='red')

    plt.tight_layout()  # adjust plot margins
    plt.savefig(file_path)  # save figure
    plt.clf()  # clear figure
    plt.close()


def scatter(x, y,
            file_path,
            colors='',
            size=20,
            title='',
            xlabel='',
            ylabel=''):
    """Create a 2D scatter plot. Many of the optional arguements
    deal with formatting the plot.

    Parameters
    ----------
    x : list|array
        container for x-axis data
    y : list|array
        container for y-axis data
    file_path : str
        path to save figure
    colors : str|list, default=''
        either single color (e.g. 'blue') or a list
    size : int|list, default=20
        int for marker size or a list of ints
    title : str, default=''
        title for plot
    xlabel : str, dfault=''
        x-axis label
    ylabel : str, default=''
        y-axis label
    """
    if colors:
        # plot user's color if supplied
        plt.scatter(x, y, c=colors, s=size)
    else:
        # otherwise rely on default color
        plt.scatter(x, y, s=size)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(file_path)
    plt.clf()  # clear figure
    plt.close()


def errorbars(x, y, err,
              save_path='',
              title='',
              xlabel='',
              ylabel='',
              label=''):
    if label:
        plt.errorbar(x, y, yerr=err, label=label, fmt='-o')
        plt.legend(loc='best')
    else:
        plt.errorbar(x, y, yerr=err, fmt='-o')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if save_path:
        plt.savefig(save_path)
        plt.close()


def correlation_plot(x, y,
                     save_path,
                     title,
                     xlabel, ylabel):
    plt.scatter(x, y)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    line_x = np.arange(x.min(), x.max())
    line_y = slope*line_x + intercept
    plt.plot(line_x, line_y,
             label='$%.2fx + %.2f$, $R^2=%.2f$' % (slope, intercept, r_value**2))
    plt.legend(loc='best')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.clf()  # clear figure
    plt.close()


def boxplot(df, by,
            column,
            save_path,
            xlabel,
            ylabel,
            title):
    # make box plot
    if column:
        axes = df.boxplot(column=column, by=by)
    else:
        axes = df.boxplot(by=by)

    if type(axes) is np.ndarray:
        # multiple boxplot case
        multi_box_flag = True
    else:
        # only one facet for boxplot
        multi_box_flag = False

    # format plot
    if multi_box_flag:
        # plot with multiple box plot facets
        for ax in axes:
            ax.set_xlabel(xlabel)
        axes[0].set_ylabel(ylabel)
        fig = axes[0].get_figure()
        fig.suptitle('')
        plt.tight_layout()
        plt.savefig(save_path)
        plt.clf()
        plt.close()
    else:
        # plot when just 1 box plot facet
        fig = axes.get_figure()
        fig.suptitle('')  # remove auto title from pandas
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.tight_layout()
        plt.savefig(save_path)
        plt.clf()  # clear figure
        plt.close()
