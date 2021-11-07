"""
Figure properties for visualizing time-course simulations and sensitivity analysis results.

Example code: Nakakuki_Cell_2010_

.. _Nakakuki_Cell_2010: https://github.com/biomass-dev/biomass/blob/master/biomass/models/Nakakuki_Cell_2010/viz.py
"""

from typing import List, Optional, Union

from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

__all__ = ["SingleObservable", "MultipleObservables", "SensitivityOptions"]


class SingleObservable(object):
    """
    Visualization options for time-course simulation (single-observable).

    Attributes
    ----------
    divided_by : int or float (default: 1)
        Convert time unit. (e.g. sec -> min).

    figsize : tuple (default: (4, 3))
        Width, height in inches.

    xlim : tuple
        Set the x limits of the current axes.

    xticks : list (default: :obj:`None`)
        Set the current tick locations of the x-axis.

    xlabel : str (default: 'Time')
        Set the label for the x-axis.

    ylim : tuple
        Set the y limits of the current axes.

    yticks : list (default: :obj:`None`)
        Set the current tick locations of the y-axis.

    ylabel : str (default: obs_name.replace('_', ' '))
        Set the label for the y-axis.

    exp_data : bool (default: True)
        if False, experimental data will not be shown.

    legend_kws : dict, optional
        Keyword arguments to pass to ``matplotlib.pyplot.legend()``.

    cmap : list or tuple
        Set colormap.

    shape : list or tuple of strings (default: ``matplotlib.lines.Line2D.filled_markers``)
        Set markers.

    dont_show : list of strings
        Set conditions you don't want to plot.
    """

    def __init__(self, cm: ListedColormap, obs_name: str) -> None:
        self.divided_by: int = 1
        self.figsize: tuple = (4, 3)
        self.xlim: tuple = ()
        self.xticks: Optional[list] = None
        self.xlabel: str = "Time"
        self.ylim: tuple = ()
        self.yticks: Optional[list] = None
        self.ylabel: str = obs_name.replace("__", "\n").replace("_", " ")
        self.exp_data: bool = True
        self.legend_kws: Optional[dict] = None
        self.cmap: list = [cm.colors[j] for j in range(10)]
        self.shape: Union[list, tuple] = Line2D.filled_markers
        self.dont_show: List[str] = []


class MultipleObservables(object):
    """
    Visualization options for time-course simulation (multi-observables).

    Attributes
    ----------
    fname : str
        Output file name.

    figsize : tuple (default: (4, 3))
        Figure size.

    observables : list of strings
        Specify which observables to plot.

    condition : str, optional
        Specify which simulation condition to be plotted.

    xlim : tuple
        Set the x limits of the current axes.

    xticks : list (default: :obj:`None`)
        Set the current tick locations of the x-axis.

    xlabel : str (default: 'Time')
        Set the label for the x-axis.

    ylim : tuple
        Set the y limits of the current axes.

    yticks : list (default: :obj:`None`)
        Set the current tick locations of the y-axis.

    ylabel : str (default: '')
        Set the label for the y-axis.

    legend_kws : dict, optional
        Keyword arguments to pass to ``matplotlib.pyplot.legend()``.

    cmap : list or tuple
        Set colormap.

    shape : list or tuple of strings (default: ``matplotlib.lines.Line2D.filled_markers``)
        Set markers.
    """

    def __init__(self, cm: ListedColormap) -> None:
        self.fname: str = "multi_observables"
        self.figsize: tuple = (4, 3)
        self.observables: List[str] = []
        self.condition: Optional[str] = None
        self.xlim: tuple = ()
        self.xticks: Optional[list] = None
        self.xlabel: str = "Time"
        self.ylim: tuple = ()
        self.yticks: Optional[list] = None
        self.ylabel: str = ""
        self.legend_kws: Optional[dict] = dict(
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            borderaxespad=0,
            labelspacing=1.25,
            frameon=False,
            fontsize=12,
        )
        self.cmap: list = [cm.colors[j] for j in range(10)]
        self.shape: Union[list, tuple] = Line2D.filled_markers


class SensitivityOptions(object):
    """
    Visualization options for sensitivity analysis results.

    Attributes
    ----------
    figsize : tuple
        Figure size.

    width : float
        The width(s) of the bars when choosing `style` == 'barplot'.

    legend_kws : dict, optional
        Keyword arguments to pass to ``matplotlib.pyplot.legend()``.

    cmap : list or tuple
        Set colormap.
    """

    def __init__(self, cm: ListedColormap) -> None:
        self.figsize: tuple = (12, 5)
        self.width: float = 0.3
        self.legend_kws: Optional[dict] = dict(loc="upper left", frameon=False)
        self.cmap: list = [cm.colors[j] for j in range(10)]
