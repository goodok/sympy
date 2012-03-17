"""
Maintain converting from one type of series to anather.


This is ceseparated file to avoid conflicts bitween modules.


        seriesexpr.py (common)
      /           \
     /             \
taylor.py       power.py
    \             /
     \           /
      \         /
       \       /
        \     /
         \   /
           |
       this.file.
"""

from taylor import TaylorSeries, TaylorSeriesExpr
from power import PowerSeries, PowerSeriesExpr


def to_power_series(self):
    return PowerSeries(self.x, sequence = self.sequence.factorialize())

TaylorSeriesExpr.to_power_series = to_power_series

def to_taylor_series(self):
    return TaylorSeries(self.x, sequence = self.sequence.unfactorialize())

PowerSeriesExpr.to_taylor_series = to_taylor_series


