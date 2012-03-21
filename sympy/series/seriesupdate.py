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

from power import PowerSeries, PowerSeriesExpr
from power_e import PowerESeries, PowerESeriesExpr

def to_power_series(self):
    return PowerSeries(self.x, sequence = self.sequence.unfactorialize())

PowerESeriesExpr.to_power_series = to_power_series

def to_power_e_series(self):
    return PowerESeries(self.x, sequence = self.sequence.factorialize())

PowerSeriesExpr.to_power_e_series = to_power_e_series


