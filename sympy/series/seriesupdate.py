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

from power_0 import PowerSeries0, PowerSeries0Expr
from power_0e import PowerSeries_0E, PowerSeries_0E_Expr

def to_power_series(self):
    return PowerSeries0(self.x, sequence = self.sequence.unfactorialize())

PowerSeries_0E_Expr.to_power_series = to_power_series

def to_power_e_series(self):
    return PowerSeries_0E(self.x, sequence = self.sequence.factorialize())

PowerSeries0Expr.to_power_e_series = to_power_e_series


