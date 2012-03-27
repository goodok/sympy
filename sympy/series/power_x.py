# -*- coding: utf-8 -*-
"""
Formal power series near some point.

It is based on:
    a_n*x**n

"""

from sympy.core import (Basic, Expr, Add, Mul, Pow, sympify)
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.singleton import (Singleton, S)
from sympy.core import (Pow)
from sympy.functions.combinatorial.factorials import factorial
from sympy.core.sets import Interval
from sympy.core.cache import cacheit

from sympy.sequences import Sequence
from sympy.sequences.expr import SeqCauchyMul, SeqCauchyPow, FaDeBruno
from seriesexpr import SeriesExpr, SeriesSliced, SeriesAdd, SeriesMul, SeriesCoeffMul, SeriesAtom, SeriesNested

from power import PowerSeriesExpr, PowerSeries, PowerSeriesMul, PowerSeriesPow, PowerSeriesPow, PowerSeriesNested

class PowerSeriesXExprOp(PowerSeriesExpr):
    is_PowerXSeries = True

    _type_must = "PowerSeriesX"
    _type_is = "PowerSeriesX"

    @property
    def _SeriesAdd(self): return PowerSeriesXAdd

    @property
    def _SeriesMul(self): return PowerSeriesXMul

    @property
    def _SeriesPow(self): return PowerSeriesXPow

    @property
    def _SeriesSliced(self): return PowerSeriesXSliced

    @property
    def _SeriesExpr(self): return PowerSeriesXExpr

    @property
    def _SeriesNested(self): return PowerSeriesXNested

    @property
    def _Reverse(self): return Reverse

    @property
    def _SeriesCoeffMul(self): return PowerSeriesXCoeffMul

    @classmethod
    def _cls_SeriesCoeffMul(cls): return PowerSeriesXCoeffMul

    @classmethod
    def _cls_SeriesMul(cls): return PowerSeriesXMul


class PowerSeriesXSliced(PowerSeriesXExprOp, SeriesSliced):
    # to maintain is_PowerSeriesXX
    pass

class PowerSeriesXExpr(PowerSeriesXExprOp):

    @cacheit
    def getitem_index(self, i):
        c = self.sequence
        return c[i]*Pow(self.x - self.point, i)

    # abstract
    def to_power_e_series(self):
        pass


class PowerSeriesX(PowerSeriesXExpr, SeriesAtom):
    """
    Examples:

    >>> from sympy import S, oo
    >>> from sympy.abc import x, k, a
    >>> from sympy.sequences import Sequence
    >>> from sympy.series.power_x import PowerSeriesX

    >>> seq = Sequence((1, oo), formula = (k, S(1)/k))
    >>> seq
    SeqFormula([1, oo), k, 1/k)

    >>> PowerSeriesX(x, sequence=seq, point=1)
    x - 1 + (x - 1)**2/2 + (x - 1)**3/3 + (x - 1)**4/4 + (x - 1)**5/5 + ...

    >>> PowerSeriesX(x, periodical = (1, -1), point=a)
    1 + a - x + (-a + x)**2 - (-a + x)**3 + (-a + x)**4 - (-a + x)**5 + ...


    Notes
    =====

    Defining through the sequences is similar to Generating Function definition
    and Discrete Laplace Tranform.
    """
    def __new__(self, x, sequence_name=None, **kwargs):
        if sequence_name:
            sequence = Sequence(sequence_name, **kwargs)
        else:
            sequence = kwargs.get("sequence", None)
            if sequence==None:
                sequence = Sequence(**kwargs)
        point =  kwargs.get("point", S.Zero)
        point = sympify(point)
        obj = SeriesExpr.__new__(self, x, sequence, point)
        return obj

    @classmethod
    def _from_args(cls, x, sequence, point):
        return cls.__new__(cls, x, sequence=sequence, point=point)

    @property
    def point(self):
        return self._args[2]


class PowerSeriesXAdd(PowerSeriesXExpr, SeriesAdd):
    """
    Summation of power series.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power_x import PowerSeriesX

    >>> a = PowerSeriesX(x, 'a', point=1)
    >>> b = PowerSeriesX(x, 'b', point=1)
    >>> a
    a[0] + (x - 1)*a[1] + (x - 1)**2*a[2] + (x - 1)**3*a[3] + ...
    >>> b
    b[0] + (x - 1)*b[1] + (x - 1)**2*b[2] + (x - 1)**3*b[3] + ...

    >>> a + b
    a[0] + b[0] + (x - 1)*(a[1] + b[1]) + (x - 1)**2*(a[2] + b[2]) + ...

    >>> (a + b)[3]
    (x - 1)**3*(a[3] + b[3])

    >>> (a + b).coeff(3)
    a[3] + b[3]

    """
    #TODO: consider the case of different points
    @property
    def point(self):
        return self._args[0].point


class PowerSeriesXMul(PowerSeriesXExpr, PowerSeriesMul):
    """
    A product of power series Expressions.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power_x import PowerSeriesX
    >>> from sympy.sequences import Sequence

    >>> a = PowerSeriesX(x, 'a', point=2)
    >>> b = PowerSeriesX(x, 'b', point=2)
    >>> c = a*b

    >>> c[0]
    a[0]*b[0]

    >>> c[1]
    (x - 2)*(a[0]*b[1] + a[1]*b[0])

    >>> c[2]
    (x - 2)**2*(a[0]*b[2] + a[1]*b[1] + a[2]*b[0])

    """
    #TODO: consider the case of different points
    @property
    def point(self):
        return self._args[0].point


class PowerSeriesXCoeffMul(PowerSeriesXExpr, SeriesCoeffMul):
    """
    Multiplication power series by scalar (common coefficient).

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power_x import PowerSeriesX
    >>> from sympy.sequences import Sequence

    >>> a = PowerSeriesX(x, 'a', point=3)
    >>> 2*a
    2*a[0] + 2*(x - 3)*a[1] + 2*(x - 3)**2*a[2] + 2*(x - 3)**3*a[3] + ...

    >>> (2*a).coefficient
    2

    >>> (2*a).series
    a[0] + (x - 3)*a[1] + (x - 3)**2*a[2] + (x - 3)**3*a[3] + (x - 3)**4*a[4] + ...

    >>> (2*a).coeff(3)
    2*a[3]

    >>> (2*a)[3]
    2*(x - 3)**3*a[3]

    """
    @property
    def point(self):
        return self.series.point


class PowerSeriesXPow(PowerSeriesXExpr, PowerSeriesPow):
    """
    Power of formal power series.

    Example
    =======
    >>> from sympy import oo
    >>> from sympy.abc import x, z
    >>> from sympy.series.power_x import PowerSeriesX

    >>> a = PowerSeriesX(x, 'a', point=z)
    >>> c = a**2
    >>> c
    a[0]**2 + 2*(x - z)*a[0]*a[1] + (x - z)**2*(2*a[0]*a[2] + a[1]**2) + ...

    >>> c[4]
    (x - z)**4*(2*a[0]*a[4] + 2*a[1]*a[3] + a[2]**2)

    """
    @property
    def point(self):
        return self.base.point


class PowerSeriesXNested(PowerSeriesNested, PowerSeriesX):
    #TODO: consider what's going on if points defers.
    @property
    def point(self):
        return self.f.point

class Reverse(PowerSeriesX, PowerSeries):
    """
    Reversion of power series.
    """
    @property
    def point(self):
        return self.original.point

