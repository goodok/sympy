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
from sympy.polys.polytools import Poly

from sympy.sequences import Sequence
from sympy.sequences.expr import SeqCauchyMul, SeqCauchyPow, FaDeBruno
from seriesexpr import SeriesExpr, SeriesSliced, SeriesAdd, SeriesMul, SeriesCoeffMul, SeriesAtom, SeriesNested, SeriesGen

from power_0 import PowerSeries0Expr, PowerSeries0, PowerSeries0Mul, PowerSeries0Pow, PowerSeries0Pow, PowerSeries0Nested

from sympy.functions.elementary.trigonometric import sin, cos
from sympy.functions.elementary.exponential import exp

from sympy.tensor.diff_operator import DiffOperatorExpr

class PowerSeriesExprOp(PowerSeries0Expr):
    is_PowerSeries = True
    is_PowerSeries0 = False
    _type_must = "PowerSeries"
    _type_is = "PowerSeries"

    @property
    def _SeriesAdd(self): return PowerSeriesAdd

    @property
    def _SeriesMul(self): return PowerSeriesMul

    @property
    def _SeriesPow(self): return PowerSeriesPow

    @property
    def _SeriesSliced(self): return PowerSeriesSliced

    @property
    def _SeriesExpr(self): return PowerSeriesExpr

    @property
    def _SeriesNested(self): return PowerSeriesNested

    @property
    def _Reverse(self): return Reverse

    @property
    def _SeriesCoeffMul(self): return PowerSeriesCoeffMul

    @classmethod
    def _cls_SeriesCoeffMul(cls): return PowerSeriesCoeffMul

    @classmethod
    def _cls_SeriesMul(cls): return PowerSeriesMul

class PowerSeriesExpr(PowerSeriesExprOp):

    @cacheit
    def getitem_index(self, i):
        c = self.sequence
        return c[i]*Pow(self.x - self.point, i)

    def _print_unevualated_power(self, x, i):
        return Pow(Add(x, -self.point, do_not_sort_in_printing=True, evaluate=False), i)

    @classmethod
    def _convert_from_scalar(cls, scalar, example):
        res = PowerSeries(example.x, sequence=Sequence((0, 0), finitlist=(scalar,)), point = example.point)
        return res

    # abstract
    def to_power_e_series(self):
        pass

    def shift(self, n):
        """
        >>> from sympy.series import PowerSeries
        >>> from sympy.abc import x
        >>> from sympy.printing.repr import srepr
        >>> ps = PowerSeries(x, periodical=(1, 2, 3, 4, 5, 6, 7))
        >>> ps
        1 + 2*x + 3*x**2 + 4*x**3 + 5*x**4 + 6*x**5 + 7*x**6 + x**7 + ...

        >>> ps.shift(-2)
        3 + 4*x + 5*x**2 + 6*x**3 + 7*x**4 + x**5 + 2*x**6 + 3*x**7 + 4*x**8 + ...

        """
        new_seq = self.sequence.shift(n)
        return self._from_args(self.x, new_seq, self.point)

    def __call__(self, arg):
        if isinstance(self.x, DiffOperatorExpr):
            c = self.sequence
            s = 0
            for i in range(c.start_index, 8):
                s += c[i]*self.x**i
            return s(arg)

    def _from_sample(self, sequence):
        return PowerSeries(self.x, sequence=sequence, point=self.point)

    def _apply_abstract_function(self, func_cls):
        """
        Apply abstract function as composition.

        f(PowerSeries())

        """
        f = func_cls.__class__
        seq = Sequence(function=lambda k: f(S.Zero).diff(self.x, k, evaluate=False)).unfactorialize()
        ps = self._from_sample(sequence=seq)
        return ps.compose(self)


    def _apply_function(self, func_cls):
        """
        Apply abstract function as composition.

        sin(PowerSeries())
        """
        ps = None
        name = func_cls.__name__
        point = self.point
        c0 = self.coeff(0)
        if name == "sin":
            if c0 is S.Zero:
                seq = Sequence(periodical=(0, 1, 0, -1)).unfactorialize()
                ps = self._from_sample(sequence=seq)
                return ps.compose(self)
            else:
                ps_zero_free = self[1:]
                return cos(c0)*sin(ps_zero_free) + sin(c0)*cos(ps_zero_free)
        elif name == "cos":
            if c0 is S.Zero:
                seq = Sequence(periodical=(1, 0, -1, 0)).unfactorialize()
                ps = self._from_sample(sequence=seq)
                return ps.compose(self)
            else:
                ps_zero_free = self[1:]
                return cos(c0)*cos(ps_zero_free) - sin(c0)*sin(ps_zero_free)

        elif name == "exp":
            if c0 is S.Zero:
                seq = Sequence(periodical=(1,)).unfactorialize()
                ps = self._from_sample(sequence=seq)
                return ps.compose(self)
            else:
                ps_zero_free = self[1:]
                return exp(c0)*func_cls(ps_zero_free)

class PowerSeriesGen(PowerSeriesExpr, SeriesGen):
    is_SeriesGen = True
    def __new__(cls, x, **kwargs):
        # TODO: redefine
        point = sympify(kwargs.get("point", S.Zero))
        res = PowerSeries(x, sequence=Sequence((1, 1), finitlist=(1,)), point=point)
        res.is_SeriesGen = True
        return res

class PowerSeriesSliced(PowerSeriesExpr, SeriesSliced):
    @property
    def point(self):
        return self.original.point

class PowerSeries(PowerSeriesExpr, SeriesAtom):
    """
    Examples:

    >>> from sympy import S, oo
    >>> from sympy.abc import x, k, a
    >>> from sympy.sequences import Sequence
    >>> from sympy.series.power import PowerSeries

    >>> seq = Sequence((1, oo), formula = (k, S(1)/k))
    >>> seq
    SeqFormula([1, oo), k, 1/k)

    >>> PowerSeries(x, sequence=seq, point=1)
    x - 1 + (x - 1)**2/2 + (x - 1)**3/3 + (x - 1)**4/4 + (x - 1)**5/5 + ...

    >>> PowerSeries(x, periodical = (1, -1), point=a)
     1 + a - x + (x - a)**2 - (x - a)**3 + (x - a)**4 - (x - a)**5 + ...


    Notes
    =====

    Defining through the sequences is similar to Generating Function definition
    and Discrete Laplace Tranform.
    """
    def __new__(cls, x=None, sequence_name=None, **kwargs):
        if sequence_name:
            sequence = Sequence(sequence_name, **kwargs)
        else:
            poly = kwargs.pop("poly", None)
            if poly:
                return cls._from_poly(poly, **kwargs) # recur
            else:
                sequence = kwargs.get("sequence", None)
                if sequence==None:
                    sequence = Sequence(**kwargs)
        point =  kwargs.get("point", S.Zero)
        point = sympify(point)
        assert x
        assert sequence
        obj = SeriesExpr.__new__(cls, x, sequence, point)
        return obj

    @property
    def point(self):
        return self._args[2]

    @classmethod
    def _from_args(cls, x, sequence, point):
        return cls.__new__(cls, x, sequence=sequence, point=point)

    @classmethod
    def _from_poly(cls, poly, **kwargs):
        assert poly.is_Poly
        assert poly.is_univariate
        x = poly.gen

        point =  kwargs.get("point", S.Zero)
        if point is not S.Zero:
            poly = Poly(poly.as_expr().subs(x, x + point), x)

        end = poly.degree()
        start = poly.monoms()[-1][0]
        coeffs = poly.all_coeffs()
        coeffs.reverse()
        coeffs = tuple(coeffs[start:])
        sequence = Sequence(Interval(start, end), finitlist=coeffs)
        return cls.__new__(cls, x, sequence=sequence, **kwargs)

class PowerSeriesAdd(PowerSeriesExpr, SeriesAdd):
    """
    Summation of power series.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power import PowerSeries

    >>> a = PowerSeries(x, 'a', point=1)
    >>> b = PowerSeries(x, 'b', point=1)
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


class PowerSeriesMul(PowerSeriesExpr, PowerSeries0Mul):
    """
    A product of power series Expressions.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power import PowerSeries
    >>> from sympy.sequences import Sequence

    >>> a = PowerSeries(x, 'a', point=2)
    >>> b = PowerSeries(x, 'b', point=2)
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


class PowerSeriesCoeffMul(PowerSeriesExpr, SeriesCoeffMul):
    """
    Multiplication power series by scalar (common coefficient).

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power import PowerSeries
    >>> from sympy.sequences import Sequence

    >>> a = PowerSeries(x, 'a', point=3)
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


class PowerSeriesPow(PowerSeriesExpr, PowerSeries0Pow):
    """
    Power of formal power series.

    Example
    =======
    >>> from sympy import oo
    >>> from sympy.abc import x, z
    >>> from sympy.series.power import PowerSeries

    >>> a = PowerSeries(x, 'a', point=z)
    >>> c = a**2
    >>> c
    a[0]**2 + 2*(x - z)*a[0]*a[1] + (x - z)**2*(2*a[0]*a[2] + a[1]**2) + ...

    >>> c[4]
    (x - z)**4*(2*a[0]*a[4] + 2*a[1]*a[3] + a[2]**2)

    """
    @property
    def point(self):
        return self.base.point


class PowerSeriesNested(PowerSeries0Nested, PowerSeries):
    #TODO: consider what's going on if points defers.
    @property
    def point(self):
        return self.b.point

class Reverse(PowerSeries, PowerSeries0):
    """
    Reversion of power series.
    """
    @property
    def point(self):
        return self.original.point

