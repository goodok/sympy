# -*- coding: utf-8 -*-
from sympy.core import (Basic, Expr, Add, Mul, Pow)
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.singleton import (Singleton, S)
from sympy.core import (Pow)
from sympy.functions.combinatorial.factorials import factorial
from sympy.core.sets import Interval
from sympy.core.cache import cacheit

from sympy.abc import k

from seriesexpr import SeriesExpr, SeriesSliced, SeriesAdd, SeriesMul, SeriesCoeffMul, SeriesAtom, SeriesNested
#from sympy.sequences import Sequence
from sympy.sequences.expr import SeqCauchyMul, SeqCauchyPow, FaDeBruno #SeqMulEW

class PowerSeriesExprOp(SeriesExpr):
    is_PowerSeries = True

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


class PowerSeriesSliced(PowerSeriesExprOp, SeriesSliced):
    # to maintain is_PowerSeries
    pass

class PowerSeriesExpr(PowerSeriesExprOp):

    @cacheit
    def getitem_index(self, i):
        c = self.sequence
        return c[i]*Pow(self.x, i)

    # abstract
    def to_power_e_series(self):
        pass


class PowerSeries(PowerSeriesExpr, SeriesAtom):
    """
    Examples:

    >>> from sympy import S, oo
    >>> from sympy.abc import x, k
    >>> from sympy.sequences import Sequence
    >>> from sympy.series.power import PowerSeries

    >>> seq = Sequence((1, oo), formula = (k, S(1)/k))
    >>> seq
    SeqFormula([1, oo), k, 1/k)

    >>> PowerSeries(x, sequence=seq)
    x + x**2/2 + x**3/3 + x**4/4 + x**5/5 + x**6/6 + x**7/7 + x**8/8 + ...

    >>> PowerSeries(x, periodical = (1, 0))
    1 + x**2 + x**4 + x**6 + x**8 + ...


    Notes
    =====

    Defining through the sequences is similar to Generating Function definition
    and Discrete Laplace Tranform.
    """
    pass

class PowerSeriesAdd(PowerSeriesExpr, SeriesAdd):
    """
    Summation of power series.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power import PowerSeries

    >>> a = PowerSeries(x, 'a')
    >>> b = PowerSeries(x, 'b')
    >>> a
    a[0] + x*a[1] + x**2*a[2] + x**3*a[3] + x**4*a[4] + x**5*a[5] + x**6*a[6] + ...
    >>> b
    b[0] + x*b[1] + x**2*b[2] + x**3*b[3] + x**4*b[4] + x**5*b[5] + x**6*b[6] + ...

    >>> a + b
    a[0] + b[0] + x*(a[1] + b[1]) + x**2*(a[2] + b[2]) + ...

    >>> (a + b)[3]
    x**3*(a[3] + b[3])

    >>> (a + b).coeff(3)
    a[3] + b[3]

    """
    pass


class PowerSeriesMul(PowerSeriesExpr, SeriesMul):
    """
    A product of power series Expressions.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power import PowerSeries
    >>> from sympy.sequences import Sequence

    >>> a = PowerSeries(x, 'a')
    >>> b = PowerSeries(x, 'b')
    >>> c = a*b

    >>> c[0]
    a[0]*b[0]

    >>> c[1]
    x*(a[0]*b[1] + a[1]*b[0])

    >>> c[2]
    x**2*(a[0]*b[2] + a[1]*b[1] + a[2]*b[0])

    """
    pass

    @property
    @cacheit
    def sequence(self):
        return SeqCauchyMul(*(s.sequence for s in self.args))

class PowerSeriesCoeffMul(PowerSeriesExpr, SeriesCoeffMul):
    """
    Multiplication power series by scalar (common coefficient).

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power import PowerSeries
    >>> from sympy.sequences import Sequence

    >>> a = PowerSeries(x, 'a')
    >>> 2*a
    2*a[0] + 2*x*a[1] + 2*x**2*a[2] + 2*x**3*a[3] + ...

    >>> (2*a).coefficient
    2

    >>> (2*a).series
    a[0] + x*a[1] + x**2*a[2] + x**3*a[3] + x**4*a[4] + ...

    >>> (2*a).coeff(3)
    2*a[3]

    >>> (2*a)[3]
    2*x**3*a[3]

    """
    pass

class PowerSeriesPow(PowerSeriesExpr, Pow):
    """
    Power of formal power series.

    Example
    =======
    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power import PowerSeries

    >>> a = PowerSeries(x, 'a')
    >>> c = a**2
    >>> c
    a[0]**2 + 2*x*a[0]*a[1] + x**2*(2*a[0]*a[2] + a[1]**2) + ...

    >>> c[4]
    x**4*(2*a[0]*a[4] + 2*a[1]*a[3] + a[2]**2)

    """

    @property
    def x(self):
        return self.base.x

    @property
    @cacheit
    def sequence(self):
        return SeqCauchyPow(self.base.sequence, self.exp)

class PowerSeriesNested(SeriesNested, PowerSeries):

    @property
    @cacheit
    def sequence(self):
        return FaDeBruno_powers(self.f.sequence, self.g.sequence)

class FaDeBruno_powers(FaDeBruno):
    """
    This is similar to SeqExp_FaDeBruno but for formal Power series.

        g(x)=\sum_{n=1}^\infty {b_n} x^n
        f(x)=\sum_{n=1}^\infty {a_n} x^n

    That is without factorials.

    """
    # we use very rough method:
    # substitute f and g sequences with g' = {g_n*n!} and f' = {f_n*n!}
    # and use FaDeBruno
    # then devide result by n!.

    # TODO: implement exponenentionize method for sequences, whis can be used often
    # `SeqMulEW(self.f, Sequence(formula=(k, factorial(k))))`

    @property
    def _g(self):
        return self.g.factorialize()

    @property
    def _f(self):
        return self.f.factorialize()


    @property
    @cacheit
    def sequence_result(self):
        return FaDeBruno(self._g, self._f).unfactorialize()

    @cacheit
    def __getitem__(self, i):
        return self.sequence_result[i]


class Reverse(PowerSeries):
    """
    Reversion of power series.

    References
    ==========

    .. [1] Donald E. "Knuth Art of Computer Programming, Volume 2: Seminumerical Algorithms",
    3rd ed., sec 4.7 "Manipulation of power series", p 526.
    .. [2] http://en.wikipedia.org/wiki/Lagrange_inversion_theorem
    .. [3] Fredrik Johansson, A fast algorithm for reversion of power series
    .. [4] Faà di Bruno's Formula

    """
    def __new__(cls, original):
        assert original.is_Series
        original_seq = original.sequence
        assert original_seq[0] == S.Zero
        assert original_seq[1] <> S.Zero
        obj = SeriesExpr.__new__(cls, original)
        return obj

    @property
    def x(self):
        return self.original.x

    @property
    def original(self):
        return self._args[0]

    @property
    def original_seq(self):
        return self.original.sequence

    @property
    def sequence(self):
        return self.original_seq.reverse()
