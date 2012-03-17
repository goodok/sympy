# -*- coding: utf-8 -*-
from sympy.core import (Basic, Expr, Add, Mul, Pow)
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.singleton import (Singleton, S)
from sympy.core import (Pow)
from sympy.functions.combinatorial.factorials import factorial
from sympy.core.sets import Interval
from sympy.core.cache import cacheit

from sympy.abc import k

from seriesexpr import SeriesExpr, SeriesAdd, SeriesMul, SeriesCoeffMul, SeriesAtom, SeriesNested
from sequences import Sequence
from sequencesexpr import SeqAdd, SeqCauchyMul, SeqCauchyPow, FaDeBruno, SeqMulEW

class PowerSeriesExprOp(SeriesExpr):
    is_PowerSeries = True

    def __neg__(self):
        return PowerSeriesCoeffMul(S.NegativeOne, self)
    def __abs__(self):
        raise NotImplementedError

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return PowerSeriesAdd(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return PowerSeriesAdd(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return PowerSeriesAdd(self, -other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return PowerSeriesAdd(other, -self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return PowerSeriesMul(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return PowerSeriesMul(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        # if other == -S.One: return PowerSeriesInverse(self)
        return PowerSeriesPow(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise NotImplementedError("Sequence Power not defined")
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return PowerSeriesMul(self, other**S.NegativeOne)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        raise NotImplementedError()

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

class PowerSeriesExpr(PowerSeriesExprOp):

    @cacheit
    def getitem_index(self, i):
        c = self.sequence
        return c[i]*Pow(self.x, i)

    def compose(self, other):
        return PowerSeriesNested(self, other)

    # abstract
    def to_taylor_series(self):
        pass


class PowerSeries(PowerSeriesExpr, SeriesAtom):
    """
    Examples:

    >>> from sympy import S, oo
    >>> from sympy.abc import x, k
    >>> from sympy.series.sequences import Sequence
    >>> from sympy.series.power import PowerSeries

    >>> seq = Sequence((1, oo), formula = (k, S(1)/k))
    >>> seq
    SeqFormula([1, oo), k, 1/k)

    >>> PowerSeries(x, sequence=seq)
    x + x**2/2 + x**3/3 + x**4/4 + x**5/5 + x**6/6 + x**7/7 + x**8/8 + ...

    >>> seq = Sequence((0, oo), periodical = (1, 0))
    >>> seq
    SeqPer([0, oo), (1, 0))

    >>> PowerSeries(x, sequence=seq)
    1 + x**2 + x**4 + x**6 + x**8 + ...

    """

    def __getitem__(self, i):
        if isinstance(i, slice):
            # new seqence
            sequence = self.sequence[i]
            return PowerSeries(self.x, sequence=sequence)
        else:
            return self.getitem_index(i)

class PowerSeriesAdd(PowerSeriesExpr, SeriesAdd):
    """
    Summation of power series.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power import PowerSeries
    >>> from sympy.series.sequences import Sequence

    >>> a = PowerSeries(x, sequence=Sequence((0, oo), 'a'))
    >>> b = PowerSeries(x, sequence=Sequence((0, oo), 'b'))
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
    def __new__(cls, *args):

        #TODO: is it correct, to check arg!=0? args must be Expr type
        args = [arg for arg in args if arg!=0]

        #TODO: create ScalAdd to keep scalar separatly
        if not all(arg.is_PowerSeries for arg in args):
            raise ValueError("Mix of PowerSeries and Scalar symbols")

        expr = Add.__new__(cls, *args)

        if expr.is_Mul:
            return PowerSeriesMul(*expr.args)
        return expr

    @property
    @cacheit
    def sequence(self):
        return SeqAdd(*(s.sequence for s in self.args))

class PowerSeriesMul(PowerSeriesExpr, SeriesMul):
    """
    A product of power series Expressions.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power import PowerSeries
    >>> from sympy.series.sequences import Sequence

    >>> a = PowerSeries(x, sequence=Sequence((0, oo), 'a'))
    >>> b = PowerSeries(x, sequence=Sequence((0, oo), 'b'))
    >>> c = a*b

    >>> c[0]
    a[0]*b[0]

    >>> c[1]
    x*(a[0]*b[1] + a[1]*b[0])

    >>> c[2]
    x**2*(a[0]*b[2] + a[1]*b[1] + a[2]*b[0])

    """

    def __new__(cls, *args):

        if cls.check_zero(args):
            return S.Zero

        coeff, sers = cls.carry_out_coeff(args)

        # if only one seqs then return it
        if len(sers)==1:
            res = sers[0]
        else:
            # form product
            nc, sers, order_symbol = cls.flatten(sers)
            res = Mul.__new__(cls, *sers)

        # wrap with coefficient
        if coeff == S.One:
            return res
        else:
            return PowerSeriesCoeffMul(coeff, res)

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
    >>> from sympy.series.sequences import Sequence

    >>> a = PowerSeries(x, sequence=Sequence((0, oo), 'a'))
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

    def __getitem__(self, i):
        if isinstance(i, slice):
            return PowerSeriesCoeffMul(self.coefficient, self.ts[i])
        else:
            return self.coefficient * self.series[i]

class PowerSeriesPow(PowerSeriesExpr, Pow):
    """
    Power of formal power series.

    Example
    =======
    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power import PowerSeries
    >>> from sympy.series.sequences import Sequence

    >>> a = PowerSeries(x, sequence=Sequence((0, oo), 'a'))
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

    @property
    @cacheit
    def interval(self):
        return self.sequence.interval

    def __getitem__(self, i):
        return self.getitem_dispatche(i)

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
        return self.g.unfactorialize()

    @property
    def _f(self):
        return self.f.unfactorialize()


    @property
    @cacheit
    def sequence_result(self):
        return FaDeBruno(self._g, self._f).factorialize()

    @cacheit
    def __getitem__(self, i):
        return self.sequence_result[i]



class Reverse(PowerSeriesExpr):
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
    pass
