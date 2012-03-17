
from sympy.core import (Basic, Expr, Add, Mul, Pow)
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.singleton import (Singleton, S)
from sympy.core import (Pow)
from sympy.functions.combinatorial.factorials import factorial
from sympy.core.cache import cacheit
from sympy.core.sets import Interval

from seriesexpr import SeriesExpr, SeriesAdd, SeriesMul, SeriesCoeffMul, SeriesAtom
from sequencesexpr import SeqAdd, SeqExpCauchyMul, SeqExpCauchyPow, SeqExp_FaDeBruno

class TaylorSeriesExprOp(SeriesExpr):
    is_TaylorSeries = True

    def __neg__(self):
        return TaylorSeriesCoeffMul(S.NegativeOne, self)
    def __abs__(self):
        raise NotImplementedError

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return TaylorSeriesAdd(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return TaylorSeriesAdd(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return TaylorSeriesAdd(self, -other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return TaylorSeriesAdd(other, -self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return TaylorSeriesMul(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return TaylorSeriesMul(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        #if other == -S.One:
        #    return Inverse(self)
        return TaylorSeriesPow(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise NotImplementedError("Sequence Power not defined")
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return TaylorSeriesMul(self, other**S.NegativeOne)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        raise NotImplementedError()

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

class TaylorSeriesExpr(TaylorSeriesExprOp):

    @cacheit
    def getitem_index(self, i):
        a =  self.sequence[i]
        if (a != S.Zero) and (i != 0):
            a = a / factorial(i) * Pow(self.x, i)
        return a

class TaylorSeries(TaylorSeriesExpr, SeriesAtom):
    """
    Formal Taylor series.

    Examples
    ========

    >>> from sympy.series import Sequence, TaylorSeries
    >>> from sympy import S, oo
    >>> from sympy.abc import x, k
    >>> seq = Sequence((1, oo), formula = (k, S(1)/k))
    >>> seq
    SeqFormula([1, oo), k, 1/k)

    >>> TaylorSeries(x, sequence=seq)
    x + x**2/4 + x**3/18 + x**4/96 + x**5/600 + ...

    >>> a = Sequence((0, oo), periodical = (1, 0))
    >>> b = Sequence((0, oo), periodical = (0, 1))

    Define hyperbolic cos series with the help of sequence:

    >>> tcosh = TaylorSeries(x, sequence=a)
    >>> tcosh
    1 + x**2/2 + x**4/24 + ...

    >>> tsinh = TaylorSeries(x, sequence=b)
    >>> tsinh
    x + x**3/6 + x**5/120 + x**7/5040 + ...

    >>> tcosh**2
    1 + x**2 + x**4/3 + 2*x**6/45 + x**8/315 + ...

    >>> tsinh**2
    x**2 + x**4/3 + 2*x**6/45 + x**8/315 + ...

    >>> c = tcosh**2 - tsinh**2
    >>> c
    1 + ...

    Slow:
    >>> c[100]      # doctest: +SKIP
    0

    """

    def coeff(self, i):
        return self.sequence[i]

    def __getitem__(self, i):
        return self.getitem_index(i)


class TaylorSeriesAdd(TaylorSeriesExpr, SeriesAdd):
    """    """
    def __new__(cls, *args):

        #TODO: is it correct, to check arg!=0? args must be Expr type
        args = [arg for arg in args if arg!=0]

        #TODO: create ScalAdd to keep scalar separatly
        if not all(arg.is_TaylorSeries for arg in args):
            raise ValueError("Mix of TaylorSeries and Scalar symbols")

        expr = Add.__new__(cls, *args)

        if expr.is_Mul:
            return TaylorSeriesMul(*expr.args)
        return expr

    @property
    @cacheit
    def sequence(self):
        return SeqAdd(*(s.sequence for s in self.args))

class TaylorSeriesMul(TaylorSeriesExpr, SeriesMul):
    """A Product of Sequence Expressions."""

    def __new__(cls, *args):
        # TODO: join with similar code of PowerSeriesMul, SeriesMul

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
            return TaylorSeriesCoeffMul(coeff, res)

    @property
    @cacheit
    def sequence(self):
        return SeqExpCauchyMul(*(s.sequence for s in self.args))

class TaylorSeriesCoeffMul(TaylorSeriesExpr, SeriesCoeffMul):
    # TODO: join with PowerSeries and use class method?
    def __getitem__(self, i):
        if isinstance(i, slice):
            return TaylorSeriesCoeffMul(self.coefficient, self.ts[i])
        else:
            return self.coefficient * self.series[i]

class TaylorSeriesPow(TaylorSeriesExpr, Pow):
    """
    Examples
    ========

    >>> from sympy.series import Sequence, TaylorSeries
    >>> from sympy import S, oo
    >>> from sympy.abc import x, k
    >>> seq = Sequence((1, oo), formula = (k, S(1)/k))
    >>> seq
    SeqFormula([1, oo), k, 1/k)

    >>> TaylorSeries(x, sequence=seq)
    x + x**2/4 + x**3/18 + x**4/96 + x**5/600 + ...

    >>> a = Sequence((0, oo), periodical = (1, 0))
    >>> b = Sequence((0, oo), periodical = (0, 1))

    Define hyperbolic cos series with the help of sequence:

    >>> tcosh = TaylorSeries(x, sequence=a)
    >>> tcosh
    1 + x**2/2 + x**4/24 + ...

    >>> tsinh = TaylorSeries(x, sequence=b)
    >>> tsinh
    x + x**3/6 + x**5/120 + x**7/5040 + ...

    >>> (tcosh + tsinh)**1000
    1 + 1000*x + 500000*x**2 + 500000000*x**3/3 + 125000000000*x**4/3 + ...

    """

    @property
    def x(self):
        return self.base.x

    @property
    @cacheit
    def sequence(self):
        return SeqExpCauchyPow(self.base.sequence, self.exp)

    @property
    @cacheit
    def interval(self):
        return self.sequence.interval

    def __getitem__(self, i):
        return self.getitem_dispatche(i)


class TaylorSeriesNested(TaylorSeries):
    def __new__(cls, *args):
        expr = TaylorSeriesExpr.__new__(cls, *args)
        return expr

    @property
    def g(self):
        return self.args[0]

    @property
    def f(self):
        return self.args[1]

    @property
    def x(self):
        return self.g.x

    @property
    @cacheit
    def sequence(self):
        return SeqExp_FaDeBruno(self.f.sequence, self.g.sequence)

    def __getitem__(self, i):
        return self.getitem_dispatche(i)

    @property
    @cacheit
    def interval(self):
        return self.f.interval
