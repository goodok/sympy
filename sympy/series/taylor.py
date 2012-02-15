
from sympy.core import (Basic, Expr, Add, Mul, Pow)
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.singleton import (Singleton, S)
from sympy.core import (Pow)
from sympy.functions.combinatorial.factorials import factorial
from sympy.core.cache import cacheit
from sympy.core.sets import Interval

from seriesexpr import SeriesExpr, SeriesAdd, SeriesMul, SeriesCoeffMul, SeriesAtom
from sequencesexpr import SeqExpCauchyMul, SeqExpCauchyPow

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
        if other == -S.One:
            return Inverse(self)
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
    pass

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

    def coeff(self, i):
        return self.sequence[i]

class TaylorSeriesMul(TaylorSeriesExpr, SeriesMul):
    """A Product of Sequence Expressions."""

    def __new__(cls, *args):
        if any(arg.is_zero for arg in args):
            return S.Zero

        # collect only series
        series = [arg for arg in args if arg.is_TaylorSeries]

        # collect scalar coefficients
        coeffs = [arg for arg in args if not arg.is_TaylorSeries]

        # calculate the multyplicity of coefficients
        if coeffs==[]:
            coeff = S.One
        else:
            coeff = Mul(*coeffs)

        # if only one seqs then return it
        if len(series)==1:
            if coeff == S.One:
                return series[0]
            else:
                return TaylorSeriesCoeffMul(coeff, series[0])

        # further
        expr = Mul.__new__(cls, *args)
        return expr

    @property
    def x(self):
        return self.args[0].x

    @property
    @cacheit
    def interval(self):
        start = Add(*(s.start_index for s in self.args))
        stop = Add(*(s.stop_index for s in self.args))
        res = Interval(start, stop)
        return res

    @property
    @cacheit
    def sequence(self):
        return SeqExpCauchyMul(*(s.sequence for s in self.args))

    @cacheit
    def __getitem__(self, i):
        if self.is_out_of_range(i):
            return S.Zero
        else:
            c = self.sequence
            return c[i] * Pow(self.x, i)/ factorial(i)



class TaylorSeriesCoeffMul(TaylorSeriesExpr, SeriesCoeffMul):
    def __getitem__(self, i):
        if isinstance(i, slice):
            return TaylorSeriesCoeffMul(self.coeff, self.ts[i])
        else:
            return self.coeff * self.series[i]

class TaylorSeriesPow(TaylorSeriesExpr, Pow):

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

    @cacheit
    def __getitem__(self, i):
        if self.is_out_of_range(i):
            return S.Zero
        else:
            c = self.sequence
            return c[i]*Pow(self.x, i)/factorial(i)


class TaylorSeries(TaylorSeriesExpr, SeriesAtom):
    """
    Examples:

    >>> from sympy.series import Sequence, TaylorSeries
    >>> from sympy import S, oo
    >>> from sympy.abc import x, k
    >>> seq = Sequence((1, oo), formula = (k, S(1)/k))
    >>> seq
    SeqFormula([1, oo), k, 1/k)

    >>> TaylorSeries(x, sequence=seq)
    x + x**2/4 + x**3/18 + x**4/96 + x**5/600 + ...

    >>> seq = Sequence((0, oo), periodical = (1, 0))
    >>> seq
    SeqPer([0, oo), (1, 0))

    >>> TaylorSeries(x, sequence=seq)
    1 + x**2/2 + x**4/24 + ...

    """

    def coeff(self, i):
        return self.sequence[i]

    def __getitem__(self, i):
        a =  self.sequence[i]
        if (a != S.Zero) and (i != 0):
            a = a / factorial(i) * Pow(self.x, i)
        return a

