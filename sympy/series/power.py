
from sympy.core import (Basic, Expr, Add, Mul, Pow)
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.singleton import (Singleton, S)
from sympy.core import (Pow)
from sympy.functions.combinatorial.factorials import factorial
from sympy.core.sets import Interval
from sympy.core.cache import cacheit

from seriesexpr import SeriesExpr, SeriesAdd, SeriesMul, SeriesCoeffMul, SeriesAtom
from sequencesexpr import SeqCauchyMul, SeqCauchyPow

class PowerSeriesExpr(SeriesExpr):
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
        if other == -S.One:
            return Inverse(self)
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

class PowerSeriesAdd(PowerSeriesExpr, SeriesAdd):
    """    """
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

class PowerSeriesMul(PowerSeriesExpr, SeriesMul):
    """A Product of Sequence Expressions."""

    def __new__(cls, *args):
        if any(arg.is_zero for arg in args):
            return S.Zero

        # collect only series
        series = [arg for arg in args if arg.is_PowerSeries]

        # collect scalar coefficients
        coeffs = [arg for arg in args if not arg.is_PowerSeries]

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
                return PowerSeriesCoeffMul(coeff, series[0])

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
        return SeqCauchyMul(*(s.sequence for s in self.args))
    
    @cacheit
    def __getitem__(self, i):
        if self.is_out_of_range(i):
            return S.Zero
        else:
            c = self.sequence
            return c[i]*Pow(self.x, i)


class PowerSeriesCoeffMul(PowerSeriesExpr, SeriesCoeffMul):
    def __getitem__(self, i):
        if isinstance(i, slice):
            return PowerSeriesCoeffMul(self.coeff, self.ts[i])
        else:
            return self.coeff * self.series[i]


class PowerSeriesPow(PowerSeriesExpr, Pow):

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

    @cacheit
    def __getitem__(self, i):
        if self.is_out_of_range(i):
            return S.Zero
        else:
            c = self.sequence
            return c[i]*Pow(self.x, i)

class PowerSeries(PowerSeriesExpr, SeriesAtom):
    """
    Examples:

    >>> from sympy.series import Sequence, PowerSeries
    >>> from sympy import S, oo
    >>> from sympy.abc import x, k
    >>> seq = Sequence((1, oo), formula = (k, S(1)/k))
    >>> seq
    SeqFormula([1, oo), k, 1/k)

    >>> PowerSeries(x, sequence=seq)
    x + x**2/4 + x**3/18 + x**4/96 + x**5/600 + ...

    >>> seq = Sequence((0, oo), periodical = (1, 0))
    >>> seq
    SeqPer([0, oo), (1, 0))

    >>> PowerSeries(x, sequence=seq)
    1 + x**2/2 + x**4/24 + ...

    """

    def __getitem__(self, i):
        a =  self.sequence[i]
        if (a != S.Zero) and (i != 0):
            a = a * Pow(self.x, i)
        return a

