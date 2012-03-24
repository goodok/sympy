# -*- coding: utf-8 -*-
"""
Formal power series near zero.

It is the same as PowerSeries (power.py) but with internal bases in form
    a_n*x**n/n!
"""

from sympy.core import (Basic, Expr, Add, Mul, Pow)
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.singleton import (Singleton, S)
from sympy.core import (Pow)
from sympy.functions.combinatorial.factorials import factorial
from sympy.core.cache import cacheit
from sympy.core.sets import Interval

from sympy.sequences.expr import SeqExpCauchyMul, SeqExpCauchyPow, FaDeBruno

from seriesexpr import SeriesExpr, SeriesSliced, SeriesAdd, SeriesMul, SeriesCoeffMul, SeriesAtom, SeriesNested
from power import Reverse as _Reverse


class PowerESeriesExprOp(SeriesExpr):
    is_PowerESeries = True

    _type_must = "PowerESeries"
    _type_is = "PowerESeries"

    @property
    def _SeriesAdd(self): return PowerESeriesAdd

    @property
    def _SeriesMul(self): return PowerESeriesMul

    @property
    def _SeriesPow(self): return PowerESeriesPow

    @property
    def _SeriesSliced(self): return PowerESeriesSliced

    @property
    def _SeriesExpr(self): return PowerESeriesExpr

    @property
    def _SeriesNested(self): return PowerESeriesNested

    @property
    def _Reverse(self): return Reverse

    @property
    def _SeriesMul(self): return PowerESeriesMul

    @property
    def _SeriesCoeffMul(self): return PowerESeriesCoeffMul

    @classmethod
    def _cls_SeriesCoeffMul(cls): return PowerESeriesCoeffMul

    @classmethod
    def _cls_SeriesMul(cls): return PowerESeriesMul



class PowerESeriesSliced(SeriesSliced, PowerESeriesExprOp):
    pass

class PowerESeriesExpr(PowerESeriesExprOp):

    @cacheit
    def getitem_index(self, i):
        a =  self.sequence[i]
        if (a != S.Zero) and (i != 0):
            a = a / factorial(i) * Pow(self.x, i)
        return a

    def shift(self, n):
        """
        >>> from sympy.series import PowerESeries
        >>> from sympy.abc import x
        >>> from sympy.printing.repr import srepr
        >>> ps = PowerESeries(x, periodical=(1, 2, 3, 4, 5, 6, 7))
        >>> ps
         1 + 2*x + 3*x**2/2 + 2*x**3/3 + 5*x**4/24 + x**5/20 + ...
         >>> srepr(ps)
         "PowerESeries(Symbol('x'), SeqPer(Interval(Integer(0), oo, False, True), (1, 2, 3, 4, 5, 6, 7)))"

        >>> ps << 2
        3/2 + 2*x/3 + 5*x**2/24 + x**3/20 + ...
        """
        new_seq = self.sequence.shift_exp(n)
        return self._from_args(self.x, new_seq)


    def compose(self, other):
        return PowerESeriesNested(self, other)

    # abstract (will be sertted in seriesupdate module)
    def to_power_series(self):
        pass

    def reverse(self):
        return Reverse(self)

class PowerESeries(PowerESeriesExpr, SeriesAtom):
    """
    Formal Power series with exponentional basis.

    Examples
    ========

    >>> from sympy.series import PowerESeries
    >>> from sympy.sequences import Sequence
    >>> from sympy import S, oo
    >>> from sympy.abc import x, k
    >>> seq = Sequence((1, oo), formula = (k, S(1)/k))
    >>> seq
    SeqFormula([1, oo), k, 1/k)

    >>> PowerESeries(x, sequence=seq)
    x + x**2/4 + x**3/18 + x**4/96 + x**5/600 + ...

    Define hyperbolic cos series with the help of sequence:

    >>> tcosh = PowerESeries(x, periodical = (1, 0))
    >>> tcosh
    1 + x**2/2 + x**4/24 + ...

    >>> tsinh = PowerESeries(x, periodical = (0, 1))
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
    pass

class PowerESeriesAdd(PowerESeriesExpr, SeriesAdd):
    """    """
    def __new__(cls, *args):

        #TODO: is it correct, to check arg!=0? args must be Expr type
        args = [arg for arg in args if arg!=0]

        #TODO: create ScalAdd to keep scalar separatly
        if not all(arg.is_PowerESeries for arg in args):
            raise ValueError("Mix of PowerESeries and Scalar symbols")

        expr = Add.__new__(cls, *args)

        if expr.is_Mul:
            # TODO: use _SeriesMul
            return PowerESeriesMul(*expr.args)
        return expr

class PowerESeriesMul(PowerESeriesExpr, SeriesMul):
    """A Product of series Expressions."""

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
            # TODO: use self._SeriesCoeffMul
            return PowerESeriesCoeffMul(coeff, res)

    @property
    @cacheit
    def sequence(self):
        return SeqExpCauchyMul(*(s.sequence for s in self.args))

class PowerESeriesCoeffMul(PowerESeriesExpr, SeriesCoeffMul):
    pass

class PowerESeriesPow(PowerESeriesExpr, Pow):
    """
    Examples
    ========

    >>> from sympy.series import PowerESeries
    >>> from sympy.sequences import Sequence
    >>> from sympy import S, oo
    >>> from sympy.abc import x, k
    >>> seq = Sequence((1, oo), formula = (k, S(1)/k))
    >>> seq
    SeqFormula([1, oo), k, 1/k)

    >>> PowerESeries(x, sequence=seq)
    x + x**2/4 + x**3/18 + x**4/96 + x**5/600 + ...

    Define hyperbolic cos series with the help of sequence:

    >>> tcosh = PowerESeries(x, periodical = (1, 0))
    >>> tcosh
    1 + x**2/2 + x**4/24 + ...

    >>> tsinh = PowerESeries(x, periodical = (0, 1))
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

class PowerESeriesNested(SeriesNested, PowerESeries):

    @property
    @cacheit
    def sequence(self):
        return FaDeBruno(self.g.sequence, self.f.sequence)


class Reverse(PowerESeries, _Reverse):
    # Note: we use power series algorithm
    @property
    def sequence(self):
        return self.original_seq.unfactorialize().reverse().factorialize()
