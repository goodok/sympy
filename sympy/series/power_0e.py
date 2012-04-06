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

from sympy.sequences import Sequence
from sympy.sequences.expr import SeqExpCauchyMul, SeqExpCauchyPow, FaDeBruno

from seriesexpr import SeriesExpr, SeriesSliced, SeriesAdd, SeriesMul, SeriesCoeffMul, SeriesAtom, SeriesNested
from power_0 import Reverse as _Reverse


class PowerSeries_0E_ExprOp(SeriesExpr):

    _op_priority = 14.0

    is_PowerSeries_0E = True

    _type_must = "PowerSeries_0E"
    _type_is = "PowerSeries_0E"

    @property
    def _SeriesAdd(self): return PowerSeries_0E_Add

    @property
    def _SeriesMul(self): return PowerSeries_0E_Mul

    @property
    def _SeriesPow(self): return PowerSeries_0E_Pow

    @property
    def _SeriesSliced(self): return PowerSeries_0E_Sliced

    @property
    def _SeriesExpr(self): return PowerSeries_0E_Expr

    @property
    def _SeriesNested(self): return PowerSeries_0E_Nested

    @property
    def _Reverse(self): return Reverse

    @property
    def _SeriesMul(self): return PowerSeries_0E_Mul

    @property
    def _SeriesCoeffMul(self): return PowerSeries_0E_CoeffMul

    @classmethod
    def _cls_SeriesCoeffMul(cls): return PowerSeries_0E_CoeffMul

    @classmethod
    def _cls_SeriesMul(cls): return PowerSeries_0E_Mul



class PowerSeries_0E_Expr(PowerSeries_0E_ExprOp):
    @cacheit
    def getitem_index(self, i):
        a =  self.sequence[i]
        if (a != S.Zero) and (i != 0):
            a = a / factorial(i) * Pow(self.x, i)
        return a

    def _print_unevualated_power(self, x, i):
        return S.One/factorial(i) * Pow(x, i)

    def shift(self, n):
        """
        >>> from sympy.series import PowerSeries_0E
        >>> from sympy.abc import x
        >>> from sympy.printing.repr import srepr
        >>> ps = PowerSeries_0E(x, periodical=(1, 2, 3, 4, 5, 6, 7))
        >>> ps
        1 + 2*x + 3*x**2/2 + 2*x**3/3 + 5*x**4/24 + x**5/20 + ...

        >>> ps << 2
        3/2 + 2*x/3 + 5*x**2/24 + x**3/20 + ...
        """
        new_seq = self.sequence.shift_exp(n)
        return self._from_args(self.x, new_seq)


    def compose(self, other):
        return PowerSeries_0E_Nested(self, other)

    # abstract (will be sertted in seriesupdate module)
    def to_power_series(self):
        pass

    def reverse(self):
        return Reverse(self)


class PowerSeries_0E_Sliced(SeriesSliced, PowerSeries_0E_Expr):
    pass


class PowerSeries_0E(PowerSeries_0E_Expr, SeriesAtom):
    """
    Formal Power series with exponentional basis.

    Examples
    ========

    >>> from sympy.series import PowerSeries_0E
    >>> from sympy.sequences import Sequence
    >>> from sympy import S, oo
    >>> from sympy.abc import x, k
    >>> seq = Sequence((1, oo), formula = (k, S(1)/k))
    >>> seq
    SeqFormula([1, oo), k, 1/k)

    >>> PowerSeries_0E(x, sequence=seq)
    x + x**2/4 + x**3/18 + x**4/96 + x**5/600 + ...

    Define hyperbolic cos series with the help of sequence:

    >>> tcosh = PowerSeries_0E(x, periodical = (1, 0))
    >>> tcosh
    1 + x**2/2 + x**4/24 + ...

    >>> tsinh = PowerSeries_0E(x, periodical = (0, 1))
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
    @classmethod
    def _from_poly(cls, poly, **kwargs):
        assert poly.is_Poly
        assert poly.is_univariate
        x = poly.gen
        end = poly.degree()
        start = poly.monoms()[-1][0]
        finitlist = poly.all_coeffs()[:-start]
        finitlist.reverse()
        finitlist = tuple(finitlist)
        sequence = Sequence(Interval(start, end), finitlist=finitlist)
        return cls.__new__(cls, x, sequence=sequence.factorialize())


class PowerSeries_0E_Add(SeriesAdd, PowerSeries_0E_Expr):
    """    """
    _op_priority = 14.0
    def __new__(cls, *args):

        #TODO: is it correct, to check arg!=0? args must be Expr type
        args = [arg for arg in args if arg!=0]

        #TODO: create ScalAdd to keep scalar separatly
        if not all(arg.is_PowerSeries_0E for arg in args):
            raise ValueError("Mix of PowerSeries_0E and Scalar symbols")

        expr = SeriesAdd.__new__(cls, *args)

        if expr.is_Mul:
            # TODO: use _SeriesMul
            return PowerSeries_0E_Mul(*expr.args)
        return expr

class PowerSeries_0E_Mul(PowerSeries_0E_Expr, SeriesMul):
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
            return PowerSeries_0E_CoeffMul(coeff, res)

    @property
    @cacheit
    def sequence(self):
        return SeqExpCauchyMul(*(s.sequence for s in self.args))

class PowerSeries_0E_CoeffMul(PowerSeries_0E_Expr, SeriesCoeffMul):
    pass

class PowerSeries_0E_Pow(PowerSeries_0E_Expr, Pow):
    """
    Examples
    ========

    >>> from sympy.series import PowerSeries_0E
    >>> from sympy.sequences import Sequence
    >>> from sympy import S, oo
    >>> from sympy.abc import x, k
    >>> seq = Sequence((1, oo), formula = (k, S(1)/k))
    >>> seq
    SeqFormula([1, oo), k, 1/k)

    >>> PowerSeries_0E(x, sequence=seq)
    x + x**2/4 + x**3/18 + x**4/96 + x**5/600 + ...

    Define hyperbolic cos series with the help of sequence:

    >>> tcosh = PowerSeries_0E(x, periodical = (1, 0))
    >>> tcosh
    1 + x**2/2 + x**4/24 + ...

    >>> tsinh = PowerSeries_0E(x, periodical = (0, 1))
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

class PowerSeries_0E_Nested(SeriesNested, PowerSeries_0E):

    @property
    @cacheit
    def sequence(self):
        return FaDeBruno(self.a.sequence, self.b.sequence)


class Reverse(PowerSeries_0E, _Reverse):
    # Note: we use power series algorithm
    @property
    def sequence(self):
        return self.original_seq.unfactorialize().reverse().factorialize()
