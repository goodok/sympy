# -*- coding: utf-8 -*-
"""
Formal power series near zero.

It is based on:
    a_n*x**n

"""

from sympy.core import (Basic, Expr, Add, Mul, Pow)
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.singleton import (Singleton, S)
from sympy.core import (Pow)

from sympy.core.function import AppliedUndef
from sympy.abc import k

from sympy.functions.combinatorial.factorials import factorial
from sympy.functions.elementary.trigonometric import sin, cos
from sympy.core.sets import Interval
from sympy.core.cache import cacheit

from seriesexpr import SeriesExpr, SeriesGen, SeriesSliced, SeriesAdd, SeriesMul, SeriesCoeffMul, SeriesPow, SeriesAtom, SeriesNested
from sympy.sequences import Sequence
from sympy.sequences.expr import SeqCauchyMul, SeqCauchyPow, PlainFaDeBruno

class PowerSeries0ExprOp(SeriesExpr):
    is_PowerSeries0 = True

    _type_must = "PowerSeries0"
    _type_is = "PowerSeries0"

    @property
    def _SeriesAdd(self): return PowerSeries0Add

    @property
    def _SeriesMul(self): return PowerSeries0Mul

    @property
    def _SeriesPow(self): return PowerSeries0Pow

    @property
    def _SeriesSliced(self): return PowerSeries0Sliced

    @property
    def _SeriesExpr(self): return PowerSeries0Expr

    @property
    def _SeriesNested(self): return PowerSeries0Nested

    @property
    def _Reverse(self): return Reverse

    @property
    def _SeriesCoeffMul(self): return PowerSeries0CoeffMul

    @classmethod
    def _cls_SeriesCoeffMul(cls): return PowerSeries0CoeffMul

    @classmethod
    def _cls_SeriesMul(cls): return PowerSeries0Mul


class PowerSeries0Expr(PowerSeries0ExprOp):

    @cacheit
    def getitem_index(self, i):
        c = self.sequence
        return c[i]*Pow(self.x, i)

    @classmethod
    def _convert_from_scalar(cls, scalar, example):
        res = PowerSeries0(example.x, sequence=Sequence((0, 0), finitlist=(scalar,)))
        return res


    # abstract
    def to_power_e_series(self):
        pass

    def _from_sample(self, sequence):
        return PowerSeries0(self.x, sequence=sequence)

    def _apply_abstract_function(self, func_cls):
        f = func_cls.__class__
        seq = Sequence(function=lambda k: f(S.Zero).diff(self.x, k, evaluate=False)).unfactorialize()
        ps = PowerSeries0(self.x, sequence=seq)
        return ps.compose(self)

    def _apply_function(self, func_cls):
        ps = None
        name = func_cls.__name__
        if name == "sin":
            c0 = self.coeff(0)
            if c0 is S.Zero:
                seq = Sequence(periodical=(0, 1, 0, -1)).unfactorialize()
                ps = self._from_sample(sequence=seq)
                return ps.compose(self)
            else:
                ps_zero_free = self[1:]
                return cos(c0)*sin(ps_zero_free) + sin(c0)*cos(ps_zero_free)
        elif name == "cos":
            c0 = self.coeff(0)
            if c0 is S.Zero:
                seq = Sequence(periodical=(1, 0, -1, 0)).unfactorialize()
                ps = self._from_sample(sequence=seq)
                return ps.compose(self)
            else:
                ps_zero_free = self[1:]
                return cos(c0)*cos(ps_zero_free) - sin(c0)*sin(ps_zero_free)

        elif name == "exp":
            c0 = self.coeff(0)
            if c0 is S.Zero:
                seq = Sequence(periodical=(1,)).unfactorialize()
                ps = self._from_sample(sequence=seq)
                return ps.compose(self)
            else:
                ps_zero_free = self[1:]
                return exp(c0)*func_cls(ps_zero_free)

class PowerSeries0Sliced(PowerSeries0Expr, SeriesSliced):
    # to maintain is_PowerSeries0
    pass


class PowerSeriesGen0(PowerSeries0Expr, SeriesGen):
    is_SeriesGen = True
    def __new__(cls, x, **kwargs):
        # TODO: redefine
        res = PowerSeries0(x, sequence=Sequence((1, 1), finitlist=(1,)) )
        res.is_SeriesGen = True
        return res

class PowerSeries0(PowerSeries0Expr, SeriesAtom):
    """
    Examples:

    >>> from sympy import S, oo
    >>> from sympy.abc import x, k
    >>> from sympy.sequences import Sequence
    >>> from sympy.series.power_0 import PowerSeries0

    >>> seq = Sequence((1, oo), formula = (k, S(1)/k))
    >>> seq
    SeqFormula([1, oo), k, 1/k)

    >>> PowerSeries0(x, sequence=seq)
    x + x**2/2 + x**3/3 + x**4/4 + x**5/5 + x**6/6 + x**7/7 + x**8/8 + ...

    >>> PowerSeries0(x, periodical = (1, 0))
    1 + x**2 + x**4 + x**6 + x**8 + ...


    Notes
    =====

    Defining through the sequences is similar to Generating Function definition
    and Discrete Laplace Tranform.
    """

class PowerSeries0Add(PowerSeries0Expr, SeriesAdd):
    """
    Summation of power series.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power_0 import PowerSeries0

    >>> a = PowerSeries0(x, 'a')
    >>> b = PowerSeries0(x, 'b')
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


class PowerSeries0Mul(PowerSeries0Expr, SeriesMul):
    """
    A product of power series Expressions.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power_0 import PowerSeries0
    >>> from sympy.sequences import Sequence

    >>> a = PowerSeries0(x, 'a')
    >>> b = PowerSeries0(x, 'b')
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

class PowerSeries0CoeffMul(PowerSeries0Expr, SeriesCoeffMul):
    """
    Multiplication power series by scalar (common coefficient).

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power_0 import PowerSeries0
    >>> from sympy.sequences import Sequence

    >>> a = PowerSeries0(x, 'a')
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

class PowerSeries0Pow(PowerSeries0Expr, SeriesPow):
    """
    Power of formal power series.

    Example
    =======
    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series.power_0 import PowerSeries0

    >>> a = PowerSeries0(x, 'a')
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

class PowerSeries0Nested(SeriesNested, PowerSeries0):

    @property
    @cacheit
    def sequence(self):
        return PlainFaDeBruno(self.a.sequence, self.b.sequence)

class Reverse(PowerSeries0):
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
