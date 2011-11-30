from sympy.core import (Basic, Expr, Add, Mul, Pow, S)
from sympy.core.decorators import _sympifyit, call_highest_priority

class TaylorSeriesExpr(Expr):
    """ TaylorSeries Expression Class"""

    _op_priority = 13.0

    is_TaylorSeries = True
    is_Identity = False
    is_ZeroSequence = False

    def __neg__(self):
        return TaylorSeriesMul(S.NegativeOne, self)
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


class TaylorSeriesAdd(TaylorSeriesExpr, Add):
    """    """

    def __new__(cls, *args):

        args = map(sequenceify, args)

        #TODO: is it correct, to check arg!=0? args must be Expr type
        args = [arg for arg in args if arg!=0]

        if not all(arg.is_TaylorSeries for arg in args):
            raise ValueError("Mix of Sequence and Scalar symbols")

        expr = Add.__new__(cls, *args)

        if expr.is_Mul:
            return TaylorSeriesMul(*expr.args)
        return expr

    def _hashable_content(self):
        return tuple(sorted(self._args, key=hash))

    def _sympystr(self, printer, *args):
        return printer._print_Add(self)

    def as_ordered_terms(self, order=None):
        return self.args



class TaylorSeriesMul(TaylorSeriesExpr, Mul):
    """A Product of Sequence Expressions."""

    def __new__(cls, *args):

        # Check that the shape of the args is consistent
        seqs = [arg for arg in args if arg.is_TaylorSeries]

        if any(arg.is_zero for arg in args):
            return Zero

        if expr.is_Add:
            return TaylorSeriesAdd(*expr.args)
        if expr.is_Pow:
            return TaylorSeriesPow(*expr.args)
        if not expr.is_Mul:
            return expr

        if any(arg.is_TaylorSeries and arg.is_Zero for arg in expr.args):
            return Zero

        return expr

