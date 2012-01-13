from sympy.core import (Basic, Expr, Add, Mul, Pow, S)
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.cache import cacheit

class TaylorSeriesExprOp(Expr):
    """ TaylorSeries Expression Class"""

    _op_priority = 13.0

    is_TaylorSeries = True
    is_Identity = False
    is_ZeroSequence = False

    show_method='series'
    show_n = 6

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


class TaylorSeriesExpr(TaylorSeriesExprOp):

    def _hashable_content(self):
        return tuple(self._args)

    @property
    @cacheit
    def start_index(self):
        return self.interval._inf

    @property
    @cacheit
    def stop_index(self):
        return self.interval._sup

    @property
    @cacheit
    def is_infinite(self):
        return self.stop_index == S.Infinity

    def is_out_of_range(self, i):
        if isinstance(i, Symbol):
            return False
        if i < self.start_index:
            return True
        if not self.is_infinite:
            if i > self.stop_index:
                return True
        return False

    def show(self, n=5, **kwargs):
        TaylorSeriesExpr.show_method = kwargs.get('method', 'series')
        self.show_n = n
        return self

    def _sympystr(self, printer, *args):
        if self.show_method=='series':
            l = [self[i] for i in range(self.start_index, self.show_n + 1)]
            l = [i for i in l if i != S.Zero]
            l = [printer._print(i) for i in l]
            return " + ". join(l) + " + ... "
        else:
            _args = [args.show(method="expr") for i in args]
            return printer._print_Basic(self, *_args)

class TaylorSeriesAdd(TaylorSeriesExpr, Add):
    """    """

    def __new__(cls, *args):

        #TODO: is it correct, to check arg!=0? args must be Expr type
        args = [arg for arg in args if arg!=0]

        #TODO: create ScalAdd to keep scalar separatly
        if not all(arg.is_TaylorSeries for arg in args):
            raise ValueError("Mix of Sequence and Scalar symbols")

        expr = Add.__new__(cls, *args)

        if expr.is_Mul:
            return TaylorSeriesMul(*expr.args)
        return expr

    def _hashable_content(self):
        return tuple(sorted(self._args, key=hash))

    @classmethod
    def flatten(cls, args_seq):
        return args_seq, [], None

    def as_ordered_terms(self, order=None):
        return self.args

    @property
    def interval(self):
        res = S.EmptySet
        for ts in self.args:
            res = res | ts.interval
        return res

    def __getitem__(self, i):
        if isinstance(i, slice):
            pass
        else:
            return Add(*(ts[i] for ts in self.args))

class TaylorSeriesMul(TaylorSeriesExpr, Mul):
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

        # further - element-wise multiplicity
        expr = Mul.__new__(cls, *args)

        return expr

    def _hashable_content(self):
        return tuple(sorted(self._args, key=hash))

    @classmethod
    def flatten(cls, args_seq):
        return args_seq, [], None

    @property
    def interval(self):
        res = S.EmptySet
        for ts in self.args:
            res = res | ts.interval
        return res


    def __getitem__(self, i):
        return 0

class TaylorSeriesCoeffMul(TaylorSeriesExpr, Mul):
    def __new__(cls, coeff, ts):
        expr = Mul.__new__(cls, coeff, ts)
        return expr

    @classmethod
    def flatten(cls, args_ts):
        coeff = args_ts[0]
        ts = args_ts[1]
        if isinstance(ts, TaylorSeriesCoeffMul):
            coeff *= ts.coeff
            ts = ts.ts
            args_ts = [coeff, ts]
        return args_ts, [], None

    @property
    def coeff(self):
        return self.args[0]

    @property
    def ts(self):
        return self.args[1]

    @property
    def interval(self):
        return self.ts.interval

    def __getitem__(self, i):
        if isinstance(i, slice):
            return TaylorSeriesCoeffMul(self.coeff, self.ts[i])
        else:
            return self.coeff * self.ts[i]

