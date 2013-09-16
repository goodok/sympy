"""
The module implemented diff operator.
"""
from sympy import Expr, Symbol, Eq, Mul, Add, Pow, expand, sympify, Tuple
from sympy.core.basic import Basic
from sympy.core.operations import AssocOp

from sympy.core.singleton import (Singleton, S)
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.cache import cacheit

from sympy.core.function import diff
from sympy.series.order import Order

class DiffOperatorExpr(Expr):

    _op_priority = 12.0

    is_DiffOperator = True
    is_Identity = False
    is_commutative = True

    # The following is adapted from the core Expr object

    def __neg__(self):
        return DOMul(S.NegativeOne, self)
    def __abs__(self):
        raise NotImplementedError

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return DOAdd(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return DOAdd(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return DOAdd(self, -other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return DOAdd(other, -self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return DOMul(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return DOMul(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return DOPow(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise NotImplementedError("DiffOperator Power is not defined")
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return DOMul(self, other**S.NegativeOne)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        return DOMul(other, self**S.NegativeOne)

class DiffOperatorOne(DiffOperatorExpr):
    """
    Represents the identic operator ( DiffOperator**0).
    """

    def __call__(self, e):
        return e

    def _latex(self, p):
        return "\partial^0"

    __metaclass__ = Singleton

    is_Identity = True

class DiffOperator(DiffOperatorExpr):
    """
    >>> from sympy.abc import x, a
    >>> from sympy.core.function import Function
    >>> from sympy.tensor.diff_operator import DiffOperator
    >>> f = Function("f")

    >>> Dx = DiffOperator(x)


    >>> Dx + Dx
    DiffOperator(x) + DiffOperator(x)

    >>> Dx(f(x))
    Derivative(f(x), x)

    >>> Dx + Dx**2
    DiffOperator(x) + DiffOperator(x)**2

    >>> (Dx + Dx**2)(f(x))
    Derivative(f(x), x) + Derivative(f(x), x, x)

    >>> (a*Dx + Dx**3)(x**4)
    4*a*x**3 + 24*x

    """

    def __new__(cls, variable):
        obj = DiffOperatorExpr.__new__(cls, variable)
        return obj

    @property
    def variable(self):
        return self._args[0]

    def _hashable_content(self):
        return tuple(sorted(self._args, key=hash))

    def __call__(self, e):
        return diff(e, self.variable)

    def _eval_nseries(self, x, n, logx):
        return self + Order(n)

    def _latex(self, p):
        return "\partial_{%s}" % (p._print(self.variable))

    def _needs_brackets(self):
        return False

class DOAdd(DiffOperatorExpr, Add):
    """
    A Sum of the DiffOperator expressions.
    """

    def __new__(cls, *args):

        expr = Add.__new__(cls, *args)

        if expr.is_Mul:
            return DOMul(*expr.args)
        return expr

    def __call__(self, e):
        return Add(*tuple(d(e) for d in self.args))


    @classmethod
    def flatten(cls, args):
        new_seq = []
        i = 0
        while args:
            o = args.pop()
            if o.__class__ is cls:
                args.extend(o.args)
            else:
                new_seq.append(o)
        new_seq.reverse()
        return new_seq, [], None


class DOMul(DiffOperatorExpr, Mul):
    def __new__(cls, *args):

        if any(arg.is_zero for arg in args):
            return S.Zero

        l = list(args)
        for i in range(len(l)):
            arg = l[i]
            if isinstance(arg, DOAdd):
                others = DOMul(*tuple(l[:i] + l[i+1:]))
                return DOAdd(*(others*a for a in arg._args))

        expr = Mul.__new__(cls, *args)
        return expr

    # is it needed?
    @classmethod
    def flatten(cls, args):
        new_seq = []
        i = 0
        while args:
            o = args.pop()
            if o.__class__ is cls:
                args.extend(o.args)
            else:
                new_seq.append(o)
        new_seq.reverse()
        return new_seq, [], None

      

    def as_ordered_terms(self, order=None):
        return self.args

    @property
    def variable(self):
        for d in self.args:
            if isinstance(d, DiffOperatorExpr):
                return d.variable

    def __call__(self, e):
        for d in self.args:
            if isinstance(d, DiffOperatorExpr):
                e = d(e)
            else:
                e = e*d
        return e


class DOPow(DiffOperatorExpr, Pow):
    def __new__(cls, *args):

        expr = Pow.__new__(cls, *args, evaluate=False)
        return expr

    def __call__(self, e):
        base = self.base
        if isinstance(base, DOMul):
            parts = (arg for arg in base.args)
            res = S.One
            for part in parts:
                if isinstance(part, DiffOperatorExpr):
                    res *= (part**self.exp)(e)
                else:
                    res *= (part**self.exp)
            return res
        else:
            return diff(e, self.base.variable, self.exp)

    def _latex(self, p):
        return "%s^{%s}j" % (self.base._latex(p), self.exp)


def poly_as_expr_do(p):
    """Convert a multinomial form into an expression.

    >>> from sympy.tensor import DiffOperator as d, poly_as_expr_do
        
    """
    rep = p.rep.to_sympy_dict()
    gens = p.gens

    result = []

    for monom, coeff in rep.iteritems():
        term = [coeff]

        any_is_do = False
        if isinstance(coeff, DiffOperatorExpr):
            any_is_do = True

        for g, m in zip(gens, monom):
            if isinstance(g, DiffOperatorExpr):
                term.append(DOPow(g, m))
                any_is_do = True
            else:
                term.append(Pow(g, m))
        if not any_is_do:
            term.append(DiffOperatorOne())
        result.append(DOMul(*term))

    return DOAdd(*result)

