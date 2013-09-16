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

    >>> from sympy.abc import x
    >>> from sympy.core.function import Function
    >>> from sympy.tensor.diff_operator import DiffOperator
    >>> f = Function("f")
    >>> Dx = DiffOperator(x)

    >>> Dx + Dx
    2*DiffOperator(x)

    >>> 2 + Dx
    DiffOperator(x) + 2*DiffOperatorOne()

    >>> (2 + Dx)(f(x))
    2*f(x) + Derivative(f(x), x)

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
    2*DiffOperator(x)

    >>> Dx(f(x))
    Derivative(f(x), x)

    >>> Dx + Dx**2
    DiffOperator(x)**2 + DiffOperator(x)

    >>> (Dx + Dx**2)(f(x))
    Derivative(f(x), x) + Derivative(f(x), x, x)

    >>> 2 + Dx
    DiffOperator(x) + 2*DiffOperatorOne()

    >>> (2 + Dx)(f(x))
    2*f(x) + Derivative(f(x), x)

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

    >>> from sympy.abc import x, a
    >>> from sympy.core.function import Function
    >>> from sympy.tensor.diff_operator import DiffOperator

    >>> DiffOperator(x) + (- 1*a*DiffOperator(x))
    (-a)*DiffOperator(x) + DiffOperator(x)


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
        """
        >>> from sympy.tensor.diff_operator import DiffOperator as d, DOAdd
        >>> from sympy.abc import x, y, z

        >>> DOAdd(d(x), DOAdd(d(y), d(z)))
        DiffOperator(x) + DiffOperator(y) + DiffOperator(z)
        """
        new_seq = []
        i = 0
        while args:
            o = args.pop()
            if o.__class__ is cls:
                args.extend(o.args)
            else:
                new_seq.append(o)
        new_seq.reverse()

        new_seq = cls.collec_homogeneous_terms(new_seq)
        return new_seq, [], None

    @classmethod
    def collec_homogeneous_terms(cls, seq):
        """
        Collect homogeneous terms
    
        Example:
        --------

        >>> from sympy.tensor.diff_operator import DiffOperator as d
        >>> from sympy.abc import x, y

        >>> d(x) + d(x)
        2*DiffOperator(x)

        """

        terms = {}      # term -> coeff
                        # e.g. x**2 -> 5   for ... + 5*x**2 + ...

        single_coeff = S.Zero  # standalone term (Number or zoo will always be in slot 0)
                        # e.g. 3 + ...
        for o in seq:
            if o.is_Number:
                single_coeff += o
                continue
            elif o.is_Add:
                # NB: here we assume Add is always commutative
                seq.extend(o.args)  # TODO zerocopy?
                continue
            elif o.is_Mul:
                c, s = o.as_coeff_Mul()

                # 3*...
                # unevaluated 2-arg Mul, but we always unfold it so
                # it can combine with other terms (just like is done
                # with the Pow below)
                if c.is_Number and s.is_Add:
                    seq.extend([c*a for a in s.args])
                    continue
            elif o.is_Pow:
                b, e = o.as_base_exp()
                if b.is_Number and (e.is_Integer or (e.is_Rational and e.is_negative)):
                    seq.append(Pow(b, e))
                    continue
                c, s = S.One, o
            else:
                # everything else
                c = S.One
                s = o
            # now we have:
            # o = c*s, where
            #
            # c is a Number
            # s is an expression with number factor extracted

            # let's collect terms with the same s, so e.g.
            # 2*x**2 + 3*x**2  ->  5*x**2
            if s in terms:
                terms[s] += c
            else:
                terms[s] = c

        newseq = []
        noncommutative = False
        for s,c in terms.items():
            # 0*s
            if c is S.Zero:
                continue
            # 1*s
            elif c is S.One:
                newseq.append(s)
            # c*s
            else:
                if s.is_Mul:
                    # Mul, already keeps its arguments in perfect order.
                    # so we can simply put c in slot0 and go the fast way.
                    cs = s._new_rawargs(*((c,) + s.args))
                    newseq.append(cs)

                else:
                    # alternatively we have to call all Mul's machinery (slow)
                    newseq.append(DOMul(c,s))
        if  single_coeff is not S.Zero:
            newseq.insert(0, DOMul(single_coeff, DiffOperatorOne()))
        return newseq

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
        elif isinstance(base, DiffOperatorExpr):
            return diff(e, self.base.variable, self.exp)
        else:
            return self.base**self.exp

    def _latex(self, p):
        return "%s^{%s}j" % (self.base._latex(p), self.exp)


def dofy(expr):
    """
    Recursively walks down an expression tree changing Expr's to DiffOperatorExpr's
    i.e. Add -> DOAdd
         Mul -> DOMul

    Only changes those Exprs which contain DiffOperator

    This function is useful when traditional SymPy functions which use Mul and
    Add are called on DiffOperatorExpressions. Examples flatten, expand, simplify...

    Calling matrixify after calling these functions will reset classes back to
    their matrix equivalents
    """
    class_dict = {DiffOperatorOne:DiffOperatorOne, Mul:DOMul, Add:DOAdd, DOMul:DOMul, DOAdd:DOAdd,
            Pow:DOPow, DOPow:DOPow}

    if expr.__class__ not in class_dict:
        return expr

    args = map(dofy, expr.args) # Recursively call down the tree

    new_cls = class_dict[expr.__class__]

    if not any(isinstance(arg, DiffOperatorExpr) for arg in args):
        return expr
    
    if new_cls == DiffOperatorOne:
        return DiffOperatorOne()
    elif new_cls == DOAdd:
        return DOAdd(*args)
    elif new_cls == DOMul:
        return DOMul(*args)
    elif new_cls == DOPow:
        return DOPow(*args)
    else:
        return Basic.__new__(class_dict[expr.__class__], *args)

def has_do(expr):
    res = False
    args = expr.args
    if any(isinstance(arg, DiffOperatorExpr) for arg in args):
        res = True
    else:
        res = any(has_do(arg) for arg in args)
    return res

def poly_as_do_expr(p):
    """Convert a multinomial form into an diff operator expression.

    >>> from sympy.tensor.diff_operator import DiffOperator as d, poly_as_do_expr
    >>> from sympy.abc import x, y, z
    >>> from sympy import Poly
    >>> from sympy.functions.elementary.exponential import exp
    >>> from sympy.core.function import Function
    >>> f = Function("f")
    >>> from sympy.core.singleton import (Singleton, S)
    
    >>> e = d(x)*z + d(x)**0 - S(5)/24*d(x)*d(y)*z + S(7)/24*d(x)**2*d(y)*z**2
    >>> e(f(x, y))
    7*z**2*Derivative(f(x, y), x, x, y)/24 + z*Derivative(f(x, y), x) - 5*z*Derivative(f(x, y), x, y)/24 + f(x, y)

    >>> p = Poly(e, z)
    >>> e = poly_as_do_expr(p)
    >>> e
    7*z**2*DiffOperator(x)**2*DiffOperator(y)/24 + (-5)*z*DiffOperator(x)*DiffOperator(y)/24 + z*DiffOperator(x) + DiffOperatorOne()

    >>> e(exp(x))
    z*exp(x) + exp(x)

    >>> e(f(x, y))
    7*z**2*Derivative(f(x, y), x, x, y)/24 + z*Derivative(f(x, y), x) - 5*z*Derivative(f(x, y), x, y)/24 + f(x, y)
    """
    e = p.as_expr()
    res = dofy(e)
    return res

def poly_as_do_expr_2(p):
    """Convert a multinomial form into an diff operator expression.

    >>> from sympy.tensor.diff_operator import DiffOperator as d, poly_as_do_expr        
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

