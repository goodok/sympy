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

class IEFH_Expr(Expr):

    _op_priority = 13.0

    is_IEFHExpr = True
    is_Identity = False
    is_commutative = True

    # The following is adapted from the core Expr object

    def __neg__(self):
        return IEFH_Mul(S.NegativeOne, self)
    def __abs__(self):
        raise NotImplementedError

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return IEFH_Add(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return IEFH_Add(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return IEFH_Add(self, -other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return IEFH_Add(other, -self)


    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return IEFH_Mul(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return IEFH_Mul(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return IEFH_Pow(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise NotImplementedError("IEFH Power is not defined")

class IEFH(IEFH_Expr):
    """
    >>> from sympy.solitons.IEFH_expr import IEFH
    >>> h = IEFH("H")
    >>> h
    H

    """

    def __new__(cls, name):
        name = Symbol(name)
        obj = IEFH_Expr.__new__(cls, name)
        return obj

    @property
    def name(self):
        return self._args[0]

    def _latex(self, p):
        return "{%s}" % (p._print(self.name))

    def _sympystr(self, printer, *args):
        return self.name.name

I = IEFH("I")
E = IEFH("E")
F = IEFH("F")
H = IEFH("H")

class IEFH_Add(IEFH_Expr, Add):
    """
    A Sum of the DiffOperator expressions.

    >>> from sympy.solitons.IEFH_expr import IEFH
    >>> h = IEFH("H")
    >>> e = IEFH("E")
    >>> e + h
    E + H

    >>> e + 2*h + h + 1
    E + 3*H + 1

    """

    def __new__(cls, *args):

        expr = Add.__new__(cls, *args)

        if expr.is_Mul:
            return IEFH_Mul(*expr.args)
        return expr

    def __call__(self, e):
        return Add(*tuple(d(e) for d in self.args))


    @classmethod
    def flatten(cls, args):
        """
        >>> from sympy.solitons.IEFH_expr import IEFH, IEFH_Add
        >>> h = IEFH("H")
        >>> e = IEFH("E")
        >>> f = IEFH("F")

        >>> IEFH_Add(e, IEFH_Add(h, f))
        E + F + H
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

        >>> from sympy.solitons.IEFH_expr import I, E, F, H
        >>> from sympy.abc import a

        >>> H + H
        2*H

        >>> 2*F*a- 2*F*a

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
                    newseq.append(IEFH_Mul(c,s))
        if  single_coeff is not S.Zero:
            newseq.insert(0, single_coeff)
        return newseq


class IEFH_Mul(IEFH_Expr, Mul):
    """
        Example:
        --------

        >>> from sympy.solitons.IEFH_expr import IEFH as m
        >>> h = m("H")
        >>> h*h
        I
    """
    def __new__(cls, *args):

        if any(arg.is_zero for arg in args):
            return S.Zero

        l = list(args)
        for i in range(len(l)):
            arg = l[i]
            if isinstance(arg, IEFH_Add):
                others = IEFH_Mul(*tuple(l[:i] + l[i+1:]))
                return IEFH_Add(*(others*a for a in arg._args))

        expr = Mul.__new__(cls, *args)
        return expr

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
        new_seq = cls.do_mul(new_seq)
        return new_seq, [], None

    @classmethod
    def do_mul(cls, args):
        coeff = S.One
        m = IEFH('I')
        for a in args:
            if isinstance(a, IEFH):
                coeff, m = cls.mul_mat(coeff, m, a)
            else:
                coeff = coeff*a
        if coeff == S.One:
            return [m]
        return [coeff, m]

    @classmethod
    def mul_mat(cls, coeff, a, b):
        an = a.name.name
        bn = b.name.name
        if an=="I":
            return coeff, b
        elif bn=="I":
            return coeff, a
        elif an=="E":
            if bn=="E":
                return coeff, IEFH("I")
            elif bn=="H":
                return -coeff, IEFH("F")
            elif bn=="F":
                return -coeff, IEFH("H")
            else:
                pass
        elif an=="H":
            if bn=="H":
                return coeff, IEFH("I")
            elif bn=="E":
                return coeff, IEFH("F")
            elif bn=="F":
                return coeff, IEFH("E")
            else:
                pass
        elif an=="F":
            if bn=="F":
                return -coeff, IEFH("I")
            elif bn=="E":
                return coeff, IEFH("H")
            elif bn=="H":
                return -coeff, IEFH("E")
        else:
            pass

    def as_ordered_terms(self, order=None):
        return self.args


    def as_coeff_Mul(self, rational=False):
        """Efficiently extract the coefficient of a product. """
        coeff, args = self.args[0], self.args[1:]

        if coeff.is_Number and not (rational and not coeff.is_Rational):
            if len(args) == 1:
                return coeff, args[0]
            else:
                return coeff, self._new_rawargs(*args)
        else:
            return S.One, self


class IEFH_Pow(IEFH_Expr, Pow):
    def __new__(cls, *args):
#        if isinstance(args[0], IEFH) and args[1].is_Number:
#            bn = args[0].name.name
#            e = args[1]
#            if bn =="I":
#                return IEFH("I")
#            elif (bn=="E") or (bn=="H") :
#                if e % 2:
#                    return IEFH(bn)
#                else:
#                    return IEFH("I")
#            elif bn=="F":
#                e = e % 4
#                if e==0:
#                    return IEFH("I")
#                elif e==1:
#                    return IEFH("F")
#                elif e==2:
#                    return -IEFH("I")
#                else:
#                    return -IEFH("F")
        mul_args = (args[0],)*args[1].p
        return IEFH_Mul(*mul_args)

#        expr = Pow.__new__(cls, *args, evaluate=False)
#        return expr

    def _latex(self, p):
        return "%s^{%s}j" % (self.base._latex(p), self.exp)

def iefh_fy(expr):
    """
    Recursively walks down an expression tree changing Expr's to IEFH_Expr's
    i.e. Add -> IEFH_Add
         Mul -> IEFH_Mul

    Only changes those Exprs which contain DiffOperator

    This function is useful when traditional SymPy functions which use Mul and
    Add are called on DiffOperatorExpressions. Examples flatten, expand, simplify...

    """
    class_dict = {Mul:IEFH_Mul, Add:IEFH_Add, IEFH_Mul:IEFH_Mul, IEFH_Add:IEFH_Add,
        Pow:IEFH_Pow, IEFH_Pow:IEFH_Pow}

    if expr.__class__ not in class_dict:
        return expr

    args = map(iefh_fy, expr.args) # Recursively call down the tree

    new_cls = class_dict[expr.__class__]

    if not any(isinstance(arg, IEFH_Expr) for arg in args):
        return expr
    
    if new_cls == IEFH_Add:
        return IEFH_Add(*args)
    elif new_cls == IEFH_Mul:
        return IEFH_Mul(*args)
    elif new_cls == IEFH_Pow:
        return IEFH_Pow(*args)
    else:
        return Basic.__new__(class_dict[expr.__class__], *args)
