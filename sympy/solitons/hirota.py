from sympy import Expr
from sympy.core.basic import Basic
from sympy.core.singleton import (Singleton, S)

from sympy.core.function import diff, UndefinedFunction
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow

from sympy.functions.combinatorial.factorials import binomial
from sympy.matrices import Matrix
from sympy.core.containers import Tuple
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.cache import cacheit


from sympy.printing.pretty.stringpict import prettyForm
from sympy.core.function import _coeff_isneg

from IEFH import I, E, F, H, as_symbol

class Hirota(Expr):
    """
    Direct calculation

    Examples:
    ---------
    >>> from sympy.abc import x
    >>> from sympy.core.function import Function
    >>> from sympy.solitons.hirota import Hirota

    >>> f = Function("f")
    >>> g = Function("g")
    >>> d = Hirota(f, g, (x, 1))
    >>> d.eval()
    f(x)*Derivative(g(x), x) - g(x)*Derivative(f(x), x)

    >>> Hirota(f, g, (x, 2)).eval()
    f(x)*Derivative(g(x), x, x) + g(x)*Derivative(f(x), x, x) - 2*Derivative(f(x), x)*Derivative(g(x), x)
    """
    def __new__(cls, *args, **assumptions):
        _args = args + ()
        obj = Expr.__new__(cls, *_args, **assumptions)
        obj._argset = _args
        return obj
    @property
    def f(self):
        return self._args[0]
    @property
    def g(self):
        return self._args[1]
    @property    
    def x(self):
        return self._args[2][0]
    @property
    def pow(self):
        return self._args[2][1]

    def eval(self):
        res = S.Zero    
        for k in range(self.pow+1):
            res += diff(self.f(self.x), self.x, k)* \
                    diff(self.g(self.x), self.x, self.pow-k)* \
                    (S.NegativeOne**k)*binomial(self.pow, k)
        return res

    def _latex(self, printer, *args):
        x = printer._print(self.x, *args)
        f = printer._print(self.f, *args)
        g = printer._print(self.g, *args)
        p = "^{%s}" % printer._print(self.pow, *args)
        if self.pow == S.One:
            p = ""
        return  r"D_{%s}%s %s \circ %s" % (x, p, f, g)


class HirotaExpr(Expr):

    _op_priority = 12.0

    is_Identity = False
    is_zero = False
    is_commutative = True

    # The following is adapted from the core Expr object

    def __neg__(self):
        return HirotaMul(S.NegativeOne, self)
    def __abs__(self):
        raise NotImplementedError

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return HirotaAdd(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return HirotaAdd(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return HirotaAdd(self, -other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return HirotaAdd(other, -self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return HirotaMul(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return HirotaMul(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return HirotaPow(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise NotImplementedError("Sequence Power not defined")
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return HirotaMul(self, other**S.NegativeOne)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        return HirotaMul(other, self**S.NegativeOne)


class HirotaOne(HirotaExpr):
    """
    Represents the empty Hirota Operator.
    """

    __metaclass__ = Singleton

    is_Identity = True


#######################
#  Hirota Un Applyed
#######################

class HirotaUnapplyedExpr(Expr):

    _op_priority = 12.0

    is_Identity = False
    is_zero = False
    is_commutative = True

    # The following is adapted from the core Expr object

    def __neg__(self):
        return HU_Mul(S.NegativeOne, self)
    def __abs__(self):
        raise NotImplementedError

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return HU_Add(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return HU_Add(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return HU_Add(self, -other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return HU_Add(other, -self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return HU_Mul(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return HU_Mul(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return HU_Pow(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise NotImplementedError("Power of HirotaUnapplyed is not defined")
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        raise NotImplementedError("Division of HirotaUnapplyed is not defined")
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        raise NotImplementedError("Division of HirotaUnapplyed is not defined")


class HirotaUnapplyedBase(HirotaUnapplyedExpr):
    def __call__(self, f, g):
        return HirotaApplyed(self, f, g)
    def _needs_brackets(self):
        return False

class HirotaUnapplyed(HirotaUnapplyedBase):

    #default matrix generator (default for standart Hirota operator)
    _A = -Matrix([[1, 0], [0, -1]])

    def __new__(cls, *args, **assumptions):
        _args = args + ()
        assert len(args)==1
        matrix = assumptions.pop("matrix", cls._A)
        dont_print_matrix = assumptions.pop("dont_print_matrix", False)
        _args = _args + (matrix, dont_print_matrix)
        obj = Expr.__new__(cls, *_args, **assumptions)
        obj._argset = _args
        return obj

    def _hashable_content(self):
        return tuple(sorted(self._args, key=hash))

    @property
    def x(self):
        return self._args[0]
    @property
    def A(self):
        return self._args[1]

    @property
    def dont_print_matrix(self):
        return self._args[2]

    @property
    def all_vars(self):
        return set((self.x, ))

    def eval(self, f, g, matrix, all_vars):
        function_args = (self.x, )
        if isinstance(f, UndefinedFunction):
            f = f(*tuple(all_vars))
        if isinstance(g, UndefinedFunction):
            g = g(*tuple(all_vars))
        m = Matrix([[f, f], [g, g]])
        hs = HS(m, matrix)
        hs = hs.diff(self.x)
        return hs.det()

    def set_matrix(self, matrix):
        return HirotaUnapplyed(self.x, matrix=matrix,
            dont_print_matrix=self.dont_print_matrix)


    def _latex(self, printer, *args):
        """
        >>> from sympy.abc import x
        >>> from sympy.solitons.hirota import HirotaUnapplyed
        >>> from sympy.printing.latex import latex
        >>> latex(HirotaUnapplyed(x))
        'D_{x} |_{- H}'
        >>> latex(HirotaUnapplyed(x,dont_print_matrix=True))
        'D_{x}'
        """
        res = ""
        sx = printer._print(self.x, *args)
        res += "D_{%s}" % (sx)
        if not self.dont_print_matrix and not printer._settings["dont_print_hirota_matrix"]:
            m = printer._print(as_symbol(self.A), *args)
            res = r"%s |_{%s}" % (res, m)
        return  res

    def _pretty(self, printer, *args):
        """
        >>> from sympy.abc import x
        >>> from sympy.solitons.hirota import HirotaUnapplyed
        >>> from sympy.printing.pretty.pretty import pprint
        >>> pprint(HirotaUnapplyed(x))
        Dx|  
          |-H
        >>> pprint(HirotaUnapplyed(x, dont_print_matrix=True))
        Dx
        """
        pform = printer._print_Indexed("D", (self.x, ))
        if  not self.dont_print_matrix  and not printer._settings["dont_print_hirota_matrix"]:
            pform = prettyForm(*pform.below(" "))
            pform = prettyForm(*pform.parens('', '|'))
            sm = printer._print(as_symbol(self.A))
            pm = printer._print("")
            pm = prettyForm(*pm.below(sm))
            pform = prettyForm(*pform.right(pm))
        return pform



class HU_Add(HirotaUnapplyedBase, Add):
    """
    A Sum of the DiffOperator expressions.
    """

    def __new__(cls, *args):

        expr = Add.__new__(cls, *args)

        if expr.is_Mul:
            return HU_Mul(*expr.args)
        return expr

    def _needs_brackets(self):
        return self.dont_print_matrix

    @property
    def dont_print_matrix(self):
        return any(h.dont_print_matrix for h in self._args)
    @property
    def A(self):
        return self._args[0].A

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
        return self._args

    def expand(self):
        # see 'auto expand' in HU_Mul.__new__
        return self

    @property
    def all_vars(self):
        res = set()
        for d in self.args:
            if isinstance(d, HirotaUnapplyedBase):
                all_vars = d.all_vars
                res.update(all_vars)
        return res

    def eval(self, f, g, matrix, all_vars):
        res = S.Zero
        for a in self.args:
            res += a.eval(f, g, matrix, all_vars)
        return res


    def set_matrix(self, matrix):
        new_args = ()
        for a in self._args:
            if isinstance(a, HirotaUnapplyedBase):
                new_args += (a.set_matrix(matrix), )
            else:
                new_args += (a, )
        return HU_Add(*new_args)


    def _latex(self, printer, *args):
        r"""
        >>> from sympy.abc import x, y
        >>> from sympy.solitons.hirota import HirotaUnapplyed
        >>> from sympy.printing.latex import latex
        >>> d = HirotaUnapplyed(x)**2*HirotaUnapplyed(y)*2
        >>> d = d + HirotaUnapplyed(y)**3
        >>> latex(d)
        '\\left(2 D_{x}^{2} D_{y} + D_{y}^{3}\\right) |_{- H}'

        >>> d = HirotaUnapplyed(x)**2*HirotaUnapplyed(y)*2
        >>> d = d + HirotaUnapplyed(y, dont_print_matrix=True)**3
        >>> latex(d)
        '2 D_{x}^{2} D_{y} + D_{y}^{3}'
        """
        _dpm_bak = printer._settings["dont_print_hirota_matrix"]
        printer._settings["dont_print_hirota_matrix"] = True

        tex = printer._print_Add(self, "none")

        printer._settings["dont_print_hirota_matrix"] = _dpm_bak
        if not printer._settings["dont_print_hirota_matrix"] and \
            not self.dont_print_matrix:
            m = printer._print(as_symbol(self.A), *args)
            tex = r"\left(%s\right) |_{%s}" % (tex, m)
        return tex


    def _pretty(self, printer, *args):
        """
        >>> from sympy.abc import x, y
        >>> from sympy.solitons.hirota import HirotaUnapplyed
        >>> from sympy.printing.pretty.pretty import pprint

        >>> d = HirotaUnapplyed(x)**2*HirotaUnapplyed(y)*2
        >>> d = d + HirotaUnapplyed(y)**3
        >>> pprint(d)
        /    2        3\|  
        |2*Dx *Dy + Dy ||  
        \              /|-H

        >>> d = HirotaUnapplyed(x)**2*HirotaUnapplyed(y)*2
        >>> d = d + HirotaUnapplyed(y, dont_print_matrix=True)**3
        >>> pprint(d)
            2        3
        2*Dx *Dy + Dy
        """
        _dpm_bak = printer._settings["dont_print_hirota_matrix"]
        printer._settings["dont_print_hirota_matrix"] = True

        pform = printer._print_Add(self, "none")

        printer._settings["dont_print_hirota_matrix"] = _dpm_bak

        if  not self.dont_print_matrix:
            pform = prettyForm(*pform.below(" "))
            pform = prettyForm(*pform.parens('(', ')'))
            pform = prettyForm(*pform.parens('', '|'))
            sm = printer._print(as_symbol(self.A))
            pm = printer._print("")
            pm = prettyForm(*pm.below(sm))
            pform = prettyForm(*pform.right(pm))
        return pform


class HU_Mul(HirotaUnapplyedBase, Mul):
    def __new__(cls, *args):


        if any(arg.is_zero for arg in args):
            return S.Zero

        # auto expand
        l = list(args)
        for i in range(len(l)):
            arg = l[i]
            if isinstance(arg, HU_Add):
                others = HU_Mul(*tuple(l[:i] + l[i+1:]))
                return HU_Add(*(others*a for a in arg._args))

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

        new_seq = cls.auto_collect(new_seq)

        return new_seq, [], None

    @classmethod
    def auto_collect(cls, args):
        """
        Collect similar terms and transform them to powers.
        """
        coeff = S.One       # standalone term
                            # e.g. 3 * ...
        powers = []
        for o in args:
            if not isinstance(o, HirotaUnapplyedExpr):
                coeff *= o
            elif isinstance(o, HU_Pow):
                powers.append((o.base, o.exp))
            else:
                powers.append((o, S.One))

        # gather powers
        h_bases = {}
        l_bases = []
        
        for b, e in powers:
            if not h_bases.has_key(b):
                h_bases[b] = e
                l_bases.append(b)
            else:
                h_bases[b] += e

        new_seq = []
        if coeff != S.One:
            new_seq.append(coeff)

        for b in l_bases:
            e = h_bases[b]
            if e == S.One:
                new_seq.append(b)
            else:
                new_seq.append(HU_Pow(b, e))
        return new_seq

    def as_ordered_terms(self, order=None):
        return self.args

    @property
    def A(self):
        for h in self.args:
            if isinstance(h, HirotaUnapplyedExpr):
                return h.A

    @property
    def dont_print_matrix(self):
        for h in self._args:
            if isinstance(h, HirotaUnapplyedExpr):
                if h.dont_print_matrix:
                    return True
        return False

    def expand(self):
        # see 'auto expand' in __new__
        return self

    @property
    def all_vars(self):
        res = set()
        for d in self.args:
            if isinstance(d, HirotaUnapplyedBase):
                res.update(d.all_vars)
        return res

    def set_matrix(self, matrix):
        new_args = ()
        for a in self._args:
            if isinstance(a, HirotaUnapplyedBase):
                new_args += (a.set_matrix(matrix), )
            else:
                new_args += (a, )
        return HU_Mul(*new_args)


    def eval(self, f, g, matrix, all_vars):
        coeff = S.One
        vars_powers = ()
        for o in self._args:
            if not isinstance(o, HirotaUnapplyedBase):
                coeff *= o
            elif isinstance(o, HU_Pow):
                vars_powers += ((o.base.x, o.exp), )
            else:
                vars_powers += ((o.x, S.One), )

        if isinstance(f, UndefinedFunction):
            f = f(*tuple(all_vars))
        if isinstance(g, UndefinedFunction):
            g = g(*tuple(all_vars)) 
        m = Matrix([[f, f], [g, g]])
        hs = HS(m, matrix)
        for (v, p) in vars_powers:
            hs = hs.diff_n(v, p)
        return coeff*hs.det()


    def _latex(self, printer, *args):
        """
        >>> from sympy.abc import x, y
        >>> from sympy.solitons.hirota import HirotaUnapplyed
        >>> from sympy.printing.latex import latex
        >>> latex(HirotaUnapplyed(x)**2*HirotaUnapplyed(y)*2)
        '2 D_{x}^{2} D_{y} |_{- H}'

        >>> d = HirotaUnapplyed(x)**2
        >>> d = d*HirotaUnapplyed(y, dont_print_matrix=True)*2
        >>> latex(d)
        '2 D_{x}^{2} D_{y}'
        """
        _dpm = printer._settings["dont_print_hirota_matrix"]
        printer._settings["dont_print_hirota_matrix"] = True

        tex = printer._print_Mul(self)

        printer._settings["dont_print_hirota_matrix"] = _dpm
        if not printer._settings["dont_print_hirota_matrix"] and \
                not self.dont_print_matrix:
            m = printer._print(as_symbol(self.A), *args)
            tex = r"%s |_{%s}" % (tex, m)
        return tex

    def _pretty(self, printer, *args):
        """
        >>> from sympy.abc import x, y
        >>> from sympy.solitons.hirota import HirotaUnapplyed
        >>> from sympy.printing.pretty.pretty import pprint
        >>> pprint(HirotaUnapplyed(x)**2*HirotaUnapplyed(y)*2)
            2   |  
        2*Dx *Dy|  
                |-H

        >>> d = HirotaUnapplyed(x)**2
        >>> d = d*HirotaUnapplyed(y, dont_print_matrix=True)*2
        >>> pprint(d)
            2   
        2*Dx *Dy
        """
        _dpm_bak = printer._settings["dont_print_hirota_matrix"]
        printer._settings["dont_print_hirota_matrix"] = True

        pform = printer._print_Mul(self)

        printer._settings["dont_print_hirota_matrix"] = _dpm_bak

        if not printer._settings["dont_print_hirota_matrix"] and \
                not self.dont_print_matrix:
            pform = prettyForm(*pform.below(" "))
            pform = prettyForm(*pform.parens('', '|'))
            sm = printer._print(as_symbol(self.A))
            pm = printer._print("")
            pm = prettyForm(*pm.below(sm))
            pform = prettyForm(*pform.right(pm))
        return pform



class HU_Pow(HirotaUnapplyedBase, Pow):
    def __new__(cls, *args):
        expr = Pow.__new__(cls, *args, evaluate=False)
        return expr

    @property
    def x(self):
        return self.base.x

    @property
    def A(self):
        return self.base.A

    @property
    def all_vars(self):
        return self.base.all_vars


    @property
    def dont_print_matrix(self):
        for h in self._args:
            if isinstance(h, HirotaUnapplyedExpr):
                if h.dont_print_matrix:
                    return True
        return False

    def eval(self, f, g, matrix, all_vars):
        function_args = (self.base.x, )
        if isinstance(f, UndefinedFunction):
            f = f(*tuple(all_vars))
        if isinstance(g, UndefinedFunction):
            g = g(*tuple(all_vars)) 
        m = Matrix([[f, f], [g, g]])
        hs = HS(m, matrix)
        hs = hs.diff_n(self.base.x, self.exp)
        return hs.det()

    def set_matrix(self, matrix):
        new_args = ()
        base = self.base.set_matrix(matrix)
        return HU_Pow(base, self.exp)


    def _latex(self, printer, *args):
        """
        >>> from sympy.abc import x
        >>> from sympy.solitons.hirota import HirotaUnapplyed
        >>> from sympy.printing.latex import latex
        >>> latex(HirotaUnapplyed(x)**2)
        'D_{x}^{2} |_{- H}'
        >>> latex(HirotaUnapplyed(x, dont_print_matrix=True)**2)
        'D_{x}^{2}'
        """
        _dpm = printer._settings["dont_print_hirota_matrix"]
        printer._settings["dont_print_hirota_matrix"] = True

        tex = printer._print_Pow(self)

        printer._settings["dont_print_hirota_matrix"] = _dpm
        if not printer._settings["dont_print_hirota_matrix"] and \
                not self.dont_print_matrix:
            m = printer._print(as_symbol(self.A), *args)
            tex = r"%s |_{%s}" % (tex, m)
        return tex

    def _pretty(self, printer, *args):
        """
        >>> from sympy.abc import x
        >>> from sympy.solitons.hirota import HirotaUnapplyed
        >>> from sympy.printing.pretty.pretty import pprint
        >>> pprint(HirotaUnapplyed(x)**2)
          2|  
        Dx |  
           |-H
        >>> pprint(HirotaUnapplyed(x, dont_print_matrix=True)**2)
          2
        Dx
        """
        _dpm_bak = printer._settings["dont_print_hirota_matrix"]
        printer._settings["dont_print_hirota_matrix"] = True

        pform = printer._print_Pow(self)

        printer._settings["dont_print_hirota_matrix"] = _dpm_bak

        if not printer._settings["dont_print_hirota_matrix"] and \
                not self.dont_print_matrix:
            pform = prettyForm(*pform.below(" "))
            pform = prettyForm(*pform.parens('', '|'))
            sm = printer._print(as_symbol(self.A))
            pm = printer._print("")
            pm = prettyForm(*pm.below(sm))
            pform = prettyForm(*pform.right(pm))
        return pform



#######################
# Hirota Applyed
#######################


class HirotaApplyedExpr(Expr):

    _op_priority = 12.0

    is_Identity = False
    is_zero = False
    is_commutative = True

    # The following is adapted from the core Expr object

    def __neg__(self):
        return HA_Mul(S.NegativeOne, self)
    def __abs__(self):
        raise NotImplementedError

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return HA_Add(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return HA_Add(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return HA_Add(self, -other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return HA_Add(other, -self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return HA_Mul(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return HA_Mul(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return HA_Pow(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise NotImplementedError("Power of HirotaUnapplyed is not defined")
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        raise NotImplementedError("Division of HirotaUnapplyed is not defined")
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        raise NotImplementedError("Division of HirotaUnapplyed is not defined")


class HirotaApplyedBase(HirotaApplyedExpr):
    def _needs_brackets(self):
        return False


class HirotaApplyed(HirotaApplyedBase):
    def __new__(cls, *args):
        expr = Expr.__new__(cls, *args, evaluate=False)
        return expr

    @property
    def operator(self):
        return self._args[0]

    @property
    def base_operator(self):
        return self._args[0]

    @property
    def coeff(self):
        return S.One

    @property
    def f(self):
        return self._args[1]

    @property
    def g(self):
        return self._args[2]

    @property
    def A(self):
        return self.operator.A

    def get_all_vars(self):
        return self.operator.all_vars


    def _subs(self, n, new_n):
        """
        >>> from sympy.core.function import Function
        >>> from sympy.core.symbol import Symbol
        >>> from sympy.solitons.hirota import HirotaUnapplyed as D
        >>> from sympy.abc import x, y
        >>> n = Symbol("n")
        >>> f = Function("f")
        >>> g = Function("g")

        >>> d = D(x)**2*D(y)
        >>> d2 = D(x)**n*D(y)
        >>> d2 = d2.subs(n, 2)

        >>> d(f, g).eval() == d2(f, g).eval()
        True
        """
        operator = self.operator.subs(n, new_n)
        f = self.f
        g = self.g
        if n == f:
            f = new_n
        if n == g:
            g = new_n
        return HirotaApplyed(operator, f, g)

    def set_matrix(self, matrix):
        operator = self.operator.set_matrix(matrix)
        return HirotaApplyed(operator, self.f, self.g)

    def eval(self, matrix=None, all_vars=None):
        """
        >>> from sympy.abc import x, y
        >>> from sympy.core.function import Function
        >>> from sympy.solitons.hirota import HirotaUnapplyed as D

        >>> f = Function("f")
        >>> g = Function("g")

        >>> d=D(x)(f, g)
        >>> d.eval()
        f(x)*Derivative(g(x), x) - g(x)*Derivative(f(x), x)


        >>> d=(D(x)**2)(f, g)
        >>> d.eval()
        f(x)*Derivative(g(x), x, x) + g(x)*Derivative(f(x), x, x) - 2*Derivative(f(x), x)*Derivative(g(x), x)

        >>> d=(D(x)**0)(f, g)
        >>> d.eval()
        f(x)*g(x)

        >>> d=(2*D(y)*D(x))(f, g)
        >>> d.eval()
        2*f(x, y)*Derivative(g(x, y), x, y) + 2*g(x, y)*Derivative(f(x, y), x, y) - 2*Derivative(f(x, y), x)*Derivative(g(x, y), y) - 2*Derivative(f(x, y), y)*Derivative(g(x, y), x)

        >>> d1 = D(x)**3 + D(y)
        >>> e1 = d1(f, f(x, y).diff(x))
        >>> e1.eval()
        f(x, y)*Derivative(f(x, y), x, y) + f(x, y)*Derivative(f(x, y), x, x, x, x) - Derivative(f(x, y), x)*Derivative(f(x, y), y) - 4*Derivative(f(x, y), x)*Derivative(f(x, y), x, x, x) + 3*Derivative(f(x, y), x, x)**2


        >>> d2 = (D(x)**3 + D(y))*D(x)
        >>> e2 = d2(f, f)
        >>> e2.eval()
        2*f(x, y)*Derivative(f(x, y), x, y) + 2*f(x, y)*Derivative(f(x, y), x, x, x, x) - 2*Derivative(f(x, y), x)*Derivative(f(x, y), y) - 8*Derivative(f(x, y), x)*Derivative(f(x, y), x, x, x) + 6*Derivative(f(x, y), x, x)**2

        >>> 2*e1.eval() == e2.eval()
        True

        """
        if all_vars==None:
            all_vars = self.get_all_vars()
        if matrix==None:
            matrix = self.A
        return self.operator.eval(self.f, self.g, matrix, all_vars)

    def _latex(self, printer, *args):
        r"""
        >>> from sympy.abc import x, y
        >>> from sympy.core.function import Function
        >>> from sympy.printing.latex import latex
        >>> from sympy.solitons.hirota import HirotaUnapplyed as D

        >>> f = Function("f")
        >>> g = Function("g")

        >>> d=2*D(y)*D(x)
        >>> da = d(f, g)
        >>> latex(da)
        '2 D_{x} D_{y} |_{- H} f \\circ g'

        >>> d = D(x)**2 + 2*D(y)*D(x)
        >>> da = d(f, g)
        >>> latex(da)
        '\\left(D_{x}^{2} + 2 D_{x} D_{y}\\right) |_{- H} f \\circ g'

        >>> d = D(x)**2 + 2*D(y)*D(x, dont_print_matrix=True)
        >>> da = d(f, g)
        >>> latex(da)
        '\\left( D_{x}^{2} + 2 D_{x} D_{y} \\right) f \\circ g'
        """
        so = printer._print(self.operator)
        if self.operator._needs_brackets():
            so = r"\left( %s \right)" % so
        sf = printer._print(self.f)
        sg = printer._print(self.g)
        tex = r"%s %s \circ %s" % (so, sf, sg)
        return tex


    def _pretty(self, printer, *args):
        r"""
        >>> from sympy.abc import x, y
        >>> from sympy.core.function import Function
        >>> from sympy.printing.pretty.pretty import pprint
        >>> from sympy.solitons.hirota import HirotaUnapplyed as D
        >>> from sympy.solitons.hirota import Hirota

        >>> f = Function("f")
        >>> g = Function("g")

        >>> d=2*D(y)*D(x)
        >>> da = d(f, g)
        >>> pprint(da)
        2*Dx*Dy|  f o g
               |-H

        >>> d = D(x)**2 + 2*D(y)*D(x)
        >>> da = d(f, g)
        >>> pprint(da)
        /  2          \|       
        |Dx  + 2*Dx*Dy||  f o g
        \             /|-H 

        >>> d = D(x)**2 + 2*D(y)*D(x, dont_print_matrix=True)
        >>> da = d(f, g)
        >>> pprint(da)
        /  2          \     
        \Dx  + 2*Dx*Dy/f o g
        """
        pform = printer._print(self.operator)
        if self.operator._needs_brackets():
            pform = prettyForm(*pform.parens('(', ')'))
        pf = printer._print(self.f)
        pg = printer._print(self.g)
        pform = prettyForm(*pform.right(pf))
        if printer._use_unicode:
            pform = prettyForm(*pform.right(u" \u25CB "))
        else:
            pform = prettyForm(*pform.right(" o "))
        pform = prettyForm(*pform.right(pg))
        return pform


class HA_Add(HirotaApplyedBase, Add):
    """
    A Sum of the HirotaApplyed expressions.
    """

    def __new__(cls, *args):

        expr = Add.__new__(cls, *args)

        if expr.is_Mul:
            return HA_Mul(*expr.args)
        return expr

    def _needs_brackets(self):
        return self.dont_print_matrix

    @property
    def dont_print_matrix(self):
        return any(h.dont_print_matrix for h in self._args)
    @property
    def A(self):
        return self._args[0].A

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
        return self._args

    def expand(self):
        # see 'auto expand' in __new__
        return self

    def get_all_vars(self):
        res = set()
        for d in self.args:
            if isinstance(d, HirotaApplyedBase):
                res.update(d.get_all_vars())
        return res

    def eval(self, matrix=None, all_vars=None):
        if matrix==None:
            matrix = self.A
        if all_vars==None:
            all_vars = self.get_all_vars()
        res = S.Zero
        for a in self.args:
            res += a.eval(matrix=matrix, all_vars=all_vars)
        return res


    def set_matrix(self, matrix):
        new_args = ()
        for a in self._args:
            if isinstance(a, HirotaApplyed):
                new_args += (a.set_matrix(matrix), )
            else:
                new_args += (a, )
        return HA_Add(*new_args)


    def _latex(self, printer, *args):
        r"""
        >>> from sympy.abc import x, y
        >>> from sympy.core.function import Function
        >>> from sympy.printing.latex import latex
        >>> from sympy.solitons.hirota import HirotaUnapplyed as D

        >>> f = Function("f")
        >>> g = Function("g")

        >>> d = D(x)*D(y)
        >>> da = d(f, f) + d(g, g)
        >>> latex(da)
        'D_{x} D_{y} |_{- H} \\left( f \\circ f + g \\circ g \\right)'

        >>> da = d(f, f) - d(g, g)
        >>> latex(da)
        'D_{x} D_{y} |_{- H} \\left( f \\circ f - g \\circ g \\right)'
        """
        args = self._args
        if all(ha.base_operator == args[0].base_operator for ha in args):
            tex = printer._print(args[0].base_operator)
            l = []
            for ha in args:
                sf = printer._print(ha.f)
                sg = printer._print(ha.g)
                sfg = r"%s \circ %s" % (sf, sg)
                if ha.coeff == S.One:
                    coeff = ""
                    sign = " + "
                elif ha.coeff == -S.One:
                    coeff = ""
                    sign = " - "
                else:
                    coeff = printer._print(ha.coeff)
                    if not _coeff_isneg(ha.coeff):
                        sign = " + "
                    else:
                        sign = ""
                l.append((sign, coeff, sfg))
            s = []

            (sign, coeff, sfg) = l[0]
            s.append("%s%s"% (coeff, sfg))

            for sign, coeff, sfg in l[1:]:
                s.append("%s%s%s"% (sign, coeff, sfg))
            tex = r"%s \left( %s \right)" % (tex, "".join(s))
        else:
            tex = printer._print_Add(self, "none")
        return tex

    def _pretty(self, printer, *args):
        pform = printer._print_Add(self, "none")
        return pform


class HA_Mul(HirotaApplyedBase, Mul):
    def __new__(cls, *args):


        if any(arg.is_zero for arg in args):
            return S.Zero

        # auto expand
        l = list(args)
        for i in range(len(l)):
            arg = l[i]
            if isinstance(arg, HA_Add):
                others = HA_Mul(*tuple(l[:i] + l[i+1:]))
                return HA_Add(*(others*a for a in arg._args))

        expr = Mul.__new__(cls, *args)
        return expr

    @property
    def f(self):
        applyed_operators = self.split_coeff_operators[1]
        res =  None
        if len(applyed_operators)==1:
            res = applyed_operators[0].f
        return res

    @property
    def g(self):
        applyed_operators = self.split_coeff_operators[1]
        res =  None
        if len(applyed_operators)==1:
            res = applyed_operators[0].g
        return res

    @property
    def base_operator(self):
        applyed_operators = self.split_coeff_operators[1]
        res =  None
        if len(applyed_operators)==1:
            res = applyed_operators[0].operator
        return res


    @property
    def coeff(self):
        return self.split_coeff_operators[0]

    @property
    def split_coeff_operators(self):
        coeff = S.One
        base_operators = ()
        for ha in self._args:
            if isinstance(ha, HirotaApplyedBase):
                base_operators += (ha, )
            else:
                coeff *= ha
        return (coeff, base_operators)

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

        new_seq = cls.auto_collect(new_seq)

        return new_seq, [], None

    @classmethod
    def auto_collect(cls, args):
        """
        Collect similar terms and transform them to powers.
        """
        coeff = S.One       # standalone term
                            # e.g. 3 * ...
        powers = []
        for o in args:
            if not isinstance(o, HirotaApplyedExpr):
                coeff *= o
            elif isinstance(o, HA_Pow):
                powers.append((o.base, o.exp))
            else:
                powers.append((o, S.One))

        # gather powers
        h_bases = {}
        l_bases = []
        
        for b, e in powers:
            if not h_bases.has_key(b):
                h_bases[b] = e
                l_bases.append(b)
            else:
                h_bases[b] += e

        new_seq = []
        if coeff != S.One:
            new_seq.append(coeff)

        for b in l_bases:
            e = h_bases[b]
            if e == S.One:
                new_seq.append(b)
            else:
                new_seq.append(HU_Pow(b, e))
        return new_seq

    def as_ordered_terms(self, order=None):
        return self.args

    @property
    def A(self):
        for h in self.args:
            if isinstance(h, HirotaApplyedExpr):
                return h.A

    @property
    def dont_print_matrix(self):
        for h in self._args:
            if isinstance(h, HirotaApplyedExpr):
                if h.dont_print_matrix:
                    return True
        return False

    def expand(self):
        # see 'auto expand' in __new__
        return self

    def set_matrix(self, matrix):
        new_args = ()
        for a in self._args:
            if isinstance(a, HirotaApplyedBase):
                new_args += (a.set_matrix(matrix), )
            else:
                new_args += (a, )
        return HA_Mul(*new_args)

    def get_all_vars(self):
        res = set()
        for d in self.args:
            if isinstance(d, HirotaApplyedBase):
                res.update(d.get_all_vars())
        return res

    def eval(self, matrix=None, all_vars=None):
        """
        >>> from sympy.abc import x, y
        >>> from sympy.core.function import Function
        >>> from sympy.printing.latex import latex
        >>> from sympy.solitons.hirota import HirotaUnapplyed as D

        >>> f = Function("f")
        >>> g = Function("g")

        >>> d = D(x)**4 - D(x)*D(y)**3
        >>> dh = d(f, f) - d(g, g)
        >>> dh.eval()
        2*f(x, y)*Derivative(f(x, y), x, x, x, x) - 2*f(x, y)*Derivative(f(x, y), x, y, y, y) - 2*g(x, y)*Derivative(g(x, y), x, x, x, x) + 2*g(x, y)*Derivative(g(x, y), x, y, y, y) - 8*Derivative(f(x, y), x)*Derivative(f(x, y), x, x, x) + 2*Derivative(f(x, y), x)*Derivative(f(x, y), y, y, y) + 6*Derivative(f(x, y), y)*Derivative(f(x, y), x, y, y) + 8*Derivative(g(x, y), x)*Derivative(g(x, y), x, x, x) - 2*Derivative(g(x, y), x)*Derivative(g(x, y), y, y, y) - 6*Derivative(g(x, y), y)*Derivative(g(x, y), x, y, y) + 6*Derivative(f(x, y), x, x)**2 - 6*Derivative(f(x, y), x, y)*Derivative(f(x, y), y, y) - 6*Derivative(g(x, y), x, x)**2 + 6*Derivative(g(x, y), x, y)*Derivative(g(x, y), y, y)
        """
        if matrix==None:
            matrix = self.A
        if all_vars==None:
            all_vars = self.get_all_vars()
        coeff, operators = self.split_coeff_operators
        res = coeff
        for a in operators:
            res *= a.eval(matrix=matrix, all_vars=all_vars)
        return res


    def _latex(self, printer, *args):
        tex = printer._print_Mul(self)
        return tex

    def _pretty(self, printer, *args):
        pform = printer._print_Mul(self)
        return pform

class HA_Pow(HirotaApplyedBase, Pow):
    def __new__(cls, *args):
        expr = Pow.__new__(cls, *args, evaluate=False)
        return expr




# depricated
class HirotaR(HirotaExpr):
    """
    Recoursive calculations

    Examples:
    ---------
    >>> from sympy.abc import x, y
    >>> from sympy.core.function import Function
    >>> from sympy.solitons.hirota import HirotaR
    >>> from sympy.solitons.IEFH import I, E, F, H
    >>> f = Function("f")
    >>> g = Function("g")

    >>> d = HirotaR(f, g, ((x, 1), ))
    >>> d.eval()
    f(x)*Derivative(g(x), x) - g(x)*Derivative(f(x), x)

    >>> HirotaR(f, g, ((x, 2), )).eval()
    f(x)*Derivative(g(x), x, x) + g(x)*Derivative(f(x), x, x) - 2*Derivative(f(x), x)*Derivative(g(x), x)

    >>> d = HirotaR(f, g, ((x, 0), ), matrix=-H)
    >>> d.eval()
    f(x)*g(x)


    >>> d = HirotaR(f, g, ((x, 2), (y, 1)), matrix=-H)
    >>> d.eval()
    f(x, y)*Derivative(g(x, y), x, x, y) - g(x, y)*Derivative(f(x, y), x, x, y) - 2*Derivative(f(x, y), x)*Derivative(g(x, y), x, y) - Derivative(f(x, y), y)*Derivative(g(x, y), x, x) + 2*Derivative(g(x, y), x)*Derivative(f(x, y), x, y) + Derivative(g(x, y), y)*Derivative(f(x, y), x, x)


    >>> d = HirotaR(f, g, ((x, 0), ), matrix=E)
    >>> d.eval()
    f(x)**2/2 - g(x)**2/2

    >>> d = HirotaR(f, g, ((x, 2), ), matrix=E)
    >>> d.eval()
    f(x)*Derivative(f(x), x, x) - g(x)*Derivative(g(x), x, x) - Derivative(f(x), x)**2 + Derivative(g(x), x)**2



    >>> d = HirotaR(f, g, ((x, 2), (y, 1)), matrix=-F)
    >>> d.eval()
    f(x, y)*Derivative(g(x, y), x, x, y) - g(x, y)*Derivative(f(x, y), x, x, y) - 2*Derivative(f(x, y), x)*Derivative(g(x, y), x, y) - Derivative(f(x, y), y)*Derivative(g(x, y), x, x) + 2*Derivative(g(x, y), x)*Derivative(f(x, y), x, y) + Derivative(g(x, y), y)*Derivative(f(x, y), x, x)


    >>> d = HirotaR(f, g, ((x, 1), (y, 1)), matrix=-F)
    >>> d.eval()
    -f(x, y)*Derivative(f(x, y), x, y) - g(x, y)*Derivative(g(x, y), x, y) + Derivative(f(x, y), x)*Derivative(f(x, y), y) + Derivative(g(x, y), x)*Derivative(g(x, y), y)

    >>> d = HirotaR(f, g, ((x, 3), (y, 1)), matrix=-F)
    >>> d.eval()
     f(x, y)*Derivative(f(x, y), x, x, x, y) + g(x, y)*Derivative(g(x, y), x, x, x, y) - 3*Derivative(f(x, y), x)*Derivative(f(x, y), x, x, y) - Derivative(f(x, y), y)*Derivative(f(x, y), x, x, x) - 3*Derivative(g(x, y), x)*Derivative(g(x, y), x, x, y) - Derivative(g(x, y), y)*Derivative(g(x, y), x, x, x) + 3*Derivative(f(x, y), x, x)*Derivative(f(x, y), x, y) + 3*Derivative(g(x, y), x, x)*Derivative(g(x, y), x, y)

    """
    #default matrix generator (default for standart Hirota operator)
    _A = -Matrix([[1, 0], [0, -1]])

    def __new__(cls, *args, **assumptions):
        _args = args + ()
        _matrix = assumptions.pop("matrix", cls._A)
        _args = _args + (_matrix, )
        obj = Expr.__new__(cls, *_args, **assumptions)
        obj._argset = _args
        return obj
    @property
    def f(self):
        return self._args[0]
    @property
    def g(self):
        return self._args[1]
    @property    
    def vals(self):
        return self._args[2]
    @property
    def A(self):
        return self._args[3]

    def subs(self, n, new_n):
        vals = Tuple(*self.vals)
        vals = vals.subs(n, new_n)
        return HirotaR(self.f, self.g, vals._args, self.A)

    def eval(self):
        vals_x = tuple(v[0] for v in self.vals)
        f = self.f(*vals_x)
        g = self.g(*vals_x) 
        m = Matrix([[f, f], [g, g]])
        hs = HS(m, self.A)
        for (v, p) in self.vals:
            hs = hs.diff_n(v, p)
        return hs.det()

    def _latex(self, printer, *args):
        f = printer._print(self.f, *args)
        g = printer._print(self.g, *args)
        res = ""
        for (x, p) in self.vals:
            sx = printer._print(x, *args)
            sp = "^{%s}" % printer._print(p, *args)
            if p == S.One:
                sp = " "
            res += "D_{%s}%s" % (sx, sp)
        m = printer._print(as_symbol(self.A), *args)
        return  r" %s |_{%s} %s  \circ %s" % (res, m, f, g)

class HirotaAdd(HirotaExpr, Add):
    """
    A Sum of the DiffOperator expressions.
    """

    def __new__(cls, *args):

        expr = Add.__new__(cls, *args)

        if expr.is_Mul:
            return DOMul(*expr.args)
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
        return new_seq, [], None

    def as_ordered_terms(self, order=None):
        return self._args


class HirotaMul(HirotaExpr, Mul):
    def __new__(cls, *args):

        if any(arg.is_zero for arg in args):
            return S.Zero

        expr = Mul.__new__(cls, *args)
        return expr

    # is it needed?
    @classmethod
    def flatten(cls, args_seq):
        return args_seq, [], None

    def as_ordered_terms(self, order=None):
        return self.args

    @property
    def vals(self):
        for d in self.args:
            if isinstance(d, HirotaExpr):
                return d.vals

    def __call__(self, e):
        for d in self.args:
            if isinstance(d, HirotaExpr):
                e = d(e)
            else:
                e = e*d
        return e


class HS(Expr):
    """
    Helper structure for Hirota operator calculations.
    """


    # default matrix as generator
    _A = -Matrix([[1, 0], [0, -1]])
    # F = Matrix([[0, 1], [-1, 0]])
    def __new__(cls, *args, **assumptions):
        _args = args + ()
        if len(_args)==1:
            _args = _args + (cls._A, )
        obj = Expr.__new__(cls, *_args, **assumptions)
        obj._argset = _args
        return obj

    @property
    def m(self):
        return self._args[0]

    @property
    def A(self):
        return self._args[1]

    def diff(self, x):
        m = self.m        
        v0 = self.A*m.col(0).diff(x)
        v1 = m.col(1)
        m1 = self._create_matrix(v0, v1)

        v0 = m.col(0)
        v1 = self.A*m.col(1).diff(x)
        m2 = self._create_matrix(v0, v1)

        return  HS_Add(HS(m1, self.A), HS(m2, self.A))

    def diff_n(self, x, n):
        res = self
        for i in range(n):
            res = res.diff(x)
        return res

    def det(self):
        m = self.m        
        v0 = m.col(0)
        v1 = self.A*m.col(1)
        return self._create_matrix(v0, v1).det()/2


    def _create_matrix(self, v0, v1):
        return Matrix([[v0[0], v1[0]], [v0[1], v1[1]]])


class HS_Add(Expr):
    """
    The sum of helper structures
    """
    def __new__(cls, *args, **assumptions):
        _args = args + ()
        obj = Expr.__new__(cls, *_args, **assumptions)
        return obj
    
    def diff(self, x):
        terms = ()
        for h in self._args:
            terms += h.diff(x)._args
        return HS_Add(*terms)

    def diff_n(self, x, n):
        res = self
        for i in range(n):
            res = res.diff(x)
        return res


    def det(self):
        res = S.Zero
        for h in self._args:
            res += h.det()
        return res


