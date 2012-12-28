from sympy import Expr
from sympy.core.basic import Basic
from sympy import S

from sympy.core.function import diff
from sympy.core.add import Add
from sympy.functions.combinatorial.factorials import binomial
from sympy.matrices import Matrix

class Hirota(Expr):
    """
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




class HirotaR(Expr):
    """
    Recoursive calculations

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
    def vals(self):
        return self._args[2]

    def eval(self):
        vals_x = tuple(v[0] for v in self.vals)
        f = self.f(*vals_x)
        g = self.g(*vals_x) 
        m = Matrix([[f, f], [g, g]])
        hs = HS(m)
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
        return  r" %s %s \circ %s" % (res, f, g)


class HS(Expr):

    A = -Matrix([[1, 0], [0, -1]])
    F = Matrix([[0, 1], [-1, 0]])
    def __new__(cls, *args, **assumptions):
        _args = args + ()
        obj = Expr.__new__(cls, *_args, **assumptions)
        obj._argset = _args
        return obj

    @property
    def m(self):
        return self._args[0]

    def diff(self, x):
        m = self.m        
        v0 = self.A*m.col(0).diff(x)
        v1 = m.col(1)
        m1 = self._create_matrix(v0, v1)

        v0 = m.col(0)
        v1 = self.A*m.col(1).diff(x)
        m2 = self._create_matrix(v0, v1)

        return  HS_Add(HS(m1), HS(m2))

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
 

