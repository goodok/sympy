from sympy import Expr
from sympy.core.basic import Basic
from sympy import S

from sympy.core.function import diff
from sympy.functions.combinatorial.factorials import binomial

class Hirota(Expr):
    """
    Examples:
    ---------
    >>> from sympy.abc import x
    >>> from sympy.core.function import Function
    >>> from sympy.solitons.hirota import Hirota

    >>> f = Function("f")
    >>> g = Function("g")
    >>> d = Hirota(f, g, x)
    >>> d.eval()
    f(x)*Derivative(g(x), x) - g(x)*Derivative(f(x), x)

    >>> Hirota(f, g, x, 2).eval()
    f(x)*Derivative(g(x), x, x) + g(x)*Derivative(f(x), x, x) - 2*Derivative(f(x), x)*Derivative(g(x), x)
    """
    def __new__(cls, *args, **assumptions):
        if len(args)<4:
            args = args + (S.One, )
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
        return self._args[2]
    @property
    def pow(self):
        return self._args[3]

    def eval(self):
        res = S.Zero    
        for k in range(self.pow+1):
            res += diff(self.f(self.x), self.x, k)* \
                    diff(self.g(self.x), self.x, self.pow-k)* \
                    (S.NegativeOne**k)*binomial(self.pow, k)
        return res
