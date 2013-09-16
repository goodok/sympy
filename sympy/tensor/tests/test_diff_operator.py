from sympy.core import symbols, Integer, Symbol, Tuple
from sympy.utilities.pytest import raises

from sympy.core.singleton import (Singleton, S)

from sympy.tensor.diff_operator import DiffOperator, DiffOperator as d, poly_as_do_expr, DiffOperatorOne
from sympy.abc import x, y, z
from sympy import Poly
from sympy.printing.str import sstrrepr
from sympy.functions.elementary.exponential import exp
from sympy import Function


def test_dofy():

    e = d(x)*z + d(x)**0 - S(5)*d(x)*d(y)*z/24 + d(x)**2*d(y)*z**2
    p = Poly(e, z)
    e = poly_as_do_expr(p)
    e(exp(x))
    
    #assert e == z**2*DiffOperator(x)**2*DiffOperator(y) + (-5)*z*DiffOperator(x)*DiffOperator(y)/24 + z*DiffOperator(x) + DiffOperatorOne()

def test_pow():
    f = Function("f")
    e = d(x)**2*d(y)

    from sympy.solitons.IEFH import I, E, F, H
    from sympy.solitons.vertex import VO
    vo = (VO(x, -I, z)*VO(y, -I, z))[0,0]
    p = Poly(vo, z).truncate(2)
    vo(exp(x))
    vo =  poly_as_do_expr(p)

    vo(exp(x))
    vo(f(x, y))
