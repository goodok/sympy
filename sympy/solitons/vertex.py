from sympy import Expr
from sympy.core.basic import Basic
from sympy.core.singleton import (Singleton, S)

from sympy.core.function import diff, UndefinedFunction
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow

from sympy.functions.combinatorial.factorials import binomial, factorial
from sympy.matrices import Matrix


from sympy.tensor.diff_operator import DiffOperator
from IEFH import I, E, F, H, as_symbol


def VO(x, M, zeta, n = 9):
    """
    Calculate $\exp(M\zeta\partial_{x})$
    """
    d = DiffOperator(x)
    res = M**0*d**0
    for i in range(1, n):
        coeff = M**i/factorial(i)*zeta**i
        res += coeff*d**i
    return res

def ExpM(M, zeta, n = 9):
    """
    Calculate $\exp(M\zeta)$
    """
    res = M**0*zeta**0
    for i in range(1, n):
        coeff = M**i/factorial(i)
        res += coeff*zeta**i
    return res


### Utils not used in this module

def DP_apply(M, v):
    """Apply Matrix DiffOperator to Vector"""
    if M[0,0] is S.Zero:
        m00 = S.Zero
    else:
        m00 = M[0,0](v[0])

    if M[1,0] is S.Zero:
        m10 = S.Zero
    else:
        m10 = M[0,1](v[1])

    if M[0,1] is S.Zero:
        m01 = S.Zero
    else:
        m01 = M[1,0](v[0])

    if M[1,1] is S.Zero:
        m11 = S.Zero
    else:
        m11 = M[1,1](v[1])

    m0 = m00 + m10
    m1 = m01 + m11
    return Matrix([[m0],[m1]])
