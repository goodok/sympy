
from sympy.tensor.diff_operator import DiffOperator
from sympy.abc import x, y
from sympy.polys.polytools import Poly


def test_Poly():
    dx = DiffOperator(x)
    p = Poly(dx + dx**2, (dx, ))
