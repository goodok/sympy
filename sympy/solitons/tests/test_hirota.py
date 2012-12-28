
from sympy.abc import x
from sympy.core.function import Function
from sympy.core.function import diff
from sympy.printing.latex import latex

from sympy.solitons.hirota import Hirota


def test_hirota():
    f = Function("f")
    g = Function("g")
    d = Hirota(f, g, (x, 1))
    assert d.eval() == f(x)*diff(g(x),x) - diff(f(x), x)*g(x)


def test_hirota_latex():
    f = Function("f")
    g = Function("g")
    d = Hirota(f, g, (x, 1))
    assert latex(d) == r"D_{x} f \circ g"


