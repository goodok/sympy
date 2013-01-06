
from sympy.abc import x, y, z
from sympy.core.function import Function
from sympy.core.function import diff
from sympy.printing.latex import latex

from sympy.solitons.hirota import Hirota, HirotaUnapplyed as D, HirotaR

from sympy import Derivative

f = Function("f")
g = Function("g")


F = f(x, y)
Fx = f(x, y).diff(x)
Fxx = f(x, y).diff(x, 2)
Fy = f(x, y).diff(y)
Fyy = f(x, y).diff(y, 2)

Fxy = f(x, y).diff(x).diff(y)
Fxxy = f(x, y).diff(x, 2).diff(y)
Fxyy = f(x, y).diff(x).diff(y, 2)

G = g(x, y)
Gx = g(x, y).diff(x)
Gxx = g(x, y).diff(x, 2)
Gy = g(x, y).diff(y)
Gyy = g(x, y).diff(y, 2)

Gxy = g(x, y).diff(x).diff(y)
Gxxy = g(x, y).diff(x, 2).diff(y)
gxyy = g(x, y).diff(x).diff(y, 2)


def test_hirota():
    f = Function("f")
    g = Function("g")
    d = Hirota(f, g, (x, 1))
    assert d.eval() == f(x)*diff(g(x),x) - diff(f(x), x)*g(x)


def test_hirotaR():
    d = HirotaR(f, g, ((x, 2), (y, 1)) )
    assert d.eval() == F*Gxxy - G*Fxxy - 2*Fx*Gxy - Fy*Gxx + 2*Gx*Fxy +Gy*Fxx

def test_HirotaUnapplyed():
    f = Function("f")
    g = Function("g")
    F = f(x, y)
    G = g(x, y)

    d = D(x)**2*D(y)
    da = d(f, g)
    assert da.eval() == F*Gxxy - G*Fxxy - 2*Fx*Gxy - Fy*Gxx + 2*Gx*Fxy +Gy*Fxx

    #F = f(x, y, z)
    #G = g(x, y, z)
    F = f(x, z, y)
    G = g(x, z, y)
    d = D(x)*D(y) + D(z)**2
    dh = d(f, f) + d(g, g)
    s = dh.get_all_vars()
    assert dh.eval() == 2*F*F.diff(x, y) + 2*G*G.diff(x, y) \
                      + 2*F*F.diff(z, z) + 2*G*G.diff(z, z) \
                      - 2*F.diff(x)*F.diff(y) - 2*G.diff(x)*G.diff(y) \
                      - 2*F.diff(z)**2 - 2*G.diff(z)**2




def test_hirota_latex():
    f = Function("f")
    g = Function("g")
    d = Hirota(f, g, (x, 1))
    assert latex(d) == r"D_{x} f \circ g"


