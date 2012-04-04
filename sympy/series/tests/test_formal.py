# -*- coding: utf-8 -*-

from sympy import S, oo
from sympy.abc import x, y, k
from sympy.utilities.pytest import XFAIL, SKIP, raises
from sympy.printing.pretty import pprint
from sympy.printing.pretty import pretty as xpretty
from sympy.core.sets import Interval
from sympy.core.symbol import Symbol, symbols
from sympy.core.cache import clear_cache
from sympy import Poly



def pretty(expr, order=None):
    """ASCII pretty-printing"""
    return xpretty(expr, order=order, use_unicode=False, wrap_line=False)


def upretty(expr, order=None):
    """Unicode pretty-printing"""
    return xpretty(expr, order=order, use_unicode=True, wrap_line=False)

from sympy.sequences import Sequence

from sympy.series import PowerSeries, PowerSeries_0E
from sympy.series.power_0 import PowerSeries0

def test_powerseries_constructor():
    x = Symbol('x')
    a = PowerSeries0(x, periodical=(1, 0))
    assert a[0] == 1
    assert a[1] == 0
    assert a[2] == x**2

    a = PowerSeries0(x, 'a')
    assert str(a[1]) in ['a[1]*x', 'x*a[1]']

def test_powerseries_print():

    from sympy.series.power_0 import PowerSeries0

    seq = Sequence((0, oo), periodical=(1, 2, 3))
    ps = PowerSeries0(x, sequence=seq)
    ps2 = ps[11:12]
    s = str(ps2)
    assert s == '3*x**11 + x**12 + ...'

    ps = PowerSeries0(x, periodical=(1, 2, 3))
    s = str(ps)
    assert s == '1 + 2*x + 3*x**2 + x**3 + 2*x**4 + 3*x**5 + x**6 + 2*x**7 + 3*x**8 + ...'

    assert upretty(ps) == \
u"""\
             2    3      4      5    6      7      8      \n\
1 + 2⋅x + 3⋅x  + x  + 2⋅x  + 3⋅x  + x  + 2⋅x  + 3⋅x  + ...\
"""
    assert pretty(ps) ==\
"""\
             2    3      4      5    6      7      8      \n\
1 + 2*x + 3*x  + x  + 2*x  + 3*x  + x  + 2*x  + 3*x  + ...\
"""


def test_powerseries_print_latex():
    from sympy.printing.latex import latex
    from sympy.series.power import PowerSeries0
    seq = Sequence((0, oo), periodical=(1, -2, 3))
    ps = PowerSeries0(x, sequence=seq)
    must = r"1 - 2 x + 3 x^{2} + x^{3} - 2 x^{4} + 3 x^{5} + x^{6}" + \
            " - 2 x^{7} + 3 x^{8} + \dotsb"
    assert latex(ps) == must

def test_types():
    # PowerSeries_0E
    a = PowerSeries_0E(x, periodical=(0, 1))
    b = PowerSeries_0E(x, periodical=(1, 0))

    assert a.is_PowerSeries_0E
    assert b.is_PowerSeries_0E
    assert (a + b).is_PowerSeries_0E
    assert (2*a + 3*b).is_PowerSeries_0E
    assert a[1:5].is_PowerSeries_0E
    a_ = a[1:5]
    b_ = b[2:7]
    c = a + b
    assert c.is_PowerSeries_0E

    from sympy.series.power_0e import PowerSeries_0E_Sliced
    a_ = PowerSeries_0E_Sliced(a, Interval(1, 5))
    assert a_._op_priority >= 13.0
    a_ = a[1:5]
    assert a_._op_priority >= 13.0
    assert a[1:5]._op_priority >= 13.0
    assert (a[1:5] + b[2:7]).is_PowerSeries_0E

    assert (a*b).is_PowerSeries_0E
    assert (a**3).is_PowerSeries_0E
    assert b.compose(a).is_PowerSeries_0E
    assert a.reverse().is_PowerSeries_0E

    c = PowerSeries0(x, periodical=(0, 1))
    d = PowerSeries0(x, periodical=(1, 0))

    # PowerSeries0
    assert c.is_PowerSeries0
    assert d.is_PowerSeries0
    assert (c + d).is_PowerSeries0
    assert (2*c + 3*d).is_PowerSeries0
    assert c[1:5].is_PowerSeries0
    assert (c[1:5] + d[2:7]).is_PowerSeries0

    assert (c*d).is_PowerSeries0
    assert (c**3).is_PowerSeries0
    assert d.compose(c).is_PowerSeries0
    assert c.reverse().is_PowerSeries0


def test_types_raises():
    a = PowerSeries_0E(x, periodical=(0, 1))
    b = PowerSeries_0E(x, periodical=(1, 0))
    c = PowerSeries0(x, periodical=(0, 1))
    d = PowerSeries0(x, periodical=(1, 0))

    raises(ValueError, '(a + c)')

@XFAIL
def test_types_raises_fails():
    a = PowerSeries_0E(x, periodical=(0, 1))
    b = PowerSeries_0E(x, periodical=(1, 0))
    c = PowerSeries0(x, periodical=(0, 1))
    d = PowerSeries0(x, periodical=(1, 0))
    raises(ValueError, '(a*c)')
    raises(ValueError, 'b.compose(c)')


def test_powerseries_power():

    seq = Sequence((0, oo), periodical=(1, 2, 3))
    ps = PowerSeries0(x, sequence=seq)
    c = ps**(-S.One)

    ps = PowerSeries0(x, sequence=Sequence((2, oo), 'a'))
    a = ps.sequence
    assert ps.sequence.order == 2
    c = ps**2

    assert c[0] == 0
    assert c[4] == x**4*a[2]**2

def test_PowerSeries_0E():
    seq = Sequence(Interval(3, oo), formula=(k, S(1)/k))
    ts = PowerSeries_0E(x, sequence=seq)
    assert str(ts) == 'x**3/18 + x**4/96 + x**5/600 + x**6/4320 + x**7/35280 + x**8/322560 + x**9/3265920 + x**10/36288000 + x**11/439084800 + ...'
    assert ts == PowerSeries_0E(x, sequence=seq)

    a = PowerSeries_0E(x, periodical=(0, 1))
    b = PowerSeries_0E(x, periodical=(1, 0))
    c = a + b
    d = 1 * c

    tcosh = PowerSeries_0E(x, periodical = (1, 0))
    tsinh = PowerSeries_0E(x, periodical = (0, 1))
    a = tcosh**2
    b = tsinh**2
    c = a - b
    assert c[0] == 1
    assert c[1] == 0
    assert c[2] == 0

#related with flatten and canonicalization
def test_PowerSeries_0E_coefficient():
    from sympy.series.power_0e import PowerSeries_0E_CoeffMul, PowerSeries_0E_Mul
    x, y = symbols('x, y')
    a = PowerSeries_0E(x, periodical=(0, 1))
    b = PowerSeries_0E(x, periodical=(1, 0))
    c = PowerSeries_0E(x, periodical=(1, 0))

    ts = (3*a)
    str(2*a)
    str(3*a)
    assert ts.coefficient == 3

    ts = 2*(y*(3*a))
    assert ts.coefficient == 6*y

    assert (3*a*y*2).coefficient == 6*y

    abc = a*3*b*c
    assert abc == 3*a*b*c
    assert isinstance(abc, PowerSeries_0E_CoeffMul)

    abc = a*3*b*y*c
    assert abc == 3*y*a*b*c
    assert isinstance(abc, PowerSeries_0E_CoeffMul)
    assert abc.coefficient == 3*y


    # manual construction with internal coefficient of arguments
    # note that it can be temporary, and later automatical simplified to carry out common coefficient
    bc = 3*b*c
    abc = PowerSeries_0E_Mul(2*y*a, bc)
    assert isinstance(abc, PowerSeries_0E_CoeffMul)
    assert isinstance(abc.series, PowerSeries_0E_Mul)
    abc_series = abc.series
    assert len(abc_series.args) == 3
    assert abc_series == a*b*c
    r = abc[0]
    assert r == 0


def test_powerseries_coefficient():
    from sympy.series.power_0 import PowerSeries0CoeffMul, PowerSeries0Mul
    a = PowerSeries0(x, periodical=(0, 1))
    b = PowerSeries0(x, periodical=(1, 0))
    c = PowerSeries0(x, periodical=(1, 0))

    ts = (3*a)
    str(2*a)
    str(3*a)
    assert ts.coefficient == 3

    ts = 2*(y*(3*a))
    assert ts.coefficient == 6*y

    assert (3*a*y*2).coefficient == 6*y

    abc = a*3*b*c
    assert abc == 3*a*b*c
    assert isinstance(abc, PowerSeries0CoeffMul)

    abc = a*3*b*y*c
    assert abc == 3*y*a*b*c
    assert isinstance(abc, PowerSeries0CoeffMul)
    assert abc.coefficient == 3*y


    # manual construction with internal coefficient of arguments
    # note that it can be temporary, and later automatical simplified to carry out common coefficient
    bc = 3*b*c
    abc = PowerSeries0Mul(2*y*a, bc)
    assert isinstance(abc, PowerSeries0CoeffMul)
    assert isinstance(abc.series, PowerSeries0Mul)
    abc_series = abc.series
    assert len(abc_series.args) == 3
    assert abc_series == a*b*c
    r = abc[0]
    assert r == 0

def test_shift():
    ps = PowerSeries0(x, periodical=(1, 2, 3, 4, 5, 6, 7))
    assert str(ps.shift(-2)) == "3 + 4*x + 5*x**2 + 6*x**3 + 7*x**4 + x**5 + 2*x**6 + 3*x**7 + 4*x**8 + ..."

    ps = PowerSeries(x, periodical=(1, 2, 3, 4, 5, 6, 7))
    assert str(ps.shift(-2)) == "3 + 4*x + 5*x**2 + 6*x**3 + 7*x**4 + x**5 + 2*x**6 + 3*x**7 + 4*x**8 + ..."

def test_sliced():
    from sympy.series import PowerSeries

    a = PowerSeries0(x, periodical=(5, 7))
    b = a[2:5]
    str(a[2:5])

    a = PowerSeries(x, periodical=(5, 7))
    str(a[2:5])


def test_nested():
    from sympy.functions.combinatorial.factorials import factorial

    from sympy.series.power_0 import PowerSeries0Nested
    ps = PowerSeries0(x, formula=(k, (-(-1)**k + 1)/factorial(k)/2))
    ns = PowerSeries0Nested(ps, ps)
    assert ns[7] == S(8)*x**7/315

    from sympy.series.power import PowerSeriesNested
    ps = PowerSeries(x, formula=(k, (-(-1)**k + 1)/factorial(k)/2))
    ns = PowerSeriesNested(ps, ps)
    assert ns[7] == S(8)*x**7/315


def test_frompoly():
    p = Poly(x**5 + 3*x**4 + 5*x**2, x)
    ps = PowerSeries0(poly=p)
    assert ps[5] == x**5

    ps = PowerSeries(poly=p)
    assert ps[5] == x**5

    p = Poly(x**5, x)
    ps = PowerSeries(poly=p, point=1)
    assert ps[5] == (x-1)**5
    assert ps[4] == 5*(x-1)**4

def test_compose():
    a = PowerSeries(x, 'a')
    b = PowerSeries(x, 'b')[1:]
    res = a.compose(b)
    v = res[0]
    v = res[1]

    p = Poly(x+x**2)
    pp= PowerSeries0(poly=p)
    psin = PowerSeries_0E(x, periodical=(0, 1)).to_power_series()
    res = psin.compose(pp)
    res[10]


def test_compose_2():
    a = PowerSeries0(x, 'a')
    b = PowerSeries0(x, 'b')[1:]
    c = a.compose(b)
    assert c[0] == a[0]
    sa = a.sequence
    sb = b.sequence
    assert c.coeff(1) == sa[1]*sb[1]
    assert c.coeff(2) == sa[1]*sb[2] + sa[2]*sb[1]**2
    assert c.coeff(3) == sa[1]*sb[3] + 2*sa[2]*sb[1]*sb[2] + sa[3]*sb[1]**3

