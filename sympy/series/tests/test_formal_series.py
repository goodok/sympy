# -*- coding: utf-8 -*-

from sympy import S, oo
from sympy.abc import x, y, k
from sympy.utilities.pytest import XFAIL, SKIP, raises
from sympy.printing.pretty import pprint
from sympy.printing.pretty import pretty as xpretty
from sympy.core.sets import Interval
from sympy.core.symbol import Symbol, symbols
from sympy.core.cache import clear_cache


def pretty(expr, order=None):
    """ASCII pretty-printing"""
    return xpretty(expr, order=order, use_unicode=False, wrap_line=False)


def upretty(expr, order=None):
    """Unicode pretty-printing"""
    return xpretty(expr, order=order, use_unicode=True, wrap_line=False)

from sympy.sequences import Sequence

from sympy.series import PowerESeries, PowerSeries

def test_powerseries_constructor():
    x = Symbol('x')
    a = PowerSeries(x, periodical=(1, 0))
    assert a[0] == 1
    assert a[1] == 0
    assert a[2] == x**2

    a = PowerSeries(x, 'a')
    assert str(a[1]) in ['a[1]*x', 'x*a[1]']

def test_powerseries_print():

    from sympy.series.power import PowerSeries

    seq = Sequence((0, oo), periodical=(1, 2, 3))
    ps = PowerSeries(x, sequence=seq)
    ps2 = ps[11:12]
    s = str(ps2)
    assert s == '3*x**11 + x**12 + ...'

    ps = PowerSeries(x, periodical=(1, 2, 3))
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
    from sympy.series.power import PowerSeries
    seq = Sequence((0, oo), periodical=(1, -2, 3))
    ps = PowerSeries(x, sequence=seq)
    must = r"1 - 2 x + 3 x^{2} + x^{3} - 2 x^{4} + 3 x^{5} + x^{6}" + \
            " - 2 x^{7} + 3 x^{8} + \dotsb"
    assert latex(ps) == must

def test_types():
    # PowerESeries
    a = PowerESeries(x, periodical=(0, 1))
    b = PowerESeries(x, periodical=(1, 0))

    assert a.is_PowerESeries
    assert b.is_PowerESeries
    assert (a + b).is_PowerESeries
    assert (2*a + 3*b).is_PowerESeries
    assert a[1:5].is_PowerESeries
    a_ = a[1:5]
    b_ = b[2:7]
    c = a + b
    assert c.is_PowerESeries

    from sympy.series.power_e import PowerESeriesSliced
    a_ = PowerESeriesSliced(a, Interval(1, 5))
    assert a_._op_priority >= 13.0
    a_ = a[1:5]
    assert a_._op_priority >= 13.0
    assert a[1:5]._op_priority >= 13.0
    assert (a[1:5] + b[2:7]).is_PowerESeries

    assert (a*b).is_PowerESeries
    assert (a**3).is_PowerESeries
    assert b.compose(a).is_PowerESeries
    assert a.reverse().is_PowerESeries

    c = PowerSeries(x, periodical=(0, 1))
    d = PowerSeries(x, periodical=(1, 0))

    # PowerSeries
    assert c.is_PowerSeries
    assert d.is_PowerSeries
    assert (c + d).is_PowerSeries
    assert (2*c + 3*d).is_PowerSeries
    assert c[1:5].is_PowerSeries
    assert (c[1:5] + d[2:7]).is_PowerSeries

    assert (c*d).is_PowerSeries
    assert (c**3).is_PowerSeries
    assert d.compose(c).is_PowerSeries
    assert c.reverse().is_PowerSeries


def test_types_raises():
    a = PowerESeries(x, periodical=(0, 1))
    b = PowerESeries(x, periodical=(1, 0))
    c = PowerSeries(x, periodical=(0, 1))
    d = PowerSeries(x, periodical=(1, 0))

    raises(ValueError, '(a + c)')

@XFAIL
def test_types_raises_fails():
    a = PowerESeries(x, periodical=(0, 1))
    b = PowerESeries(x, periodical=(1, 0))
    c = PowerSeries(x, periodical=(0, 1))
    d = PowerSeries(x, periodical=(1, 0))
    raises(ValueError, '(a*c)')
    raises(ValueError, 'b.compose(c)')


def test_powerseries_power():

    seq = Sequence((0, oo), periodical=(1, 2, 3))
    ps = PowerSeries(x, sequence=seq)
    c = ps**(-S.One)

    ps = PowerSeries(x, sequence=Sequence((2, oo), 'a'))
    a = ps.sequence
    assert ps.sequence.order == 2
    c = ps**2

    assert c[0] == 0
    assert c[4] == x**4*a[2]**2

def test_PowerESeries():
    seq = Sequence(Interval(3, oo), formula=(k, S(1)/k))
    ts = PowerESeries(x, sequence=seq)
    assert str(ts) == 'x**3/18 + x**4/96 + x**5/600 + x**6/4320 + x**7/35280 + x**8/322560 + x**9/3265920 + x**10/36288000 + x**11/439084800 + ...'
    assert ts == PowerESeries(x, sequence=seq)

    a = PowerESeries(x, periodical=(0, 1))
    b = PowerESeries(x, periodical=(1, 0))
    c = a + b
    d = 1 * c

    tcosh = PowerESeries(x, periodical = (1, 0))
    tsinh = PowerESeries(x, periodical = (0, 1))
    a = tcosh**2
    b = tsinh**2
    c = a - b
    assert c[0] == 1
    assert c[1] == 0
    assert c[2] == 0

#related with flatten and canonicalization
def test_PowerESeries_coefficient():
    from sympy.series.power_e import PowerESeriesCoeffMul, PowerESeriesMul
    x, y = symbols('x, y')
    a = PowerESeries(x, periodical=(0, 1))
    b = PowerESeries(x, periodical=(1, 0))
    c = PowerESeries(x, periodical=(1, 0))

    ts = (3*a)
    str(2*a)
    str(3*a)
    assert ts.coefficient == 3

    ts = 2*(y*(3*a))
    assert ts.coefficient == 6*y

    assert (3*a*y*2).coefficient == 6*y

    abc = a*3*b*c
    assert abc == 3*a*b*c
    assert isinstance(abc, PowerESeriesCoeffMul)

    abc = a*3*b*y*c
    assert abc == 3*y*a*b*c
    assert isinstance(abc, PowerESeriesCoeffMul)
    assert abc.coefficient == 3*y


    # manual construction with internal coefficient of arguments
    # note that it can be temporary, and later automatical simplified to carry out common coefficient
    bc = 3*b*c
    abc = PowerESeriesMul(2*y*a, bc)
    assert isinstance(abc, PowerESeriesCoeffMul)
    assert isinstance(abc.series, PowerESeriesMul)
    abc_series = abc.series
    assert len(abc_series.args) == 3
    assert abc_series == a*b*c
    r = abc[0]
    assert r == 0


def test_powerseries_coefficient():
    from sympy.series.power import PowerSeriesCoeffMul, PowerSeriesMul
    a = PowerSeries(x, periodical=(0, 1))
    b = PowerSeries(x, periodical=(1, 0))
    c = PowerSeries(x, periodical=(1, 0))

    ts = (3*a)
    str(2*a)
    str(3*a)
    assert ts.coefficient == 3

    ts = 2*(y*(3*a))
    assert ts.coefficient == 6*y

    assert (3*a*y*2).coefficient == 6*y

    abc = a*3*b*c
    assert abc == 3*a*b*c
    assert isinstance(abc, PowerSeriesCoeffMul)

    abc = a*3*b*y*c
    assert abc == 3*y*a*b*c
    assert isinstance(abc, PowerSeriesCoeffMul)
    assert abc.coefficient == 3*y


    # manual construction with internal coefficient of arguments
    # note that it can be temporary, and later automatical simplified to carry out common coefficient
    bc = 3*b*c
    abc = PowerSeriesMul(2*y*a, bc)
    assert isinstance(abc, PowerSeriesCoeffMul)
    assert isinstance(abc.series, PowerSeriesMul)
    abc_series = abc.series
    assert len(abc_series.args) == 3
    assert abc_series == a*b*c
    r = abc[0]
    assert r == 0
