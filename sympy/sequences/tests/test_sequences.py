
from sympy import S, oo
from sympy.abc import k
from sympy.matrices import MatrixSymbol
from sympy.utilities.pytest import XFAIL, SKIP
from sympy.printing.pretty import pprint
from sympy.printing.pretty import pretty
from sympy.core.sets import Interval
from sympy.core.symbol import Symbol, symbols
from sympy.core.cache import clear_cache

from sympy.sequences import Sequence, SequenceSymbol, abstract_sequences
from sympy.sequences.kinds import SeqPer, SeqFormula, SeqFunc
from sympy.sequences.expr import SeqAdd, SeqCoeffMul, SeqCauchyMul
from sympy.series import TaylorSeries, PowerSeries


def test_sequence_index():
    seq = Sequence(Interval(3, oo), formula=(k, S(1)/k))
    assert seq[5:7] == Sequence((5, 7), formula=(k, S(1)/k))

def test_preodiacal():
    a = Sequence(periodical=(1))
    assert a.baselist == (1,)

    a = SeqPer(Interval(2, oo), (1, 2, 3))
    assert a.baselist ==  (1, 2, 3)

    b = SeqPer(Interval(2, oo), (1, 2, 3))

    assert a==b

    assert a[0] == S.Zero
    assert a[1] == S.Zero
    assert a[2] == 1
    assert a[3] == 2

    assert a[1:5] == SeqPer(Interval(2, 5), (1, 2, 3))
    assert a[3:5] == SeqPer(Interval(3, 5), (2, 3, 1))

    assert str(a) == "SeqPer([2, oo), (1, 2, 3))"
    assert pretty(a) == "[0, ..., 1, 2, 3, 1, 2, 3, 1, ...]"

    a = SeqPer((0, oo), (0, 1))
    c = a**2

    e = c[2]
    assert e==1

def test_formula():
    from sympy.abc import k

    a = SeqFormula(Interval(2, oo), k, S(1)/k)
    b = SeqFormula(Interval(2, oo), k, S(1)/k)


    assert a==b
    assert a.k == k
    assert a.formula == S(1)/k

    assert a[0] == S.Zero
    assert a[1] == S.Zero
    assert a[2] == S(1)/2
    assert a[3] == S(1)/3

    assert a[1:5] == SeqFormula((2, 5), k, S(1)/k)
    assert a[3:5] == SeqFormula(Interval(3, 5), k, S(1)/k)

    assert str(a) == "SeqFormula([2, oo), k, 1/k)"
    assert pretty(a) == "[0, ..., 1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, ...]"


def test_function():
    from sympy.abc import k
    from sympy import Function
    f = lambda k: S(1)/k**2
    a = SeqFunc(Interval(2, oo), f)
    b = SeqFunc(Interval(2, oo), f)

    assert a==b
    assert a.function == f

    assert a[0] == S.Zero
    assert a[1] == S.Zero
    assert a[2] == S(1)/4
    assert a[3] == S(1)/9

    #assert a[1:5] == SeqFunc((2, 5), f)
    #assert a[3:5] == SeqFunc((3, 5), f)

    assert pretty(a) == "[0, ..., 1/4, 1/9, 1/16, 1/25, 1/36, 1/49, 1/64, ...]"


def test_sequance_factory():
    seq = Sequence(Interval(3, oo), formula=(k, S(1)/k))
    s = pretty(seq)
    assert seq[4:] == Sequence(Interval(4, oo), formula=(k, S(1)/k))

def test_add_empty():
    seq = SeqPer((8, 12), (0, 2))
    r = SeqAdd(seq, S.EmptySequence)
    assert r==seq
    r = seq + S.EmptySequence
    assert r==seq

def test_add():
    a = Sequence((0, oo), periodical = (1, 0))
    b = Sequence((0, oo), periodical = (0, 1))
    c = a - b
    assert c.interval == Interval(0, oo)

    a = Sequence((0, oo), periodical = (1, 0))
    b = Sequence((0, oo), periodical = (0, 2))
    c = 3*a + 7*b
    assert c.is_Sequence

    c = 3*a + 7*b
    d = 2*c
    assert d.is_Sequence

#related with flatten and canonicalization
def test_coefficient():
    a = Sequence(periodical = (1, 0))
    y = Symbol('y')
    seq = (3*a)
    str(2*a)
    str(3*a)
    assert seq.coefficient == 3

    seq = 2*(y*(3*a))
    assert seq.coefficient == 6*y

    assert (3*a*y*2).coefficient == 6*y

    # collecting of coefitients
    a, b, c = abstract_sequences('a,b,c')

    abc = a*3*b*c
    assert abc == 3*a*b*c
    assert isinstance(abc, SeqCoeffMul)
    assert abc[0] == 3*a[0]*b[0]*c[0]

    abc = a*3*b*y*c
    assert abc == 3*y*a*b*c
    assert isinstance(abc, SeqCoeffMul)
    assert abc.coefficient == 3*y
    assert abc[0] == 3*y*a[0]*b[0]*c[0]

def test_coefficient_inside_mul():

    y = Symbol('y')
    a, b, c = abstract_sequences('a:c')

    # manual construction of SeqCauchyMul with internal coefficients of arguments
    bc = 3*b*c
    abc = SeqCauchyMul(2*y*a, bc)
    assert not isinstance(abc, SeqCoeffMul)
    assert isinstance(abc, SeqCauchyMul)
    #assert abc6.coefficient == 6*y
    r = abc[0]
    assert r == 6*y*a[0]*b[0]*c[0]

def test_mul_print():
    clear_cache()
    a, b, c = abstract_sequences('a:c')

    abc = a*b*c
    s = str(abc)
    assert s == 'a*b*c'

@XFAIL
def test_coeffmul_print():
    a, b, c = abstract_sequences('a b c')

    abc = a*b*3*c
    s = str(abc)    # '3*(a*b*c)' now
    assert s == '3*a*b*c'

def test_symbol():
    a = SequenceSymbol((0, 3), 'a')
    i = Symbol('i')
    assert a[i] != S.Zero
    assert a.interval == Interval(0, 3)
    assert a[0] != a[3]
    assert a[4] == S.Zero

    from sympy.interactive.printing import init_printing
    init_printing(list_sequences=True)

    s = str(a)
    s = str(a[i])
    s = str(a[0])


    # cancelation
    c = a**2
    e1 = c[1].args[1].args[0] # a[0] - 0 is Symbol
    e2 = c[1].args[2].args[0] # a[0] - 0 is int
    assert e1 == e2

def test_sequence_CauchyPower():
    a = Sequence(periodical=(0, 1, 0, -1))
    r = a**(2)
    assert r[0] == S.Zero
    assert r[1] == S.Zero


#@SLOW
def test_sequence_CauchyPower_recurr():
    clear_cache()
    a = Sequence(periodical=(1, 0))
    c = a**2
    r = c[200]

def test_sequence_reverse():
    a = Sequence(periodical=(0, 1, 0, -1))
    ssin = a.unfactorialize()
    r = ssin.reverse()
    assert r[0] == S.Zero
    assert r[1] == S.One
    assert r[5] == S(3)/40
    assert r[7] == S(5)/112

    ssin.reverse().compose(ssin)



def test_powerseries_constructor():
    x = Symbol('x')
    a = PowerSeries(x, periodical=(1, 0))
    assert a[0] == 1
    assert a[1] == 0
    assert a[2] == x**2

    a = PowerSeries(x, 'a')
    assert str(a[1]) in ['a[1]*x', 'x*a[1]']

def test_powerseries_print():
    from sympy import S, oo
    from sympy.abc import x, k
    from sympy.series.power import PowerSeries

    seq = Sequence((0, oo), periodical=(1, 2, 3))
    ps = PowerSeries(x, sequence=seq)
    ps2 = ps[11:12]
    s = str(ps2)
    assert s == '3*x**11 + x**12 + ...'

def test_powerseries_power():
    from sympy import S, oo
    from sympy.abc import x, k
    from sympy.series.power import PowerSeries

    seq = Sequence((0, oo), periodical=(1, 2, 3))
    ps = PowerSeries(x, sequence=seq)
    c = ps**(-S.One)

    ps = PowerSeries(x, sequence=Sequence((2, oo), 'a'))
    a = ps.sequence
    assert ps.sequence.first_nonzero_n == 2
    c = ps**2

    assert c[0] == 0
    assert c[4] == x**4*a[2]**2

def test_taylorseries():
    seq = Sequence(Interval(3, oo), formula=(k, S(1)/k))
    from sympy.abc import x
    ts = TaylorSeries(x, sequence=seq)
    assert str(ts) == 'x**3/18 + x**4/96 + x**5/600 + x**6/4320 + x**7/35280 + x**8/322560 + x**9/3265920 + x**10/36288000 + x**11/439084800 + ...'
    assert ts == TaylorSeries(x, sequence=seq)

    a = TaylorSeries(x, sequence=SeqPer((0, oo), (0, 1)))
    b = TaylorSeries(x, sequence=SeqPer((0, oo), (1, 0)))
    c = a + b
    d = 1 * c

#related with flatten and canonicalization
def test_taylorseries_coefficient():
    from sympy.series.taylor import TaylorSeriesCoeffMul, TaylorSeriesMul
    x, y = symbols('x, y')
    a = TaylorSeries(x, sequence=SeqPer((0, oo), (0, 1)))
    b = TaylorSeries(x, sequence=SeqPer((0, oo), (1, 0)))
    c = TaylorSeries(x, sequence=SeqPer((0, oo), (1, 0)))

    ts = (3*a)
    str(2*a)
    str(3*a)
    assert ts.coefficient == 3

    ts = 2*(y*(3*a))
    assert ts.coefficient == 6*y

    assert (3*a*y*2).coefficient == 6*y

    abc = a*3*b*c
    assert abc == 3*a*b*c
    assert isinstance(abc, TaylorSeriesCoeffMul)

    abc = a*3*b*y*c
    assert abc == 3*y*a*b*c
    assert isinstance(abc, TaylorSeriesCoeffMul)
    assert abc.coefficient == 3*y


    # manual construction with internal coefficient of arguments
    # note that it can be temporary, and later automatical simplified to carry out common coefficient
    bc = 3*b*c
    abc = TaylorSeriesMul(2*y*a, bc)
    assert isinstance(abc, TaylorSeriesCoeffMul)
    assert isinstance(abc.series, TaylorSeriesMul)
    abc_series = abc.series
    assert len(abc_series.args) == 3
    assert abc_series == a*b*c
    r = abc[0]
    assert r == 0


def test_powerseries_coefficient():
    from sympy.series.power import PowerSeriesCoeffMul, PowerSeriesMul
    x, y = symbols('x, y')
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

def test_sequences_expression_slicing():
    a, b, c = abstract_sequences('a:c')
    abc = a*b*c
    r = abc[2:5]
    assert r.interval == Interval(2, 5)

def test_taylorseries_expression_slicing():
    x = Symbol("x")
    a = PowerSeries(x, periodical=(0, 2))
    b = PowerSeries(x, periodical=(2, 0))
    c = PowerSeries(x, periodical=(3, 4))
    r = (a*b*c)[2:5]
    assert r.interval == Interval(2, 5)

def test_taylorseries_compose():
    x = Symbol("x")
    tcos = TaylorSeries(x, periodical=(1, 0, -1, 0))
    tsin = TaylorSeries(x, periodical=(0, 1, 0, -1))
    r = tcos.compose(tsin)
    assert r.coeff(4) == S(5)
    assert r.coeff(6) == -S(37)

def test_series_reverse():
    x = Symbol("x")
    tsin = TaylorSeries(x, periodical=(0, 1, 0, -1))
    psin = tsin.to_power_series()
    psin.sequence

    t1 = psin.reverse()
    r = t1[1]


@XFAIL
def test_series_print_finite():
    ps = PowerSeries(x, periodical=(5, 7))
    r = str(ps[2:4])
    assert r== '5*x**2 + 7*x**3 + 5*x**4'  # now '+ ...' trayler.

def test_taylorseries_print():
    from sympy.abc import x
    from sympy.interactive.printing import init_printing
    a = TaylorSeries(x, sequence=SeqPer((0, oo), (0, 1)))
    b = TaylorSeries(x, sequence=SeqPer((0, oo), (1, 0)))
    c = a + b
    assert str(c) == '1 + x + x**2/2 + x**3/6 + x**4/24 + x**5/120 + x**6/720 + x**7/5040 + x**8/40320 + ...'

    a = TaylorSeries(x, sequence=SeqPer((0, oo), (0, 1, 0, -1)))
    s = str(a)
    assert s == 'x - x**3/6 + x**5/120 - x**7/5040 + ...'

