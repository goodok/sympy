
from sympy import S, oo
from sympy.abc import k
from sympy.matrices import MatrixSymbol
from sympy.utilities.pytest import XFAIL, SKIP
from sympy.printing.pretty import pprint
from sympy.printing.pretty import pretty
from sympy.core.sets import Interval
from sympy.core.symbol import Symbol

from sympy.series.sequences import Sequence, SeqPer, SeqFormula, SeqFunc, SequenceSymbol
from sympy.series.sequencesexpr import SeqAdd
from sympy.series.taylor import TaylorSeries

def test_sequence_index():
    seq = Sequence(Interval(3, oo), formula=(k, S(1)/k))
    assert seq[5:7] == Sequence((5, 7), formula=(k, S(1)/k))

def test_preodiacal():
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
    assert pretty(a) == "[0, ..., 1, 2, 3, 1, 2, ...]"

    a = SeqPer((0, oo), (0, 1))
    import pudb; pudb.set_trace()
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
    assert pretty(a) == "[0, ..., 1/2, 1/3, 1/4, 1/5, 1/6, ...]"


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

    assert pretty(a) == "[0, ..., 1/4, 1/9, 1/16, 1/25, 1/36, ...]"


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


def test_taylorseries():
    seq = Sequence(Interval(3, oo), formula=(k, S(1)/k))
    from sympy.abc import x
    ts = TaylorSeries(x, sequence=seq)
    assert str(ts) == 'x**3/18 + x**4/96 + x**5/600 + x**6/4320 + ...'
    assert ts == TaylorSeries(x, sequence=seq)

    a = TaylorSeries(x, sequence=SeqPer((0, oo), (0, 1)))
    b = TaylorSeries(x, sequence=SeqPer((0, oo), (1, 0)))
    c = a + b
    d = 1 * c

def test_series_print():
    from sympy.abc import x
    from sympy.interactive.printing import init_printing
    a = TaylorSeries(x, sequence=SeqPer((0, oo), (0, 1)))
    b = TaylorSeries(x, sequence=SeqPer((0, oo), (1, 0)))
    c = a + b
    assert str(c) == '1 + x + x**2/2 + x**3/6 + x**4/24 + x**5/120 + x**6/720 + ...'

