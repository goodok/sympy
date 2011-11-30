
from sympy import S, oo
from sympy.abc import k
from sympy.matrices import MatrixSymbol
from sympy.utilities.pytest import XFAIL, SKIP
from sympy.printing.pretty import pprint
from sympy.printing.pretty import pretty
from sympy.core.sets import Interval

from sympy.series.sequences import Sequence, SeqPer, SeqFormula, SeqFunc
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


def test_taylorseries():
    seq = Sequence(Interval(3, oo), formula=(k, S(1)/k))
    from sympy.abc import x
    ts = TaylorSeries(x, sequence=seq)
    assert str(ts) == 'x**3/18 + x**4/96 + x**5/600 + ... '
    assert ts == TaylorSeries(x, sequence=seq)

