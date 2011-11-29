
from sympy.series.sequences import Sequence, SeqPer, SeqFormula, SeqFunc
from sympy import S, oo
from sympy.abc import k
from sympy.matrices import MatrixSymbol
from sympy.utilities.pytest import XFAIL, SKIP
from sympy.printing.pretty import pprint
from sympy.printing.pretty import pretty
from sympy.core.sets import Interval


@SKIP("")
def test_Add():
    from sympy.abc import x, y, z
    from sympy.core import Add

    e = Add(y, z, x)
    s = str(e)
    e = x + z + y
    s = str(e)

@SKIP("")
def test_sequence_index():
    seq = Sequence(Interval(3, oo), formula=(k, S(1)/k))
    assert seq[5:7] == Sequence((5, 7), formula=(k, S(1)/k))

@SKIP("")
def test_Matrix():
    A = MatrixSymbol('A', 3, 5)
    B = MatrixSymbol('B', 3, 5)
    #from pudb import set_trace; set_trace()
    c = A + B

@SKIP("")
def test_add():
    a = Sequence(Interval(0, oo), baselist=(1, 0), kind="periodical")
    b = Sequence(Interval(0, oo), baselist=(0, 1), kind="periodical")
    c = Sequence(Interval(0, oo), baselist=(0, 1), kind="periodical")
    from pudb import set_trace; set_trace()
    e = a + b
    s = str(e)
    # b = Sequence((0, oo), baselist = (0, 2), kind="periodical")
    # e = a - b
    # s = str(e)
    # from pudb import set_trace; set_trace()
    # s = a + b + c


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

    #print str(a)
    #assert str(a) == "SeqFunc(Interval(2, oo), 1/k**2)"
    assert pretty(a) == "[0, ..., 1/4, 1/9, 1/16, 1/25, 1/36, ...]"


def test_sequance_factory():
    #import pudb; pudb.set_trace()
    seq = Sequence(Interval(3, oo), formula=(k, S(1)/k))
    s = pretty(seq)

    assert seq[4:] == Sequence(Interval(4, oo), formula=(k, S(1)/k))

