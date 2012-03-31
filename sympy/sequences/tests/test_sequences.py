
from sympy import S, oo
from sympy.abc import k
from sympy.utilities.pytest import XFAIL, SKIP
from sympy.printing.pretty import pprint
from sympy.printing.pretty import pretty
from sympy.core.sets import Interval
from sympy.core.symbol import Symbol, symbols
from sympy.core.cache import clear_cache

from sympy.sequences import Sequence, SequenceSymbol, abstract_sequences
from sympy.sequences.kinds import SeqPer, SeqFormula, SeqFunc
from sympy.sequences.expr import SeqAdd, SeqCoeffMul, SeqCauchyMul

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
    assert s == '{a}*{b}*{c}'

@XFAIL
def test_coeffmul_print():
    a, b, c = abstract_sequences('a b c')

    abc = a*b*3*c
    s = str(abc)    # '3*(a*b*c)' now
    assert s == '3*{a}*{b}*{c}'

def test_symbol():
    a = SequenceSymbol((0, 3), 'a')
    i = Symbol('i')
    assert a[i] != S.Zero
    assert a.interval == Interval(0, 3)
    assert a[0] != a[3]
    assert a[4] == S.Zero

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



def test_sequences_expression_slicing():
    a, b, c = abstract_sequences('a:c')
    abc = a*b*c
    r = abc[2:5]
    assert r.interval == Interval(2, 5)


