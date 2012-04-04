# -*- coding: utf-8 -*-
from sympy import S, oo
from sympy.abc import k
from sympy.utilities.pytest import XFAIL, SKIP
from sympy.printing.pretty import pprint
from sympy.printing.pretty import pretty as xpretty
from sympy.printing.latex import latex
from sympy.core.sets import Interval
from sympy.core.symbol import Symbol, symbols
from sympy.core.cache import clear_cache

from sympy.sequences import Sequence, SequenceSymbol, abstract_sequences
from sympy.sequences.kinds import SeqPer, SeqFormula, SeqFunc
from sympy.sequences.expr import SeqAdd, SeqCoeffMul, SeqCauchyMul

def pretty(expr, order=None):
    """ASCII pretty-printing"""
    return xpretty(expr, order=order, use_unicode=False, wrap_line=False)


def upretty(expr, order=None):
    """Unicode pretty-printing"""
    return xpretty(expr, order=order, use_unicode=True, wrap_line=False)

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

def test_symbol_print():
    a = Sequence('a')
    e = a[0]*a[1]**2*a[2]
    assert str(e) == 'a[0]*a[1]**2*a[2]'
    assert pretty(e) == """\
     2   \n\
a0*a1 *a2\
"""
    assert upretty(e) == u"""\
     2   \n\
a₀⋅a₁ ⋅a₂\
"""
    assert latex(e) == r'a_{0} a_{1}^{2} a_{2}'

    e2 = (a[0]+a[3])*a[1]**2*a[2]
    assert latex(e2) == r'\left(a_{0} + a_{3}\right) a_{1}^{2} a_{2}'


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
    #r = c[200]

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

def test_ordered_partitions():
    from sympy.sequences.utilities import OrderedPartitions
    assert OrderedPartitions(4, 1) == ((4,),)
    assert OrderedPartitions(5, 3) == ((1, 1, 3), (1, 2, 2), (1, 3, 1),
            (2, 1, 2), (2, 2, 1), (3, 1, 1))


def test_plain_bell_poly():
    from sympy.sequences.utilities import PlainBellPoly
    B = symbols('b:10')
    B_ = B[1:]
    C = tuple([S(i) for i in range(100)])[1:]

    r = PlainBellPoly(5, 2, B_).value()
    assert r == 2*B[1]*B[4] + 2*B[2]*B[3]

    r = PlainBellPoly(7, 4, B_).value()
    assert r == 4*B[1]**3*B[4] + 12*B[1]**2*B[2]*B[3] + 4*B[1]*B[2]**3

    r = PlainBellPoly(7, 4, C).value()
    assert r == 120

def test_bell():
    from sympy.sequences.utilities import BellPoly
    B = symbols('b:10')
    B_ = B[1:]
    C = tuple([S(i) for i in range(100)])
    C = C[1:]

    r = BellPoly(5, 2, B_).value()
    assert r == 5*B[1]*B[4] + 10*B[2]*B[3]

    r = BellPoly(7, 4, B_).value()
    assert r == 35*B[1]**3*B[4] + 210*B[1]**2*B[2]*B[3] + 105*B[1]*B[2]**3

    r = BellPoly(7, 4, C).value()
    assert r == 2240


def test_equevalence_to_bell():
    from sympy.sequences.utilities import PlainBellPoly as bell_plain
    from sympy.functions.combinatorial.factorials import factorial
    from sympy.functions.combinatorial.numbers import bell
    B = symbols('b:10')
    D = tuple( [B[i]/factorial(i) for i in range(0, len(B))])
    B = B[1:]
    D = D[1:]

    assert bell(4, 1, B) == 24*bell_plain(4, 1, D).value()  # 4!
    assert bell(4, 2, B) == 12*bell_plain(4, 2, D).value()  # 4!/2
    assert bell(4, 3, B) ==  4*bell_plain(4, 3, D).value()  # 4!/3!
    assert bell(4, 4, B) ==    bell_plain(4, 4, D).value()  # 4!/4!

    assert bell(5, 1, B) == 120*bell_plain(5, 1, D).value()  # 5!
    assert bell(5, 2, B) ==  60*bell_plain(5, 2, D).value()  # 5!/2!
    assert bell(5, 3, B) ==  20*bell_plain(5, 3, D).value()  # 5!/3!
    assert bell(5, 4, B) ==   5*bell_plain(5, 4, D).value()  # 5!/4!
