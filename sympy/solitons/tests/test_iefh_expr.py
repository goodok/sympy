from sympy.solitons.IEFH_expr import I, E, F, H
from sympy import Symbol, S


def test_mul():
    assert I*I == I
    assert E*I == E
    assert F*I == F
    assert H*I == H
    assert I*E == E
    assert I*F == F
    assert I*H == H

    assert H*H == I
    assert E*E == I
    assert F*F == -I

    assert E*H == -F
    assert E*F == -H
    assert H*E == F
    assert H*F == E
    assert F*E == H
    assert F*H == -E

def test_collec_homogeneous_terms():
    a = Symbol("a")
    e = 2*F*a- 2*F*a
    assert e == S.Zero
        
