
from sympy import Sequence
from sympy import S, oo
from sympy.abc import k



def test_sequence_index():
    seq = Sequence((3, oo), formula=(k, S(1)/k))
    assert seq[5:7] == Sequence((5, 7), formula=(k, S(1)/k))

