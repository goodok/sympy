
from sympy import Expr
from sympy.core.cache import cacheit, cacheit_recurr
from sympy.core.sets import Interval
from sympy.core.singleton import (Singleton, S)
from sympy.functions.combinatorial.factorials import factorial

from expr import SeqExpr, SeqExprMethods, SeqCauchyPow, FaDeBruno

class _SeqExprMethods(object):
    """
    Short names for unitary operations constructors.
    """
    def shiftleft(self, n):
        """
        Returns this sequence shifted to the left.
        """
        return SeqShiftLeft(self, n)

    def shiftright(self, n):
        """
        Returns this sequence shifted to the right.
        """
        return SeqShiftRight(self, n)

    def shift(self, n):
        """
        Returns this sequence shifted to the right.

        If n is negative, sequence will be shifted to the left, and elemets below n
        will be discarded.

        Does not change this sequence.
        """
        if (n < S.Zero):
            return SeqShiftLeft(self, -n)
        elif n == S.Zero:
            return self
        return SeqShiftRight(self, n)

    def shift_exp(self, n):
        if (n < S.Zero):
            return SeqShiftLeftExp(self, -n)
        elif n == S.Zero:
            return self
        return SeqShiftRightExp(self, n)

    def shiftleft_exp(self, n):
        return SeqShiftLeftExp(self, n)

    def shiftright_exp(self, n):
        return SeqShiftRightExp(self, n)

    def factorialize(self):
        return Factorialize(self)

    def unfactorialize(self):
        return UnFactorialize(self)

    def compose(self, other):
        return FaDeBruno(self.factorialize(), other.factorialize()).unfactorialize()

    def reverse(self):
        return ReverseLangrange(self)


# TODO: think about hierahy to avoid this.
# Now
#  SeqExp is common object
#  Methods must present in any SeqExpr object
#  But methods uses another SeqExpr objects.
#  so those two modules (expr.py and methods.py) uses each other.

SeqExprMethods.shift = _SeqExprMethods.shift.im_func
SeqExprMethods.shiftleft = _SeqExprMethods.shiftleft.im_func
SeqExprMethods.shiftright = _SeqExprMethods.shiftright.im_func
SeqExprMethods.shift_exp = _SeqExprMethods.shift_exp.im_func
SeqExprMethods.shiftleft_exp = _SeqExprMethods.shiftleft_exp.im_func
SeqExprMethods.shiftright_exp = _SeqExprMethods.shiftright_exp.im_func

SeqExprMethods.factorialize = _SeqExprMethods.factorialize.im_func
SeqExprMethods.unfactorialize = _SeqExprMethods.unfactorialize.im_func
SeqExprMethods.compose = _SeqExprMethods.compose.im_func
SeqExprMethods.reverse = _SeqExprMethods.reverse.im_func


class SeqShiftLeft(SeqExpr):
    """
    Shift sequence to the left.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.sequences import Sequence, SeqPer
    >>> from sympy.sequences.methods import SeqShiftLeft
    >>> from sympy.printing.pretty.pretty import pprint

    >>> a = SeqPer((0, oo), (1, 2, 3))
    >>> pprint(a)
    [1, 2, 3, 1, 2, 3, 1, ...]

    >>> pprint(SeqShiftLeft(a, 1))
    [2, 3, 1, 2, 3, 1, 2, ...]


    >>> a = Sequence((0, oo), 'a')
    >>> pprint(a)
    [a[0], a[1], a[2], a[3], a[4], a[5], a[6], ...]

    >>> pprint(SeqShiftLeft(a, 2))
    [a[2], a[3], a[4], a[5], a[6], a[7], a[8], ...]


    You can also use "shiftleft" function, or "shift" with negative argument:

    >>> pprint(a.shiftleft(2))
    [a[2], a[3], a[4], a[5], a[6], a[7], a[8], ...]

    >>> pprint(a.shift(-2))
    [a[2], a[3], a[4], a[5], a[6], a[7], a[8], ...]

    >>> pprint(a.shiftright(2))
    [0, ..., a[0], a[1], a[2], a[3], a[4], ...]

    You can also use "<<" operators:

    >>> pprint(a << 2)
    [a[2], a[3], a[4], a[5], a[6], a[7], a[8], ...]

    >>> pprint(a >> 2)
    [0, ..., a[0], a[1], a[2], a[3], a[4], ...]

    See Also
    ========
    sympy.sequences.expr.SeqShiftRight, sympy.sequences.expr.SeqShiftLeftExp

    """
    # TODO:
    # SeqShiftLeft(SeqPer([0, oo), (1, 2, 3)), 2) ---> SeqPer([0, oo), (3, 2, 1)
    # seq.shiftright(2).shiftleft(2) --> seq
    # seq.left(2).shiftright(2) --> seq, but with the first two items removed

    def __new__(cls, *args):
        if (args[1]==0): return args[0]
        expr = Expr.__new__(cls, *args)
        return expr

    @property
    def offset(self):
        return self.args[1]

    def getitem_index(self, i):
        n = i + self.offset
        return self.args[0][n]

    @property
    @cacheit
    def interval(self):
        old = self.args[0].interval
        left = old.inf - self.offset
        right = old.sup - self.offset
        # check that the new interval has no negative indexes.
        new = Interval(left, right) & Interval(S.Zero, S.Infinity)
        return new


class SeqShiftRight(SeqExpr):

    """
    Shift sequence to the right.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.sequences import Sequence
    >>> from sympy.printing.pretty.pretty import pprint

    >>> a = Sequence('a')
    >>> pprint(a.shiftright(2))
    [0, ..., a[0], a[1], a[2], a[3], a[4], ...]

    See Also
    ========
    sympy.sequences.expr.SeqShiftLeft, SeqShiftRightExp

    """

    def __new__(cls, *args):
        if (args[1]==0): return args[0]
        expr = Expr.__new__(cls, *args)
        return expr

    @property
    def offset(self):
        return self.args[1]

    def getitem_index(self, i):
        n = i - self.offset
        if n >= 0:
            return self.args[0][n]
        else:
            return S.Zero

    @property
    def interval(self):
        old = self.args[0].interval
        return Interval(old.inf + self.offset, old.sup + self.offset)



class Factorialize(SeqExpr):
    """
    Return sequence with elemetes multiplied to n!.

    Laplace transformatioin.
    It is analog of `ogf()` function power series , which returns the ordinary
    generating function associated to series.
    And analog of "serlaplace" in PARI/GP.
    """
    #TODO:  name-token is laplace trasformation ?
    #TODO: Factorialize(Factorialize(seq)) --> seq

    @property
    def interval(self):
        return self.original.interval

    @property
    def original(self):
        return self._args[0]

    @cacheit
    def getitem_index(self, i):
        return self.args[0][i]*factorial(i)

class UnFactorialize(Factorialize):
    """
    Return sequence with elemetes devided by n!.

    It is analiog of `egf` function of power series, which returns the
    exponential generating function associated to series.

    This can also be computed as convol(f,exp(t)) in PARI/GP ???
    Or as SeqMulEW(self, Sequence(formula(x, 1/x!)))

    """
    @cacheit
    def getitem_index(self, i):
        return self.original[i]/factorial(i)


################################################################################
#           Operations (for exponential generating sequences)                  #
################################################################################

class SeqShiftLeftExp(SeqShiftLeft):
    """
    Shift sequence to the left exponentially.

    This helper implimetation for the PowerESeries.

    Define sequence:
        1 + x**2/2 + x**4/4! + x**6/6! + x**8/8! + ...
    with coefficients
        [1, 0, 1, 0, ...]

    Remove 1 as the first element deleted while shifting:
        x**2/2 + x**4/24 + x**6/720 + ...

    Then carry out the x**2 term:

        x**2/2 + x**4/24 + x**6/720 + ...  = x**2/(1/2 + (2!/4!)*x**2/2!
            + (4!/6!)*x**4/4!) + ...

    So we obtain another sequence:

        [1/2, 0, 1/12, 0, 1/30, 0, 1/56, ...]


    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.sequences import Sequence, SeqPer
    >>> from sympy.printing.pretty.pretty import pprint
    >>> from sympy.sequences.methods import SeqShiftLeftExp, SeqShiftRightExp

    >>> a = SeqPer((0, oo), (1, 0))

    >>> al = SeqShiftLeftExp(a, 2)
    >>> pprint(al)
    [1/2, 0, 1/12, 0, 1/30, 0, 1/56, ...]

    Or we can use methods instead of class:

    >>> pprint(a.shiftleft_exp(2))
    [1/2, 0, 1/12, 0, 1/30, 0, 1/56, ...]


    See Also
    ========

    sympy.sequences.expr.SeqShiftRightExp

    """
    # TODO:
    #  SeqShiftLeftExp(SeqShiftRightExp(seq, 2), 2) --> seq

    def getitem_index(self, i):
        offset = self.offset
        n = i + offset
        bc = factorial(i)/factorial(n) # i < n
        return self.args[0][n]*bc

class SeqShiftRightExp(SeqShiftRight):
    """
    Shift sequence to the right exponentially.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.sequences import Sequence, SeqPer
    >>> from sympy.printing.pretty.pretty import pprint
    >>> from sympy.sequences.methods import SeqShiftLeftExp, SeqShiftRightExp

    >>> a = SeqPer((0, oo), (1, 0))

    >>> ar = SeqShiftRightExp(a, 2)
    >>> pprint(ar)
    [0, ..., 2, 0, 12, 0, 30, ...]

    Or we can use methods instead of class:

    >>> pprint(a.shiftright_exp(2))
    [0, ..., 2, 0, 12, 0, 30, ...]

    Revert:

    >>> pprint(SeqShiftLeftExp(ar, 2))
    [1, 0, 1, 0, 1, 0, 1, ...]

    >>> pprint(ar.shiftleft_exp(2))
    [1, 0, 1, 0, 1, 0, 1, ...]


    See Also
    ========

    sympy.sequences.expr import SeqShiftLeftExp

    """

    def getitem_index(self, i):
        offset = self.offset
        n = i - offset
        if n < 0:
            return S.Zero
        bc = factorial(i)/factorial(n)  # i > n
        return self.args[0][n]*bc

class ReverseLangrange(SeqExpr):
    """
    Composition inversion of power series
    """
    def __new__(cls, *args):
        assert args[0].is_Sequence
        assert len(args)==1
        expr = Expr.__new__(cls, *args)
        return expr

    @property
    def original(self):
        return self._args[0]

    @property
    @cacheit
    def original_prepared(self):
        """
        Return prepeared original series.

        if original sequence is {0, a1, a2, a3...}
        then prepared sequence is {1, 0, a2, a3...}
        """
        #TODO: more convenient way to prepare
        r = self.original.shiftleft(1)
        return r

    @property
    def interval(self):
        return Interval(S.Zero, S.Infinity)

    #@cacheit_recurr(0)
    @cacheit
    def getitem_index(self, i):
        if i == S.Zero:
            return S.Zero
        s_n = self.original_prepared**(-i)  # TODO: very rough
        res = s_n[i-1]/i
        return res
