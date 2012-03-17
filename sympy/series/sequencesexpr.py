# -*- coding: utf-8 -*-
from sympy import Expr, Symbol, Eq, Mul, Add, Pow, expand, sympify, Tuple
from sympy.core.basic import Basic
from sympy.core.operations import AssocOp

from sympy.core.singleton import (Singleton, S)
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.cache import cacheit, cacheit_recurr
from sympy.core.sets import Interval
from sympy.functions.combinatorial.factorials import factorial, binomial
from sympy.functions.combinatorial.numbers import bell


################################################################################
#     Interfaces: expression operations, interval, printing, main sequence     #
################################################################################

class SeqExprOp(Expr):
    """ Sequence Expression Class
    Sequence Expressions subclass SymPy Expr's so that
    SeqAdd inherits from Add
    SeqMul inherits from Mul

    They use _op_priority to gain control with binary operations (+, *, -, **)
    are used

    They implement operations specific to Sequence Algebra.
    """

    _op_priority = 12.0

    is_Sequence = True
    is_SequenceAtom = False
    is_Identity = False
    is_EmptySequence = False
    # is_commutative = False ???

    # The following is adapted from the core Expr object

    def __neg__(self):
        return SeqMul(S.NegativeOne, self)
    def __abs__(self):
        raise NotImplementedError

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return SeqAdd(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return SeqAdd(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return SeqAdd(self, -other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return SeqAdd(other, -self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return SeqMul(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return SeqMul(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        if other == -S.One:
            return Inverse(self)
        return SeqCauchyPow(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise NotImplementedError("Sequence Power not defined")
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return SeqMul(self, other**S.NegativeOne)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        raise NotImplementedError()

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

class SeqExprInterval(object):
    @property
    def interval(self):
        # abstract property
        # implemented in SeqAdd, SeqMul and other expressions.
        """
        The interval of the sequence's index.

        Examples
        ========
        >>> from sympy import oo
        >>> from sympy.series.sequences import Sequence, SeqPer
        >>> from sympy.series.sequencesexpr import SeqShiftLeft

        >>> a = Sequence((0, 5), finitlist=(1, 2, 3, 4, 5, 6))
        >>> b = Sequence((8, oo), periodical = (0, 2))
        >>> c = a + b

        >>> c.interval
        [0, 5] U [8, oo)

        """
        return None

    @property
    @cacheit
    def start_index(self):
        return self.interval._inf

    @property
    @cacheit
    def stop_index(self):
        return self.interval._sup

    @property
    @cacheit
    def is_infinite(self):
        return self.stop_index == S.Infinity

    @property
    def length(self):
        return self.stop_index - self.start_index + 1

    def is_out_of_range(self, i):
        if isinstance(i, Symbol):
            return None
        if i < self.start_index:
            return True

        if not self.is_infinite:
            if i > self.stop_index:
                return True
        return False

    def calc_interval_from_slice(self, slc):
        slc_start = slc.start
        if slc_start == None:
            slc_start = S.Zero
        slc_stop  = slc.stop
        if slc_stop == None:
            slc_stop = S.Infinity
        return self.interval & Interval(slc_start, slc_stop)

    def getitem_dispatche(self, i):
        if isinstance(i, slice):
            return self.getitem_slicing(i)
        elif self.is_out_of_range(i):
            return S.Zero
        else:
            return self.getitem_index(i)

    def getitem_slicing(self, i):
        mask_interval = self.calc_interval_from_slice(i)
        return SeqSliced(self, mask_interval)

class SeqExprPrint(object):
    show_n = 7
    def _sympystr(self, printer, *args):
        if printer._settings["list_sequences"]:
            printset =  self._get_printset()
            if self.is_infinite or (self.length > self.show_n):
                printset.append("...")

            def _print_item(item):
                if item == '...': return "..."
                else: return printer._print(item)

            l = [_print_item(item) for item in printset]
            return '[' + ', '.join(l) + ']'
        else:
            return printer._print_Basic(self, *args)

    def _pretty(self,  printer, *args):
        printset =  self._get_printset()
        if self.is_infinite or (self.length > self.show_n):
            printset.append("...")
        return printer._print_seq(printset, '[', ']', ', ')

    def _get_printset(self):
        printset = []
        if self.start_index > 1:
          printset.append(S.Zero)
          printset.append("...")
        elif self.start_index ==1:
          printset.append(S.Zero)

        count = self.show_n
        if not self.is_infinite:
            count = min(count, self.length)
        printset.extend([self[i] for i in range(self.start_index, self.start_index + count)])
        return printset

class SeqExprMain(object):
    """
    Interface for aligned to non zero sequence.
    """
    @property
    @cacheit
    def first_nonzero_n(self):
        # or maybe main_offset
        """
        get first non zero index
        """
        # TODO:
        # This is rough, because the zero-test problem, and cycle.
        # - Embed it to the primitive, and analize its interval
        # - Embed it to the operation Add (minimum), Mul (summation)

        max_atempts  = 10
        i = 0
        start_index = self.start_index
        while (i < max_atempts):
            if not self[i + start_index] == S.Zero:
                return i + start_index
            i += 1
        raise StopIteration

    @property
    @cacheit
    def main(self):
        return SeqShiftLeft(self, self.first_nonzero_n)

    @property
    @cacheit
    def mainexp(self):
        return SeqShiftLeftExp(self, self.first_nonzero_n)


class SeqExprShifting(object):
    """
    Short methods for unitary operations.
    """
    def shiftleft(self, n):
        return SeqShiftLeft(self, n)

    def shifright(self, n):
        return SeqShiftRight(self, n)

class SeqExpr(SeqExprOp, SeqExprInterval, SeqExprPrint, SeqExprShifting, SeqExprMain):

    def _hashable_content(self):
        return tuple(sorted(self._args, key=hash))


class EmptySequence(SeqExpr):
    """
    Represents the empty sequence.
    """

    __metaclass__ = Singleton

    is_EmptySequence = True
    def __getitem__(self, i):
        return S.Zero


###########################
#       Operations
###########################

class SeqShiftLeft(SeqExpr):
    """
    Shift sequence to the left.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.series.sequences import Sequence, SeqPer
    >>> from sympy.series.sequencesexpr import SeqShiftLeft
    >>> from sympy.printing.pretty.pretty import pprint

    >>> a = SeqPer((0, oo), (1, 2, 3))
    >>> pprint(a)
    [1, 2, 3, 1, 2, 3, 1, ...]

    >>> pprint(SeqShiftLeft(a, 1))
    [2, 3, 1, 2, 3, 1, 2, ...]

    >>> a.shiftleft(2)
    SeqShiftLeft(SeqPer([0, oo), (1, 2, 3)), 2)

    >>> pprint(a.shiftleft(2))
    [3, 1, 2, 3, 1, 2, 3, ...]


    >>> a = Sequence((0, oo), 'a')
    >>> pprint(a)
    [a[0], a[1], a[2], a[3], a[4], a[5], a[6], ...]

    >>> pprint(a.shiftleft(2))
    [a[2], a[3], a[4], a[5], a[6], a[7], a[8], ...]

    >>> pprint(a.shifright(2))
    [0, 0, a[0], a[1], a[2], a[3], a[4], ...]

    See Also
    ========
    sympy.series.sequencesexpr.SeqShiftRight, SeqShiftLeftExp

    """
    # TODO:
    # SeqShiftLeft(SeqPer([0, oo), (1, 2, 3)), 2) ---> SeqPer([0, oo), (3, 2, 1)
    # seq.shifright(2).shiftleft(2) --> seq
    # seq.left(2).shifright(2) --> seq, but with the first two items removed

    def __new__(cls, *args):
        if (args[1]==0): return args[0]
        expr = Expr.__new__(cls, *args)
        return expr

    @property
    def offset(self):
        return self.args[1]

    def __getitem__(self, i):
        n = i + self.offset
        return self.args[0][n]

    @property
    def interval(self):
        # TODO: calculate
        return self.args[0].interval

class SeqShiftRight(SeqExpr):

    """
    Shift sequence to the right.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.series.sequences import Sequence
    >>> from sympy.printing.pretty.pretty import pprint

    >>> a = Sequence((0, oo), 'a')
    >>> pprint(a.shifright(2))
    [0, 0, a[0], a[1], a[2], a[3], a[4], ...]

    See Also
    ========
    sympy.series.sequencesexpr.SeqShiftLeft, SeqShiftRightExp

    """

    def __new__(cls, *args):
        if (args[1]==0): return args[0]
        expr = Expr.__new__(cls, *args)
        return expr

    @property
    def offset(self):
        return self.args[1]

    def __getitem__(self, i):
        n = i - self.offset
        if n >= 0:
            return self.args[0][n]
        else:
            return S.Zero

    @property
    def interval(self):
        # TODO: calculate
        return self.args[0].interval

class SeqSliced(SeqExpr):
    """
    Return sliced sequence or expression of sequence.

    When sequence is simple, then trivially we can recreated and change interval.
    But when expression of sequence is complex (like multiplication) then we
    wrap this expression and mask original interval of it.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.series.sequences import Sequence, abstract_sequences
    >>> from sympy.printing.pretty.pretty import pprint

    >>> a = Sequence(periodical=(5, 7))
    >>> a[2:5]
    SeqPer([2, 5], (5, 7))


    >>> b = Sequence(periodical=(1, 1))
    >>> c = a*b

    >>> c
    SeqPer([0, oo), (1, 1))*SeqPer([0, oo), (5, 7))
    >>> c.interval
    [0, oo)
    >>> pprint(c)
    [5, 12, 17, 24, 29, 36, 41, ...]

    >>> c[2:5]
    SeqSliced(SeqPer([0, oo), (1, 1))*SeqPer([0, oo), (5, 7)), [2, 5])
    >>> c[2:5].interval
    [2, 5]
    >>> pprint(c[2:5])
    [0, ..., 17, 24, 29, 36]

    """
    def __new__(cls, original, mask_interval):
        assert original.is_Sequence
        obj = SeqExpr.__new__(cls, original, mask_interval)
        return obj

    @property
    def original(self):
        return self._args[0]

    @property
    def mask_interval(self):
        return self._args[1]

    @property
    def interval(self):
        return self.mask_interval & self.original.interval

    def __getitem__(self, i):
        if self.is_out_of_range(i):
            return S.Zero
        else:
            return  self.original[i]

class SeqAdd(SeqExpr, Add):
    """
    A Sum of the Sequences expressions.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.series.sequences import Sequence
    >>> from sympy.printing.pretty.pretty import pprint

    >>> a = Sequence((0, oo), periodical = (1, 0))
    >>> b = Sequence((0, oo), periodical = (0, 2))

    >>> a + b
    SeqPer([0, oo), (1, 0)) + SeqPer([0, oo), (0, 2))

    >>> pprint(a + b)
    [1, 2, 1, 2, 1, 2, 1, ...]

    >>> a - b
    SeqPer([0, oo), (1, 0)) - SeqPer([0, oo), (0, 2))

    >>> pprint(a - b)
    [1, -2, 1, -2, 1, -2, 1, ...]

    >>> pprint(3*a - b)
    [3, -2, 3, -2, 3, -2, 3, ...]

    """

    def __new__(cls, *args):

        if not all(arg.is_Sequence for arg in args):
            raise ValueError("Mix of Sequence and Scalar symbols")

        # remove EmptySequence
        args = tuple([arg for arg in args if not arg.is_EmptySequence])
        if len(args)==0:
            return S.EmptySequence

        expr = Add.__new__(cls, *args)

        if expr.is_Mul:
            return SeqMul(*expr.args)
        return expr

    @classmethod
    def flatten(cls, args_seq):
        return args_seq, [], None

    def as_ordered_terms(self, order=None):
        return self.args

    @property
    def interval(self):
        res = S.EmptySet
        for seq in self.args:
            res = res | seq.interval
        return res

    def __getitem__(self, i):
        if isinstance(i, slice):
            return SeqAdd(*(seq[i] for seq in self.args))
        else:
            return Add(*(seq[i] for seq in self.args))

    def _sympystr(self, printer, *args):
        if printer._settings["list_sequences"]:
            return SeqExpr._sympystr(self, printer, *args)
        else:
            return printer._print_Add(self)

class SeqCoeffMul(SeqExpr, Mul):
    """
   Multiplication of sequence by scalar.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.series.sequences import Sequence
    >>> from sympy.printing.pretty.pretty import pprint

    >>> a = Sequence((0, oo), periodical = (1, 0))
    >>> 2*a
    2*SeqPer([0, oo), (1, 0))

    >>> type(2*a)
    <class 'sympy.series.sequencesexpr.SeqCoeffMul'>

    >>> type(-a)
    <class 'sympy.series.sequencesexpr.SeqCoeffMul'>

    >>> pprint(2*a)
    [2, 0, 2, 0, 2, 0, 2, ...]

    >>> b = Sequence((0, oo), periodical = (0, 2))
    >>> pprint(3*a - b)
    [3, -2, 3, -2, 3, -2, 3, ...]

    """

    def __new__(cls, coeff, seq):
        expr = Mul.__new__(cls, coeff, seq)
        return expr

    def _hashable_content(self):
        return self._args

    @classmethod
    def flatten(cls, args_seq):
        coeff = args_seq[0]
        seq = args_seq[1]
        if isinstance(seq, SeqCoeffMul):
            coeff *= seq.coefficient
            seq = seq.seq
            args_seq = [coeff, seq]
        return args_seq, [], None

    @property
    def coefficient(self):
        return self.args[0]

    @property
    def seq(self):
        return self.args[1]

    @property
    def interval(self):
        return self.seq.interval

    def __getitem__(self, i):
        if isinstance(i, slice):
            return SeqCoeffMul(self.coefficient, self.seq[i])
        else:
            return self.coefficient * self.seq[i]


    def _sympystr(self, printer, *args):
        if printer._settings["list_sequences"]:
            return SeqExpr._sympystr(self, printer, *args)
        else:
            #s = printer._print(self.coefficient, *args)
            #s += "*"
            #s += self.seq._sympystr(printer, *args)
            return printer._print_Mul(self)


class SeqMul(SeqExpr, Mul):
    """
    Dispatcher of multiplication of sequences.
    """

    def __new__(cls, *args):

        if any(arg.is_zero for arg in args):
            return S.Zero

        # if at least one sequence is empty then result is EmptySequence
        if any(seq.is_EmptySequence for seq in args if seq.is_Sequence):
            return S.EmptySequence

        # collect scalar coefficients
        # and collect sequenses (without coefficients)
        # Here we preserve order for noncommutative cases of coefficients, but
        # consider only simple cases: when sequences are commutative itself
        # TODO: use generators
        coeffs = []
        seqs = []
        for arg in args:
            if not arg.is_Sequence:
                coeffs.append(arg)
            elif isinstance(arg, SeqCoeffMul):
                coeffs.append(arg.coefficient)
                seqs.append(arg.seq)
            else:
                seqs.append(arg)

        # calculate the multiplicity of coefficients
        if coeffs==[]:
            coeff = S.One
        else:
            coeff = Mul(*coeffs)

        # if only one seqs then return it
        if len(seqs)==1:
            res = seqs[0]
        else:
            # form Cauchy product
            res = SeqCauchyMul.__new__(SeqCauchyMul, *seqs)

        # wrap with coefficient
        if coeff == S.One:
            return res
        else:
            return SeqCoeffMul(coeff, res)

    # is it needed?
    @classmethod
    def flatten(cls, args_seq):
        return args_seq, [], None

    def as_ordered_terms(self, order=None):
        return self.args

    @property
    def interval(self):
        start = Add(*(s.start_index for s in self.args))
        stop = Add(*(s.stop_index for s in self.args))
        res = Interval(start, stop)
        return res

    def _sympystr(self, printer, *args):
        if printer._settings["list_sequences"]:
            return SeqExpr._sympystr(self, printer, *args)
        else:
            return printer._print_Mul(self)


class SeqMulEW(SeqExpr, Expr):
    """
    Element-wise multiplications of sequences.

    Examples
    ========
    >>> from sympy import oo
    >>> from sympy.series.sequences import Sequence
    >>> from sympy.series.sequencesexpr import SeqMulEW
    >>> from sympy.printing.pretty.pretty import pprint

    >>> a = Sequence((0, oo), periodical = (1, 2))
    >>> b = Sequence((0, oo), periodical = (3, 4, 5))
    >>> c = SeqMulEW(a, b)

    >>> c
    SeqMulEW(SeqPer([0, oo), (1, 2)), SeqPer([0, oo), (3, 4, 5)))

    >>> pprint(c)
    [3, 8, 5, 6, 4, 10, 3, ...]

    See Also
    ========

    sympy.series.sequencesexpr.SeqCauchyMul

    """
    @property
    def interval(self):
        res = self.args[0].interval
        for seq in self.args:
            res = res & seq.interval
        return res

    def __getitem__(self, i):
        return self.getitem_dispatche(i)

    # TODO: use @cacheit_recurr
    @cacheit
    def getitem_index(self, i):
        return self.args[0][i]*self.args[1][i]

    def _sympystr(self, printer, *args):
        if printer._settings["list_sequences"]:
            return SeqExpr._sympystr(self, printer, *args)
        else:
            return printer._print_Basic(self, *args)



class SeqCauchyMul(SeqExpr, Mul):
    """
    Cauchy product of sequences.

        c_n=\sum_{k=0}^n a_k b_{n-k}.


    Examples
    ========
    >>> from sympy import oo
    >>> from sympy.series.sequences import Sequence
    >>> from sympy.printing.pretty.pretty import pprint

    >>> a = Sequence((0, oo), 'a')
    >>> b = Sequence((0, oo), 'b')
    >>> c = a*b

    >>> type(c)
    <class 'sympy.series.sequencesexpr.SeqCauchyMul'>

    >>> c[0]
    a[0]*b[0]

    >>> c[2]
    a[0]*b[2] + a[1]*b[1] + a[2]*b[0]

    >>> c[3]
    a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0]

    See Also
    ========

    sympy.series.sequencesexpr.SeqCauchyMul, SeqExpCauchyMul


    References
    ==========
    .. [1] http://en.wikipedia.org/wiki/Cauchy_product
    .. [2] Donald E. "Knuth Art of Computer Programming, Volume 2: Seminumerical Algorithms",
    3rd ed., sec 4.7 "Manipulation of power series", p 525.
    .. [3] http://en.wikipedia.org/wiki/Convolution#Fast_convolution_algorithms

    """
    def __new__(cls, *args):
        # we do not carry out and collect scalar coefficients here
        # as we do it in SeqMul constructor (when parsing '*' expression like '2*a*3*b'.
        # The process of auto simplification must be logically separated in any case
        # and the algorithms must be indepebed whether this simplification applyed or not
        # see test_coefficient_inside_mul()
        expr = Mul.__new__(cls, *args)
        return expr

    @classmethod
    def flatten(cls, args):
        #TODO: how to use AssocOp.flatten(cls, args)
        new_seq = []
        while args:
            o = args.pop()
            if o.__class__ is cls: # classes must match exactly
                args.extend(o.args)
            else:
                new_seq.append(o)
        # c_part, nc_part, order_symbols
        return [], new_seq, None


    # this is used in core in printer system, (sort_key), roots order.
    def as_ordered_terms(self, order=None):
        return self.args

    @property
    @cacheit
    def interval(self):
        """
        >>> from sympy import oo
        >>> from sympy.series.sequences import Sequence
        >>> from sympy.printing.pretty.pretty import pprint

        >>> a = Sequence((0, oo), 'a')
        >>> b = Sequence((0, oo), 'b')
        >>> c = a*b
        >>> c
        a*b
        >>> type(c)
        <class 'sympy.series.sequencesexpr.SeqCauchyMul'>

        >>> c.interval
        [0, oo)

        >>> c = (Sequence((3, 4), 'a')*Sequence((5, 7), 'b'))
        >>> c.interval
        [8, 11]

        """
        start = Add(*(s.start_index for s in self.args))
        stop = Add(*(s.stop_index for s in self.args))
        res = Interval(start, stop)
        return res

    def __getitem__(self, i):
        return self.getitem_dispatche(i)

    # TODO: use @cacheit_recurr
    @cacheit
    def getitem_index(self, i):
        c = []
        if len(self.args)==2:
            a = self.args[0]
            b = self.args[1]
        else:
            a = self.args[0]
            # recurrsion
            b = SeqCauchyMul(*self.args[1:])
        # TODO: optimize the range (if a.start_index > 0)
        # TODO: optimize k is integer or Expression
        for k in range(0, i+1):
            k = S(k)
            c.append(a[k]*b[i-k])
        return Add(*tuple(c))

    def _sympystr(self, printer, *args):
        if printer._settings["list_sequences"]:
            return SeqExpr._sympystr(self, printer, *args)
        else:
            return printer._print_Mul(self)


class SeqCauchyPow(SeqExpr, Pow):
    """

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.series.sequences import Sequence, SeqPer
    >>> from sympy.printing.pretty.pretty import pprint

    >>> a = Sequence((0, oo), 'a')
    >>> c = a**2

    >>> c[0]
    a[0]**2
    >>> c[1]
    2*a[0]*a[1]
    >>> c[2]
    2*a[0]*a[2] + a[1]**2

    >>> c[6]
    2*a[0]*a[6] + 2*a[1]*a[5] + 2*a[2]*a[4] + a[3]**2

    >>> a = SeqPer((0, oo), (1, 0))
    >>> pprint(a**2)
    [1, 0, 2, 0, 3, ...]

    >>> pprint(a**3)
    [1, 0, 3, 0, 6, 0, 10, ...]


    >>> a = SeqPer((0, oo), (0, 1))
    >>> pprint(a**2)
    [0, 0, 1, 0, 2, 0, 3, ...]

    References
    ==========

    .. [1] Donald E. "Knuth Art of Computer Programming, Volume 2: Seminumerical Algorithms",
    3rd ed., sec 4.7 "Manipulation of power series", p 526.
    .. [2] Faà di Bruno's Formula

    """

    def __new__(cls, *args):
        expr = Pow.__new__(cls, *args)
        return expr

    @property
    @cacheit
    def interval(self):
        if self.exp > 0:
            start = self.base.start_index * self.exp
            stop = self.base.stop_index * self.exp
            res = Interval(start, stop)
        else:
            res = Interval(S.Zero, S.Infinity)
        return res

    def __getitem__(self, i):
        return self.getitem_dispatche(i)

    @cacheit_recurr(0)
    def getitem_index(self, i):
        # TODO: implement generator
        base = self.base
        exp = self.exp
        if i == S.Zero:
            return Pow(base[i], exp)
        else:
            w = base.main
            # TODO: optimize k is integer or Expression
            iw = i - S(self.first_nonzero_n)
            if iw < 0:
                return S.Zero
            elif iw == 0:
                return Pow(w[0], exp)

            c = []
            for k in range(1, iw+1):
                c.append((k*exp - iw + k)*w[k]*self[S(i)-k]) # recursion
            # TODO: optimize cancel
            return (Add(*tuple(c))/iw/w[S.Zero]).cancel()

    @property
    def first_nonzero_n(self):
        return self.base.first_nonzero_n * self.exp;

    def _sympystr(self, printer, *args):
        return printer._print_Pow(self)


################################################################################
#           Operations (for exponential generating sequences)                  #
################################################################################

class SeqShiftLeftExp(SeqShiftLeft):
    """
    Shift sequence to the left exponentially.

    This helper implimetation for the TaylorSeries.

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
    >>> from sympy.series.sequences import Sequence, SeqPer
    >>> from sympy.printing.pretty.pretty import pprint
    >>> from sympy.series.sequencesexpr import SeqShiftLeftExp, SeqShiftRightExp

    >>> a = SeqPer((0, oo), (1, 0))


    >>> al = SeqShiftLeftExp(a, 2)
    >>> pprint(al)
    [1/2, 0, 1/12, 0, 1/30, 0, 1/56, ...]


    See Also
    ========

    sympy.series.sequencesexpr import SeqShiftRightExp

    """
    # TODO:
    #  SeqShiftLeftExp(SeqShiftRightExp(seq, 2), 2) --> seq

    def __getitem__(self, i):
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
    >>> from sympy.series.sequences import Sequence, SeqPer
    >>> from sympy.printing.pretty.pretty import pprint
    >>> from sympy.series.sequencesexpr import SeqShiftLeftExp, SeqShiftRightExp

    >>> a = SeqPer((0, oo), (1, 0))

    >>> ar = SeqShiftRightExp(a, 2)
    >>> pprint(ar)
    [0, 0, 2, 0, 12, 0, 30, ...]

    Revert:

    >>> pprint(SeqShiftLeftExp(ar, 2))
    [1, 0, 1, 0, 1, 0, 1, ...]


    See Also
    ========

    sympy.series.sequencesexpr import SeqShiftLeftExp

    """

    def __getitem__(self, i):
        offset = self.offset
        n = i - offset
        if n < 0:
            return S.Zero
        bc = factorial(i)/factorial(n)  # i > n
        return self.args[0][n]*bc

class SeqExpCauchyMul(SeqCauchyMul, Mul):
    """
    Product of Exponential Generation sequences.

    This is similar to SeqCauchyMul, but for exponentional coefficient we use
    binomial coefficients:

        c_n=\sum_{k=0}^n \binom{n}{k} a_k b_{n-k} .


    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.series.sequences import Sequence, SeqPer
    >>> from sympy.printing.pretty.pretty import pprint
    >>> from sympy.series.sequencesexpr import SeqExpCauchyMul

    >>> a = Sequence((0, oo), 'a')
    >>> b = Sequence((0, oo), 'b')
    >>> c = SeqExpCauchyMul(a, b)
    >>> c
    a*b
    >>> c[0]
    a[0]*b[0]
    >>> c[1]
    a[0]*b[1] + a[1]*b[0]
    >>> c[2]
    a[0]*b[2] + 2*a[1]*b[1] + a[2]*b[0]
    >>> c[3]
    a[0]*b[3] + 3*a[1]*b[2] + 3*a[2]*b[1] + a[3]*b[0]

    See Also
    ========

    sympy.series.sequencesexpr.SeqCauchyMul, sympy.series.taylor

    """
    def __getitem__(self, i):
        return self.getitem_dispatche(i)

    #@cacheit_recurr(0)
    @cacheit
    def getitem_index(self, i):
        if i == S.Zero:
            return Mul(*(seq[i] for seq in self.args))
        else:
            c = []
            if len(self.args)==2:
                a = self.args[0]
                b = self.args[1]
            else:
                a = self.args[0]
                # recurrsion
                b = SeqExpCauchyMul(*self.args[1:])
            # TODO: optimize the range (if a.start_index > 0)
            for k in range(0, i+1):
                c.append(a[k]*b[i-k]*binomial(i, k))
            return Add(*tuple(c))

class SeqExpCauchyPow(SeqCauchyPow):
    """

    Power of sequences (exponential).

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.series.sequences import Sequence, SeqPer
    >>> from sympy.printing.pretty.pretty import pprint
    >>> from sympy.series.sequencesexpr import SeqExpCauchyPow


    Consider coefficient of Taylor series of `sin(x)**2'
        x**2 + x**4/3 + 2*x**6/45  --> [0, 0, 2, 0, 8, 0, 32]

    >>> a = SeqPer((0, oo), (0, 1))
    >>> c = SeqExpCauchyPow(a,2)
    >>> pprint(c)
    [0, 0, 2, 0, 8, 0, 32, ...]

    References
    ==========

    .. [1] Donald E. "Knuth Art of Computer Programming, Volume 2: Seminumerical Algorithms",
    3rd ed., sec 4.7 "Manipulation of power series", p 526.
    .. [2] Faà di Bruno's Formula
    """

    # TODO: implement for various kind of Sequences.
    # TODO: use SeqCauchyPow: a_n = b_n/n!
    def __getitem__(self, i):
        return self.getitem_dispatche(i)

    @cacheit_recurr(0)
    def getitem_index(self, i):
        base = self.base
        exp = self.exp
        if i == S.Zero:
            return Pow(base[i], exp)
        else:
            return self.mainright[i]

    @property
    @cacheit
    def main(self):
        w = self.base.mainexp  # shifting left of self.base
        return SeqExpCauchyPow_Main(w, self.exp)

    @property
    @cacheit
    def mainright(self):
        w = SeqShiftRightExp(self.main, self.first_nonzero_n)
        return w

class SeqExpCauchyPow_Main(SeqCauchyPow):
    """
    Power of sequences (exponential) for main sequence.
    """
    @cacheit_recurr(0)      #TODO: use fist_cached_index.
    def __getitem__(self, i):
        if self.is_out_of_range(i):
            return S.Zero
        else:
            base = self.base
            exp = self.exp
            if i == S.Zero:
                return Pow(base[i], exp)
            else:
                # TODO: optimize k is integer or Expression
                assert self.base.first_nonzero_n == 0
                #w = base.mainexp
                w = base
                if i < 0:
                    return S.Zero
                elif i == 0:
                    return Pow(w[0], exp)

                c = []
                for k in range(1, i + 1):
                    bc = S.One/factorial(i - k)/factorial(k)
                    wik = self[S(i)-k]
                    c.append((k*exp - i + k)*w[k]*wik*bc) # recursion
                # TODO: optimize cancel
                r = (Add(*tuple(c))/i/w[S.Zero]).cancel()
                r = r*factorial(i)
                return r

class SeqExp_FaDeBruno(SeqExpr):
    """
    Calculate g(f(x)) series.

    Let
        g(x)=\sum_{n=1}^\infty {b_n \over n!} x^n
        f(x)=\sum_{n=1}^\infty {a_n \over n!} x^n

    Note, that a_0 == 0, b_n == 0

    Then
        g(f(x)) = \sum_{n=1}^\infty
{\sum_{k=1}^{n} b_k B_{n,k}(a_1,\dots,a_{n-k+1}) \over n!} x^n,
    """
    def __new__(cls, *args):
        assert len(args)==2
        assert all(arg.is_Sequence for arg in args)
        assert args[0][0] == S.Zero
        assert args[1][0] == S.Zero
        expr = Expr.__new__(cls, *args)
        return expr

    def _hashable_content(self):
        return self._args

    @property
    def g(self):
        return self.args[0]

    @property
    def f(self):
        return self.args[1]

    @property
    def interval(self):
        return self.f.interval

    @cacheit
    def __getitem__(self, i):
        s = S.Zero
        for k in xrange(1, i+1):
            s += self.g[k] * bell(i, k, self.f)
        return s

