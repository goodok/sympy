# -*- coding: utf-8 -*-
from sympy import Expr, Symbol, Eq, Mul, Add, Pow, expand, sympify, Tuple
from sympy.core.basic import Basic
from sympy.core.singleton import (Singleton, S)
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.cache import cacheit
#from sympy.functions.elementary.miscellaneous import Min, Max
from sympy.core.numbers import (Infinity, Zero)
from sympy.core.sets import Interval
from sympy.functions.combinatorial.factorials import factorial, binomial


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
    Interface for aligned to non zero sequenece.
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


class SeqExpr(SeqExprOp, SeqExprInterval, SeqExprPrint, SeqExprMain):

    def _hashable_content(self):
        return tuple(sorted(self._args, key=hash))


class EmptySequence(SeqExpr):
    """Represents the empty sequence."""

    __metaclass__ = Singleton

    is_EmptySequence = True
    def __getitem__(self, i):
        return S.Zero


class SeqShiftLeft(SeqExpr):

    def __new__(cls, *args):
        # if (args[1]==0): return args[0]
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

    def __new__(cls, *args):
        # if (args[1]==0): return args[0]
        expr = Expr.__new__(cls, *args)
        return expr

    @property
    def offset(self):
        return self.args[1]

    def __getitem__(self, i):
        n = i - self.offset
        if n > 0:
            return self.args[0][n]
        else:
            return S.Zero

    @property
    def interval(self):
        # TODO: calculate
        return self.args[0].interval

class SeqShiftLeftExp(SeqShiftLeft):

    def __getitem__(self, i):
        offset = self.offset
        n = i + offset
        bc = factorial(i)/factorial(n) # i < n
        return self.args[0][n]*bc

class SeqShiftRightExp(SeqShiftRight):

    def __getitem__(self, i):
        offset = self.offset
        n = i - offset
        if n < 0:
            return S.Zero
        bc = factorial(i)/factorial(n)  # i > n
        return self.args[0][n]*bc


class SeqAdd(SeqExpr, Add):
    """A Sum of Sequence Expressions."""

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

class SeqMul(SeqExpr, Mul):
    """A Product of Sequence Expressions (element-wise)."""

    def __new__(cls, *args):

        if any(arg.is_zero for arg in args):
            return S.Zero

        # collect only sequenses
        seqs = [arg for arg in args if arg.is_Sequence]

        # if at least one sequence is empty then result is EmptySequence
        if any(arg.is_EmptySequence for arg in seqs):
            return S.EmptySequence

        # collect scalar coefficients
        coeffs = [arg for arg in args if not arg.is_Sequence]

        # calculate the multiplicity of coefficients
        if coeffs==[]:
            coeff = S.One
        else:
            coeff = Mul(*coeffs)

        # if only one seqs then return it
        if len(seqs)==1:
            if coeff == S.One:
                return seqs[0]
            else:
                return SeqCoeffMul(coeff, seqs[0])

        # Cauchy product
        #expr = SeqCauchyMul.__new__(SeqCauchyMul, *seqs)
        #return expr
        return SeqCauchyMul(*seqs)
        # element-wise multiplicity
        raise NotImplemented

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

    def _sympystr(self, printer, *args):
        if printer._settings["list_sequences"]:
            return SeqExpr._sympystr(self, printer, *args)
        else:
            return printer._print_Mul(self)

class SeqCoeffMul(SeqExpr, Mul):
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
            coeff *= seq.coeff
            seq = seq.seq
            args_seq = [coeff, seq]
        return args_seq, [], None

    @property
    def coeff(self):
        return self.args[0]

    @property
    def seq(self):
        return self.args[1]

    @property
    def interval(self):
        return self.seq.interval

    def __getitem__(self, i):
        if isinstance(i, slice):
            return SeqCoeffMul(self.coeff, self.seq[i])
        else:
            return self.coeff * self.seq[i]


    def _sympystr(self, printer, *args):
        if printer._settings["list_sequences"]:
            return SeqExpr._sympystr(self, printer, *args)
        else:
            return printer._print_Mul(self)

class SeqCauchyMul(SeqExpr, Mul):
    """Cauchy product of sequences.
    
    Discrete convolution of the two sequences.
    """
    def __new__(cls, *args):
        expr = Mul.__new__(cls, *args)
        return expr

    @property
    @cacheit
    def interval(self):
        start = Add(*(s.start_index for s in self.args))
        stop = Add(*(s.stop_index for s in self.args))
        res = Interval(start, stop)
        return res

    @cacheit
    def __getitem__(self, i):
        if self.is_out_of_range(i):
            return S.Zero
        else:
            # TODO: args > 2
            if len(self.args)==2:
                c = []
                a = self.args[0]
                b = self.args[1]
                # TODO: optimize the range (if a.start_index > 0)
                # TODO: optimize k is integer or Expression
                for k in range(0, i+1):
                    k = S(k)
                    c.append(a[k]*b[i-k])
                return Add(*tuple(c))
            return S.Zero

    def _sympystr(self, printer, *args):
        if printer._settings["list_sequences"]:
            return SeqExpr._sympystr(self, printer, *args)
        else:
            return printer._print_Mul(self)


class SeqCauchyPow(SeqExpr, Pow):

    """
    Faà di Bruno's Formula
    
    Donald E. "Knuth Art of Computer Programming, Volume 2: Seminumerical Algorithms",
    3rd ed., sec 4.7 "Manipulation of power series", p 526.
    """

    def __new__(cls, *args):
        expr = Pow.__new__(cls, *args)
        return expr

    @property
    @cacheit
    def interval(self):
        start = self.base.start_index * self.exp
        stop = self.base.stop_index * self.exp
        res = Interval(start, stop)
        return res

    @cacheit
    def __getitem__(self, i):
        # TODO: implement generator
        if self.is_out_of_range(i):
            return S.Zero
        else:
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

class SeqExpCauchyMul(SeqCauchyMul, Mul):
    """Product of Exponential Generation sequences.
    """
    @cacheit
    def __getitem__(self, i):
        if self.is_out_of_range(i):
            return S.Zero
        else:
            if i == S.Zero:
                return Mul(*(seq[i] for seq in self.args))
            else:
                # TODO: args > 2
                if len(self.args)==2:
                    c = []
                    a = self.args[0]
                    b = self.args[1]
                    # TODO: optimize the range (if a.start_index > 0)
                    for k in range(0, i+1):
                        c.append(a[k]*b[i-k]*binomial(i, k))
                    return Add(*tuple(c))
            return S.Zero


class SeqExpCauchyPow_Main(SeqCauchyPow):
    """
    Only for main sequence
    """

    @cacheit
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
                    wik = self[i-k]
                    c.append((k*exp - i + k)*w[k]*wik*bc) # recursion
                # TODO: optimize cancel
                r = (Add(*tuple(c))/i/w[S.Zero]).cancel()
                r = r*factorial(i)
                return r



class SeqExpCauchyPow(SeqCauchyPow):

    # TODO: implement for various kind of Sequences.
    # TODO: use SeqCauchyPow: a_n = b_n/n!
    @cacheit
    def __getitem__(self, i):
        if self.is_out_of_range(i):
            return S.Zero
        else:
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
