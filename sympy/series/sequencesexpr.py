from sympy import Expr, Symbol, Eq, Mul, Add, Pow, expand, sympify, Tuple
from sympy.core.basic import Basic
from sympy.core.singleton import (Singleton, S)
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.cache import cacheit
#from sympy.functions.elementary.miscellaneous import Min, Max
from sympy.core.numbers import (Infinity, Zero)
from sympy.core.sets import Interval


class SeqExpr(Expr):
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
        return SeqPow(self, other)
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

class EmptySequence(SeqExpr):
    """Represents the empty sequence."""

    __metaclass__ = Singleton

    is_EmptySequence = True

    def __getitem__(self, i):
        return S.Zero


class SeqAdd(SeqExpr, Add):
    """A Sum of Sequence Expressions."""

    def __new__(cls, *args):

        args = map(sequenceify, args)

        args = tuple([arg for arg in args if not arg.is_EmptySequence])

        if not all(arg.is_Sequence for arg in args):
            raise ValueError("Mix of Sequence and Scalar symbols")

        expr = Add.__new__(cls, *args)
        expr = sequenceify(expr)

        if expr.is_Mul:
            return SeqMul(*expr.args)
        return expr

    def _hashable_content(self):
        return tuple(sorted(self._args, key=hash))

    def _sympystr(self, printer, *args):
        return printer._print_Add(self)

    def as_ordered_terms(self, order=None):
        return self.args

    @property
    def interval(self):
        res = S.EmptySet
        for seq in self.args:
            res = res | seq.interval
        return res

    @property
    @cacheit
    def start_index(self):
        return self.interval._inf

    @property
    @cacheit
    def stop_index(self):
        return self.interval._sup

    def is_out_of_range(self, i):
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

    def __getitem__(self, i):
        if isinstance(i, slice):
            return SeqAdd(*(seq[i] for seq in self.args))
        else:
            return Add(*(seq[i] for seq in self.args))


class SeqMul(SeqExpr, Mul):
    """A Product of Sequence Expressions."""

    def __new__(cls, *args):

        # Check that the shape of the args is consistent
        seqs = [arg for arg in args if arg.is_Sequence]

        if any(arg.is_Sequence and arg.is_EmptySequence for arg in args):
            return S.EmptySequence

        #expr = sequenceify(Mul.__new__(cls, *args))
        expr = sequenceify(Mul.__new__(cls, *args))
        if expr.is_Add:
            return SeqAdd(*expr.args)
        if expr.is_Pow:
            return SeqPow(*expr.args)
        if not expr.is_Mul:
            return expr

        if any(arg.is_Sequence and arg.is_EmptySequence for arg in expr.args):
            return S.EmptySequence

        return expr


def sequences_only(expr):
    #return [sym for sym in expr.free_symbols if sym.is_Sequence]
    #return [sym for sym in expr.args if sym.is_Sequence]
    return []

def sequenceify(expr):
    """
    Recursively walks down an expression tree changing Expr's to SeqExpr's
    i.e. Add -> SeqAdd
         Mul -> SeqMul

    Only changes those Exprs which contain SequenceSymbols

    This function is useful when traditional SymPy functions which use Mul and
    Add are called on SeqExpressions. Examples flatten, expand, simplify...

    Calling sequenceify after calling these functions will reset classes back to
    their Sequence equivalents

    For internal use
    """
    if len(sequences_only(expr))==0: # No Sequence symbols present
        return expr

    class_dict = {Mul:SeqMul, Add:SeqAdd, SeqMul:SeqMul, SeqAdd:SeqAdd,
            Pow:SeqPow, SeqPow:SeqPow}

    if expr.__class__ not in class_dict.keys():
        return expr

    args = map(sequenceify, expr.args) # Recursively call down the tree

    return Basic.__new__(class_dict[expr.__class__], *args)

