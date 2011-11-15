from sympy import Expr, Symbol, Eq, Mul, Add, Pow, expand, sympify, Tuple
from sympy.core.basic import Basic
from sympy.core.singleton import S
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.cache import cacheit


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
    is_ZeroSequence = False
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


class SeqAdd(SeqExpr, Add):
    """A Sum of Sequence Expressions

    SeqAdd inherits from and operates like SymPy Add

    >>> from sympy import SeqAdd, SequenceSymbol
    >>> A = SequenceSymbol('A', 5, 5)
    >>> B = SequenceSymbol('B', 5, 5)
    >>> C = SequenceSymbol('C', 5, 5)
    >>> SeqAdd(A, B, C)
    A + B + C
    """


    def __new__(cls, *args):

        args = map(sequenceify, args)

        #TODO: is it correct, to check arg!=0? args must be Expr type
        args = [arg for arg in args if arg!=0]

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




class SeqMul(SeqExpr, Mul):
    """A Product of Sequence Expressions

    SeqMul inherits from and operates like SymPy Mul

    >>> from sympy import SeqMul, SequenceSymbol
    >>> A = SequenceSymbol('A', 5, 4)
    >>> B = SequenceSymbol('B', 4, 3)
    >>> C = SequenceSymbol('C', 3, 6)
    >>> SeqMul(A, B, C)
    A*B*C
    """

    def __new__(cls, *args):

        # Check that the shape of the args is consistent
        seqs = [arg for arg in args if arg.is_Sequence]

        if any(arg.is_zero for arg in args):
            return ZeroSequence()

        #expr = sequenceify(Mul.__new__(cls, *args))
        expr = sequenceify(Mul.__new__(cls, *args))
        if expr.is_Add:
            return SeqAdd(*expr.args)
        if expr.is_Pow:
            return SeqPow(*expr.args)
        if not expr.is_Mul:
            return expr

        if any(arg.is_Sequence and arg.is_ZeroSequence for arg in expr.args):
            return ZeroSequence()

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

def linear_factors(expr, *syms):
    """Reduce a Sequence Expression to a sum of linear factors

    Given symbols and a Sequence expression linear in those symbols return a
    dict mapping symbol to the linear factor

    >>> from sympy import SequenceSymbol, linear_factors, symbols
    >>> n, m, l = symbols('n m l')
    >>> A = SequenceSymbol('A', n, m)
    >>> B = SequenceSymbol('B', m, l)
    >>> C = SequenceSymbol('C', n, l)
    >>> linear_factors(2*A*B + C, B, C)
    {B: 2*A, C: I}
    """

    expr = sequenceify(expand(expr))
    d = {}
    if expr.is_Sequence and expr.is_Symbol:
        if expr in syms:
            d[expr] = Identity(expr.n)

    if expr.is_Add:
        for sym in syms:
            total_factor = 0
            for arg in expr.args:
                factor = arg.coeff(sym)
                if not factor:
                    # .coeff fails when powers are in the expression
                    if sym in arg.free_symbols:
                        raise ValueError("Expression not linear in symbols")
                    else:
                        factor = 0
                factor = sympify(factor)
                total_factor += factor
            d[sym] = total_factor
    elif expr.is_Mul:
        for sym in syms:
            factor = expr.coeff(sym)
            if not factor:
                # .coeff fails when powers are in the expression
                if sym in expr.free_symbols:
                    raise ValueError("Expression not linear in symbols")
                else:
                    factor = 0
            factor = sympify(factor)
            d[sym] = factor

    if any(sym in Sequence_symbols(Tuple(*d.values())) for sym in syms):
        raise ValueError("Expression not linear in symbols")

    return d

