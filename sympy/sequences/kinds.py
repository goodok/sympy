from sympy.core import (Basic, Expr)
from sympy.core.singleton import (Singleton, S)
from sympy.core.power import Pow
from sympy.functions import factorial
from sympy.solvers.recurr import rsolve
from sympy.core.numbers import ilcm
from sympy.core.sets import Interval
from sympy.core.symbol import Symbol, symbols

from expr import (SeqExpr, EmptySequence)

class Sequence(SeqExpr):
    """Represents an Sequence.

    Helper fabric class for sequences constructions.

    Examples
    ========

    >>> from sympy.sequences import Sequence, SeqFormula
    >>> from sympy import oo, S
    >>> from sympy.abc import k
    >>> from sympy.printing.pretty.pretty import pprint


    >>> seq = Sequence((3, oo), formula=(k, S(1)/k))
    >>> seq
    SeqFormula([3, oo), k, 1/k)
    >>> pprint(seq)
    [0, ..., 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 1/9, ...]

    >>> Sequence((0, oo), periodical = (1, 0))
    SeqPer([0, oo), (1, 0))

    >>> Sequence((0, 2), finitlist = (1, 2, 3))
    SeqList([0, 2], (1, 2, 3))

    >>> seq = Sequence((1, oo), formula=(k, S(1)/k))
    >>> seq[2]
    1/2
    >>> seq[3:]
    SeqFormula([3, oo), k, 1/k)
    >>> seq[:4]
    SeqFormula([1, 4], k, 1/k)
    >>> seq[3:6]
    SeqFormula([3, 6], k, 1/k)

    >>> seq = SeqFormula((3, 6), k, 1/k)
    >>> seq[8:10]
    EmptySequence()

    See also
    ========

    sympy.concrete.summataions.Sum

    """

    __slots__ = ['rep', 'gens']

    is_SequenceAtom = True

    show_n = 5
    def __new__(cls, interval=None, name=None, **kwargs):
        """Create a new Sequence instance out of something useful. """

        if interval==None:
            interval = Interval(S.Zero, S.Infinity)
        elif isinstance(interval, basestring):
            name = interval
            interval = Interval(S.Zero, S.Infinity)
        elif type(interval)== tuple:
            interval = Interval(interval[0], interval[1])

        if name is not None:
            return SequenceSymbol._from_args(interval, name, **kwargs)

        baselist = kwargs.pop("periodical", None)
        if baselist:
            return SeqPer._from_args(interval, baselist, **kwargs)

        baselist = kwargs.pop("finitlist", None)
        if baselist:
            return SeqList._from_args(interval, baselist, **kwargs)

        formula = kwargs.pop("formula", None)
        if formula:
            k = formula[0]
            formula = formula[1]
            return SeqFormula._from_args(interval, k, formula, **kwargs)

        function = kwargs.pop("function", None)
        if function:
            return SeqFunc._from_args(interval, function, **kwargs)

        recurr = kwargs.pop("recurr", None)
        if recurr:
            return SeqRecurr._from_args(interval, recurr, **kwargs)

        raise ValueError(kwargs)

class SequenceBase(SeqExpr):
    """
    Base (atomic) Sequence class.

    All other kinds of sequences (SeqPer, SeqFunc, ) inherited from it.
    This class rather is an Abstract class.
    """
    is_SequenceAtom = True

    def __new__(cls, *args, **kwarg):
        if len(args) and (type(args[0])== tuple):
            _args = list(args)
            _args[0] = Interval(args[0][0], args[0][1])
            args = tuple(_args)
        obj = SeqExpr.__new__(cls, *args)
        return obj

    def _hashable_content(self):
        return tuple(self._args)

    @property
    def interval(self):
        return self._args[0]

    @property
    def is_direct_calculated(self):
        return False

class SequenceSymbol(SequenceBase, Symbol):
    """Symbolic representation of a sequence.

    Examples
    ========

    >>> from sympy import Symbol
    >>> from sympy.sequences import SequenceSymbol
    >>> from sympy.printing.pretty.pretty import pprint

    >>> a = SequenceSymbol((0, 10), 'a')
    >>> a
    {a}
    >>> pprint(a)
    [a0, a1, a2, a3, a4, a5, a6, ...]

    >>> i = Symbol('i')
    >>> a[i]
    a[i]

    >>> b = SequenceSymbol((0, 3), 'b')
    >>> a + b
    {a} + {b}
    >>> pprint(a + b)
    [a0 + b0, a1 + b1, a2 + b2, a3 + b3, a4, a5, a6, ...]

    """

    is_commutative = True

    def __new__(cls, interval, name):
        obj = SequenceBase.__new__(cls, interval, name)
        return obj

    @classmethod
    def _from_args(cls, interval, name, **opt):
        """Construct a SequenceSymbol"""
        obj = SequenceSymbol.__new__(cls, interval, name)
        return obj

    @property
    def name(self):
        return self.args[1]

    def __getitem__(self, i, **kw_args):
        if isinstance(i, slice):
            new_interval = self.calc_interval_from_slice(i)
            if new_interval == S.EmptySet:
                return S.EmptySequence
            return self._from_args(new_interval, self.name)

        else:
            if self.is_out_of_range(i):
                return S.Zero
            if isinstance(i, int):
                i = S(i)
            return IndexedSequenceSymbol(self, i, **kw_args)

    def _eval_subs(self, old, new):
        if self==old:
            return new
        else:
            raise NotImplemented

    def __call__(self, *args):
        raise TypeError( "%s object is not callable"%self.__class__ )

    def _sympystr(self, printer):
        if False: # printer._settings["list_sequences"]:
            return SeqExpr._sympystr(self, printer)
        else:
            return "{" + printer._print_Symbol(self) + "}"

    def _latex(self, p):
        base = p._print_Symbol(self)
        return r"\left\{ %s \right\}" % (base)

class IndexedSequenceSymbol(Expr):

    is_commutative = True
    is_Atom = True  # for latex printer

    def __new__(cls, base, *args, **kw_args):
        return Expr.__new__(cls, base, *args, **kw_args)

    @property
    def base(self):
        return self.args[0]

    @property
    def indices(self):
        return self.args[1:]

    def _sympystr(self, p):
        indices = map(p.doprint, self.indices)
        l = ", ".join(indices)
        base = p._print_Symbol(self.base)
        return "%s[%s]" % (base, l)

    def _pretty(self, p):
        return p._print_Indexed(self.base.name, self.indices)

    def _latex(self, p):
        """
        >>> from sympy.sequences import Sequence
        >>> from sympy.printing.latex import print_latex
        >>> a = Sequence('a')
        >>> print_latex(a[7], mode="plain")
        a_{7}

        """
        indices = map(p._print, self.indices)
        l = ", ".join(indices)
        base = p._print_Symbol(self.base)
        return "%s_{%s}" % (base, l)

def abstract_sequences(names):
    """
    Transform strings into instances of :class:`SequenceSymbol` class.


        >>> from sympy.sequences import abstract_sequences

        >>> x, y, z = abstract_sequences('x,y,z')
        >>> a, b, c = abstract_sequences('a b c')
        >>> a, b, c = abstract_sequences('a:c')
        >>> (a*b*c)[0]
        a[0]*b[0]*c[0]
    """
    return symbols(names, cls=Sequence)

class SeqPer(SequenceBase):
    """
    Periodical sequence.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.sequences import Sequence, SeqPer
    >>> from sympy.printing.pretty.pretty import pprint

    >>> seq = SeqPer((0, oo), (1, 2))
    >>> pprint(seq)
    [1, 2, 1, 2, 1, 2, 1, ...]

    >>> seq = Sequence((0, oo), periodical = (1, 2, 3))
    >>> pprint(seq)
    [1, 2, 3, 1, 2, 3, 1, ...]

"""


    def __new__(cls, interval, baselist = None, **kwargs):
        """
        Create a new periodical sequence SeqPer instance out of something useful.
        """
        if not isinstance(baselist, tuple) and baselist is not None:
            baselist = tuple((baselist,))
        baselist = tuple(S(i) for i in baselist)
        obj = SequenceBase.__new__(cls, interval, baselist)
        return obj

    @classmethod
    def _from_args(cls, interval, baselist, **opt):
        """Construct a sequence from a ``function``. """
        obj = SeqPer.__new__(cls, interval, baselist)
        return obj

    @property
    def baselist(self):
        return self._args[1]

    @property
    def period(self):
        #TODO: handle the case when baselist== (1)
        return len(self.baselist)

    def __getitem__(self, i):
        if isinstance(i, slice):
            new_interval = self.calc_interval_from_slice(i)
            if new_interval == S.EmptySet:
                return S.EmptySequence

            # shift left circular the base list if its start index has changed
            new_baselist = self.baselist
            baselist = self.baselist
            shift = (new_interval.left - self.start_index) % self.period
            if shift != 0:
                new_baselist = baselist[shift:] + baselist[:shift]

            return self._from_args(new_interval, new_baselist)

        else:
            if self.is_out_of_range(i):
                return S.Zero
            i = (i - self.start_index) % self.period
            return self.baselist[i]

class SeqList(SequenceBase):
    """
    Finite list sequence.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.sequences import Sequence, SeqList
    >>> from sympy.printing.pretty.pretty import pprint

    >>> seq = Sequence((1, 4), finitlist=(1, 2, 3, 4))
    >>> pprint(seq)
    [0, 1, 2, 3, 4]

    >>> len(seq.baselist)
    4

    >>> seq = Sequence((0, 8), finitlist=(1, 2, 3, 4, 5, 6, 7, 8, 9))
    >>> pprint(seq)
    [1, 2, 3, 4, 5, 6, 7, ...]

    """

    def __new__(cls, interval, baselist = None, **kwargs):

        """Create a new finite list seuqence SeqList instance out of something useful. """
        assert len(baselist) == interval.sup - interval.inf + 1
        obj = SequenceBase.__new__(cls, interval, baselist)
        return obj

    @classmethod
    def _from_args(cls, interval, baselist, **opt):
        """Construct a sequence from a ``function``. """
        obj = SeqList.__new__(cls, interval, baselist)
        return obj

    @property
    def baselist(self):
        return self._args[1]

    def __getitem__(self, i):
        if isinstance(i, slice):
            new_interval = self.calc_interval_from_slice(i)
            if new_interval == S.EmptySet:
                return S.EmptySequence
            a = new_interval.inf - self.start_index
            b = new_interval.sup - self.start_index + 1
            new_baselist = self.baselist[a:b]
            return self._from_args(new_interval, new_baselist)

        else:
            if self.is_out_of_range(i):
                return S.Zero
            return self.baselist[i-self.start_index]

class SeqFormula(SequenceBase):

    # It is may bee depricated, since SeqFunction is present.

    """
    Sequence defined by formula.

    Examples
    ========

    >>> from sympy.sequences import Sequence, SeqFormula
    >>> from sympy import oo, S
    >>> from sympy.abc import k
    >>> from sympy.printing.pretty.pretty import pprint

    >>> seq = Sequence((3, oo), formula=(k, S(1)/k))
    >>> seq
    SeqFormula([3, oo), k, 1/k)

    >>> pprint(seq)
    [0, ..., 1/3, 1/4, 1/5, 1/6, 1/7, ...]

    """

    def __new__(cls, interval, k, formula, **kwargs):
        obj = SequenceBase.__new__(cls, interval, k, formula)
        return obj

    @classmethod
    def _from_args(cls, interval, k, formula, **opt):
        """Construct a sequence from a ``function``. """
        obj = SeqFormula.__new__(cls, interval, k, formula)
        return obj

    @property
    def k(self):
        return self._args[1]

    @property
    def formula(self):
        return self._args[2]

    def __getitem__(self, i):
        if isinstance(i, slice):
            new_interval = self.calc_interval_from_slice(i)
            if new_interval == S.EmptySet:
                return S.EmptySequence
            return self._from_args(new_interval, self.k, self.formula)

        elif self.is_out_of_range(i):
            return S.Zero
        else:
            return self.formula.subs(self.k, i)

class SeqFunc(SequenceBase):
    """
    Sequence defined by function.

    Examples
    ========

    >>> from sympy.sequences import Sequence, SeqFunc
    >>> from sympy import oo, S
    >>> from sympy.abc import k
    >>> from sympy.printing.pretty.pretty import pprint

    >>> f = lambda k: S(1)/k
    >>> seq = Sequence((1, oo), function = f)
    >>> pprint(seq)
    [0, 1, 1/2, 1/3, 1/4, 1/5, ...]

    """
    def __new__(cls, interval, function, **kwargs):
        obj = SequenceBase.__new__(cls, interval, function)
        return obj

    @classmethod
    def _from_args(cls, interval, function, **opt):
        """Construct a sequence from a ``function``. """
        obj = SeqFunc.__new__(cls, interval, function)
        return obj

    @property
    def function(self):
        return self._args[1]

    def __getitem__(self, i):
        if isinstance(i, slice):
            new_interval = self.calc_interval_from_slice(i)
            if new_interval == S.EmptySet:
                return S.EmptySequence
            return self._from_args(new_interval, self.function)

        elif self.is_out_of_range(i):
            return S.Zero
        else:
            return self.function(i)

class SeqRecurr(SequenceBase):
    """
    Sequence defined by recurrence formula.

    Examples
    ========

    >>> from sympy.sequences import Sequence, SeqRecurr
    >>> from sympy import oo, S, Function
    >>> from sympy.abc import n
    >>> from sympy.printing.pretty.pretty import pprint

    >>> y = Function("y")
    >>> eq = (n-1)*y(n+2) - (n**2+3*n-2)*y(n+1) + 2*n*(n+1)*y(n)

    >>> seq = Sequence((2, oo), recurr=(eq, y(n), { y(0):0, y(1):3 }) )
    >>> pprint(seq)
    [0, ..., 6, 6, -24, -264, -1968, -14736, -120192, ...]
    >>> seq.formula
    3*2**n - 3*n!

    See Also
    ========

    sympy.solvers.recurr.rsolve

    """
    def __new__(cls, interval, recurr, **kwargs):

        """Create a new periodical sequence SeqPer instance out of something useful. """
        k = recurr[1].args[0]
        formula = rsolve(*recurr)
        obj = SequenceBase.__new__(cls, interval, recurr, k, formula)
        return obj

    def _hashable_content(self):
        _args = self._args
        return (_args[0], _args[2], _args[3])

    @classmethod
    def _from_args(cls, interval, recurr, **opt):
        """Construct a sequence from a ``function``. """
        obj = SeqRecurr.__new__(cls, interval, recurr)
        return obj

    @property
    def recurr(self):
        return self._args[1]

    @property
    def k(self):
        return self._args[2]

    @property
    def formula(self):
        return self._args[3]


    def __getitem__(self, i):
        if isinstance(i, slice):
            new_interval = self.calc_interval_from_slice(i)
            if new_interval == S.EmptySet:
                return S.EmptySequence
            return self._from_args(new_size, self.recurr)

        elif self.is_out_of_range(i):
            return S.Zero
        else:
            return self.formula.subs(self.k, i)


