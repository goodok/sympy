from sympy.core import (Basic, Expr)
from sympy.core.singleton import (Singleton, S)
from sympy.core.power import Pow
from sympy.functions import factorial
from sympy.solvers.recurr import rsolve
from sympy.core.numbers import ilcm
from sympy.core.sets import Interval
from sympy.core.symbol import Symbol

from sequencesexpr import (SeqExpr, EmptySequence)

class SequenceBase(SeqExpr):
    """
    Base (atomic) Sequence class.

    All other kindes of sequences (SeqPer, SeqFunc, ) inherited from it.
    This class rather is Abstract
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
    """Symbolic representation of a Sequence object"""

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
            return IndexedSequenceSymbol(self, i, **kw_args)

    def _eval_subs(self, old, new):
        if self==old:
            return new
        else:
            raise NotImplemented

    def __call__(self, *args):
        raise TypeError( "%s object is not callable"%self.__class__ )

class IndexedSequenceSymbol(Expr):

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
        return "%s[%s]" % (p.doprint(self.base), ", ".join(indices))


class SeqPer(SequenceBase):
    """
    Periodical sequence.
    """

    def __new__(cls, interval, baselist = None, **kwargs):

        """Create a new periodical sequence SeqPer instance out of something useful. """

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
    Periodical sequence.
    """

    def __new__(cls, interval, baselist = None, **kwargs):

        """Create a new periodical seuqence SeqPer instance out of something useful. """
        assert len(baselist) == interval.right - interval.left + 1
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
            a = new_interval.left - self.start_index
            b = new_interval.right - self.start_index + 1
            new_baselist = self.baselist[a:b]
            return self._from_args(new_interval, new_baselist)

        else:
            if self.is_out_of_range(i):
                return S.Zero
            return self.baselist[i-self.start_index]

class SeqFormula(SequenceBase):
    """
    Sequence defined by formula.

    It is may bee depricated, since SeqFunction is.
    """

    def __new__(cls, interval, k, formula, **kwargs):

        """Create a new periodical seuqence SeqPer instance out of something useful. """

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
    """
    def __new__(cls, interval, function, **kwargs):

        """Create a new periodical seuqence SeqPer instance out of something useful. """

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
    Sequence defined by function.
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

class Sequence(SeqExpr):
    """Represents an Sequence.

    Helper fabric class for sequences constructions.

    See also:
        solvers.recurr.rsolve()
        concrete.summataions.Sum()

    """

    __slots__ = ['rep', 'gens']

    is_SequenceAtom = True

    show_n = 5
    def __new__(cls, interval, name=None, **kwargs):
        """Create a new Sequence instance out of something useful. """

        if type(interval)== tuple:
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

