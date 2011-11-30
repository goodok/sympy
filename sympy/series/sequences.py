from sympy.core import (Basic, Expr)
from sympy.core.singleton import (Singleton, S)
from sympy.core.power import Pow
from sympy.functions import factorial
from sympy.solvers.recurr import rsolve
from sympy.core.numbers import ilcm
from sympy.core.sets import Interval


from sequencesexpr import (SeqExpr, EmptySequence)

class SequenceBase(SeqExpr):
    """
    Base (atomic) Sequence class.

    All other kindes of sequences (SeqPer, SeqFunc, ) inherited from it.
    This class rather is Abstract
    """
    is_SequenceAtom = True
    show_n = 5

    def __new__(cls, *args, **kwarg):
        if len(args) and (type(args[0])== tuple):
            _args = list(args)
            _args[0] = Interval(args[0][0], args[0][1])
            args = tuple(_args)
        obj = SeqExpr.__new__(cls, *args)
        return obj

    @property
    def interval(self):
        return self._args[0]

    @property
    def start_index(self):
        return self.interval.left

    @property
    def stop_index(self):
        return self.interval.right

    @property
    def length(self):
        return self.stop_index - self.start_index + 1

    @property
    def is_infinite(self):
        return self.stop_index == S.Infinity

    @property
    def is_direct_calculated(self):
        return False

    def calc_interval_from_slice(self, slc):
        slc_start = slc.start
        if slc_start == None:
            slc_start = S.Zero
        slc_stop  = slc.stop
        if slc_stop == None:
            slc_stop = S.Infinity
        return self.interval & Interval(slc_start, slc_stop)

    def is_out_of_range(self, i):
        if i < self.start_index:
            return True

        if not self.is_infinite:
            if i > self.stop_index:
                return True

        return False


    def _pretty(self,  printer, *args):
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

        if self.is_infinite or (self.length > self.show_n):
            printset.append("...")
        return printer._print_seq(printset, '[', ']', ', ' )


class SeqPer(SequenceBase):
    """
    Periodical sequence.
    """

    def __new__(cls, interval, baselist = None, **kwargs):

        """Create a new periodical seuqence SeqPer instance out of something useful. """

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
            a = new_interval.left - self.start
            b = new_interval.rigth - self.start
            new_baselist = baselist[a:b]
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
    def __new__(cls, interval, **kwargs):
        """Create a new Sequence instance out of something useful. """

        if type(interval)== tuple:
            interval = Interval(interval[0], interval[1])

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

