from sympy.core import (Basic, Expr, S)
from sympy.core.numbers import (Infinity, Zero)
from sympy.core.power import Pow
from sympy.functions import factorial
from sympy.solvers.recurr import rsolve
from sympy.core.numbers import ilcm


from sequencesexpr import SeqExpr

class SequenceBase(SeqExpr):
    """
    Base (atomic) Sequence class.

    All other kindes of sequences (SeqPer, SeqFunc, ) inherited from it.
    This class rather is Abstract
    """
    is_SequenceAtom = True
    show_n = 5

    def __new__(cls, *args, **kwarg):
        obj = SeqExpr.__new__(cls, *args)
        return obj

    @property
    def size(self):
        return self._args[0]

    @property
    def start_index(self):
        return self.size[0]

    @property
    def stop_index(self):
        return self.size[1]

    @property
    def length(self):
        return self.stop_index - self.start_index + 1

    @property
    def is_infinite(self):
        return isinstance(self.stop_index, Infinity)

    @property
    def is_direct_calculated(self):
        return False

    def calc_size(self, slc):
        start_index = self.start_index
        stop_index = self.stop_index
        if slc.start == None:
            slc_start = 0
        else:
            slc_start = slc.start
        start = max(slc_start, start_index)
        if self.is_infinite:
            if slc.stop == None:
                stop = stop_index
            else:
                stop = slc.stop
        else:
            stop = min(slc.stop, stop_index)
        return (start, stop)

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

    def __new__(cls, size, baselist = None, **kwargs):

        """Create a new periodical seuqence SeqPer instance out of something useful. """

        obj = SequenceBase.__new__(cls, size, baselist)
        return obj

    @classmethod
    def _from_args(cls, size, baselist, **opt):
        """Construct a sequence from a ``function``. """
        obj = SeqPer.__new__(cls, size, baselist)
        return obj

    @property
    def baselist(self):
        return self._args[1]

    @property
    def period(self):
        return len(self.baselist)

    def __getitem__(self, i):
        if isinstance(i, slice):
            new_size = self.calc_size(i)

            # shift left circular the base list if its start index has changed
            new_baselist = self.baselist
            baselist = self.baselist
            shift = (new_size[0] - self.start_index) % self.period
            if shift != 0:
                new_baselist = baselist[shift:] + baselist[:shift]

            return self._from_args(new_size, new_baselist)

        else:
            if self.is_out_of_range(i):
                return S.Zero
            i = (i - self.start_index) % self.period
            return self.baselist[i]

class SeqList(SequenceBase):
    """
    Periodical sequence.
    """

    def __new__(cls, size, baselist = None, **kwargs):

        """Create a new periodical seuqence SeqPer instance out of something useful. """
        assert len(baselist) == size[1] - size[0] + 1
        obj = SequenceBase.__new__(cls, size, baselist)
        return obj

    @classmethod
    def _from_args(cls, size, baselist, **opt):
        """Construct a sequence from a ``function``. """
        obj = SeqList.__new__(cls, size, baselist)
        return obj

    @property
    def baselist(self):
        return self._args[1]

    def __getitem__(self, i):
        if isinstance(i, slice):
            new_size = self.calc_size(i)
            a = new_size[0] - self.start
            b = new_size[1] - self.start
            new_baselist = baselist[a:b]
            return self._from_args(new_size, new_baselist)

        else:
            if self.is_out_of_range(i):
                return S.Zero
            return self.baselist[i-self.start_index]

class SeqFormula(SequenceBase):
    """
    Sequence defined by formula.

    It is may bee depricated, since SeqFunction is.
    """

    def __new__(cls, size, k, formula, **kwargs):

        """Create a new periodical seuqence SeqPer instance out of something useful. """

        obj = SequenceBase.__new__(cls, size, k, formula)
        return obj

    @classmethod
    def _from_args(cls, size, k, formula, **opt):
        """Construct a sequence from a ``function``. """
        obj = SeqFormula.__new__(cls, size, k, formula)
        return obj

    @property
    def k(self):
        return self._args[1]

    @property
    def formula(self):
        return self._args[2]

    def __getitem__(self, i):
        if isinstance(i, slice):
            new_size = self.calc_size(i)
            return self._from_args(new_size, self.k, self.formula)

        elif self.is_out_of_range(i):
            return S.Zero
        else:
            return self.formula.subs(self.k, i)

class SeqFunc(SequenceBase):
    """
    Sequence defined by function.
    """
    def __new__(cls, size, function, **kwargs):

        """Create a new periodical seuqence SeqPer instance out of something useful. """

        obj = SequenceBase.__new__(cls, size, function)
        return obj

    @classmethod
    def _from_args(cls, size, function, **opt):
        """Construct a sequence from a ``function``. """
        obj = SeqFunc.__new__(cls, size, function)
        return obj

    @property
    def function(self):
        return self._args[1]

    def __getitem__(self, i):
        if isinstance(i, slice):
            new_size = self.calc_size(i)
            return self._from_args(new_size, self.function)

        elif self.is_out_of_range(i):
            return S.Zero
        else:
            return self.function(i)

class SeqRecurr(SequenceBase):
    """
    Sequence defined by function.
    """
    def __new__(cls, size, recurr, **kwargs):

        """Create a new periodical sequence SeqPer instance out of something useful. """
        k = recurr[1].args[0]
        formula = rsolve(*recurr)
        obj = SequenceBase.__new__(cls, size, recurr, k, formula)
        return obj

    @classmethod
    def _from_args(cls, size, recurr, **opt):
        """Construct a sequence from a ``function``. """
        obj = SeqRecurr.__new__(cls, size, recurr)
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
            new_size = self.calc_size(i)
            return self._from_args(new_size, self.recurr)

        elif self.is_out_of_range(i):
            return S.Zero
        else:
            return self.formula.subs(self.k, i)

class Sequence(SeqExpr):
    """Represents an Sequence.

    Constructions:

    Through formula (finite or infinite sequence):

    >>> from sympy import Sequence
    >>> from sympy import S, oo
    >>> from sympy.abc import k

    >>> seq = Sequence((3, oo), formula=(k, S(1)/k))
    >>> seq
    [0, ..., 1/3, 1/4, 1/5, 1/6, 1/7, ...]
    >>> seq.is_infinite
    True
    >>> seq.start_index
    3
    >>> seq.stop_index
    oo
    >>> seq.info
    Sequence((3, oo), formula=(k, 1/k))

    >>> seq[4:6]
    [0, ..., 1/4, 1/5, 1/6]


    >>> f = AbstractFunction("f")          #doctest: +SKIP
    >>> Sequence((1, oo), function=f)      #doctest: +SKIP
    [f(1), f(2), f(3), ...]

    >>> f = lambda k: S(1)/k
    >>> seq = Sequence((1, oo), function = f)
    >>> seq
    [0, 1, 1/2, 1/3, 1/4, 1/5, ...]


    By recurrent formula or with the help of .rsolve() (finite or infinite sequence)

    >>> Sequence((1, oo), recurr = (a(0)=1, a(i+1) = a(i) + 1 ))    #doctest: +SKIP
    [1, 2, 3, 4, 5, ...]

    Through periodical list (finite or infinite sequence)

    >>> Sequence((0, oo), baselist=[1, 2, 3, 4], kind="periodical")
    [1, 2, 3, 4, 1, ...]

    Through the list (only finite sequence)

    >>> Sequence((2, 5), [1, 2, 3, 4])      #doctest: +SKIP
    [..., 1, 2, 3, 4]

    Representations:

    By default Sequence printed as lists (if it possible) for human reading:

    >>> seq = Sequence((1, oo), formula=(k, S(1)/k))
    >>> seq
    [0, 1, 1/2, 1/3, 1/4, 1/5, ...]

    But internaly it is still Sequence object

    >>> type(seq)
    <class 'sympy.series.sequences.Sequence'>


    To see information about object

    >>> seq.info
    Sequence((1, oo), formula=(k, 1/k))

    it is has constructor form.


    See also:
        solvers.recurr.rsolve()
        concrete.summataions.Sum()

    """

    __slots__ = ['rep', 'gens']

    is_SequenceAtom = True

    show_n = 5
    def __new__(cls, size, **kwargs):
        """Create a new Sequence instance out of something useful. """

        baselist = kwargs.pop("periodical", None)
        if baselist:
            return SeqPer._from_args(size, baselist, **kwargs)

        baselist = kwargs.pop("finitlist", None)
        if baselist:
            return SeqList._from_args(size, baselist, **kwargs)

        formula = kwargs.pop("formula", None)
        if formula:
            k = formula[0]
            formula = formula[1]
            return SeqFormula._from_args(size, k, formula, **kwargs)

        function = kwargs.pop("function", None)
        if function:
            return SeqFunc._from_args(size, function, **kwargs)

        recurr = kwargs.pop("recurr", None)
        if recurr:
            return SeqRecurr._from_args(size, recurr, **kwargs)

        raise ValueError(kwargs)

class TaylorSeries(Basic):
    """
    Examples:

    >>> from sympy import Sequence, TaylorSeries
    >>> from sympy import S, oo
    >>> from sympy.abc import x, k
    >>> seq = Sequence((1, oo), formula = (k, S(1)/k))
    >>> seq
    [0, 1, 1/2, 1/3, 1/4, 1/5, ...]

    >>> TaylorSeries(x, sequence=seq)
    x + x**2/4 + x**3/18 + x**4/96 + x**5/600 + ...

    >>> seq = Sequence((0, oo), baselist = [1, 0], kind="periodical")
    >>> seq
    [1, 0, 1, 0, ...]

    >>> TaylorSeries(x, sequence=seq)
    1 + x**2/2 + x**4/24 + ...

    """

    show_n = 5
    def __new__(self, x, **kwargs):

        sequence = kwargs.get("sequence", None)

        obj = Basic.__new__(self, x, sequence)
        return obj

    @property
    def x(self):
        return self._args[0]

    @property
    def sequence(self):
        return self._args[1]

    def __getitem__(self,key):
        a =  self.sequence[key]
        if not isinstance(a, Zero) and key <> 0:
            a = a / factorial(key) * Pow(self.x, key)
        return a


    def show(self, n=5):
        self.show_n = n
        return self

    def _sympystr(self, printer, *args):
        s = self.sequence
        l = [self[i] for i in range(s.size[0], self.show_n + 1)]
        l = [i for i in l if not isinstance(i, Zero)]
        l = [printer._print(i) for i in l]
        return " + ". join(l) + " + ... "

    def __mul__(self, other):
        assert self.x == other.x

