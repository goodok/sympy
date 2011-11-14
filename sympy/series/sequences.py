from sympy.core import (Basic, Expr, S)
from sympy.core.numbers import Infinity, Zero
from sympy.core.power import Pow
from sympy.functions import factorial
from sympy.solvers.recurr import rsolve

class Sequence(Expr):
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
    >>> seq.end_index
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

    is_Sequence = True

    show_n = 5
    def __new__(cls, size, baselist = [], **kwargs):
        """Create a new Sequence instance out of something useful. """
        formula = kwargs.pop("formula", None)
        if formula:
            return cls._from_formula(size, formula, **kwargs)

        if baselist:
            return cls._from_baselist(size, baselist, **kwargs)

        function = kwargs.pop("function", None)
        if function:
            return cls._from_function(size, function, **kwargs)

        recurr = kwargs.pop("recurr", None)
        if recurr:
            return cls._from_recurr(size, recurr, **kwargs)


    @classmethod
    def new(cls, size):
        """Construct :class:`Sequence` instance from raw representation. """
        obj = Basic.__new__(cls)

        obj.size = size

        return obj

    @classmethod
    def from_formula(cls, size, formula, **kwargs):
        """Construct a sequence from a ``formula``. """
        return cls._from_formula(size, formula, **kwargs)

    @classmethod
    def _from_formula(cls, size, formula, **opt):
        """Construct a sequence from a ``formula``. """
        obj = Expr.__new__(cls, **opt)

        kind = "formula"
        baselist = []
        period = 0

        arglist = [size, kind, baselist, period, formula[0], formula[1]]
        obj._args = tuple(arglist)
        return obj

    @classmethod
    def from_function(cls, size, function, **kwargs):
        """Construct a sequence from a ``function``. """
        return cls._from_function(size, function, **kwargs)

    @classmethod
    def _from_function(cls, size, function, **opt):
        """Construct a sequence from a ``function``. """
        obj = Expr.__new__(cls, **opt)

        kind = "function"

        arglist = [size, kind, function]
        obj._args = tuple(arglist)
        return obj

    @classmethod
    def from_baselist(cls, size, baselist, **kwargs):
        """Construct a sequence from a ``baselist``. """
        return cls._from_baselist(size, baselist, **kwargs)

    @classmethod
    def _from_baselist(cls, size, baselist, **opt):
        """Construct a sequence from a ``baselist``. """
        obj = Expr.__new__(cls, **opt)

        baselist = baselist
        kind = opt.pop("kind", "finitlist")
        if kind=="periodical":
            period = len(baselist)
        else:
            period = 0
            assert size[1]- size[0] + 1 == len(baselist)

        arglist = [size, kind, baselist, period, None]
        obj._args = tuple(arglist)
        return obj

    @classmethod
    def from_recurr(cls, size, recurr, **kwargs):
        """Construct a sequence from a ``recurr``. """
        return cls._from_recurr(size, recurr, **kwargs)

    @classmethod
    def _from_recurr(cls, size, recurr, **opt):
        """Construct a sequence from a ``recurr``. """
        obj = Expr.__new__(cls, **opt)

        kind = "recurr"
        k = recurr[1].args[0]
        formula = rsolve(*recurr)

        arglist = [size, kind, recurr, None, k, formula]
        obj._args = tuple(arglist)
        return obj


    @property
    def kind(self):
        return self._args[1]

    @property
    def size(self):
        return self._args[0]

    @property
    def start_index(self):
        return self.size[0]

    @property
    def length(self):
        return self.end_index - self.start_index + 1

    @property
    def end_index(self):
        return self.size[1]

    @property
    def is_infinite(self):
        return isinstance(self.end_index, Infinity)

    @property
    def baselist(self):
        return self._args[2]

    @property
    def period(self):
        assert self.kind == "periodical"
        assert len(self.baselist) == self._args[3]
        return len(self.baselist)

    @property
    def k(self):
        assert (self.kind == "formula") or (self.kind == "recurr")
        return self._args[4]

    @property
    def formula(self):
        assert (self.kind == "formula") or (self.kind == "recurr")
        return self._args[5]

    @property
    def function(self):
        return self._args[2]

    @property
    def recurr(self):
        return self._args[2]

    def __getitem__(self, i):
        kind = self.kind
        if isinstance(i, slice):
            slc = i
            start = max(slc.start, self.start_index)
            if self.is_infinite:
                stop = slc.stop
            else:
                stop = min(slc.stop, self.stop_index)
            new_size = (start, stop)

            if (kind == "formula"):
                return self.from_formula(new_size, (self.k, self.formula))
        else:
            if i < self.start_index:
                return S.Zero

            if not self.is_infinite:
                if i > self.end_index:
                    return S.Zero

            if (kind == "formula") or (kind == "recurr"):
                return self.formula.subs(self.k, i)
            if kind == "function":
                return self.function(i)
            elif kind == "periodical":
                i = (i - self.start_index) % self.period
                return self.baselist[i]
            elif kind == "finitlist":
                i = (i - self.start_index)
                return self.baselist[i]

    def show(self, n=5):
        self.show_n = n
        return self

    def _sympystr(self, printer, *args):
        l = []

        if self.start_index > 1:
          l.append(printer._print(S.Zero))
          l.append("...")
        elif self.start_index ==1:
          l.append(printer._print(S.Zero))

        count = self.show_n
        if not self.is_infinite:
            count = min(count, self.length)

        l.extend([printer._print(self[i]) for i in range(self.start_index, self.start_index + count)])

        if self.is_infinite or (self.length > self.show_n):
            l.append("...")
        return "[" + ", ". join(l) + "]"


    @property
    def info(self):
        # TODO: connect with printer system properly
        class SequenceInfo(object):
            def __init__(self, sequence):
                self.s = sequence._info()
            def __str__(self):
                return self.s
        return SequenceInfo(self)

    def _info(self):
        r = self.__class__.__name__ + "((%s, %s)" % self.size
        kind = self.kind
        if kind == "formula":
            r += ", formula=(%s, %s)" % (self.k, self.formula)
        elif kind == "periodical":
             r += ', baselist=%s, kind="%s"' % (self.baselist, kind)
        elif kind == "finitlist":
            r += ', baselist=%s' % (self.baselist)
        elif kind == "recurr":
            r += ', recurr=(%s, %s, %s)' % (self.recurr)
        r += ")"
        return r


    @property
    def domain(self):
        """Get the ground domain of ``self``. """
        return self.get_domain()

    def get_domain(f):
        """Get the ground domain of ``f``. """
        return f.rep.dom


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

