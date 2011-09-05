from sympy import Basic
from sympy import S
from sympy.core.numbers import Infinity, Zero
from sympy.core.power import Pow
from sympy.functions import factorial


class Sequence(Basic):
    """Represents an unevaluated limit.

    Examples:

    
    through formula (finite or infinite sequence):
 
    >>> from sympy import Sequence
    >>> from sympy import S, oo
    >>> from sympy.abc import k
    >>> Sequence((1, oo), k, S(1)/k)
    [1, 1/2, 1/3, ...]


    by recurrent formula or with the help of .rsolve() (finite or infinite sequence)

    >> Sequence((1, oo), rsolve(a(0)=1, a(i+1) = a(i) + 1 ) (see rsolve interface)
    [1, 2, 3, 4, ...]

    through periodical list (finite or infinite sequence)
 
    >>> Sequence((0, oo), baselist = [1, 2, 3, 4], kind="periodical")
    [1, 2, 3, 4, 1, ...]
    
    through list (only finite sequence)
 
    >> Sequence((2, 5), [1, 2, 3, 4])
    [..., 1, 2, 3, 4]

    """
    show_n = 5
    def __new__(self, size, k=None, func=None, **kwargs):
        
        kind = kwargs.get("kind", "normal")
        baselist = kwargs.get("baselist", [])

        obj = Basic.__new__(self, size, kind, baselist, k, func)
        return obj

    @property
    def size(self):
        return self._args[0]

    @property
    def kind(self):
        return self._args[1]

    @property
    def baselist(self):
        return self._args[2]

    @property
    def period(self):
        return len(self.baselist)
    @property
    def k(self):
        return self._args[3]

    @property
    def func(self):
        return self._args[4]


    def _sympystr(self, printer, *args):
        l = [printer._print(self[i]) for i in range(self.size[0], self.show_n + 1)]
        return "[" + ", ". join(l) + ", ... ]"
            
    def __getitem__(self,key):
        if key < self.size[0]:
            return S.Zero
        if not isinstance(self.size[1], Infinity):
            if key > self.size[1]:
                return S.Zero
        if self.kind == "normal":
            return self.func.subs(self.k, key)
        elif self.kind == "periodical":
            i = (key - self.size[0]) % self.period
            return self.baselist[i]
            
    def show(self, n=5):
        self.show_n = n
        return self


class TaylorSeries(Basic):

    """

    Examples:

    >>> from sympy import Sequence, TaylorSeries
    >>> from sympy import S, oo
    >>> from sympy.abc import x, k
    >>> s = Sequence((1, oo), k, S(1)/k)
    >>> s
    [1, 1/2, 1/3, ...]

    >>> TaylorSeries(x, sequence=s)
    x + x**2/4 + x**3/18 + x**4/96 + x**5/600 + ... 

    >>> s = Sequence((0, oo), baselist = [1, 0], kind="periodical")
    >>> s
    [1, 0, 1, 0, ...]

    >>> TaylorSeries(x, sequence=s)
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
        
                
