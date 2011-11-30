
from sympy.core.singleton import (Singleton, S)
from sympy.core import (Pow)
from sympy.functions.combinatorial.factorials import factorial

from taylorexpr import TaylorSeriesExpr

class TaylorSeries(TaylorSeriesExpr):
    """
    Examples:

    >>> from sympy.series import Sequence, TaylorSeries
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
    show_method ='series'

    def __new__(self, x, **kwargs):

        sequence = kwargs.get("sequence", None)

        obj = TaylorSeriesExpr.__new__(self, x, sequence)
        return obj

    @property
    def x(self):
        return self._args[0]

    @property
    def sequence(self):
        return self._args[1]

    def __getitem__(self, i):
        a =  self.sequence[i]
        if (a != S.Zero) and (i != 0):
            a = a / factorial(i) * Pow(self.x, i)
        return a


    def show(self, n=5, **kwargs):
        self.show_method = kwargs.get('method', 'series')
        self.show_n = n
        return self

    def _sympystr(self, printer, *args):
        if self.show_method=='series':
            s = self.sequence
            l = [self[i] for i in range(s.start_index, self.show_n + 1)]
            l = [i for i in l if i != S.Zero]
            l = [printer._print(i) for i in l]
            return " + ". join(l) + " + ... "
        else:
            return printer._print_Basic(self, *args)

    def __mul__(self, other):
        assert self.x == other.x

