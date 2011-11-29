
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

