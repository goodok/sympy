from sympy.core import (Basic, Expr, Add, Mul, Pow, S)
from sympy.core.symbol import Symbol
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.cache import cacheit

from sympy.series.sequences import SequenceSymbol
from sympy.series.sequencesexpr import SeqCoeffMul


"""
The helper module for Formal power series and Formal Taylor series.

It contains common classes and metrhods for them.
"""

################################################################################
#     Interfaces: expression operations, interval, printing,                   #
################################################################################


class SeriesExprOp(Expr):
    """
    Embedding operations to the core
    """

    _op_priority = 13.0

    is_Series = True
    is_Identity = False
    is_TaylorSeries = False

    show_n = 8

    def __neg__(self):
        return SeriesCoeffMul(S.NegativeOne, self)
    def __abs__(self):
        raise NotImplementedError

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return SeriesAdd(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return SeriesAdd(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return SeriesAdd(self, -other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return SeriesAdd(other, -self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return SeriesMul(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return SeriesMul(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        if other == -S.One:
            return Inverse(self)
        return SeriesPow(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise NotImplementedError("Sequence Power not defined")
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return SeriesMul(self, other**S.NegativeOne)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        raise NotImplementedError()

    __truediv__ = __div__
    __rtruediv__ = __rdiv__


class SeriesExprInterval(object):

    @property
    @cacheit
    def start_index(self):
        return self.interval._inf

    @property
    @cacheit
    def stop_index(self):
        return self.interval._sup

    @property
    @cacheit
    def is_infinite(self):
        return self.stop_index == S.Infinity

    def is_out_of_range(self, i):
        if isinstance(i, Symbol):
            return False
        if i < self.start_index:
            return True
        if not self.is_infinite:
            if i > self.stop_index:
                return True
        return False


class SeriesExprPrint(object):
    """
    Interface with methods for printers.
    """

    def _sympystr(self, printer, *args):
        if printer._settings["list_series"]:
            l = [self[i] for i in range(self.start_index, self.show_n + 1)]
            terms = [i for i in l if i != S.Zero]

            from sympy.printing.precedence import precedence
            PREC = 40
            l = []
            for term in terms:
                t = printer._print(term)
                if t.startswith('-'):
                    sign = "-"
                    t = t[1:]
                else:
                    sign = "+"
                if precedence(term) < PREC:
                    l.extend([sign, "(%s)"%t])
                else:
                    l.extend([sign, t])
            sign = l.pop(0)
            if sign=='+':
                sign = ""
            return sign + ' '.join(l) + " + ..."
        else:
            return printer._print_Basic(self, *args)

    def _pretty(self, printer, *args):
        from sympy.printing.pretty.stringpict import prettyForm, stringPict
        if printer._settings["list_series"]:
            # see pretty._print_Add()
            l = [self[i] for i in range(self.start_index, self.show_n + 1)]
            terms = [i for i in l if i != S.Zero]

            def pretty_negative(pform, index):
                """Prepend a minus sign to a pretty form. """
                if index == 0:
                    if pform.height() > 1:
                        pform_neg = '- '
                    else:
                        pform_neg = '-'
                else:
                    pform_neg = ' - '
                pform = stringPict.next(pform_neg, pform)
                return prettyForm(binding=prettyForm.NEG, *pform)

            pforms, indices = [], []
            for i, term in enumerate(terms):
                if term.is_Mul and term.as_coeff_mul()[0] < 0:
                    pform = printer._print(-term)
                    pforms.append(pretty_negative(pform, i))
                elif term.is_Rational and term.q > 1:
                    pforms.append(None)
                    indices.append(i)
                elif term.is_Number and term < 0:
                    pform = printer._print(-term)
                    pforms.append(pretty_negative(pform, i))
                else:
                    pforms.append(printer._print(term))

            if indices:
                large = True
                for pform in pforms:
                    if pform is not None and pform.height() > 1:
                        break
                else:
                    large = False

                for i in indices:
                    term, negative = terms[i], False
                    if term < 0:
                        term, negative = -term, True
                    if large:
                        pform = prettyForm(str(term.p))/prettyForm(str(term.q))
                    else:
                        pform = printer._print(term)
                    if negative:
                        pform = pretty_negative(pform, i)
                    pforms[i] = pform

            pforms.append(prettyForm("..."))
            return prettyForm.__add__(*pforms)

        else:
            return printer._print_Basic(self, *args)


class SeriesExpr(SeriesExprOp, SeriesExprInterval, SeriesExprPrint):
    def _hashable_content(self):
        return tuple(self._args)

    def coeff(self, i):
        return self.sequence[i]


class SeriesAtom(SeriesExpr):
    def __new__(self, x, **kwargs):
        sequence = kwargs.get("sequence", None)
        obj = SeriesExpr.__new__(self, x, sequence)
        return obj

    @property
    def x(self):
        return self._args[0]

    @property
    def sequence(self):
        return self._args[1]

    @property
    def interval(self):
        return self.sequence.interval

    def coeff(self, i):
        return self.sequence[i]

    def _sympystr(self, printer, *args):
        if printer._settings["list_series"]:
            s = self.sequence
            l = [self[i] for i in range(s.start_index, s.start_index + self.show_n + 1)]
            l = [i for i in l if i != S.Zero]
            l = [printer._print(i) for i in l]
            return " + ". join(l) + " + ..."
        else:
            if isinstance(self.sequence, SequenceSymbol):
                l = [printer._print_Symbol(self.x), printer._print_Basic(self.sequence)]
                r = self.__class__.__name__ + "(%s)" % ", ".join(l)
            else:
                r = printer._print_Basic(self, *args)
            return r


################################################################################
#                             Operations                                       #
################################################################################

class SeriesAdd(SeriesExpr, Add):

    def _hashable_content(self):
        return tuple(sorted(self._args, key=hash))

    @property
    def x(self):
        return self._args[0].x

    @classmethod
    def flatten(cls, args_seq):
        return args_seq, [], None

    def as_ordered_terms(self, order=None):
        return self.args

    @property
    def interval(self):
        res = S.EmptySet
        for ts in self.args:
            res = res | ts.interval
        return res

    def __getitem__(self, i):
        if isinstance(i, slice):
            pass
        else:
            return Add(*(ts[i] for ts in self.args))

    def _sympystr(self, printer, *args):
        if printer._settings["list_series"]:
            return SeriesExprPrint._sympystr(self, printer, *args)
        else:
            return printer._print_Add(self)


class SeriesMul(SeriesExpr, Mul):
    def _hashable_content(self):
        return tuple(sorted(self._args, key=hash))

    @classmethod
    def flatten(cls, args_seq):
        return args_seq, [], None

    def _sympystr(self, printer, *args):
        if printer._settings["list_series"]:
            return SeriesExprPrint._sympystr(self, printer, *args)
        else:
            return printer._print_Mul(self)

class SeriesCoeffMul(SeriesExpr, Mul):
    def __new__(cls, coeff, ts):
        expr = Mul.__new__(cls, coeff, ts)
        return expr

    @classmethod
    def flatten(cls, args_ts):
        coeff = args_ts[0]
        subseries = args_ts[1]
        if isinstance(subseries, SeriesCoeffMul):
            coeff *= subseries.coeffitient
            series = subseries.series
            args_ts = [coeff, series]
        return args_ts, [], None

    @property
    def coeffitient(self):
        return self.args[0]

    @property
    def series(self):
        return self.args[1]

    @property
    def interval(self):
        return self.series.interval

    @property
    @cacheit
    def sequence(self):
        return SeqCoeffMul(self.coeffitient, self.series.sequence)



    def _sympystr(self, printer, *args):
        if printer._settings["list_series"]:
            return SeriesExprPrint._sympystr(self, printer, *args)
        else:
            return printer._print_Mul(self)




