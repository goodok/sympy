from sympy.core import (Basic, Expr, Add, Mul, Pow, S)
from sympy.core.symbol import Symbol
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.cache import cacheit
from sympy.core.sets import Interval

from sympy.sequences import Sequence, SequenceSymbol
from sympy.sequences.expr import SeqCoeffMul, SeqSliced


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

    #TODO: join dublicated code with SeqExprInterval
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

    def calc_interval_from_slice(self, slc):
        slc_start = slc.start
        if slc_start == None:
            slc_start = S.Zero
        slc_stop  = slc.stop
        if slc_stop == None:
            slc_stop = S.Infinity
        return self.interval & Interval(slc_start, slc_stop)

    def getitem_dispatche(self, i):
        if isinstance(i, slice):
            return self.getitem_slicing(i)
        elif self.is_out_of_range(i):
            return S.Zero
        else:
            return self.getitem_index(i)

    def getitem_slicing(self, i):
        mask_interval = self.calc_interval_from_slice(i)
        return SeriesSliced(self, mask_interval)


class SeriesExprPrint(object):
    """
    Interface with methods for printers.
    """

    def _sympystr(self, printer, *args):
        if printer._settings["list_series"]:
            l = [self[i] for i in xrange(self.start_index, self.start_index + self.show_n + 1)]
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
            if isinstance(self.sequence, SequenceSymbol):
                l = [printer._print_Symbol(self.x), printer._print_Basic(self.sequence)]
                r = self.__class__.__name__ + "(%s)" % ", ".join(l)
            else:
                r = printer._print_Basic(self, *args)
            return r

    def _pretty(self, printer, *args):
        from sympy.printing.pretty.stringpict import prettyForm, stringPict
        if printer._settings["list_series"]:
            # see pretty._print_Add()
            l = [self[i] for i in xrange(self.start_index, self.show_n + 1)]
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
    def __new__(self, x, sequence_name=None, **kwargs):
        if sequence_name:
            sequence = Sequence(sequence_name, **kwargs)
        else:
            sequence = kwargs.get("sequence", None)
            if sequence==None:
                sequence = Sequence(**kwargs)
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


class SeriesSliced(SeriesExpr):
    """
    Return sliced series.

    When series is simple, then trivially we can recreated and change interval.
    But when expression of series is complex (like multiplication) then we
    wrap this expression and mask original interval of it.

    Examples
    ========

    >>> from sympy import oo
    >>> from sympy.abc import x
    >>> from sympy.series import PowerSeries
    >>> from sympy.printing.pretty.pretty import pprint

    >>> a = PowerSeries(x, periodical=(5, 7))
    >>> a[2:5]
    5*x**2 + 7*x**3 + 5*x**4 + 7*x**5 + ...

    >>> b = PowerSeries(x, periodical=(1, 1))
    >>> c = a*b

    >>> c
    5 + 12*x + 17*x**2 + 24*x**3 + 29*x**4 + ...
    >>> c.interval
    [0, oo)

    >>> c[2:5]
    17*x**2 + 24*x**3 + 29*x**4 + 36*x**5 + ...
    >>> c[2:5].interval
    [2, 5]

    See Also
    ========

    sympy.series.seqencesexpr.SeqSliced

    """
    # TODO: join with dublicated code of SeqSliced
    def __new__(cls, original, mask_interval):
        assert original.is_Series
        obj = SeriesExpr.__new__(cls, original, mask_interval)
        return obj

    @property
    def original(self):
        return self._args[0]

    @property
    def mask_interval(self):
        return self._args[1]

    @property
    def interval(self):
        return self.mask_interval & self.original.interval

    def __getitem__(self, i):
        if self.is_out_of_range(i):
            return S.Zero
        else:
            return  self.original[i]

    @property
    def sequnence(self):
        return SeqSliced(self.original.sequence, self.mask_interval)


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
        return self.getitem_dispatche(i)

    def getitem_index(self, i):
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
    def check_zero(cls, args):

        if any(arg.is_zero for arg in args):
            return True

        # if at least one sequence is empty then result is EmptySequence
        # TODO: test it
        if any(ser.sequence.is_EmptySequence for ser in args if ser.is_Series):
            return True

    @classmethod
    def carry_out_coeff(cls, args):

        # collect scalar coefficients
        # and collect sequenses (without coefficients)
        # Here we preserve order for noncommutative cases of coefficients, but
        # consider only simple cases: when sequences are commutative itself
        # TODO: use generators

        coeffs = []
        sers = []
        for arg in args:
            if not arg.is_Series:
                coeffs.append(arg)
            elif isinstance(arg, SeriesCoeffMul):
                coeffs.append(arg.coefficient)
                sers.append(arg.series)
            else:
                sers.append(arg)

        # calculate the multiplicity of coefficients
        if coeffs==[]:
            coeff = S.One
        else:
            coeff = Mul(*coeffs)
        return coeff, sers


    @classmethod
    def flatten(cls, args):
        #TODO: how to use AssocOp.flatten(cls, args)
        new_seq = []
        while args:
            o = args.pop()
            if o.__class__ is cls: # classes must match exactly
                args.extend(o.args)
            else:
                new_seq.append(o)
        # c_part, nc_part, order_symbols
        return [], new_seq, None

    @property
    def x(self):
        return self.args[0].x

    @property
    @cacheit
    def interval(self):
        start = Add(*(s.start_index for s in self.args))
        stop = Add(*(s.stop_index for s in self.args))
        res = Interval(start, stop)
        return res


    def __getitem__(self, i):
        return self.getitem_dispatche(i)

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
            coeff *= subseries.coefficient
            series = subseries.series
            args_ts = [coeff, series]
        return args_ts, [], None

    @property
    def coefficient(self):
        return self.args[0]

    @property
    def series(self):
        return self.args[1]

    @property
    def x(self):
        return self.series.x

    @property
    def interval(self):
        return self.series.interval

    @property
    @cacheit
    def sequence(self):
        return SeqCoeffMul(self.coefficient, self.series.sequence)

    def _sympystr(self, printer, *args):
        if printer._settings["list_series"]:
            return SeriesExprPrint._sympystr(self, printer, *args)
        else:
            return printer._print_Mul(self)

class SeriesNested(SeriesExpr):
    def __new__(cls, *args):
        expr = SeriesExpr.__new__(cls, *args)
        return expr

    @property
    def g(self):
        return self.args[0]

    @property
    def f(self):
        return self.args[1]

    @property
    def x(self):
        return self.g.x

    @property
    @cacheit
    def sequence(self):
        # abstract
        return None

    def __getitem__(self, i):
        return self.getitem_dispatche(i)

    @property
    @cacheit
    def interval(self):
        return self.f.interval


