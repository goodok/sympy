from sympy.core import (Basic, Expr, Add, Mul, Pow, S)
from sympy.core.symbol import Symbol
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.cache import cacheit
from sympy.core.sets import Interval, EmptySet

from sympy.sequences import Sequence, SequenceSymbol
from sympy.sequences.expr import (SeqAdd, SeqExprInterval, SeqExprPrint,
                        SeqSliced, SeqCoeffMul)


"""
The helper module for formal power series.

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
    is_SeriesAtom = False
    is_SeriesGen = False
    is_Identity = False

    is_PowerSeries = False
    is_PowerSeries_0E = False
    is_PowerSeries0 = False


    _type_must = "Series"
    _type_is = "Series"


    show_n = 8

    def __neg__(self):
        return self._SeriesCoeffMul(S.NegativeOne, self)
    def __abs__(self):
        raise NotImplementedError

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return self._SeriesAdd(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return self._SeriesAdd(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return self._SeriesAdd(self, -other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return self._SeriesAdd(other, -self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return self._SeriesMul(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return self._SeriesMul(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        # if other == -S.One: return PowerSeriesInverse(self)
        return self._SeriesPow(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise NotImplementedError("Series Power not defined")
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return self._SeriesMul(self, other**S.NegativeOne)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        raise NotImplementedError()

    def __lshift__(self, other):
        "Overloading for <<"
        return self.shift(-other)

    def __rshift__(self, other):
        "Overloading for >>"
        return self.shift(other)

    __truediv__ = __div__
    __rtruediv__ = __rdiv__


    @property
    def _SeriesAdd(self): return SeriesAdd

    @property
    def _SeriesMul(self): return SeriesMul

    @property
    def _SeriesPow(self): return SeriesPow

    @property
    def _SeriesSliced(self): return SeriesSliced

    @property
    def _SeriesExpr(self): return SeriesExpr

    @property
    def _SeriesNested(self): return SeriesNested

    @property
    def _SeriesMul(self): return SeriesMul

    @property
    def _SeriesCoeffMul(self): return SeriesCoeffMul

    @classmethod
    def _cls_SeriesCoeffMul(cls): return SeriesCoeffMul

    @classmethod
    def _cls_SeriesMul(cls): return SeriesMul

    def __repr__(self):
        from sympy.printing import srepr
        return srepr(self)


class SeriesExprInterval(SeqExprInterval):

    #TODO: join dublicated code with SeqExprInterval

    @property
    @cacheit
    def interval(self):
        return self.sequence.interval

class SeriesExprPrint(SeqExprPrint):
    """
    Interface with methods for printers.
    """

    def _print_unevualated_power(self, x, i):
        return Pow(x, i)

    def _sympystr(self, printer, *args):
        #if printer._settings["list_series"]:
        if True:

            s = self.sequence
            l = [(s[i], i) for i in xrange(self.start_index, self.start_index + self.show_n + 1)]
            l = [(c, i) for (c, i) in l if c != S.Zero]
            terms = [c* self._print_unevualated_power(self.x, i) for (c, i) in l]

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

    def _sympyrepr(self, printer, *args):
        r = printer.emptyPrinter(self)
        return r

    def _pretty(self, printer, *args):
        from sympy.printing.pretty.stringpict import prettyForm, stringPict
        #if printer._settings["list_series"]:
        if True:
            # see pretty._print_Add()

            s = self.sequence
            l = [(s[i], i) for i in xrange(self.start_index, self.start_index + self.show_n + 1)]
            l = [(c, i) for (c, i) in l if c != S.Zero]
            terms = [c* self._print_unevualated_power(self.x, i) for (c, i) in l]

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

    def _latex(self, printer):
        # TODO: work out with sorting in better way

        s = self.sequence
        l = [(s[i], i) for i in xrange(self.start_index, self.start_index + self.show_n + 1)]
        l = [(c, i) for (c, i) in l if c != S.Zero]

        (c, i) = l[0]
        term = c*self._print_unevualated_power(self.x, i)
        tex = printer._print(term)
        for (c, i) in l[1:]:
            term = c*self._print_unevualated_power(self.x, i)
            next = printer._print(term)
            if next[0]=='-':
                pass
            else:
                tex += " +"
            tex += " " + next
        tex += " + \dotsb"
        return tex

class SeriesExpr(SeriesExprOp, SeriesExprInterval, SeriesExprPrint):
    def _hashable_content(self):
        return tuple(self._args)

    def getitem_slicing(self, i):
        # abstract
        pass
        mask = self.calc_interval_from_slice(i)
        return self._SeriesSliced(self, mask)

    @classmethod
    def _convert_from_scalar(cls, scalar, example):
        res = SeriesAtom(example.x, sequence=Sequence((0, 0), finitlist=(scalar,)))
        return res

    #TODO: if uncommented then exception
    #@cacheit
    #def getitem_index(self, i):
        # abstract
        #c = self.sequence
        #return c[i]*Pow(self.x, i)

    def coeff(self, i):
        # TODO: the name-token override `Expr.coeff`
        """
        Return the i-th coefficient, that is near the term of i power.
        """
        return self.sequence[i]

    def compose(self, other):
        return self._SeriesNested(self, other)

    def reverse(self):
        return self._Reverse(self)

    def shift(self, n):
        """
        >>> from sympy.series import PowerSeries
        >>> from sympy.abc import x
        >>> from sympy.printing.repr import srepr
        >>> ps = PowerSeries(x, periodical=(1, 2, 3, 4, 5, 6, 7))
        >>> ps
        1 + 2*x + 3*x**2 + 4*x**3 + 5*x**4 + 6*x**5 + 7*x**6 + x**7 + ...

        >>> ps.shift(-2)
        3 + 4*x + 5*x**2 + 6*x**3 + 7*x**4 + x**5 + 2*x**6 + 3*x**7 + 4*x**8 + ...

        """
        new_seq = self.sequence.shift(n)
        return self._from_args(self.x, new_seq)

class SeriesGen(SeriesExpr):
    is_SeriesGen = True
    # TODO: redefine
    def __new__(cls, x, **kwargs):
        return SeriesAtom(x, sequence=Sequence((1, 1), finitlist=(1,)) )


class SeriesAtom(SeriesExpr):
    is_SeriesAtom = True
    def __new__(cls, x=None, sequence_name=None, **kwargs):
        if sequence_name:
            sequence = Sequence(sequence_name, **kwargs)
        else:
            poly = kwargs.get("poly", None)
            if poly:
                return cls._from_poly(poly) # recur
            else:
                sequence = kwargs.get("sequence", None)
                if sequence==None:
                    sequence = Sequence(**kwargs)
        assert x
        assert sequence
        obj = SeriesExpr.__new__(cls, x, sequence)
        return obj

    @property
    def x(self):
        return self._args[0]

    @property
    def sequence(self):
        return self._args[1]

    @classmethod
    def _from_args(cls, x, sequence, **kwargs):
        return cls.__new__(cls, x, sequence=sequence, **kwargs)

    @classmethod
    def _from_poly(cls, poly, **kwargs):
        assert poly.is_Poly
        assert poly.is_univariate
        x = poly.gen
        end = poly.degree()
        start = poly.monoms()[-1][0]
        coeffs = poly.all_coeffs()
        coeffs.reverse()
        coeffs = tuple(coeffs[start:])
        sequence = Sequence(Interval(start, end), finitlist=coeffs)
        return cls.__new__(cls, x, sequence=sequence, **kwargs)

class SeriesSliced(SeriesExpr, SeqSliced):
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
    def __new__(cls, original, mask):
        assert original.is_Series

        if isinstance(original, SeriesSliced):
            # Absorb nested
            mask = mask & original.mask
            original = original.original

        if isinstance(mask, EmptySet):
            return S.Zero

        obj = SeriesExpr.__new__(cls, original, mask)
        return obj

    @property
    def x(self):
        return self.original.x


    @property
    def sequence(self):
        return SeqSliced(self.original.sequence, self.mask)

################################################################################
#                             Operations                                       #
################################################################################

class SeriesAdd(SeriesExpr, Add):
    def __new__(cls, *args):

        #TODO: is it correct, to check arg!=0? args must be Expr type
        args = [arg for arg in args if arg!=0]

        #TODO: create ScalAdd to keep scalar separatly
        scalars = (arg for arg in args if not arg.is_Series)
        series = [arg for arg in args if arg.is_Series]
        if not all(arg._type_is == cls._type_must for arg in series):
            raise ValueError("Mix of Series and Scalar symbols")

        scalar = Add(*tuple(scalars))
        if scalar is not S.Zero:
            scalar_series = cls._convert_from_scalar(scalar, series[0])
            series[0:0] = [scalar_series]

        expr = Add.__new__(cls, *tuple(series))

        if expr.is_Mul:
            return cls._cls_SeriesMul()(*expr.args)
        return expr

    def _hashable_content(self):
        return tuple(sorted(self._args, key=hash))

    @property
    def x(self):
        return self._args[0].x

    @classmethod
    def flatten(cls, args):
        new_seq = []
        while args:
            o = args.pop()
            if o.__class__ is cls: # classes must match exactly
                args.extend(o.args)
            else:
                new_seq.append(o)
        new_seq.reverse()
        return new_seq, [], None

    def as_ordered_terms(self, order=None):
        return self.args

    @property
    def interval(self):
        res = S.EmptySet
        for ts in self.args:
            res = res | ts.interval
        return res

    def getitem_index(self, i):
        return Add(*(ts[i] for ts in self.args))

    @property
    @cacheit
    def sequence(self):
        return SeqAdd(*(s.sequence for s in self.args))

    def _sympystr(self, printer, *args):
        if True: #printer._settings["list_series"]:
            return SeriesExprPrint._sympystr(self, printer, *args)
        else:
            return printer._print_Add(self)


class SeriesMul(SeriesExpr, Mul):

    def __new__(cls, *args):
        if cls.check_zero(args):
            return S.Zero

        coeff, sers = cls.carry_out_coeff(args)

        # if only one seqs then return it
        if len(sers)==1:
            res = sers[0]
        else:
            # form product
            nc, sers, order_symbol = cls.flatten(sers)
            res = Mul.__new__(cls, *sers)

        # wrap with coefficient
        if coeff == S.One:
            return res
        else:
            return cls._cls_SeriesCoeffMul()(coeff, res)

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
        new_seq.reverse()
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

    def _sympystr(self, printer, *args):
        if True: #printer._settings["list_series"]:
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
        if True: # printer._settings["list_series"]:
            return SeriesExprPrint._sympystr(self, printer, *args)
        else:
            return printer._print_Mul(self)

class SeriesPow(SeriesExpr, Pow):
    def __new__(cls, b, e, evaluate=True):
        obj = SeriesExpr.__new__(cls, b, e)
        return obj

    @property
    def sequence(self):
        SeqCauchyPow(self.base.sequence, self.exp)

class SeriesNested(SeriesExpr):
    def __new__(cls, *args):
        if args[1].is_SeriesGen:
            return args[0]
        expr = SeriesExpr.__new__(cls, *args)
        return expr

    @property
    def a(self):
        return self.args[0]

    @property
    def b(self):
        return self.args[1]

    @property
    def x(self):
        return self.a.x

    @property
    @cacheit
    def sequence(self):
        # abstract
        return None

    @property
    @cacheit
    def interval(self):
        ai = self.a.interval
        bi = self.b.interval
        start = ai._inf * bi._inf
        stop = ai._sup * bi._sup
        return Interval(start, stop)
