from sympy import Expr
from sympy.core.basic import Basic
from sympy import S
from sympy.matrices import Matrix

from sympy.core.function import diff

I = Matrix([[1, 0], [0, 1]])
E = Matrix([[0, 1], [1, 0]])
F = Matrix([[0, 1], [-1, 0]])
H = Matrix([[1, 0], [0, -1]])



class IEFH(Expr):
    def __init__(self, i, e, f, h):
        self.i = i
        self.e = e
        self.f = f
        self.h = h
    def __add__(self, other):
        i = self.i + other.i
        e = self.e + other.e
        f = self.f + other.f
        h = self.h + other.h
        return self.__class__(i, e, f, h)
    def __sub__(self, other):
        return self + other*(-1)
    def __mul__(self, other):
        s = self
        o = other
        if isinstance(o, self.__class__):
            i = s.i*o.i + s.e*o.e - s.f*o.f + s.h*o.h
            e = s.i*o.e + s.e*o.i - s.f*o.h + s.h*o.f
            f = s.i*o.f - s.e*o.h + s.f*o.i + s.h*o.e
            h = s.i*o.h - s.e*o.f + s.f*o.e + s.h*o.i
        else:
            i = s.i*o
            e = s.e*o
            f = s.f*o
            h = s.h*o
        return self.__class__(i, e, f, h)
    
    def __div__(self, other):
        s = self
        o = other
        if isinstance(o, self.__class__):
            raise Exception("Not Implemented yet")
        else:
            i = s.i/o
            e = s.e/o
            f = s.f/o
            h = s.h/o
        return self.__class__(i, e, f, h)

    def expand(self, deep=True, **hints):
        sargs, terms = (self.i, self.e, self.f, self.h), []
        for term in sargs:
            newterm = term.expand(**hints)
            terms.append(newterm)
        terms = tuple(terms)
        return self.__class__(*terms)


    def __str__(self):
        return "(%s, %s, %s, %s)" % (self.i, self.e, self.f, self.h)
    
    def _pretty(self,  printer, *args):
        printset = (self.i, self.e, self.f, self.h)
        return printer._print_seq(printset, '(', ')', ', ')
    
    def _latex(self, printer, *args):
        texargs = [r"%s" % printer._print(symbol, *args) for symbol in (self.i, self.e, self.f, self.h)]
        return r"\left(%s\right)" % ", ".join(texargs)






class AKNS_Iterator(object):
    def __init__(self, Q0, Q1, x):
        self.Q = [Q0, Q1]
        self.n = 1
        self.x = x
        
    def next(self):
        n = self.n
        Q = self.Q
        h = self.Bn().i /2
        e = diff(Q[n].f, self.x)/2 + Q[1].e*Q[n].h
        f = diff(Q[n].e, self.x)/2 + Q[1].f*Q[n].h
        Q_new = IEFH(S.Zero, e, f, h)*2
        Q.append(Q_new)
        self.n += 1
        return Q_new.expand()
        
    def Bn(self):
        return self.B(self.n)
    def B(self, n):
        B = IEFH(S.Zero, S.Zero, S.Zero, S.Zero)
        for m in range(1, n+1):
            B += self.Q[m]*self.Q[self.n-m+1]
        return B
