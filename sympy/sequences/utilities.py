
from sympy import S, symbols
from sympy.core.cache import cacheit
from sympy.core.function import expand_mul
from sympy.functions.combinatorial.factorials import factorial

Zero = S.Zero
One = S.One

@cacheit
def OrderedPartitions(n, k):
    """
    Generate ordered partitiens.

    This is the helper function.

      k:1      2          3             4
    n 
    1: (1)
    2: (2)  (1, 1)

    3: (3)  (1, 2)    (1, 1, 1)
            (2, 1)

    4: (4)  (1, 3)    (1, 1, 2)    (1, 1, 1, 1)
            (2, 2)    (1, 2, 1)
            (3, 1)    (2, 1, 1)


    (a1, a2, ..., an) --\        (b1, b2, ..., b3)
        |                \           |
        |                 \          |
        |                  \         |
        |                   \        |
    (1 + a1, a2, ..., an)    \  (1+b1, b2, ..., b3)
                              \ (1, a1, a2, ..., an)
    """
    if n<=0:
        return ()
    if k==1: return ((n,), )
    if k>n: return ()
    res = ()
    L = OrderedPartitions(n-1, k-1)
    for l in L:
        # extending
        res += ( ((1,) + l), )
    M = OrderedPartitions(n-1, k)
    for l in M:
        # increment
        res += ( (l[0]+1,) +  l[1:], )
    return res


#############################
###     Splitted Monom    ###
#############################


# instead of keeping trailer_powers we evaluate this as expression
# Advantages:
#   * keep one rational expression instead of tuple
#   * hash only by leading index
#   * it can helps when  b[i] is Zero or One
#   * substitution b[i] to the final result is inline, while after calculation
#     with ordinary Splitted Monom we must substitut result at the end.
# Diasadvantages:
#   * calculation

# version 2

# Without momomole class (just leading index and value): 1.84 sec instead of 2.2

class PlainBellPoly_Splitted(object):
    """
      k:1      2          3             4
    n 
    1: (1)
    2: (2)  (1, 1)
        b2   b1**2

    3: (3)   (1, 2)          (1, 1, 1)        monomial representation
             (2, 1)

        b3   b1*b2           b1**3
             b2*b1

       (3)   (1, b2)        (1, b1**2)        splitted and valued representation
             (2, b1)

    4: (4)   (1, 3)         (1, 1, 2)               (1, 1, 1, 1)
             (2, 2)         (1, 2, 1)
             (3, 1)         (2, 1, 1)

        b4   b1*b3          b1*(2*b1*b2)            b1*(b1**3)
             b2*b2          b2*(b1**2)
             b3*b1

        (4)  (1, b3)        (1, 2*b1*b)            (1, (b1**3))
             (2, b2)        (2, b1**2)
             (3, b1)
    """
    def __init__(self, values, hli={}):
        self.hli = hli
        self.values = values

    def inc(self):
        hli = {}
        for i, value in self.hli.items():
            hli[i+1] = value
        return PlainBellPoly_Splitted(self.values, hli)


    def extend(self, leading_index):
        v = self.value()
        if v is not Zero:
            hli = {leading_index: v}
        return PlainBellPoly_Splitted(self.values, hli)

    def add(self, other):
        hli = self.hli
        for k, value in other.hli.items():
            if hli.has_key(k):
                av = hli[k] + value
                if av is Zero:
                    hli.pop(k)
                else:
                    hli[k] = av
            else:
                hli[k] = value
        return PlainBellPoly_Splitted(self.values, hli)

    def __len__(self):
        return len(self.hli.keys())

    def __str__(self):
        l = []
        for k, value in self.hli.items():
            l.append('%s: %s' % (k, str(value)))
        return "{%s}" % ', '.join(l)

    def value(self):
        v = Zero
        for k, value in self.hli.items():
            av = value*self.values[k]
            if not av.is_Number:
                av = expand_mul(av)   # for QQ is not needed
            v += av
        return v

@cacheit
def PlainBellPoly(n, k, values, delta_index=0):

# TODO: use specialized memorize cache to reduce memory usage.
#       * keep only previouse line of table
#       * chacheit hashs "values", remove it
    if n<=0:
        return ()
    if k==1:
        p = PlainBellPoly_Splitted(values, {n-1 + delta_index: One})
        return p
    if k>n: return ()
    L = PlainBellPoly(n-1, k-1, values)
    L = L.extend(delta_index)
    M = PlainBellPoly(n-1, k, values)
    if M is not ():
        M = M.inc()
        L = L.add(M)
    return L



######################################
#   Bell Polynome of the second kind #
######################################

class BellPoly_Splitted(object):
    def __init__(self, values, hli={}, n = 1):
        self.hli = hli
        self.values = values
        self.n = n

    def inc(self):
        hli = {}
        for i, (v, a, b) in self.hli.items():
            a_new = a + 1
            b_new = b + 1
            v_new = (a_new/b_new)*v
            hli[i+1] = (v_new, a_new, b_new)
        return BellPoly_Splitted(self.values, hli, self.n + 1)


    def extend(self, leading_index):
        # leading_index is always 0
        v_new = self.value()
        a_new = self.n + 1
        b_new = Zero
        monom = (v_new, a_new, b_new)
        if v_new is not Zero:
            hli = {leading_index: monom}
        else:
            hli = {}
        return BellPoly_Splitted(self.values, hli, self.n + 1)

    def add(self, other):
        hli = self.hli
        for k, monom in other.hli.items():
            if hli.has_key(k):
                monom_old = hli[k]
                av = monom_old[0] + monom[0]
                if av is Zero:
                    hli.pop(k)
                else:
                    # each monom(k) has the same "a" and "b"!
                    hli[k] = (av, monom[1], monom[2])
            else:
                hli[k] = monom
        return BellPoly_Splitted(self.values, hli, self.n)

    def value(self):
        value = Zero
        for k, (v, a, b) in self.hli.items():
            av = self.values[k]*v
            if not av.is_Number:
                av = expand_mul(av)   # for QQ is not needed
            value += av
        return value

    def __len__(self):
        return len(self.hli.keys())

    def __str__(self):
        l = []
        for k, monom in self.hli.items():
            l.append('%s: %s' % (k, str(monom)))
        return "n=%s, {%s}" % (self.n, ', '.join(l))



@cacheit
def BellPoly(n, k, values, delta_index=0):
    if n<=0:
        return ()
    if k==1:
        p = BellPoly_Splitted(values, {n-1 + delta_index: (One, n-1, n-1)}, n-1)
        return p
    if k>n: return ()
    L = BellPoly(n-1, k-1, values)
    L = L.extend(delta_index)
    M = BellPoly(n-1, k, values)
    if M is not ():
        M = M.inc()
        L = L.add(M)
    return L

