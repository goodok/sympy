from sympy.core.symbol import symbols
from sympy.matrices.matrices import Matrix
from sympy.functions.elementary.exponential import exp

from sympy.solitons.IEFH import I, E, F, H
from sympy.tensor.diff_operator import DiffOperator
from sympy.solitons.hirota import HirotaUnapplyed as D


a = symbols("alpha:4")
z = symbols("z:4")

X = symbols("x:8")
d1 = DiffOperator(X[1])
d2 = DiffOperator(X[2])
d3 = DiffOperator(X[3])
d4 = DiffOperator(X[4])
d5 = DiffOperator(X[5])
d6 = DiffOperator(X[6])
d7 = DiffOperator(X[7])

args = (X[1], X[2], X[3], X[4], X[5], X[6], X[7])
(x1, x2, x3, x4, x5, x6, x7) = args

Billig = []

# Vacuum zero solution
f = 1
g = 0
S0 = Matrix([[f],[g]])


# One soliton
f = 1
g = a[1]*exp(z[1]*x1 + z[1]**3*x3)

S1 = Matrix([[f],[g]])

# Two solitones

A12 = (z[2]-z[1])**2/(z[2]+z[1])**2
f = 1 - a[1]*a[2]*exp(z[2]*x1 + z[2]**3*x3)*exp(z[1]*x1 + z[1]**3*x3)*A12
g = a[1]*exp(z[1]*x1 + z[1]**3*x3) + a[2]*exp(z[2]*x1 + z[2]**3*x3)
f.simplify()

S2 = Matrix([[f],[g]])

Billig = [S0, S1, S2]
