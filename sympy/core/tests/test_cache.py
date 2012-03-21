from sympy.core.cache import cacheit, clear_cache
from sympy.utilities.pytest import XFAIL
from sympy import Symbol, symbols
from sympy.core.add import Add
from sympy.logic.boolalg import And, Or

def test_cacheit_doc():
    @cacheit
    def testfn():
        "test docstring"
        pass

    assert testfn.__doc__ == "test docstring"
    assert testfn.__name__ == "testfn"

def test_Add():
    clear_cache()
    x, y = symbols("x, y")

    a = Add(x, y)
    b = Add(x, y)
    c = Add(y, x)

    assert a == b
    assert a == c

    assert a is b
    assert a is c

def test_Add_collect():
    clear_cache()
    from sympy import S, collect, expand, factor, Wild
    a, b, c, x = symbols("a, b, c, x")

    r = collect(a*x**2 + b*x**2 + a*x - b*x + c, x)
    assert r == c + x**2*(a + b) + x*(a - b)

def test_assumtion_Add_assumtion():
    clear_cache()
    x, y = symbols("x, y")
    a = Add(x, y)

    x = Symbol("x", real=True)
    b = Add(x, y)

    assert a is not b


# logic.boolalg

def test_And():
    clear_cache()

    x, y = symbols("x, y")
    a = And(x, y)
    b = And(x, y)
    c = And(y, x)

    assert a == b
    assert a == c

    #assert a is b
    #assert a is c

@XFAIL
def test_And_failed_1():
    clear_cache()
    x, y = symbols("x, y")
    assert And(x, y) is And(x, y)

@XFAIL
def test_And_failed_2():
    clear_cache()
    x, y = symbols("x, y")
    assert And(x, y) is And(y, x)
