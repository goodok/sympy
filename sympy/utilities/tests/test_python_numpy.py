from sympy import (Rational, Symbol, list2numpy, sin, Float, Matrix, lambdify,
        symarray, symbols, array, mpmath)
from sympy.utilities.pytest import raises
from sympy.utilities.python_numpy import Array as ndarray, empty
import sympy

mpmath.mp.dps = 16
sin02 = mpmath.mpf("0.198669330795061215459412627")

# Test the basics:

def test_array():
    a = array([1, 2, 3], use_sympy_array=True)
    assert a[0] == 1
    assert a[1] == 2
    assert a[2] == 3

def test_empty():
    a = empty(2)
    assert a.shape == (2, )
    a = empty((2,))
    assert a.shape == (2,)
    a = empty((1, 2))
    assert a.shape == (1, 2)
    a = empty((1, 2, 3))
    assert a.shape == (1, 2, 3)
    a = empty((4, 1, 2, 3))
    assert a.shape == (4, 1, 2, 3)

def test_array_2d():
    a = array([[1, 2, 3], [4, 5, 6]], use_sympy_array=True)
    assert a.shape == (2, 3)
    assert a[0, 0] == 1
    assert a[0, 1] == 2
    assert a[0, 2] == 3
    assert a[1, 0] == 4
    assert a[1, 1] == 5
    assert a[1, 2] == 6
    raises(IndexError, "a[0, 3]")
    raises(IndexError, "a[1, 3]")
    raises(IndexError, "a[2, 0]")
    raises(IndexError, "a[-1, 0]")

    a = array([[1, 2], [4, 5]], use_sympy_array=True)
    assert a.shape == (2, 2)
    assert a[0, 0] == 1
    assert a[0, 1] == 2
    assert a[1, 0] == 4
    assert a[1, 1] == 5
    raises(IndexError, "a[1, 2]")
    raises(IndexError, "a[0, 2]")

    a = array([[1], [2]], use_sympy_array=True)
    assert a.shape == (2, 1)
    assert a[0, 0] == 1
    assert a[1, 0] == 2
    raises(IndexError, "a[1, 1]")
    raises(IndexError, "a[0, 1]")

def test_array2():
    a = array([1, 2, 3], use_sympy_array=True)
    b = array([1, 2, 3], use_sympy_array=True)
    c = array([1, 2, 4], use_sympy_array=True)
    assert (a == b).all()
    assert (a == b).any()
    assert not (a != b).all()
    assert not (a != b).any()

    assert (a == c).any()
    assert not (a == c).all()
    assert not (a != c).all()
    assert (a != c).any()

    raises(ValueError, "if a == b: pass")

def test_array3():
    a = array([1, 2, 3], use_sympy_array=True)
    b = array([1, 2], use_sympy_array=True)
    c = array([[1, 2, 3], [4, 5, 6]], use_sympy_array=True)
    d = array([[1, 2, 3], [4, 5, 6]], use_sympy_array=True)
    e = array([[2, 2, 3], [4, 5, 6]], use_sympy_array=True)
    assert a != b
    assert not (a == b)
    assert a != c
    assert not (a == c)
    assert (c==d).all()
    assert (c!=e).any()
    assert (c==e).any()

def test_add1():
    a = array([1, 2, 3], use_sympy_array=True)
    b = array([1+5, 2+5, 3+5], use_sympy_array=True)
    assert (a + 5 == b).all()
    assert (5 + a == b).all()
    assert (a+a == array([2, 4, 6], use_sympy_array=True)).all()

def test_add2():
    a = array([[1, 2, 3], [4, 5, 6]], use_sympy_array=True)
    b = array([[1+5, 2+5, 3+5], [4+5, 5+5, 6+5]], use_sympy_array=True)
    assert (a + 5 == b).all()
    assert (5 + a == b).all()
    assert (a+a == array([[2, 4, 6], [8, 10, 12]], use_sympy_array=True)).all()

def test_tolist():
    a = array([1, 2, 3], use_sympy_array=True)
    assert a.tolist() == [1, 2, 3]
    a = array([[1, 2, 3], [4, 5, 6]], use_sympy_array=True)
    assert a.tolist() == [[1, 2, 3], [4, 5, 6]]

# The rest of this file is copied from test_numpy.py, that tests real numpy.
# Here we use the same tests, but test SymPy's array implementation. This makes
# sure, that things behave in the same way.

# The numpy.matrix test were deleted, as SymPy already has a Matrix class.


# first, systematically check, that all operations are implemented and don't
# raise and exception

def test_systematic_basic():
    def s(sympy_object, numpy_array):
        x = sympy_object + numpy_array
        x = numpy_array + sympy_object
        x = sympy_object - numpy_array
        x = numpy_array - sympy_object
        x = sympy_object * numpy_array
        x = numpy_array * sympy_object
        x = sympy_object / numpy_array
        x = numpy_array / sympy_object
        x = sympy_object ** numpy_array
        x = numpy_array ** sympy_object
    x = Symbol("x")
    y = Symbol("y")
    sympy_objs = [
            Rational(2),
            Float("1.3"),
            x,
            y,
            pow(x,y)*y,
            5,
            5.5,
            ]
    numpy_objs = [
            array([1], use_sympy_array=True),
            array([3, 8, 1], use_sympy_array=True),
            array([x, x**2, Rational(5)], use_sympy_array=True),
            array([x/y*sin(y), 5, Rational(5)], use_sympy_array=True),
            ]
    for x in sympy_objs:
        for y in numpy_objs:
            s(x,y)


# now some random tests, that test particular problems and that also
# check that the results of the operations are correct

def test_basics():
    one = Rational(1)
    zero = Rational(0)
    x = Symbol("x")
    assert array(1, use_sympy_array=True) == array(one, use_sympy_array=True)
    assert array([one], use_sympy_array=True) == array([one], use_sympy_array=True)
    assert array([x], use_sympy_array=True) == array([x], use_sympy_array=True)
    assert array(x, use_sympy_array=True) == array(Symbol("x"), use_sympy_array=True)
    assert array(one+x, use_sympy_array=True) == array(1+x, use_sympy_array=True)

    X = array([one, zero, zero], use_sympy_array=True)
    assert (X == array([one, zero, zero], use_sympy_array=True)).all()
    assert (X == array([one, 0, 0], use_sympy_array=True)).all()

def test_arrays():
    one = Rational(1)
    zero = Rational(0)
    X = array([one, zero, zero], use_sympy_array=True)
    Y = one*X
    X = array([Symbol("a")+Rational(1,2)], use_sympy_array=True)
    Y = X+X
    assert Y == array([1+2*Symbol("a")], use_sympy_array=True)
    assert not Y != array([1+2*Symbol("a")], use_sympy_array=True)
    Y = Y + 1
    assert Y == array([2+2*Symbol("a")], use_sympy_array=True)
    assert not Y != array([2+2*Symbol("a")], use_sympy_array=True)
    Y = X-X
    assert Y == array([0], use_sympy_array=True)
    assert not Y != array([0], use_sympy_array=True)

def test_conversion1():
    x = Symbol("x")
    a = list2numpy([x**2, x], use_sympy_array=True)
    #looks like an array?
    assert isinstance(a, ndarray)
    assert a[0] == x**2
    assert a[1] == x
    assert len(a) == 2
    #yes, it's the array

def test_conversion2():
    x = Symbol("x")
    a = 2*list2numpy([x**2, x], use_sympy_array=True)
    b = list2numpy([2*x**2, 2*x], use_sympy_array=True)
    assert (a == b).all()

    one = Rational(1)
    zero = Rational(0)
    X = list2numpy([one, zero, zero], use_sympy_array=True)
    Y = one*X
    X = list2numpy([Symbol("a")+Rational(1,2)], use_sympy_array=True)
    Y = X+X
    assert Y == array([1+2*Symbol("a")], use_sympy_array=True)
    Y = Y + 1
    assert Y == array([2+2*Symbol("a")], use_sympy_array=True)
    Y = X-X
    assert Y == array([0], use_sympy_array=True)

def test_list2numpy():
    x = Symbol("x")
    assert (array([x**2, x], use_sympy_array=True) == list2numpy([x**2, x],
        use_sympy_array=True)).all()

def test_Matrix1():
    x = Symbol("x")
    m = Matrix([[x, x**2], [5, 2/x]])
    assert (array(m.subs(x, 2), use_sympy_array=True) == array([[2, 4],[5, 1]],
        use_sympy_array=True)).all()
    m = Matrix([[sin(x), x**2], [5, 2/x]])
    assert (array(m.subs(x, 2), use_sympy_array=True) == array([[sin(2), 4],[5,
        1]], use_sympy_array=True)).all()

def test_Matrix3():
    x = Symbol("x")
    a = array([[2, 4],[5, 1]], use_sympy_array=True)
    assert Matrix(a) == Matrix([[2, 4], [5, 1]])
    assert Matrix(a) != Matrix([[2, 4], [5, 2]])
    a = array([[sin(2), 4], [5, 1]], use_sympy_array=True)
    assert Matrix(a) == Matrix([[sin(2), 4],[5, 1]])
    assert Matrix(a) != Matrix([[sin(0), 4],[5, 1]])

def test_Matrix_sum():
    x, y, z = Symbol('x'), Symbol('y'), Symbol('z')
    M = Matrix([[1,2,3],[x,y,x],[2*y,-50,z*x]])
    m = array([[2,3,4],[x,5,6],[x,y,z**2]], use_sympy_array=True)
    assert M+m == Matrix([[3,5,7],[2*x,y+5,x+6],[2*y+x,y-50,z*x+z**2]])
    # This:
    #assert m+M == Matrix([[3,5,7],[2*x,y+5,x+6],[2*y+x,y-50,z*x+z**2]])
    # Was changed to:
    assert (m+M == array([[3,5,7],[2*x,y+5,x+6],[2*y+x,y-50,z*x+z**2]],
        use_sympy_array=True)).all()
    assert M+m == M.add(m)

def _test_Matrix_mul():
    x, y, z = Symbol('x'), Symbol('y'), Symbol('z')
    M = Matrix([[1,2,3],[x,y,x]])
    m = matrix([[2,4],[x,6],[x,z**2]])
    assert M*m == Matrix([
                         [         2 + 5*x,        16 + 3*z**2],
                         [2*x + x*y + x**2, 4*x + 6*y + x*z**2],
                         ])

    assert m*M == Matrix([
                         [   2 + 4*x,      4 + 4*y,      6 + 4*x],
                         [       7*x,    2*x + 6*y,          9*x],
                         [x + x*z**2, 2*x + y*z**2, 3*x + x*z**2],
                         ])
    a = array([2], use_sympy_array=True)
    assert a[0] * M == 2 * M
    assert M * a[0] == 2 * M

def _test_Matrix_array():
    class matarray(object):
        def __array__(self):
            return array([[1,2,3],[4,5,6],[7,8,9]], use_sympy_array=True)
    matarr = matarray()
    assert Matrix(matarr) == Matrix([[1,2,3],[4,5,6],[7,8,9]])

def _test_issue629():
    x = Symbol("x")
    assert (Rational(1,2)*array([2*x, 0], use_sympy_array=True) == array([x,
        0], use_sympy_array=True)).all()
    assert (Rational(1,2)+array([2*x, 0], use_sympy_array=True) ==
            array([2*x+Rational(1,2), Rational(1,2)], use_sympy_array=True)).all()
    assert (Float("0.5")*array([2*x, 0], use_sympy_array=True) ==
            array([Float("1.0")*x, 0], use_sympy_array=True)).all()
    assert (Float("0.5")+array([2*x, 0], use_sympy_array=True) ==
            array([2*x+Float("0.5"), Float("0.5")], use_sympy_array=True)).all()

def _test_lambdify():
    x = Symbol("x")
    f = lambdify(x, sin(x), "numpy")
    prec = 1e-15
    assert -prec < f(0.2) - sin02 < prec
    try:
        f(x) # if this succeeds, it can't be a numpy function
        assert False
    except AttributeError:
        pass

def _test_lambdify_matrix():
    x = Symbol("x")
    f = lambdify(x, Matrix([[x, 2*x],[1, 2]]), "numpy")
    assert (f(1) == matrix([[1,2],[1,2]])).all()

def _test_lambdify_matrix_multi_input():
    x,y,z=sympy.symbols('x,y,z')
    M=sympy.Matrix([[x**2, x*y, x*z],
                    [y*x, y**2, y*z],
                    [z*x, z*y, z**2]])
    f = lambdify((x,y,z), M, "numpy")

    xh,yh,zh = 1.0, 2.0, 3.0
    expected = matrix([[xh**2, xh*yh, xh*zh],
                       [yh*xh, yh**2, yh*zh],
                       [zh*xh, zh*yh, zh**2]])
    actual = f(xh,yh,zh)
    assert numpy.allclose(actual,expected)

def _test_lambdify_matrix_vec_input():
    X=sympy.DeferredVector('X')
    M=Matrix([[X[0]**2, X[0]*X[1], X[0]*X[2]],
              [X[1]*X[0], X[1]**2, X[1]*X[2]],
              [X[2]*X[0], X[2]*X[1], X[2]**2]])
    f = lambdify(X, M, "numpy")

    Xh = array([1.0, 2.0, 3.0], use_sympy_array=True)
    expected = matrix([[Xh[0]**2, Xh[0]*Xh[1], Xh[0]*Xh[2]],
                       [Xh[1]*Xh[0], Xh[1]**2, Xh[1]*Xh[2]],
                       [Xh[2]*Xh[0], Xh[2]*Xh[1], Xh[2]**2]])
    actual = f(Xh)
    assert numpy.allclose(actual,expected)

def _test_lambdify_transl():
    from sympy.utilities.lambdify import NUMPY_TRANSLATIONS
    for sym, mat in NUMPY_TRANSLATIONS.iteritems():
        assert sym in sympy.__dict__
        assert mat in numpy.__dict__

def _test_symarray():
    """Test creation of numpy arrays of sympy symbols."""

    import numpy as np
    import numpy.testing as npt

    syms = symbols('_0 _1 _2')
    s1 = symarray("", 3)
    s2 = symarray("", 3)
    npt.assert_array_equal (s1, np.array(syms, dtype=object))
    assert s1[0] is s2[0]

    a = symarray('a', 3)
    b = symarray('b', 3)
    assert not(a[0] is b[0])

    asyms = symbols('a_0 a_1 a_2')
    npt.assert_array_equal (a, np.array(asyms, dtype=object))

    # Multidimensional checks
    a2d = symarray('a', (2,3))
    assert a2d.shape == (2,3)
    a00, a12 = symbols('a_0_0, a_1_2')
    assert a2d[0,0] is a00
    assert a2d[1,2] is a12

    a3d = symarray('a', (2,3,2))
    assert a3d.shape == (2,3,2)
    a000, a120, a121 = symbols('a_0_0_0, a_1_2_0 a_1_2_1')
    assert a3d[0,0,0] is a000
    assert a3d[1,2,0] is a120
    assert a3d[1,2,1] is a121