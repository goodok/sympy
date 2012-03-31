from sympy.core.cache import cacheit, cacheit_recurr

def test_cacheit_doc():
    @cacheit
    def testfn():
        "test docstring"
        pass

    assert testfn.__doc__ == "test docstring"
    assert testfn.__name__ == "testfn"


def test_cacheit_recurr():
    class A(object):
        @cacheit_recurr(None)
        #@cacheit
        def __getitem__(self, i):
            if i<=1:
                return 1
            return self[i-1] + self[i-2]

        def fist_cached_index(self):
            return 5

    class B(object):
        @cacheit_recurr(6)
        #@cacheit
        def __getitem__(self, i):
            if i<=1:
                return 1
            return self[i-1] + self[i-2]

    a = A()
    assert a[3] == 3
    assert a[4] == 5
    assert a[5] == 8
    assert a[6] == 13
    assert a.__getitem__._cache_it_cache.items()[0][1] == [8, 13]

    b = B()
    assert b[3] == 3
    assert b[4] == 5
    assert b[5] == 8
    assert b[6] == 13
    assert b.__getitem__._cache_it_cache.items()[0][1] == [13]
