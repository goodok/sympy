""" Caching facility for SymPy """

# TODO: refactor CACHE & friends into class?

# global cache registry:
CACHE = []  # [] of
            #    (item, {} or tuple of {})

from sympy.core.decorators import wraps

def print_cache():
    """print cache content"""

    for item, cache in CACHE:
        item = str(item)
        head = '='*len(item)

        print head
        print item
        print head

        if not isinstance(cache, tuple):
            cache = (cache,)
            shown = False
        else:
            shown = True

        for i, kv in enumerate(cache):
            if shown:
                print '\n*** %i ***\n' % i

            for k, v in kv.iteritems():
                print '  %s :\t%s' % (k, v)

def clear_cache():
    """clear cache content"""
    for item, cache in CACHE:
        if not isinstance(cache, tuple):
            cache = (cache,)

        for kv in cache:
            kv.clear()

########################################

def __cacheit_nocache(func):
    return func

def __cacheit(func):
    """caching decorator.

       important: the result of cached function must be *immutable*


       Examples
       ========

       >>> from sympy.core.cache import cacheit
       >>> @cacheit
       ... def f(a,b):
       ...    return a+b

       >>> @cacheit
       ... def f(a,b):
       ...    return [a,b] # <-- WRONG, returns mutable object

       to force cacheit to check returned results mutability and consistency,
       set environment variable SYMPY_USE_CACHE to 'debug'
    """

    func._cache_it_cache = func_cache_it_cache = {}
    CACHE.append((func, func_cache_it_cache))

    @wraps(func)
    def wrapper(*args, **kw_args):
        """
        Assemble the args and kw_args to compute the hash.
        """
        if kw_args:
            keys = kw_args.keys()
            keys.sort()
            items = [(k+'=', kw_args[k]) for k in keys]
            k = args + tuple(items)
        else:
            k = args
        k = k + tuple(map(lambda x: type(x), k))
        try:
            return func_cache_it_cache[k]
        except KeyError:
            pass
        func_cache_it_cache[k] = r = func(*args, **kw_args)
        return r
    return wrapper

def __cacheit_debug(func):
    """cacheit + code to check cache consistency"""
    cfunc = __cacheit(func)

    @wraps(func)
    def wrapper(*args, **kw_args):
        # always call function itself and compare it with cached version
        r1 = func (*args, **kw_args)
        r2 = cfunc(*args, **kw_args)

        # try to see if the result is immutable
        #
        # this works because:
        #
        # hash([1,2,3])         -> raise TypeError
        # hash({'a':1, 'b':2})  -> raise TypeError
        # hash((1,[2,3]))       -> raise TypeError
        #
        # hash((1,2,3))         -> just computes the hash
        hash(r1), hash(r2)

        # also see if returned values are the same
        assert r1 == r2

        return r1
    return wrapper

def __cacheit_recurr_nocache(func):
    return func

def __cacheit_recurr(fist_cached_index):
    """caching decorator for recurrsions,

    To envelope the recurrsion to the cycle and avoid the stack overflow.

    See analog: sympy.utilities.memorization
    """
    def decorator(func):
        func._cache_it_cache = func_cache_it_cache = {}
        CACHE.append((func, func_cache_it_cache))

        @wraps(func)
        def wrapper(*args, **kw_args):
            """
            Assemble the args and kw_args to compute the hash.
            """
            if kw_args:
                keys = kw_args.keys()
                keys.sort()
                items = [(k+'=', kw_args[k]) for k in keys]
                k = args + tuple(items)
            else:
                k = args
            # k = (self, i)
            n = k[1]
            _self = k[0]

            fci = fist_cached_index
            if fci == None:
                fci = _self.fist_cached_index()
            if n < fci:
                return func(_self, n)

            k = k[:-1] + tuple(map(lambda x: type(x), k)) # k = (self, type(self), type(1))
            if func_cache_it_cache.has_key(k):
                cache_list = func_cache_it_cache[k]
            else:
                func_cache_it_cache[k] = cache_list = []
            L = len(cache_list)
            if n - fci <= L-1:
                return cache_list[n - fci]
            for i in xrange(L + fci, n+1):
                cache_list.append(func(_self, i))
            return cache_list[-1]
        return wrapper
    return decorator

def __cacheit_recurr_debug(func):
    """cacheit + code to check cache consistency"""
    cfunc = __cacheit_recurr(func)

    @wraps(func)
    def wrapper(*args, **kw_args):
        # always call function itself and compare it with cached version
        r1 = func (*args, **kw_args)
        r2 = cfunc(*args, **kw_args)

        # try to see if the result is immutable
        #
        # this works because:
        #
        # hash([1,2,3])         -> raise TypeError
        # hash({'a':1, 'b':2})  -> raise TypeError
        # hash((1,[2,3]))       -> raise TypeError
        #
        # hash((1,2,3))         -> just computes the hash
        hash(r1), hash(r2)

        # also see if returned values are the same
        assert r1 == r2

        return r1
    return wrapper

def _getenv(key, default=None):
    from os import getenv
    return getenv(key, default)

# SYMPY_USE_CACHE=yes/no/debug
USE_CACHE = _getenv('SYMPY_USE_CACHE', 'yes').lower()

if USE_CACHE == 'no':
    cacheit  = __cacheit_nocache
    cacheit_recurr = __cacheit_recurr_nocache
elif USE_CACHE == 'yes':
    cacheit  = __cacheit
    cacheit_recurr = __cacheit_recurr
elif USE_CACHE == 'debug':
    # a lot slower
    cacheit  = __cacheit_debug
    cacheit_recurr = __cacheit_recurr_debug
else:
    raise RuntimeError('unrecognized value for SYMPY_USE_CACHE: %s' % USE_CACHE)
