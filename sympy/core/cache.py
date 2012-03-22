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
            kw_args_sorted = tuple(items)
        else:
            kw_args_sorted = ()

        # try get from cache without unflatten (that is row args)
        k = args + kw_args_sorted
        k = k + tuple(map(lambda x: type(x), k))
        synonym = k
        try:
            return func_cache_it_cache[k]
        except KeyError:
            pass

        # try flatten args, that is caconical order of args
        cls = args[0]
        # we check existens of method of "hashable_content_before_creation"
        # course this method can accept not only SymPy expression, but and
        # e.g. pythonic ints.
        # we check func.__name__ != ... because we allow to use "cacheit" and
        # for this function itself, but must avoid infinity recoursive.
        if hasattr(cls, "_hashable_content_before_creation") and \
                    func.__name__ != "_hashable_content_before_creation":
            seq = args[1:]
            flatten_args = cls._hashable_content_before_creation(seq)
            flatten_args = (cls,) + flatten_args

            k = flatten_args + kw_args_sorted
            k = k + tuple(map(lambda x: type(x), k))
            try:
                r = func_cache_it_cache[k]
                func_cache_it_cache[synonym] = r
                return r
            except KeyError:
                pass

        r = func(*args, **kw_args)
        func_cache_it_cache[k] = r
        func_cache_it_cache[synonym] = r
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

def _getenv(key, default=None):
    from os import getenv
    return getenv(key, default)

# SYMPY_USE_CACHE=yes/no/debug
USE_CACHE = _getenv('SYMPY_USE_CACHE', 'yes').lower()

if USE_CACHE == 'no':
    cacheit  = __cacheit_nocache
elif USE_CACHE == 'yes':
    cacheit  = __cacheit
elif USE_CACHE == 'debug':
    cacheit  = __cacheit_debug   # a lot slower
else:
    raise RuntimeError('unrecognized value for SYMPY_USE_CACHE: %s' % USE_CACHE)
