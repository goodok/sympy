from sympy.simplify.hyperexpand import (ShiftA, ShiftB, UnShiftA, UnShiftB,
                       MeijerShiftA, MeijerShiftB, MeijerShiftC, MeijerShiftD,
                       MeijerUnShiftA, MeijerUnShiftB, MeijerUnShiftC,
                       MeijerUnShiftD,
                       ReduceOrder, reduce_order, apply_operators,
                       devise_plan, make_derivative_operator, Formula,
                       hyperexpand, IndexPair, IndexQuadruple,
                       reduce_order_meijer,
                       build_hypergeometric_formula)
from sympy import hyper, I, S, meijerg, Piecewise, exp_polar
from sympy.utilities.pytest import raises
from sympy.abc import z, a, b, c
from sympy.utilities.randtest import test_numerically as tn
from sympy.utilities.pytest import XFAIL, skip, slow
from random import randrange

from sympy import (cos, sin, log, exp, asin, lowergamma, atanh, besseli,
                   gamma, sqrt, pi, erf, exp_polar)

def test_branch_bug():
    assert hyperexpand(hyper((-S(1)/3, S(1)/2), (S(2)/3, S(3)/2), -z)) == \
           -z**S('1/3')*lowergamma(exp_polar(I*pi)/3, z)/5 \
           + sqrt(pi)*erf(sqrt(z))/(5*sqrt(z))
    assert hyperexpand(meijerg([S(7)/6, 1], [], [S(2)/3], [S(1)/6, 0], z)) == \
           2*z**S('2/3')*(2*sqrt(pi)*erf(sqrt(z))/sqrt(z) - 2*lowergamma(S(2)/3, z)/z**S('2/3'))*gamma(S(2)/3)/gamma(S(5)/3)

def test_hyperexpand():
    # Luke, Y. L. (1969), The Special Functions and Their Approximations,
    # Volume 1, section 6.2

    assert hyperexpand(hyper([], [], z)) == exp(z)
    assert hyperexpand(hyper([1, 1], [2], -z)*z) == log(1 + z)
    assert hyperexpand(hyper([], [S.Half], -z**2/4)) == cos(z)
    assert hyperexpand(z*hyper([], [S('3/2')], -z**2/4)) == sin(z)
    assert hyperexpand(hyper([S('1/2'), S('1/2')], [S('3/2')], z**2)*z) \
           == asin(z)

def can_do(ap, bq, numerical=True, div=1, lowerplane=False):
    from sympy import exp_polar, exp
    r = hyperexpand(hyper(ap, bq, z))
    if r.has(hyper):
        return False

    if not numerical:
        return True

    repl = {}
    for n, a in enumerate(r.free_symbols - set([z])):
        repl[a] = randcplx(n)/div
    [a, b, c, d] = [2, -1, 3, 1]
    if lowerplane:
        [a, b, c, d] = [2, -2, 3, -1]
    return tn(hyper(ap, bq, z).subs(repl), r.replace(exp_polar, exp).subs(repl), z,
              a=a, b=b, c=c, d=d)

def test_roach():
    # Kelly B. Roach.  Meijer G Function Representations.
    # Section "Gallery"
    assert can_do([S(1)/2], [S(9)/2])
    assert can_do([], [1, S(5)/2, 4])
    assert can_do([-S.Half, 1, 2], [3, 4])
    assert can_do([S(1)/3], [-S(2)/3, -S(1)/2, S(1)/2, 1])
    assert can_do([-S(3)/2, -S(1)/2], [-S(5)/2, 1])
    assert can_do([-S(3)/2,], [-S(1)/2, S(1)/2]) # shine-integral

@XFAIL
def test_roach_fail():
    assert can_do([-S(3)/2, -S(1)/2], [2]) # elliptic integrals
    assert can_do([-S(1)/2, 1], [S(1)/4, S(1)/2, S(3)/4]) # PFDD
    assert can_do([S(3)/2], [S(5)/2, 5]) # struve function
    assert can_do([-S(1)/2, S(1)/2, 1], [S(3)/2, S(5)/2]) # polylog, pfdd
    assert can_do([1, 2, 3], [S(1)/2, 4]) # XXX ?
    assert can_do([S(1)/2], [-S(1)/3, -S(1)/2, -S(2)/3]) # PFDD ?

# For the long table tests, see end of file

def test_polynomial():
    from sympy import oo
    assert hyperexpand(hyper([], [-1], z)) == oo
    assert hyperexpand(hyper([-2], [-1], z)) == oo
    assert hyperexpand(hyper([0, 0], [-1], z)) == 1
    assert can_do([-5, -2, randcplx(), randcplx()], [-10, randcplx()])

def test_hyperexpand_bases():
    assert hyperexpand(hyper([2], [a], z)) == \
  a + z**(-a + 1)*(-a**2 + 3*a + z*(a - 1) - 2)*exp(z)*lowergamma(a - 1, z) - 1
    # TODO [a+1, a-S.Half], [2*a]
    assert hyperexpand(hyper([1, 2], [3], z)) == -2/z - 2*log(-z + 1)/z**2
    assert hyperexpand(hyper([S.Half, 2], [S(3)/2], z)) == \
      -1/(2*z - 2) + atanh(sqrt(z))/sqrt(z)/2
    assert hyperexpand(hyper([S(1)/2, S(1)/2], [S(5)/2], z)) == \
               (-3*z + 3)/4/(z*sqrt(-z + 1)) \
               + (6*z - 3)*asin(sqrt(z))/(4*z**(S(3)/2))
    assert hyperexpand(hyper([1, 2], [S(3)/2], z)) == -1/(2*z - 2) \
            - asin(sqrt(z))/(sqrt(z)*(2*z - 2)*sqrt(-z + 1))
    assert hyperexpand(hyper([-S.Half - 1, 1, 2], [S.Half, 3], z)) == \
             sqrt(z)*(6*z/7 - S(6)/5)*atanh(sqrt(z)) \
           + (-30*z**2 + 32*z - 6)/35/z - 6*log(-z + 1)/(35*z**2)
    assert hyperexpand(hyper([1+S.Half, 1, 1], [2, 2], z)) == \
           -4*log(sqrt(-z + 1)/2 + S(1)/2)/z
    # TODO hyperexpand(hyper([a], [2*a + 1], z))
    # TODO [S.Half, a], [S(3)/2, a+1]
    assert hyperexpand(hyper([2], [b, 1], z)) == \
             z**(-b/2 + S(1)/2)*besseli(b - 1, 2*sqrt(z))*gamma(b) \
           + z**(-b/2 + 1)*besseli(b, 2*sqrt(z))*gamma(b)
    # TODO [a], [a - S.Half, 2*a]

def test_hyperexpand_parametric():
    assert hyperexpand(hyper([a, S(1)/2 + a], [S(1)/2], z)) \
        == (1 + sqrt(z))**(-2*a)/2 + (1 - sqrt(z))**(-2*a)/2
    assert hyperexpand(hyper([a, -S(1)/2 + a], [2*a], z)) \
        == 2**(2*a - 1)*((-z + 1)**(S(1)/2) + 1)**(-2*a + 1)

def test_shifted_sum():
    from sympy import simplify
    assert simplify(hyperexpand(z**4*hyper([2], [3, S('3/2')], -z**2))) \
           == -S(1)/2 + cos(2*z)/2 + z*sin(2*z) - z**2*cos(2*z)

def randrat():
    """ Steer clear of integers. """
    return S(randrange(25) + 10)/50

def randcplx(offset=-1):
    """ Polys is not good with real coefficients. """
    return randrat() + I*randrat() + I*(1 + offset)

def test_formulae():
    from sympy.simplify.hyperexpand import FormulaCollection
    formulae = FormulaCollection().formulae
    for formula in formulae:
        h = hyper(formula.indices.ap, formula.indices.bq, formula.z)
        rep = {}
        for n, sym in enumerate(formula.symbols):
            rep[sym] = randcplx(n)

        # NOTE hyperexpand returns truly branched functions. We know we are
        #      on the main sheet, but numerical evaluation can still go wrong
        #      (e.g. if exp_polar cannot be evalf'd).
        #      Just replace all exp_polar by exp, this usually works.

        # first test if the closed-form is actually correct
        h = h.subs(rep)
        closed_form = formula.closed_form.subs(rep).rewrite('nonrepsmall')
        z = formula.z
        assert tn(h, closed_form.replace(exp_polar, exp), z)

        # now test the computed matrix
        cl = (formula.C * formula.B)[0].subs(rep).rewrite('nonrepsmall')
        assert tn(closed_form.replace(exp_polar, exp), cl.replace(exp_polar, exp), z)
        deriv1 = z*formula.B.applyfunc(lambda t: t.rewrite('nonrepsmall')).diff(z)
        deriv2 = formula.M * formula.B
        for d1, d2 in zip(deriv1, deriv2):
            assert tn(d1.subs(rep).replace(exp_polar, exp),
                      d2.subs(rep).rewrite('nonrepsmall').replace(exp_polar, exp), z)

def test_meijerg_formulae():
    from sympy.simplify.hyperexpand import MeijerFormulaCollection
    formulae = MeijerFormulaCollection().formulae
    for sig in formulae:
        for formula in formulae[sig]:
          g = meijerg(formula.indices.an, formula.indices.ap,
                      formula.indices.bm, formula.indices.bq,
                      formula.z)
          rep = {}
          for sym in formula.symbols:
              rep[sym] = randcplx()

          # first test if the closed-form is actually correct
          g = g.subs(rep)
          closed_form = formula.closed_form.subs(rep)
          z = formula.z
          assert tn(g, closed_form, z)
          #print closed_form

          # now test the computed matrix
          cl = (formula.C * formula.B)[0].subs(rep)
          assert tn(closed_form, cl, z)
          deriv1 = z*formula.B.diff(z)
          deriv2 = formula.M * formula.B
          for d1, d2 in zip(deriv1, deriv2):
              assert tn(d1.subs(rep), d2.subs(rep), z)

def op(f): return z*f.diff(z)

def test_plan():
    assert devise_plan(IndexPair([0], ()), IndexPair([0], ()), z) == []
    raises(ValueError, 'devise_plan(IndexPair([1], ()), IndexPair((), ()), z)')
    raises(ValueError, 'devise_plan(IndexPair([2], [1]), IndexPair([2], [2]), z)')
    raises(KeyError,
           'devise_plan(IndexPair([2], []), IndexPair([S("1/2")], []), z)')

    # We cannot use pi/(10000 + n) because polys is insanely slow.
    a1, a2, b1 = map(lambda n: randcplx(n), range(3))
    b1 += 2*I
    h = hyper([a1, a2], [b1], z)

    h2 = hyper((a1 + 1, a2), [b1], z)
    assert tn(apply_operators(h, devise_plan(IndexPair((a1 + 1, a2), [b1]),
                                      IndexPair((a1, a2), [b1]), z), op),
       h2, z)

    h2 = hyper((a1 + 1, a2 - 1), [b1], z)
    assert tn(apply_operators(h, devise_plan(IndexPair((a1 + 1, a2 - 1), [b1]),
                                      IndexPair((a1, a2), [b1]), z), op),
       h2, z)

def test_plan_derivatives():
    import random
    random.setstate((3, (2316134452, 2939535594, 1180604944, 805005976, 687247053, 1551873565, 3596896848, 349442968, 7165242, 2847966214, 279753078, 1067881182, 81681645, 2057975529, 4107681706, 2805876213, 895569731, 4293487064, 1451048427, 980163346, 2409305917, 773389990, 2679063466, 1236957122, 485652260, 1902098688, 2277499113, 794111716, 3227749256, 2622871864, 1672450452, 1588656375, 1728624350, 1505957928, 1070464002, 1862345757, 996167638, 87258197, 2163175040, 439719245, 3446407317, 4183808605, 2009673949, 2887501351, 1853753786, 2615619868, 887075592, 229355813, 1851672136, 3028339444, 3200999588, 1457012466, 1729326915, 3599423730, 3056214251, 724148542, 2466042451, 2567015796, 2611023122, 780673977, 2979556354, 4050415155, 406688071, 951497796, 2141558973, 428053397, 1724777538, 3993291586, 1551854956, 1577299338, 347720801, 740164880, 1540793363, 518992846, 90907102, 1076230674, 3713067252, 2012195112, 529181255, 1871210103, 2417135810, 530274576, 662370732, 531129842, 2272101152, 2691862973, 1828821237, 274937281, 3355554499, 2091217095, 3829529556, 1196247116, 56330095, 2362918206, 348610938, 1913734075, 2438364900, 1312450427, 3676645734, 4233748562, 737956889, 2797193457, 4132564605, 3026500601, 3032503617, 3281909376, 2132227573, 2867027873, 2991836870, 893056435, 1968251915, 1307205091, 552579683, 1250524595, 52154542, 3369811210, 3065294489, 256280912, 2634403447, 1586774283, 3368223471, 3372752397, 2327884349, 2462985767, 900243879, 831604262, 2023390132, 2828078423, 1930978898, 2734783207, 803493449, 2256223058, 2018163374, 495539891, 2752447597, 2867698104, 2736162515, 1691071182, 1938994128, 2077096464, 222818990, 1495156010, 3063602609, 3278268988, 1243055691, 2097748691, 947057608, 238533549, 3487049459, 543871179, 2714885714, 913862633, 58133197, 3202647005, 2331197475, 2499691729, 3851396247, 2487849682, 1119661484, 3744471528, 1553121367, 2007201668, 3307407771, 3172832305, 1296240854, 2929620022, 1268098341, 3422025327, 185809165, 3882199091, 811890386, 3124596579, 3937063770, 3666428829, 156722548, 1625189854, 3551597506, 40988505, 448578359, 1629172646, 1419939713, 1947666069, 2662684530, 1542187404, 1957938959, 743996794, 3873091265, 336312820, 284703669, 3668190370, 3558443320, 185436310, 2144542652, 3083294211, 1461040621, 1018415209, 3309803046, 1703567371, 3685561911, 2631099231, 4222659517, 606555516, 1067344535, 781835130, 1822614546, 3865820860, 1617230309, 1141104383, 3461880669, 2212351565, 2608582405, 4284435701, 434910988, 1653091042, 3180601581, 4185691584, 3663037907, 2288293804, 3778516578, 2137138832, 2644627823, 234364794, 3662586565, 3100738146, 683999299, 3945996813, 630864, 971223469, 1709295542, 2132910301, 1580622664, 132404934, 3079075995, 4224907559, 92257524, 1608692248, 145961315, 2314423878, 2037556850, 702024500, 4120617109, 717603922, 923013309, 4060053535, 3285548625, 676632915, 3209979902, 3582353832, 3554477510, 1164846175, 3444588245, 2372927676, 4171537015, 3578808515, 1528153351, 3533734580, 643044367, 4046986669, 1369431315, 1492411458, 2254878000, 3740302077, 4190206929, 1116386841, 1648147371, 2928680994, 1475984097, 2256391526, 1519952303, 2572362454, 4005047751, 3138725238, 205787315, 2953574678, 3065290588, 437404544, 3413149831, 611245023, 2498823126, 3195989240, 341370669, 3554157976, 2735719448, 3130220910, 1600145902, 1426718745, 869415012, 773069322, 2170142077, 2025051272, 124084053, 3150237911, 3973987382, 3535112337, 4255024619, 3642677101, 3247314272, 2374384516, 410466860, 4249681313, 702157370, 3415563290, 3141670979, 1194119419, 4147806055, 173234826, 3092030617, 3379723909, 1629207910, 1548799047, 4121007354, 3567234236, 2432279533, 57999133, 2374446893, 944370999, 2240495416, 899839007, 578071280, 834059854, 1775692916, 1887587877, 1595980965, 310669373, 4293724667, 2241783245, 97781244, 1856138337, 1128858864, 3582429303, 326642084, 2370621057, 1506611015, 184312478, 4031483705, 3490684065, 3522739925, 2486836714, 2034368036, 900585578, 3480555139, 3306256779, 2758343659, 4027857202, 3884369916, 312664992, 151265907, 1806565054, 4066082281, 968625487, 2844758154, 1678065553, 1568071994, 1170020712, 3130389177, 2751082176, 2138567465, 1932404322, 1415599620, 3809827932, 2034165800, 3082205537, 740006189, 1741898081, 3259711565, 2353983939, 496666559, 1667393751, 1316281406, 469792389, 1304962287, 400438561, 1836367015, 1143926770, 3944238743, 3096828587, 827697911, 2460071428, 3359226460, 1225571537, 3182418372, 1875183977, 2457074373, 374787659, 3093873746, 3689306864, 1689610990, 1760693481, 1852743591, 3185082822, 173217405, 3961773525, 344540660, 310762262, 863703567, 556907037, 1216965241, 170228730, 338297387, 3701948240, 719918379, 1067024735, 3162521554, 3074430859, 2237582600, 449626666, 195169909, 1658584982, 2254967368, 2265531073, 3119627674, 3377262836, 2786987950, 2636869204, 245162880, 3190597842, 3669240139, 1817639969, 3429889808, 1417481196, 2145830939, 2954866785, 1856107107, 1042058327, 1232926931, 1841014554, 1933490219, 1788891993, 3300860199, 2757295918, 1956986147, 616367099, 452872962, 4155648820, 732194387, 704531372, 3820152218, 2411297884, 759112178, 2831496738, 2060454645, 1253116388, 2327740455, 1437865471, 1610926959, 3213885230, 2545745604, 2476951086, 769839511, 3613057110, 3403894216, 3844820275, 2706854805, 1448972211, 4132866683, 151884579, 376187291, 3055895279, 1496824206, 2930641185, 3224514851, 1481839062, 508172661, 1622176994, 1448043711, 1884960342, 125551660, 2432510766, 2679416537, 2708202264, 1772936405, 2006211631, 761467388, 3789828996, 4076523416, 776479291, 548567206, 1930417463, 3237938438, 1253629661, 2164518038, 1614985910, 3387348572, 1364715535, 3673367649, 3877456507, 1698313752, 1078076859, 207888577, 1504728569, 58702769, 2254214484, 3176625406, 4134186199, 2008262668, 3588692839, 3543829738, 3598014515, 2749030941, 54110657, 3177007956, 1192423610, 2825391280, 1735375896, 1679028242, 1432814255, 1581558122, 2069462238, 1351006453, 4280652467, 1209759921, 4222492642, 3628550236, 3955057428, 413525627, 3883981353, 2894088655, 838319235, 1903274254, 1098275482, 320558285, 1304611655, 595004992, 4150138991, 3457250254, 3113355644, 2246976179, 3168136198, 307586044, 662915776, 1431857819, 1518683146, 3942753446, 3633034834, 3410518253, 1492179093, 4108722463, 1338766228, 1528716302, 545624105, 3336676993, 2903310170, 420163928, 3082942681, 725461532, 3911656311, 14448456, 1142136591, 2775738263, 983709490, 840541491, 1387171814, 4142296030, 3132638378, 3942466163, 2806681501, 1428129098, 1576009588, 3201671054, 2492910204, 3566692026, 1526388087, 1485008314, 1849441138, 1680485417, 737566886, 2337044107, 4194901934, 1919014967, 4181849878, 383886011, 2840624566, 2070434158, 1733311891, 873781217, 2852372091, 2918464282, 2083925300, 3418657739, 3005873265, 935869311, 2182104567, 565212310, 3362034969, 1502963523, 2295424916, 1275550970, 4233566270, 4197862855, 661569690, 954705687, 2963716479, 1740655690, 1308468633, 409163461, 2210242570, 2317605731, 2144326061, 2694620653, 787921023, 2961592158, 3736326421, 273759413, 1270445013, 3828208415, 1279974758, 731066837, 3338906460, 487482633, 2460889926, 1324790961, 4153290165, 1278268978, 3651853866, 2917371678, 205231226, 1968794984, 3354262318, 589427050, 19966254, 3833703531, 1078974679, 740399374, 227532387, 490), None))
    a1, a2, a3 = 1, 2, S('1/2')
    b1, b2 = 3, S('5/2')
    h = hyper((a1, a2, a3), (b1, b2), z)
    h2 = hyper((a1 + 1, a2 + 1, a3 + 2), (b1 + 1, b2 + 1), z)
    ops = devise_plan(IndexPair((a1 + 1, a2 + 1, a3 + 2), (b1 + 1, b2 + 1)),
                      IndexPair((a1, a2, a3), (b1, b2)), z)
    f = Formula((a1, a2, a3), (b1, b2), z, h, [])
    deriv = make_derivative_operator(f.M, z)
    assert tn((apply_operators(f.C, ops, deriv)*f.B)[0], h2, z)

    h2 = hyper((a1, a2 - 1, a3 - 2), (b1 - 1, b2 - 1), z)
    ops = devise_plan(IndexPair((a1, a2 - 1, a3 - 2), (b1 - 1, b2 - 1)),
                      IndexPair((a1, a2, a3), (b1, b2)), z)
    assert tn((apply_operators(f.C, ops, deriv)*f.B)[0], h2, z)

def test_reduction_operators():
    a1, a2, b1 = map(lambda n: randcplx(n), range(3))
    h = hyper([a1], [b1], z)

    assert ReduceOrder(2, 0) is None
    assert ReduceOrder(2, -1) is None
    assert ReduceOrder(1, S('1/2')) is None

    h2 = hyper((a1, a2), (b1, a2), z)
    assert tn(ReduceOrder(a2, a2).apply(h, op), h2, z)

    h2 = hyper((a1, a2 + 1), (b1, a2), z)
    assert tn(ReduceOrder(a2 + 1, a2).apply(h, op), h2, z)

    h2 = hyper((a2 + 4, a1), (b1, a2), z)
    assert tn(ReduceOrder(a2 + 4, a2).apply(h, op), h2, z)

    # test several step order reduction
    ap = (a2 + 4, a1, b1 + 1)
    bq = (a2, b1, b1)
    nip, ops = reduce_order(IndexPair(ap, bq))
    assert nip.ap == (a1,)
    assert nip.bq == (b1,)
    assert tn(apply_operators(h, ops, op), hyper(ap, bq, z), z)

def test_shift_operators():
    a1, a2, b1, b2, b3 = map(lambda n: randcplx(n), range(5))
    h = hyper((a1, a2), (b1, b2, b3), z)

    raises(ValueError, 'ShiftA(0)')
    raises(ValueError, 'ShiftB(1)')

    assert tn(ShiftA(a1).apply(h, op), hyper((a1 + 1, a2), (b1, b2, b3), z), z)
    assert tn(ShiftA(a2).apply(h, op), hyper((a1, a2 + 1), (b1, b2, b3), z), z)
    assert tn(ShiftB(b1).apply(h, op), hyper((a1, a2), (b1 - 1, b2, b3), z), z)
    assert tn(ShiftB(b2).apply(h, op), hyper((a1, a2), (b1, b2 - 1, b3), z), z)
    assert tn(ShiftB(b3).apply(h, op), hyper((a1, a2), (b1, b2, b3 - 1), z), z)

def test_ushift_operators():
    a1, a2, b1, b2, b3 = map(lambda n: randcplx(n), range(5))
    h = hyper((a1, a2), (b1, b2, b3), z)

    raises(ValueError, 'UnShiftA((1,), (), 0, z)')
    raises(ValueError, 'UnShiftB((), (-1,), 0, z)')
    raises(ValueError, 'UnShiftA((1,), (0, -1, 1), 0, z)')
    raises(ValueError, 'UnShiftB((0, 1), (1,), 0, z)')

    s = UnShiftA((a1, a2), (b1, b2, b3), 0, z)
    assert tn(s.apply(h, op), hyper((a1 - 1, a2), (b1, b2, b3), z), z)
    s = UnShiftA((a1, a2), (b1, b2, b3), 1, z)
    assert tn(s.apply(h, op), hyper((a1, a2 - 1), (b1, b2, b3), z), z)

    s = UnShiftB((a1, a2), (b1, b2, b3), 0, z)
    assert tn(s.apply(h, op), hyper((a1, a2), (b1 + 1, b2, b3), z), z)
    s = UnShiftB((a1, a2), (b1, b2, b3), 1, z)
    assert tn(s.apply(h, op), hyper((a1, a2), (b1, b2 + 1, b3), z), z)
    s = UnShiftB((a1, a2), (b1, b2, b3), 2, z)
    assert tn(s.apply(h, op), hyper((a1, a2), (b1, b2, b3 + 1), z), z)


def can_do_meijer(a1, a2, b1, b2, numeric=True):
    """
    This helper function tries to hyperexpand() the meijer g-function
    corresponding to the parameters a1, a2, b1, b2.
    It returns False if this expansion still contains g-functions.
    If numeric is True, it also tests the so-obtained formula numerically
    (at random values) and returns False if the test fails.
    Else it returns True.
    """
    from sympy import unpolarify, expand
    r = hyperexpand(meijerg(a1, a2, b1, b2, z))
    if r.has(meijerg):
        return False
    # NOTE hyperexpand() returns a truly branched function, whereas numerical
    #      evaluation only works on the main branch. Since we are evaluating on
    #      the main branch, this should not be a problem, but expressions like
    #      exp_polar(I*pi/2*x)**a are evaluated incorrectly. We thus have to get
    #      rid of them. The expand heuristically does this...
    r = unpolarify(expand(r, force=True, power_base=True, power_exp=False,
                          mul=False, log=False, multinomial=False, basic=False))

    if not numeric:
        return True

    repl = {}
    for n, a in enumerate(meijerg(a1, a2, b1, b2, z).free_symbols - set([z])):
        repl[a] = randcplx(n)
    return tn(meijerg(a1, a2, b1, b2, z).subs(repl), r.subs(repl), z)

def test_meijerg_expand():
    from sympy import combsimp, simplify
    # from mpmath docs
    assert hyperexpand(meijerg([[],[]], [[0],[]], -z)) == exp(z)

    assert hyperexpand(meijerg([[1,1],[]], [[1],[0]], z)) == \
        log(z + 1)
    assert hyperexpand(meijerg([[1,1],[]], [[1],[1]], z)) == \
        z/(z + 1)
    assert hyperexpand(meijerg([[],[]], [[S(1)/2],[0]], (z/2)**2)) \
           == sin(z)/sqrt(pi)
    assert hyperexpand(meijerg([[],[]], [[0], [S(1)/2]], (z/2)**2)) \
           == cos(z)/sqrt(pi)
    assert can_do_meijer([], [a], [a-1, a-S.Half], [])
    assert can_do_meijer([], [], [a/2], [-a/2], False) # branches...
    assert can_do_meijer([a], [b], [a], [b, a - 1])

    # wikipedia
    assert hyperexpand(meijerg([1], [], [], [0], z)) == \
       Piecewise((0, abs(z) < 1), (1, abs(1/z) < 1),
                 (meijerg([1], [], [], [0], z), True))
    assert hyperexpand(meijerg([], [1], [0], [], z)) == \
       Piecewise((1, abs(z) < 1), (0, abs(1/z) < 1),
                 (meijerg([], [1], [0], [], z), True))

    # The Special Functions and their Approximations
    assert can_do_meijer([], [], [a + b/2], [a, a - b/2, a + S.Half])
    assert can_do_meijer([], [], [a], [b], False) # branches only agree for small z
    assert can_do_meijer([], [S.Half], [a], [-a])
    assert can_do_meijer([], [], [a, b], [])
    assert can_do_meijer([], [], [a, b], [])
    assert can_do_meijer([], [], [a, a+S.Half], [b, b+S.Half])
    assert can_do_meijer([], [], [a, -a], [0, S.Half], False) # dito
    assert can_do_meijer([], [], [a, a+S.Half, b, b+S.Half], [])
    assert can_do_meijer([S.Half], [], [0], [a, -a])
    assert can_do_meijer([S.Half], [], [a], [0, -a], False) # dito
    assert can_do_meijer([], [a - S.Half], [a, b], [a - S.Half], False)
    assert can_do_meijer([], [a+S.Half], [a+b, a-b, a], [], False)
    assert can_do_meijer([a+S.Half], [], [b, 2*a-b, a], [], False)

    # This for example is actually zero.
    assert can_do_meijer([], [], [], [a, b])

    # Testing a bug:
    assert hyperexpand(meijerg([0, 2], [], [], [-1, 1], z)) == \
        Piecewise((0, abs(z) < 1),
                  (z*(1 - 1/z**2)/2, abs(1/z) < 1),
                  (meijerg([0, 2], [], [], [-1, 1], z), True))

    # Test that the simplest possible answer is returned:
    assert combsimp(simplify(hyperexpand(meijerg([1], [1-a], [-a/2, -a/2 + S(1)/2],
                                                 [], 1/z)))) == \
           -2*sqrt(pi)*(sqrt(z + 1) + 1)**a/a

    # Test that hyper is returned
    assert hyperexpand(meijerg([1], [], [a], [0, 0], z)) == \
           z**a*gamma(a)*hyper((a,), (a + 1, a + 1), z*exp_polar(I*pi))/gamma(a + 1)**2

def test_meijerg_lookup():
    from sympy import uppergamma
    assert hyperexpand(meijerg([a], [], [b, a], [], z)) == \
           z**b*exp(z)*gamma(-a + b + 1)*uppergamma(a - b, z)
    assert hyperexpand(meijerg([0], [], [0, 0], [], z)) == \
           exp(z)*uppergamma(0, z)
    assert can_do_meijer([a], [], [b, a+1], [])
    assert can_do_meijer([a], [], [b+2, a], [])
    assert can_do_meijer([a], [], [b-2, a], [])

@XFAIL
def test_meijerg_expand_fail():
    # These basically test hyper([], [1/2 - a, 1/2 + 1, 1/2], z),
    # which is *very* messy. But since the meijer g actually yields a
    # sum of bessel functions, things can sometimes be simplified a lot and
    # are then put into tables...
    assert can_do_meijer([], [], [a + S.Half], [a, a - b/2, a + b/2])
    assert can_do_meijer([], [], [0, S.Half], [a, -a])
    assert can_do_meijer([], [], [3*a - S.Half, a, -a - S.Half], [a - S.Half])
    assert can_do_meijer([], [], [0, a - S.Half, -a - S.Half], [S.Half])
    assert can_do_meijer([], [], [a, b + S(1)/2, b], [2*b - a])
    assert can_do_meijer([], [], [a, b + S(1)/2, b, 2*b - a])
    assert can_do_meijer([S.Half], [], [-a, a], [0])

def test_meijerg():
    # carefully set up the parameters.
    # NOTE: this used to fail sometimes. I believe it is fixed, but if you
    #       hit an inexplicable test failure here, please let me know the seed.
    a1, a2 = map(lambda n: randcplx() - 5*I - n*I, range(2))
    b1, b2 = map(lambda n: randcplx() + 5*I + n*I, range(2))
    b3, b4, b5, a3, a4, a5 = map(lambda n: randcplx(), range(6))
    g = meijerg([a1], [a3, a4], [b1], [b3, b4], z)

    assert ReduceOrder.meijer_minus(3, 4) is None
    assert ReduceOrder.meijer_plus(4, 3) is None

    g2 = meijerg([a1, a2], [a3, a4], [b1], [b3, b4, a2], z)
    assert tn(ReduceOrder.meijer_plus(a2, a2).apply(g, op), g2, z)

    g2 = meijerg([a1, a2], [a3, a4], [b1], [b3, b4, a2 + 1], z)
    assert tn(ReduceOrder.meijer_plus(a2, a2 + 1).apply(g, op), g2, z)

    g2 = meijerg([a1, a2 - 1], [a3, a4], [b1], [b3, b4, a2 + 2], z)
    assert tn(ReduceOrder.meijer_plus(a2 - 1, a2 + 2).apply(g, op), g2, z)

    g2 = meijerg([a1], [a3, a4, b2 - 1], [b1, b2 + 2], [b3, b4], z)
    assert tn(ReduceOrder.meijer_minus(b2 + 2, b2 - 1).apply(g, op), g2, z, tol=1e-6)

    # test several-step reduction
    an = [a1, a2]
    bq = [b3, b4, a2 + 1]
    ap = [a3, a4, b2 - 1]
    bm = [b1, b2 + 1]
    niq, ops = reduce_order_meijer(IndexQuadruple(an, ap, bm, bq))
    assert niq.an == (a1,)
    assert set(niq.ap) == set([a3, a4])
    assert niq.bm == (b1,)
    assert set(niq.bq) == set([b3, b4])
    assert tn(apply_operators(g, ops, op), meijerg(an, ap, bm, bq, z), z)

def test_meijerg_shift_operators():
    # carefully set up the parameters. XXX this still fails sometimes
    a1, a2, a3, a4, a5, b1, b2, b3, b4, b5 = \
        map(lambda n: randcplx(n), range(10))
    g = meijerg([a1], [a3, a4], [b1], [b3, b4], z)

    assert tn(MeijerShiftA(b1).apply(g, op),
              meijerg([a1], [a3, a4], [b1 + 1], [b3, b4], z), z)
    assert tn(MeijerShiftB(a1).apply(g, op),
              meijerg([a1 - 1], [a3, a4], [b1], [b3, b4], z), z)
    assert tn(MeijerShiftC(b3).apply(g, op),
              meijerg([a1], [a3, a4], [b1], [b3 + 1, b4], z), z)
    assert tn(MeijerShiftD(a3).apply(g, op),
              meijerg([a1], [a3 - 1, a4], [b1], [b3, b4], z), z)

    s = MeijerUnShiftA([a1], [a3, a4], [b1], [b3, b4], 0, z)
    assert tn(s.apply(g, op), meijerg([a1], [a3, a4], [b1 - 1], [b3, b4], z), z)

    s = MeijerUnShiftC([a1], [a3, a4], [b1], [b3, b4], 0, z)
    assert tn(s.apply(g, op), meijerg([a1], [a3, a4], [b1], [b3 - 1, b4], z), z)

    s = MeijerUnShiftB([a1], [a3, a4], [b1], [b3, b4], 0, z)
    assert tn(s.apply(g, op), meijerg([a1 + 1], [a3, a4], [b1], [b3, b4], z), z)

    s = MeijerUnShiftD([a1], [a3, a4], [b1], [b3, b4], 0, z)
    assert tn(s.apply(g, op), meijerg([a1], [a3 + 1, a4], [b1], [b3, b4], z), z)

def test_meijerg_confluence():
    def t(m, a, b):
        from sympy import sympify, Piecewise
        a, b = sympify([a, b])
        m_ = m
        m = hyperexpand(m)
        if not m == Piecewise((a, abs(z) < 1), (b, abs(1/z) < 1), (m_, True)):
            return False
        if not (m.args[0].args[0] == a and m.args[1].args[0] == b):
            return False
        z0 = randcplx()/10
        if abs(m.subs(z, z0).n() - a.subs(z, z0).n()).n() > 1e-10:
            return False
        if abs(m.subs(z, 1/z0).n() - b.subs(z, 1/z0).n()).n() > 1e-10:
            return False
        return True

    assert t(meijerg([], [1, 1], [0, 0], [], z), -log(z), 0)
    assert t(meijerg([], [3, 1], [0, 0], [], z), -z**2/4 + z - log(z)/2 - S(3)/4, 0)
    assert t(meijerg([], [3, 1], [-1, 0], [], z),
             z**2/12 -z/2 + log(z)/2 + S(1)/4 + 1/(6*z), 0)
    assert t(meijerg([], [1, 1, 1, 1], [0, 0, 0, 0], [], z), -log(z)**3/6, 0)
    assert t(meijerg([1, 1], [], [], [0, 0], z), 0, -log(1/z))
    assert t(meijerg([1, 1], [2, 2], [1, 1], [0, 0], z),
             -z*log(z) + 2*z, -log(1/z) + 2)
    assert t(meijerg([S(1)/2], [1, 1], [0, 0], [S(3)/2], z), log(z)/2 - 1, 0)

    def u(an, ap, bm, bq):
        m = meijerg(an, ap, bm, bq, z)
        m2 = hyperexpand(m, allow_hyper=True)
        if m2.has(meijerg) and not (m2.is_Piecewise and len(m2.args) == 3):
            return False
        return tn(m, m2, z)
    assert u([], [1], [0, 0], [])
    assert u([1, 1], [], [], [0])
    assert u([1, 1], [2, 2, 5], [1, 1, 6], [0, 0])
    assert u([1, 1], [2, 2, 5], [1, 1, 6], [0])

def test_lerchphi():
    from sympy import combsimp, exp_polar, polylog, log, lerchphi
    assert hyperexpand(hyper([1, a], [a + 1], z)/a) == lerchphi(z, 1, a)
    assert hyperexpand(hyper([1, a, a], [a + 1, a + 1], z)/a**2) == lerchphi(z, 2, a)
    assert hyperexpand(hyper([1, a, a, a], [a + 1, a + 1, a + 1], z)/a**3) == \
           lerchphi(z, 3, a)
    assert hyperexpand(hyper([1] + [a]*10, [a + 1]*10, z)/a**10) \
           == lerchphi(z, 10, a)
    assert combsimp(hyperexpand(meijerg([0, 1-a], [], [0], [-a],
                    exp_polar(-I*pi)*z))) == \
           lerchphi(z, 1, a)
    assert combsimp(hyperexpand(meijerg([0, 1-a, 1-a], [], [0], [-a, -a],
                    exp_polar(-I*pi)*z))) == \
           lerchphi(z, 2, a)
    assert combsimp(hyperexpand(meijerg([0, 1-a, 1-a, 1-a], [], [0], [-a, -a, -a],
                    exp_polar(-I*pi)*z))) == \
           lerchphi(z, 3, a)

    assert hyperexpand(z*hyper([1, 1], [2], z)) == -log(1 + -z)
    assert hyperexpand(z*hyper([1, 1, 1], [2, 2], z)) == polylog(2, z)
    assert hyperexpand(z*hyper([1, 1, 1, 1], [2, 2, 2], z)) == polylog(3, z)

    assert hyperexpand(hyper([1, a, 1 + S(1)/2], [a + 1, S(1)/2], z)) == \
           -2*a/(z - 1) + (-2*a**2 + a)*lerchphi(z, 1, a)

    # Now numerical tests. These make sure reductions etc are carried out
    # correctly

    # a rational function (polylog at negative integer order)
    assert can_do([2, 2, 2], [1, 1])

    # NOTE these contain log(1-x) etc ... better make sure we have |z| < 1
    # reduction of order for polylog
    assert can_do([1, 1, 1, b + 5], [2, 2, b], div=10)

    # reduction of order for lerchphi
    # XXX lerchphi in mpmath is flaky
    assert can_do([1, a, a, a, b + 5], [a + 1, a + 1, a + 1, b], numerical=False)

    # test a bug
    assert hyperexpand(hyper([S(1)/2, S(1)/2, S(1)/2, 1],
                             [S(3)/2, S(3)/2, S(3)/2], S(1)/4)) == \
           -polylog(3, exp_polar(I*pi)/2) + polylog(3, S(1)/2)

def test_partial_simp():
    # First test that hypergeometric function formulae work.
    a, b, c, d, e = map(lambda _: randcplx(), range(5))
    for idxp in [IndexPair([a, b, c], [d, e]), IndexPair([], [a, b, c, d, e])]:
        f = build_hypergeometric_formula(idxp)
        z = f.z
        assert f.closed_form == hyper(idxp.ap, idxp.bq, z)
        deriv1 = f.B.diff(z)*z
        deriv2 = f.M*f.B
        for func1, func2 in zip(deriv1, deriv2):
            assert tn(func1, func2, z)

    # Now test that formulae are partially simplified.
    from sympy.abc import a, b, z
    assert hyperexpand(hyper([3, a], [1, b], z)) == \
           (-a*b/2 + a*z/2 + 2*a)*hyper([a + 1], [b], z) \
         + (a*b/2 - 2*a + 1)*hyper([a], [b], z)
    assert tn(hyperexpand(hyper([3, d], [1, e], z)), hyper([3, d], [1, e], z), z)
    assert hyperexpand(hyper([3], [1, a, b], z)) == \
           hyper((), (a, b), z) \
           + z*hyper((), (a + 1, b), z)/(2*a) \
           - z*(b - 4)*hyper((), (a + 1, b + 1), z)/(2*a*b)
    assert tn(hyperexpand(hyper([3], [1, d, e], z)), hyper([3], [1, d, e], z), z)

def test_hyperexpand_special():
    assert hyperexpand(hyper([a, b], [c], 1)) == \
           gamma(c)*gamma(c - a - b)/gamma(c - a)/gamma(c - b)
    assert hyperexpand(hyper([a, b], [1 + a - b], -1)) == \
           gamma(1 + a/2)*gamma(1 + a - b)/gamma(1 + a)/gamma(1 + a/2 - b)
    assert hyperexpand(hyper([a, b], [1 + b - a], -1)) == \
           gamma(1 + b/2)*gamma(1 + b - a)/gamma(1 + b)/gamma(1 + b/2 - a)
    assert hyperexpand(meijerg([1 - z - a/2], [1 - z + a/2], [b/2], [-b/2], 1)) == \
           gamma(1 - 2*z)*gamma(z + a/2 + b/2)/gamma(1 - z + a/2 - b/2) \
           /gamma(1 - z - a/2 + b/2)/gamma(1 - z + a/2 + b/2)

def test_Mod1_behavior():
    from sympy import Symbol, simplify, lowergamma
    n = Symbol('n', integer=True)
    # Note: this should not hang.
    assert simplify(hyperexpand(meijerg([1], [], [n + 1], [0], z))) == \
           lowergamma(n + 1, z)

@slow
def test_prudnikov_misc():
    assert can_do([1, (3 + I)/2, (3 - I)/2], [S(3)/2, 2])
    assert can_do([S.Half, a - 1], [S(3)/2, a + 1], lowerplane=True)
    assert can_do([], [b + 1])
    assert can_do([a], [a - 1, b + 1])

    assert can_do([a], [a - S.Half, 2*a])
    assert can_do([a], [a - S.Half, 2*a + 1])
    assert can_do([a], [a - S.Half, 2*a - 1])
    assert can_do([a], [a + S.Half, 2*a])
    assert can_do([a], [a + S.Half, 2*a + 1])
    assert can_do([a], [a + S.Half, 2*a - 1])
    assert can_do([S.Half], [b, 2-b])
    assert can_do([S.Half], [b, 3-b])
    assert can_do([1], [2, b])

    assert can_do([a, a+S.Half], [2*a, b, 2*a - b + 1])
    assert can_do([a, a+S.Half], [S.Half, 2*a, 2*a + S.Half])
    assert can_do([a], [a+1], lowerplane=True) # lowergamma

@slow
def test_prudnikov_1():
    # A. P. Prudnikov, Yu. A. Brychkov and O. I. Marichev (1990).
    # Integrals and Series: More Special Functions, Vol. 3,.
    # Gordon and Breach Science Publisher

    # 7.3.1
    assert can_do([a, -a], [S.Half])
    assert can_do([a, 1 - a], [S.Half])
    assert can_do([a, 1 - a], [S(3)/2])
    assert can_do([a, 2 - a], [S.Half])
    assert can_do([a, 2 - a], [S(3)/2])
    assert can_do([a, 2 - a], [S(3)/2])
    assert can_do([a, a + S(1)/2], [2*a - 1])
    assert can_do([a, a + S(1)/2], [2*a])
    assert can_do([a, a + S(1)/2], [2*a + 1])
    assert can_do([a, a + S(1)/2], [S(1)/2])
    assert can_do([a, a + S(1)/2], [S(3)/2])
    assert can_do([a, a/2 + 1], [a/2])
    assert can_do([1, b], [2])
    assert can_do([1, b], [b + 1], numerical=False) # Lerch Phi
             # NOTE: branches are complicated for |z| > 1

    assert can_do([a], [2*a])
    assert can_do([a], [2*a + 1])
    assert can_do([a], [2*a - 1])

@slow
def test_prudnikov_2():
    h = S.Half
    assert can_do([-h, -h], [h])
    assert can_do([-h, h], [3*h])
    assert can_do([-h, h], [5*h])
    assert can_do([-h, h], [7*h])
    assert can_do([-h, 1], [h])

    for p in [-h, h]:
      for n in [-h, h, 1, 3*h, 2, 5*h, 3, 7*h, 4]:
          for m in [-h, h, 3*h, 5*h, 7*h]:
              assert can_do([p, n], [m])
      for n in [1, 2, 3, 4]:
          for m in [1, 2, 3, 4]:
              assert can_do([p, n], [m])

@slow
def test_prudnikov_3():
    h = S.Half
    assert can_do([S(1)/4, S(3)/4], [h])
    assert can_do([S(1)/4, S(3)/4], [3*h])
    assert can_do([S(1)/3, S(2)/3], [3*h])
    assert can_do([S(3)/4, S(5)/4], [h])
    assert can_do([S(3)/4, S(5)/4], [3*h])

    for p in [1, 2, 3, 4]:
      for n in [-h, h, 1, 3*h, 2, 5*h, 3, 7*h, 4, 9*h]:
          for m in [1, 3*h, 2, 5*h, 3, 7*h, 4]:
              assert can_do([p, m], [n])


@slow
def test_prudnikov_4():
    h = S.Half
    for p in [3*h, 5*h, 7*h]:
      for n in [-h, h, 3*h, 5*h, 7*h]:
          for m in [3*h, 2, 5*h, 3, 7*h, 4]:
              assert can_do([p, m], [n])
      for n in [1, 2, 3, 4]:
          for m in [2, 3, 4]:
              assert can_do([p, m], [n])

@slow
def test_prudnikov_5():
    h = S.Half

    for p in [1, 2, 3]:
        for q in range(p, 4):
            for r in [1, 2, 3]:
                for s in range(r, 4):
                    assert can_do([-h, p, q], [r, s])

    for p in [h, 1, 3*h, 2, 5*h, 3]:
        for q in [h, 3*h, 5*h]:
            for r in [h, 3*h, 5*h]:
                for s in [h, 3*h, 5*h]:
                    if s <= q and s <= r:
                        assert can_do([-h, p, q], [r, s])

    for p in [h, 1, 3*h, 2, 5*h, 3]:
        for q in [1, 2, 3]:
            for r in [h, 3*h, 5*h]:
                for s in [1, 2, 3]:
                    assert can_do([-h, p, q], [r, s])

@slow
def test_prudnikov_6():
    h = S.Half

    for m in [3*h, 5*h]:
        for n in [1, 2, 3]:
            for q in [h, 1, 2]:
                for p in [1, 2, 3]:
                     assert can_do([h, q, p], [m, n])
            for q in [1, 2, 3]:
                for p in [3*h, 5*h]:
                     assert can_do([h, q, p], [m, n])

    for q in [1, 2]:
      for p in [1, 2, 3]:
         for m in [1, 2, 3]:
             for n in [1, 2, 3]:
                 assert can_do([h, q, p], [m, n])

    assert can_do([h, h, 5*h], [3*h, 3*h])
    assert can_do([h, 1, 5*h], [3*h, 3*h])
    assert can_do([h, 2, 2], [1, 3])

    # pages 435 to 457 contain more PFDD and stuff like this

@slow
def test_prudnikov_7():
    assert can_do([3], [6])

    h = S.Half
    for n in [h, 3*h, 5*h, 7*h]:
        assert can_do([-h], [n])
    for m in [-h, h, 1, 3*h, 2, 5*h, 3, 7*h, 4]: # HERE
        for n in [-h, h, 3*h, 5*h, 7*h, 1, 2, 3, 4]:
            assert can_do([m], [n])

@slow
def test_prudnikov_8():
    h = S.Half

    # 7.12.2
    for a in [1, 2, 3]:
        for b in [1, 2, 3]:
            for c in range(1, a+1):
                for d in [h, 1, 3*h, 2, 5*h, 3]:
                    assert can_do([a, b], [c, d])
        for b in [3*h, 5*h]:
            for c in [h, 1, 3*h, 2, 5*h, 3]:
                for d in [1, 2, 3]:
                    assert can_do([a, b], [c, d])

    for a in [-h, h, 3*h, 5*h]:
        for b in [1, 2, 3]:
            for c in [h, 1, 3*h, 2, 5*h, 3]:
                for d in [1, 2, 3]:
                    assert can_do([a, b], [c, d])
        for b in [h, 3*h, 5*h]:
            for c in [h, 3*h, 5*h, 3]:
                for d in [h, 1, 3*h, 2, 5*h, 3]:
                    if c <= b:
                        assert can_do([a, b], [c, d])

@slow
def test_prudnikov_9():
    # 7.13.1 [we have a general formula ... so this is a bit pointless]
    for i in range(9):
        assert can_do([], [(S(i) + 1)/2])
    for i in range(5):
        assert can_do([], [-(2*S(i) + 1)/2])

@slow
def test_prudnikov_10():
    # 7.14.2
    h = S.Half
    for p in [-h, h, 1, 3*h, 2, 5*h, 3, 7*h, 4]:
      for m in [1, 2, 3, 4]:
          for n in range(m, 5):
              assert can_do([p], [m, n])

    for p in [1, 2, 3, 4]:
      for n in [h, 3*h, 5*h, 7*h]:
          for m in [1, 2, 3, 4]:
              assert can_do([p], [n, m])

    for p in [3*h, 5*h, 7*h]:
      for m in [h, 1, 2, 5*h, 3, 7*h, 4]:
         assert can_do([p], [h, m])
         assert can_do([p], [3*h, m])

    for m in [h, 1, 2, 5*h, 3, 7*h, 4]:
       assert can_do([7*h], [5*h, m])

    assert can_do([-S(1)/2], [S(1)/2, S(1)/2]) # shine-integral shi

@slow
def test_prudnikov_11():
    # 7.15
    assert can_do([a, a+S.Half], [2*a, b, 2*a - b])
    assert can_do([a, a+S.Half], [S(3)/2, 2*a, 2*a - S(1)/2])

    assert can_do([S(1)/4, S(3)/4], [S(1)/2, S(1)/2, 1])
    assert can_do([S(5)/4, S(3)/4], [S(3)/2, S(1)/2, 2])
    assert can_do([S(5)/4, S(3)/4], [S(3)/2, S(3)/2, 1])
    assert can_do([S(5)/4, S(7)/4], [S(3)/2, S(5)/2, 2])

    assert can_do([1, 1], [S(3)/2, 2, 2]) # cosh-integral chi

@slow
def test_prudnikov_12():
    # 7.16
    assert can_do([], [a, a + S.Half, 2*a], False) # branches only agree for some z!
    assert can_do([], [a, a + S.Half, 2*a+1], False) # dito
    assert can_do([], [S.Half, a, a+S.Half])
    assert can_do([], [S(3)/2, a, a+S.Half])

    assert can_do([], [S(1)/4, S(1)/2, S(3)/4])
    assert can_do([], [S(1)/2, S(1)/2, 1])
    assert can_do([], [S(1)/2, S(3)/2, 1])
    assert can_do([], [S(3)/4, S(3)/2, S(5)/4])
    assert can_do([], [1, 1, S(3)/2])
    assert can_do([], [1, 2, S(3)/2])
    assert can_do([], [1, S(3)/2, S(3)/2])
    assert can_do([], [S(5)/4, S(3)/2, S(7)/4])
    assert can_do([], [2, S(3)/2, S(3)/2])

@XFAIL
def test_prudnikov_fail_2F1():
    assert can_do([a, b], [b + 1]) # incomplete beta function
    assert can_do([-1, b], [c])    # Poly. also -2, -3 etc

    # TODO polys

    # Legendre functions:
    assert can_do([a, b], [a + b + S.Half])
    assert can_do([a, b], [a + b - S.Half])
    assert can_do([a, b], [a + b + S(3)/2])
    assert can_do([a, b], [(a + b + 1)/2])
    assert can_do([a, b], [(a + b)/2 + 1])
    assert can_do([a, b], [a - b + 1])
    assert can_do([a, b], [a - b + 2])
    assert can_do([a, b], [2*b])
    assert can_do([a, b], [S.Half])
    assert can_do([a, b], [S(3)/2])
    assert can_do([a, 1 - a], [c])
    assert can_do([a, 2 - a], [c])
    assert can_do([a, 3 - a], [c])
    assert can_do([a, a + S(1)/2], [c])
    assert can_do([1, b], [c])
    assert can_do([1, b], [S(3)/2])

    h = S.Half
    # Elliptic integrals
    for p in [-h, h]:
        for m in [h, 3*h, 5*h, 7*h]:
            for n in [1, 2, 3, 4]:
                assert can_do([p, m], [n])

    assert can_do([S(1)/4, S(3)/4], [1])

    # PFDD
    o = S(1)
    assert can_do([o/8, 1], [o/8*9])
    assert can_do([o/6, 1], [o/6*7])
    assert can_do([o/6, 1], [o/6*13])
    assert can_do([o/5, 1], [o/5*6])
    assert can_do([o/5, 1], [o/5*11])
    assert can_do([o/4, 1], [o/4*5])
    assert can_do([o/4, 1], [o/4*9])
    assert can_do([o/3, 1], [o/3*4])
    assert can_do([o/3, 1], [o/3*7])
    assert can_do([o/8*3, 1], [o/8*11])
    assert can_do([o/5*2, 1], [o/5*7])
    assert can_do([o/5*2, 1], [o/5*12])
    assert can_do([o/5*3, 1], [o/5*8])
    assert can_do([o/5*3, 1], [o/5*13])
    assert can_do([o/8*5, 1], [o/8*13])
    assert can_do([o/4*3, 1], [o/4*7])
    assert can_do([o/4*3, 1], [o/4*11])
    assert can_do([o/3*2, 1], [o/3*5])
    assert can_do([o/3*2, 1], [o/3*8])
    assert can_do([o/5*4, 1], [o/5*9])
    assert can_do([o/5*4, 1], [o/5*14])
    assert can_do([o/6*5, 1], [o/6*11])
    assert can_do([o/6*5, 1], [o/6*17])
    assert can_do([o/8*7, 1], [o/8*15])

@XFAIL
def test_prudnikov_fail_3F2():
    assert can_do([a, a + S(1)/3, a + S(2)/3], [S(1)/3, S(2)/3])
    assert can_do([a, a + S(1)/3, a + S(2)/3], [S(2)/3, S(4)/3])
    assert can_do([a, a + S(1)/3, a + S(2)/3], [S(4)/3, S(5)/3])

    # page 421
    assert can_do([a, a + S(1)/3, a + S(2)/3], [3*a/2, (3*a+1)/2])

    # pages 422 ...
    assert can_do([-S.Half, S.Half, S.Half], [1, 1]) # elliptic integrals
    assert can_do([-S.Half, S.Half, 1], [S(3)/2, S(3)/2])
    # TODO LOTS more

    # PFDD
    assert can_do([S(1)/8, S(3)/8, 1], [S(9)/8, S(11)/8])
    assert can_do([S(1)/8, S(5)/8, 1], [S(9)/8, S(13)/8])
    assert can_do([S(1)/8, S(7)/8, 1], [S(9)/8, S(15)/8])
    assert can_do([S(1)/6, S(1)/3, 1], [S(7)/6, S(4)/3])
    assert can_do([S(1)/6, S(2)/3, 1], [S(7)/6, S(5)/3])
    assert can_do([S(1)/6, S(2)/3, 1], [S(5)/3, S(13)/6])
    assert can_do([S.Half, 1, 1], [S(1)/4, S(3)/4])
    # LOTS more


@XFAIL
def test_prudnikov_fail_other():
    # 7.11.2

    # 7.12.1
    assert can_do([1, a], [b, 1 - 2*a + b]) # ???

    # 7.14.2
    assert can_do([-S(1)/2], [S(1)/2, 1]) # struve
    assert can_do([1], [S(1)/2, S(1)/2])  # struve
    assert can_do([S(1)/4], [S(1)/2, S(5)/4]) # PFDD
    assert can_do([S(3)/4], [S(3)/2, S(7)/4]) # PFDD
    assert can_do([1], [S(1)/4, S(3)/4]) # PFDD
    assert can_do([1], [S(3)/4, S(5)/4]) # PFDD
    assert can_do([1], [S(5)/4, S(7)/4]) # PFDD
    # TODO LOTS more

    # 7.15.2
    assert can_do([S(1)/2, 1], [S(3)/4, S(5)/4, S(3)/2]) # PFDD
    assert can_do([S(1)/2, 1], [S(7)/4, S(5)/4, S(3)/2]) # PFDD

    # 7.16.1
    assert can_do([], [S(1)/3, S(2/3)]) # PFDD
    assert can_do([], [S(2)/3, S(4/3)]) # PFDD
    assert can_do([], [S(5)/3, S(4/3)]) # PFDD

    # XXX this does not *evaluate* right??
    assert can_do([], [a, a + S.Half, 2*a-1])
