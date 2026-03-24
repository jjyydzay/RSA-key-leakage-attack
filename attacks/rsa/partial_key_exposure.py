import logging
import os
import sys
import time
from math import ceil, floor
from math import gcd
from math import sqrt
from sympy import gcdex
import multiprocessing
from multiprocessing import Pool, cpu_count

from sage.all import sqrt
from sage.all import var
from sage.rings.polynomial.polynomial_modn_dense_ntl import small_roots
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.all import QQ
from sage.all import RR
from sage.all import ZZ
from sage.all import Zmod
from sage.all import is_prime, power_mod
from sage.all import parallel
from sage.rings.integer import Integer
from shared.small_roots import jochemsz_may_modular
from shared.small_roots import jochemsz_may_integer
from sage.all import ZZ
from sage.all import continued_fraction
from attacks.factorization import known_phi
from shared.small_roots import herrmann_may

path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(os.path.abspath(__file__)))))
if sys.path[1] != path:
    sys.path.insert(1, path)

from attacks.factorization import known_phi
from attacks.rsa import small_exponent
from shared.hensel import hensel_roots
from shared.small_roots import blomer_may
from shared.small_roots import ernst
from shared.small_roots import howgrave_graham
# from shared.Shiho_small_root import coron
sys.set_int_max_str_digits(0)


def Wiener(N, e):
    """
    Recovers the prime factors of a modulus and the private exponent if the private exponent is too small.
    :param N: the modulus
    :param e: the public exponent
    :return: a tuple containing the prime factors and the private exponent, or None if the private exponent was not found
    """
    convergents = continued_fraction(ZZ(e) / ZZ(N)).convergents()
    for c in convergents:
        k = c.numerator()
        d = c.denominator()
        if pow(pow(2, e, N), d, N) != 2:
            continue

        phi = (e * d - 1) // k
        factors = known_phi.factorize(N, phi)
        if factors:
            return *factors, int(d)

def BonehDurfee(N, e, factor_bit_length, partial_p=None, delta=0.25, m=1, t=None, roots_method="groebner"):
    """
    Recovers the prime factors if the private exponent is too small.
    This implementation exploits knowledge of least significant bits of prime factors, if available.
    More information: Boneh D., Durfee G., "Cryptanalysis of RSA with Private Key d Less than N^0.292"
    :param N: the modulus
    :param e: the public exponent
    :param factor_bit_length: the bit length of the prime factors
    :param partial_p: the partial prime factor p (PartialInteger) (default: None)
    :param delta: a predicted bound on the private exponent (d < N^delta) (default: 0.25)
    :param m: the m value to use for the small roots method (default: 1)
    :param t: the t value to use for the small roots method (default: automatically computed using m)
    :return: a tuple containing the prime factors, or None if the factors were not found
    """
    # Use additional information about factors to speed up Boneh-Durfee.
    p_lsb, p_lsb_bit_length = (0, 0) if partial_p is None else partial_p.get_known_lsb()
    q_lsb = (pow(p_lsb, -1, 2 ** p_lsb_bit_length) * N) % (2 ** p_lsb_bit_length)
    A = ((N >> p_lsb_bit_length) + pow(2, -p_lsb_bit_length, e) * (p_lsb * q_lsb - p_lsb - q_lsb + 1))

    x, y = ZZ["x", "y"].gens()
    f = x * (A + y) + pow(2, -p_lsb_bit_length, e)
    X = int(RR(e) ** delta)
    Y = int(2 ** (factor_bit_length - p_lsb_bit_length + 1))
    # print("x, y = ZZ[\"x\", \"y\"].gens()")
    # print("pol =", f)
    # print("XX =", X)
    # print("YY =", Y)
    t = int((1 - 2 * delta) * m) if t is None else t
    logging.info(f"Trying {m = }, {t = }...")
    for x0, y0 in herrmann_may.modular_bivariate(f, e, m, t, X, Y, roots_method=roots_method):
        z = int(f(x0, y0))
        if z % e == 0:
            k = pow(x0, -1, e)
            s = (N + 1 + k) % e
            phi = N - s + 1
            factors = known_phi.factorize(N, phi)
            print(factors)
            p, q = factors
            phi_n = (p - 1) * (q - 1)
            d, _, _ = gcdex(e, phi_n)
            d = d % phi_n
            if factors:
                return *factors, int(d)

    return None


def _bdf_corollary_1_prime(e, f, N, m, t, X):
    logging.debug(f"Solving f wiht root x0=0")
    p = abs(int(f(0)))
    if 1 < p < N and N % p == 0:
        q = N // p
        phi = (p - 1) * (q - 1)
        yield p, q, pow(e, -1, phi)
    for x0, in howgrave_graham.modular_univariate(f, N, m, t, X):
        p = abs(int(f(x0)))
        if 1 < p < N and N % p == 0:
            q = N // p
            phi = (p - 1) * (q - 1)
            yield p, q, pow(e, -1, phi)    
            
def _bdf_theorem_3_3(N, e, d_bit_length, d1, d1_bit_length, m, t, thetaLogN=2, k=None, known_p=None):
    logging.info(f"Trying {m = }, {t = }...")
    p = Zmod(e)["p"].gen()
    x = Zmod(N)["x"].gen()
    X = int(2**(thetaLogN)*RR(N)/(2**d1_bit_length))
    if k is None:
        d0 = d1 << (d_bit_length - d1_bit_length)
        k_ = (e * d0 - 1) // N
        denominator = int((1 - 2**(-100)) * N)
        k_prime = (e * (d0 + 2 ** (d_bit_length - d1_bit_length)) - 1) // denominator
        print("k bound is:", k_, "and", k_prime)
        logging.info("Generating solutions for k candidates...")
        for k in range(k_ - 1, k_prime + 2):
            S = int(N+1-(e*d0-1)//k)
            D = int(sqrt(S**2-4*N))
            pApproximation = int((S+D)//2)
            print(f"Solving RSA with k={k}")
            logging.debug(f"Solving RSA with k={k}")
            if known_p is not None:
                delta_p = abs(known_p-pApproximation)
                logging.debug(f"Solving univariate equation with {int(delta_p).bit_length()} bits")
                print(f"Solving univariate equation with {int(delta_p).bit_length()} bits")
            logging.debug(f"Set coppersmith bound {int(X).bit_length()} bits")
            print(f"Set coppersmith bound {int(X).bit_length()} bits")
            f = x + pApproximation
            for p_, q_, d_ in _bdf_corollary_1_prime(e, f, N, m, t, X):
                return p_, q_, d_

    return None


def _bdf_corollary_1(e, f, N, m, t, X):

    for x0, in howgrave_graham.modular_univariate(f, N, m, t, X):
        p = int(f(x0))
        if 1 < p < N and N % p == 0:
            q = N // p
            phi = (p - 1) * (q - 1)
            yield p, q, pow(e, -1, phi)


def _bdf_theorem_6(N, e, d_bit_length, d1, d1_bit_length):
    d0 = d1 << (d_bit_length - d1_bit_length)
    k_ = (e * d0 - 1) // N
    logging.info("Generating solutions for k candidates...")
    for k in range(k_ - 40, k_ + 40):
        yield d0, k


def _bdf_3(N, e, d_bit_length, d0, d0_bit_length, r, m, t):
    n = N.bit_length()
    logging.info(f"Trying {m = }, {t = }...")
    p = ZZ["p"].gen()
    x = Zmod(N)["x"].gen()
    X = int(2 * RR(N) ** (1 / 2) / r)  # Equivalent to 2^(n / 2 + 1) / r
    logging.info("Generating solutions for k candidates...")
    for k in range(1, e):
        f = k * p ** 2 + (e * d0 - (1 + k * (N + 1))) * p + k * N
        for p0 in hensel_roots(f, 2, d0_bit_length):
            f = x * r + p0
            for p_, q_, d_ in _bdf_corollary_1(e, f, N, m, t, X):
                return p_, q_, d_

    return None


def _bdf_4_1(N, e, d_bit_length, d1, d1_bit_length, m, t):
    logging.info(f"Trying {m = }, {t = }...")
    p = Zmod(e)["p"].gen()
    x = Zmod(N)["x"].gen()
    X = int(2 * RR(N) ** (1 / 2) / e)  # Equivalent to 2^(n / 2 + 1) / e
    for _, k in _bdf_theorem_6(N, e, d_bit_length, d1, d1_bit_length):
        f = k * p ** 2 - (1 + k * (N + 1)) * p + k * N
        for p0 in f.roots(multiplicities=False):
            f = x * e + int(p0)
            for p_, q_, d_ in _bdf_corollary_1(e, f, N, m, t, X):
                return p_, q_, d_

    return None


def _bdf_4_2(N, e, d_bit_length, d1, d1_bit_length):
    for d0, k in _bdf_theorem_6(N, e, d_bit_length, d1, d1_bit_length):
        if gcd(e, k) != 1:
            continue

        d1 = pow(e, -1, k)
        for v in range(ceil(e / k) + 1):
            d2 = int(QQ(d0) / k + v - QQ(d1) / k)
            d = k * d2 + d1
            if pow(pow(2, e, N), d, N) == 2:
                phi = (e * d - 1) // k
                factors = known_phi.factorize(N, phi)
                if factors:
                    return *factors, d

    return None


def _bdf_4_3(N, e, d_bit_length, d0, d0_bit_length, d1, d1_bit_length, r, m, t):
    logging.info(f"Trying {m = }, {t = }...")
    p = ZZ["p"].gen()
    x = Zmod(N)["x"].gen()
    X = int(2 * RR(N) ** (1 / 2) / r)  # Equivalent to 2^(n / 2 + 1) / r

    for _, k in _bdf_theorem_6(N, e, d_bit_length, d1, d1_bit_length):
        f = k * p ** 2 + (e * d0 - (1 + k * (N + 1))) * p + k * N
        # print("f =", f)
        for p0 in hensel_roots(f, 2, d0_bit_length):
            f = x * r + p0  
            for p_, q_, d_ in _bdf_corollary_1(e, f, N, m, t, X):
                return p_, q_, d_

    return None


def _bm_4(N, e, d_bit_length, d1, d1_bit_length, m, t):
    d_ = d1 << (d_bit_length - d1_bit_length)
    k_ = (e * d_ - 1) // (N + 1)

    x, y, z = ZZ["x", "y", "z"].gens()
    # x1, x2, x3 = ZZ["x1", "x2", "x3"].gens()
    f = e * x + (k_ + y) * z + e * d_ - 1
    # f = e * x1 + (k_ + x3) * x2 + e * d_ - 1

    X = 2 ** (d_bit_length - d1_bit_length)  # Equivalent to N^delta
    Y = int(4 * e / RR(N) ** (1 / 2))  # Equivalent to 4N^(alpha - 1 / 2)
    Z = int(3 * RR(N) ** (1 / 2))
    print(f"Trying {m = }, {t = }...")
    for x0, y0, z0 in blomer_may.modular_trivariate(f, N, m, t, X, Y, Z):  #groebner
        print("x0 =", x0, "y0 =", y0, "z0 =", z0)
        d = d_ + x0
        phi = N - z0
        if pow(pow(2, e, N), d, N) == 2:
        # if power_mod(power_mod(2, e, N), d, N) == 2:
            factors = known_phi.factorize(N, phi)
            if factors:
                return *factors, d

    return None


def _bm_6(N, e, d_bit_length, d0, d0_bit_length, M, m, t):
    # x1, x2 = ZZ["x1", "x2"].gens()
    y, z = ZZ["y", "z"].gens()
    f = y * (N - z) - e * d0 + 1
    Y = e  # Equivalent to N^alpha
    Z = int(3 * RR(N) ** (1 / 2))
    start_time = time.time()
    blomer_may.modular_bivariate(f, e * M, m, t, Y, Z, roots_method='julian')
    end_time = time.time()
    print("modular_bivariate:", end_time - start_time)
    for y0, z0 in blomer_may.modular_bivariate(f, e * M, m, t, Y, Z, roots_method='julian'):
        print("y0 =", y0, "z0 =", z0)
        phi = N - z0
        d = pow(e, -1, phi)
        d = Integer(d)
        if pow(pow(2, e, N), d, N) == 2:
            factors = known_phi.factorize(N, phi)
            if factors:
                return *factors, d

    return None


def _ernst_4_1_1(N, e, d_bit_length, d1, d1_bit_length, m, t):
    d_ = d1 << (d_bit_length - d1_bit_length)
    R = e * d_ - 1

    x, y, z = ZZ["x", "y", "z"].gens()
    f = e * x - N * y + y * z + R
    X = 2 ** (d_bit_length - d1_bit_length)  # Equivalent to N^delta
    Y = 2 ** d_bit_length  # Equivalent to N^beta
    Z = int(3 * RR(N) ** (1 / 2))
    W = N * Y
    logging.info(f"Trying {m = }, {t = }...")
    # for x0, y0, z0 in ernst.integer_trivariate_1(f, m, t, W, X, Y, Z, roots_method='julian'):
    print("[WARNING] We are using jochemsz_may_integer to substitute ernst-1.")
    strategy = jochemsz_may_integer.Ernst1Strategy(t)
    for x0, y0, z0 in jochemsz_may_integer.integer_multivariate(f, m, W, [X, Y, Z], strategy, roots_method='julian'):
        d = d_ + x0
        print("we found d:", d)
        phi = N - z0
        if pow(pow(2, e, N), d, N) == 2:
            factors = known_phi.factorize(N, phi)
            if factors:
                return *factors, d

    return None

def _ernst_Appendix_C(N, e, d_bit_length, d0, d0_bit_length, d1, d1_bit_length, m, t):
    M1 = 2 ** d0_bit_length
    M2 = 2 ** (d_bit_length - d1_bit_length)
    R = e * d0  + e * M2 * d1 - 1
    x, y, z = ZZ["x", "y", "z"].gens()
    f = e * M1 * x - N * y + y * z + R
    X = 2 ** (d_bit_length - d1_bit_length - d0_bit_length)  # Equivalent to N^delta
    Y = 2 ** d_bit_length  # Equivalent to N^beta
    Z = int(3 * RR(N) ** (1 / 2))
    W = max(e*M1*X, N*Y, Y*Z, R)
    # print(f)
    # print("X=",X,",Y=",Y,",Z=",Z)
    # tau = t / m
    # print("RR(X) ** (1 + 3 * tau) * RR(Y) ** (2 + 3 * tau) * RR(Z) ** (1 + 3 * tau + 3 * tau ** 2):\n",
    #     RR(X) ** (1 + 3 * tau) * RR(Y) ** (2 + 3 * tau) * RR(Z) ** (1 + 3 * tau + 3 * tau ** 2))
    # print("W**(1 + 3 * tau):\n",RR(W) ** (1 + 3 * tau))
    logging.info(f"Trying {m = }, {t = }...")
    print("Trying m=",m,",t=",t)
    # print("[WARNING] We are using jochemsz_may_integer to substitute ernst-1.")
    # strategy = jochemsz_may_integer.Ernst1Strategy([t, 0, 0])
    # for x0, y0, z0 in jochemsz_may_integer.integer_multivariate(f, m, W, [X, Y, Z], strategy):
    for x0, y0, z0 in ernst.integer_trivariate_1(f, m, t, W, X, Y, Z, check_bounds=False, roots_method='julian'):
        d = d0 + M1 * x0 + M2 * d1
        print("we found d:", d)
        phi = N - z0
        if pow(pow(2, e, N), d, N) == 2:
            factors = known_phi.factorize(N, phi)
            if factors:
                return *factors, d

    return None

def _ernst_4_1_2(N, e, d_bit_length, d1, d1_bit_length, m, t):
    d_ = d1 << (d_bit_length - d1_bit_length)
    k_ = (e * d_ - 1) // N
    R = e * d_ - 1 - k_ * N

    x, y, z = ZZ["x", "y", "z"].gens()
    f = e * x - N * y + y * z + k_ * z + R
    X = 2 ** (d_bit_length - d1_bit_length)  # Equivalent to N^delta
    Y = 4 * int(max(2 ** (d_bit_length - d1_bit_length), 2 ** d_bit_length / RR(N) ** (1 / 2)))  # Equivalent to 4N^max(delta, beta - 1 / 2)
    Z = int(3 * RR(N) ** (1 / 2))
    W = N * Y
    logging.info(f"Trying {m = }, {t = }...")
    # for x0, y0, z0 in ernst.integer_trivariate_2(f, m, t, W, X, Y, Z, roots_method='julian'):
    print("[WARNING] We are using jochemsz_may_integer to substitute ernst-2.")
    strategy = jochemsz_may_integer.Ernst2Strategy(t)
    for x0, y0, z0 in jochemsz_may_integer.integer_multivariate(f, m, W, [X, Y, Z], strategy, roots_method='julian'):
        d = d_ + x0
        phi = N - z0
        if pow(pow(2, e, N), d, N) == 2:
            factors = known_phi.factorize(N, phi)
            if factors:
                return *factors, d

    return None

def _ernst_Appendix_C_1(N, e, d_bit_length, d0, d0_bit_length, d1, d1_bit_length, m, t):
    M1 = 2 ** d0_bit_length
    M2 = 2 ** (d_bit_length - d1_bit_length)
    d_ = d1 << (d_bit_length - d1_bit_length)
    k_ = (e * d_ - 1) // N
    R = e * d0  + e * M2 * d1 - 1 - k_ * N

    x, y, z = ZZ["x", "y", "z"].gens()
    f = e * M1 * x - N * y + y * z + k_ * z + R
    X = 2 ** (d_bit_length - d1_bit_length - d0_bit_length)  # Equivalent to N^delta
    Y = 4 * int(max(2 ** (d_bit_length - d1_bit_length), 2 ** d_bit_length / RR(N) ** (1 / 2)))  # Equivalent to 4N^max(delta, beta - 1 / 2)
    Z = int(3 * RR(N) ** (1 / 2))
    W = max(e*M1*X, N*Y, Y*Z, k_*Z, R)
    logging.info(f"Trying {m = }, {t = }...")
    print("Trying m=",m,",t=",t)
    for x0, y0, z0 in ernst.integer_trivariate_2(f, m, t, W, X, Y, Z, check_bounds=False, roots_method='julian'):
        d = d0 + x0 * M1 + d1 * M2
        print("we found d:", d)
        phi = N - z0
        if pow(pow(2, e, N), d, N) == 2:
            factors = known_phi.factorize(N, phi)
            if factors:
                return *factors, d

    return None

def _ernst_4_2(N, e, d_bit_length, d1, d1_bit_length, m, t):
    d_ = d1 << (d_bit_length - d1_bit_length)
    k_ = (e * d_ - 1) // N
    R = e * d_ - 1 - k_ * N

    x, y, z = ZZ["x", "y", "z"].gens()
    f = e * x - N * y + y * z + k_ * z + R
    X = 2 ** (d_bit_length - d1_bit_length)  # Equivalent to N^delta
    Y = 4 * int(max((e * 2 ** (d_bit_length - d1_bit_length)) / N, e / RR(N) ** (1 / 2)))  # Equivalent to 4N^max(alpha + delta - 1, alpha - 1 / 2)
    Z = int(3 * RR(N) ** (1 / 2))
    W = N * Y
    # print(f)
    # print("X=",X,",Y=",Y,",Z=",Z)
    logging.info(f"Trying {m = }, {t = }...")
    # for x0, y0, z0 in ernst.integer_trivariate_2(f, m, t, W, X, Y, Z, roots_method='julian'):
    print("[WARNING] We are using jochemsz_may_integer to substitute ernst-2.")
    strategy = jochemsz_may_integer.Ernst2Strategy(t)
    for x0, y0, z0 in jochemsz_may_integer.integer_multivariate(f, m, W, [X, Y, Z], strategy, roots_method='julian'):
        d = d_ + x0
        # print("d=", bin(d))
        phi = N - z0
        if pow(pow(2, e, N), d, N) == 2:
            factors = known_phi.factorize(N, phi)
            if factors:
                return *factors, d

    return None

def _ernst_Appendix_C_2(N, e, d_bit_length, d0, d0_bit_length, d1, d1_bit_length, m, t):
    M1 = 2 ** d0_bit_length
    M2 = 2 ** (d_bit_length - d1_bit_length)
    d_ = d1 << (d_bit_length - d1_bit_length)
    k_ = (e * d_ - 1) // N
    R = e * d0  + e * M2 * d1 - 1 - k_ * N

    x, y, z = ZZ["x", "y", "z"].gens()
    f = e * M1 * x - N * y + y * z + k_ * z + R
    X = 2 ** (d_bit_length - d1_bit_length - d0_bit_length)  # Equivalent to N^delta
    Y = 4 * int(max((e * 2 ** (d_bit_length - d1_bit_length)) / N, e / RR(N) ** (1 / 2)))  # Equivalent to 4N^max(alpha + delta - 1, alpha - 1 / 2)
    Z = int(3 * RR(N) ** (1 / 2))
    W = max(e*M1*X, N*Y, Y*Z, k_*Z, R)
    logging.info(f"Trying {m = }, {t = }...")
    print("Trying m=",m,",t=",t)
    for x0, y0, z0 in ernst.integer_trivariate_2(f, m, t, W, X, Y, Z, check_bounds=False, roots_method='julian'):
        d = d0 + x0 * M1 + d1 * M2
        phi = N - z0
        if pow(pow(2, e, N), d, N) == 2:
            factors = known_phi.factorize(N, phi)
            if factors:
                return *factors, d

    return None    

def _ernst_4_3(N, e, d_bit_length, d0, d0_bit_length, M, m, t):
    R = e * d0 - 1

    x, y, z = ZZ["x", "y", "z"].gens()
    f = e * M * x - N * y + y * z + R
    X = 2 ** (d_bit_length - d0_bit_length)  # Equivalent to N^delta
    Y = 2 ** d_bit_length  # Equivalent to N^beta
    Z = int(3 * RR(N) ** (1 / 2))
    W = N * Y
    logging.info(f"Trying {m = }, {t = }...")
    
    # for x0, y0, z0 in ernst.integer_trivariate_1(f, m, t, W, X, Y, Z, roots_method='julian'):
    print("[WARNING] We are using jochemsz_may_integer to substitute ernst-1.")
    strategy = jochemsz_may_integer.Ernst1Strategy(t)
    for x0, y0, z0 in jochemsz_may_integer.integer_multivariate(f, m, W, [X, Y, Z], strategy, roots_method='julian'):
        d = x0 * M + d0
        phi = N - z0
        if pow(pow(2, e, N), d, N) == 2:
            factors = known_phi.factorize(N, phi)
            if factors:
                return *factors, d

    return None

def pre_attack(N, e, partial_d, factor_e=True, m=1, t=None, traversal_bit_length=10, option=None):
    d_length = partial_d.bit_length
    dl, dl_length = partial_d.get_known_lsb()
    dm, dm_length = partial_d.get_known_msb()
    print("===test===")
    if(dl_length > 0 and dm_length > 0):
        print("Known Both LSBs and MSBs Attack")
    if(dl_length > 0 and dm_length == 0):
        print("Known LSBs Attack")
    if(dl_length == 0 and dm_length > 0):
        print("Known MSBs Attack")
    if(dl_length == 0 and dm_length == 0):
        print("Small Private Exponent Attack")
    print("N =", N)
    print("e =", e)
    print("d_length =", d_length)
    try:
        print("e length =", e.bit_length())
    except:
        print("e length =", e.nbits())
    print("dl =", dl)
    print("dl_length =", dl_length)
    print("dm =", dm)
    print("dm_length =", dm_length)
    print("m =", m)
    print("t =", t)
    print("MSBs =", dm_length / d_length)
    print("LSBs =", dl_length / d_length)

    print("===test===")


    return attack(N, e, d_length, dl, dl_length, dm, dm_length, factor_e=factor_e, m=m, t=t, traversal_bit_length=traversal_bit_length, option=option)
    # return 0,0,0
'''
def attack_single_lsbit(args):
    N, e, d_length, dl, dl_length, dm, dm_length, factor_e, m, t, lsbits, traversal_bit_length = args
    ddll = lsbits * 2 ** dl_length + dl
    print("lsbits", lsbits)
    ddll_length = dl_length + traversal_bit_length
    ddmm = dm
    ddmm_length = dm_length
    return attack0(N, e, d_length, ddll, ddll_length, ddmm, ddmm_length, factor_e, m, t)

def attack_single_msbit(args):
    N, e, d_length, dl, dl_length, dm, dm_length, factor_e, m, t, msbits, traversal_bit_length = args
    ddmm = msbits + dm * 2 ** (traversal_bit_length)
    print("msbits", msbits)
    ddmm_length = dm_length + traversal_bit_length
    ddll = dl
    ddll_length = dl_length
    return attack0(N, e, d_length, ddll, ddll_length, ddmm, ddmm_length, factor_e, m, t)

def attack(N, e, d_length, dl, dl_length, dm, dm_length, factor_e=True, m=1, t=None, traversal_bit_length=10, option=None):
    # d_length = partial_d.bit_length
    # dl, dl_length = partial_d.get_known_lsb()
    # dm, dm_length = partial_d.get_known_msb()
    assert dl_length >= 0 or dm_length >= 0, "At least some lsb or msb of d must be known."

    if option == 'LSB':
        print("Traversal extra lsb")
        with Pool(processes=cpu_count()-2) as pool:
            # Prepare arguments for each lsbits
            args = [(N, e, d_length, dl, dl_length, dm, dm_length, factor_e, m, t, lsbits, traversal_bit_length)
                    for lsbits in range(1, 2 ** traversal_bit_length)]
            results = pool.map(attack_single_lsbit, args)
        
        # Filter out None results and return the first valid result
        for res in results:
            if res is not None:
                return res
        return None

    elif option == 'MSB':
        print("Traversal extra msb")
        with Pool(processes=cpu_count()-2) as pool:
            # Prepare arguments for each msbits
            args = [(N, e, d_length, dl, dl_length, dm, dm_length, factor_e, m, t, msbits, traversal_bit_length)
                    for msbits in range(1, 2 ** traversal_bit_length)]
            results = pool.map(attack_single_msbit, args)
        
        # Filter out None results and return the first valid result
        for res in results:
            if res is not None:
                return res
        return None

    elif option is None:
        print("No traversal extra LSB or MSB.")
        return attack0(N, e, d_length, dl, dl_length, dm, dm_length, factor_e, m, t)
    
    else:
        print("You must choose using extra traversal or not! (LSB, MSB or None (default))")
'''

class AttackResultFound(Exception):
    def __init__(self, result):
        self.result = result

def attack_single_lsbit(args):
    N, e, d_length, dl, dl_length, dm, dm_length, factor_e, m, t, lsbits, traversal_bit_length = args
    ddll = lsbits * 2 ** dl_length + dl
    print("lsbits", lsbits)
    ddll_length = dl_length + traversal_bit_length
    ddmm = dm
    ddmm_length = dm_length
    
    result = attack0(N, e, d_length, ddll, ddll_length, ddmm, ddmm_length, factor_e, m, t)
    
    if result is not None:
        raise AttackResultFound(result)

def attack_single_msbit(args):
    N, e, d_length, dl, dl_length, dm, dm_length, factor_e, m, t, msbits, traversal_bit_length = args
    ddmm = msbits + dm * 2 ** (traversal_bit_length)
    print("msbits", msbits)
    ddmm_length = dm_length + traversal_bit_length
    ddll = dl
    ddll_length = dl_length
    
    result = attack0(N, e, d_length, ddll, ddll_length, ddmm, ddmm_length, factor_e, m, t)
    
    if result is not None:
        raise AttackResultFound(result)

def attack(N, e, d_length, dl, dl_length, dm, dm_length, factor_e=True, m=1, t=None, traversal_bit_length=10, option=None):
    # d_length = partial_d.bit_length
    # dl, dl_length = partial_d.get_known_lsb()
    # dm, dm_length = partial_d.get_known_msb()
    assert dl_length >= 0 or dm_length >= 0, "At least some lsb or msb of d must be known."

    if option == 'LSB':
        print("Traversal extra lsb")
        with Pool(processes=cpu_count()-2) as pool:
            # Prepare arguments for each lsbits
            args = [(N, e, d_length, dl, dl_length, dm, dm_length, factor_e, m, t, lsbits, traversal_bit_length)
                    for lsbits in range(1, 2 ** traversal_bit_length)]
            try:
                results = pool.map(attack_single_lsbit, args)
            except AttackResultFound as e:
                return e.result
        
        return None

    elif option == 'MSB':
        print("Traversal extra msb")
        with Pool(processes=cpu_count()-2) as pool:
            # Prepare arguments for each msbits
            args = [(N, e, d_length, dl, dl_length, dm, dm_length, factor_e, m, t, msbits, traversal_bit_length)
                    for msbits in range(1, 2 ** traversal_bit_length)]
            try:
                results = pool.map(attack_single_msbit, args)
            except AttackResultFound as e:
                return e.result

        return None

    elif option is None:
        print("No traversal extra LSB or MSB.")
        return attack0(N, e, d_length, dl, dl_length, dm, dm_length, factor_e, m, t)
    
    else:
        print("You must choose using extra traversal or not! (LSB, MSB or None (default))")


def attack0(N, e, d_bit_length, d0, d0_bit_length, d1, d1_bit_length, factor_e, m, t):
    """
    Recovers the prime factors of a modulus and the private exponent if part of the private exponent is known.
    More information: Boneh D., Durfee G., Frankel Y., "An Attack on RSA Given a Small Fraction of the Private Key Bits"
    More information: Blomer J., May A., "New Partial Key Exposure Attacks on RSA"
    More information: Ernst M. et al., "Partial Key Exposure Attacks on RSA Up to Full Size Exponents"
    :param N: the modulus
    :param e: the public exponent
    :param partial_d: the partial private exponent d (PartialInteger)
    :param factor_e: whether we should attempt to factor e (for BDF) if it is not prime (default: True)
    :param m: the m value to use for the small roots method (default: 1)
    :param t: the t value to use for the small roots method (default: automatically computed using m)
    :return: a tuple containing the prime factors and the private exponent, or None if the private exponent was not found
    """
    # d_bit_length = partial_d.bit_length
    # d0, d0_bit_length = partial_d.get_known_lsb()
    # d1, d1_bit_length = partial_d.get_known_msb()
    # assert d0_bit_length > 0 or d1_bit_length > 0, "At least some lsb or msb of d must be known."

    n = N.bit_length()
    # Subtract one here, because 2^t < e < 2^(t + 1).
    # t_ = e.bit_length() - 1
    # t_ = e.nbits() - 1
    try:
        t_ = e.bit_length() - 1
    except AttributeError:
        t_ = e.nbits() - 1
    alpha = t_ / n
    beta = d_bit_length / n
    kappa = d0_bit_length / n

    flag = False

    if (beta <= 0.25 and d0_bit_length < 1 and d1_bit_length < 1):
        print("Using Wiener's continued fractions approach")
        p_, q_, d_ = Wiener(N, e)
        # print("p_, q_, d_:", p_, q_, d_)
        return (p_, q_, d_)

    if (beta > 0.25 and beta <= 0.27 and d0_bit_length < 1 and d1_bit_length < 1):
        print("Using original Boneh-Durfee attack")
        delta = beta + 0.0001
        p_, q_, d_ = BonehDurfee(N, e, 512, delta=delta, m=m, t=t, roots_method='groebner') # groebner
        # print("p_, q_:", p_, q_)
        return (p_, q_, d_)

    if (beta > 0.27 and beta < 0.293 and d0_bit_length < 1 and d1_bit_length < 1):
        print("Using Boneh-Durfee attack Variant in [LZQ23]")
        parallel_num = multiprocessing.cpu_count() // 2
        delta = beta
        print("delta =", RR(delta))
        if (d_bit_length - 275 < 5):
            s = 5
        else:
            s = d_bit_length - 275
        g = ceil(s/3)
        print("s =", s)
        print("g =", g)
        p_, q_, d_ = small_exponent.BonehDurfee_LZQ23(N, e, delta, m, t, s, g, parallel_num)

        return (p_, q_, d_)
    
    assert beta >= 0.292, "Use Wiener's or the Boneh-Durfee attack if d is very small."
    print("Attacking...Known LSB length is:",d0_bit_length,", Known MSB length is:",d1_bit_length,", d length is:",d_bit_length,", e length is:",t_)
    
    if d0_bit_length > 0 and d1_bit_length > 0:
        # Known lsbs and msbs.
        M = 2 ** d0_bit_length
        if 1 <= t_ <= (n / 2) and d0_bit_length >= n / 4 and d1_bit_length >= t_:
            logging.info("Using Boneh-Durfee-Frankel (Section 4.3)...")
            print("Using Boneh-Durfee-Frankel (Section 4.3)...")

            assert t is not None, "t can not be None for Boneh-Durfee-Frankel small roots."
            return _bdf_4_3(N, e, d_bit_length, d0, d0_bit_length, d1, d1_bit_length, M, m, t)

        delta = (d_bit_length - d0_bit_length - d1_bit_length) / n

        epsilon = 0.05 # Sometimes we need to modify this parameter!

        if (alpha >= 0.75):
            if delta < 5 / 6 - 1 / 3 * sqrt(1 + 6 * beta) - epsilon:
                # When e is full size and d is not.
                logging.info("Using Ernst (Appendix C)...")
                print("Using Ernst (Appendix C)...")
                t = int((1 / 2 - delta) * m) if t is None else t
                # return _ernst_Appendix_C(N, e, d_bit_length, d0, d0_bit_length, d1, d1_bit_length, m, t)
                result = _ernst_Appendix_C(N, e, d_bit_length, d0, d0_bit_length, d1, d1_bit_length, m, t)
                if result is None:
                    logging.error(f"Failed to execute _ernst_Appendix_C: {result}")
                    print("Failed to execute _ernst_Appendix_C")
                else:
                    return result
            if (delta <= (3 - 4 * kappa - 4 * kappa ** 2) / (16 + 16 * kappa) - epsilon and beta <= (11 + 4 * kappa - 4 * kappa ** 2) / (16 + 16 * kappa)) or (delta <= 1 / 3 + 1 / 3 * beta - 1 / 3 * sqrt(4 * beta ** 2 + 2 * beta - 2) - epsilon and beta >= (11 + 4 * kappa - 4 * kappa ** 2) / (16 + 16 * kappa)):
                # When e is full size and d is not, and using approximation of k.
                logging.info("Using Ernst (Appendix C_1)...")
                print("Using Ernst (Appendix C_1)...")
                gamma = max(delta, beta - 1 / 2)
                t = int(((1 / 2 - delta - gamma) / (2 * gamma)) * m) if t is None else t
                # return _ernst_Appendix_C_1(N, e, d_bit_length, d0, d0_bit_length, d1, d1_bit_length, m, t)
                result = _ernst_Appendix_C_1(N, e, d_bit_length, d0, d0_bit_length, d1, d1_bit_length, m, t)
                if result is None:
                    logging.error(f"Failed to execute _ernst_Appendix_C_1: {result}")
                    print("Failed to execute _ernst_Appendix_C_1")
                else:
                    return result

        if (alpha > 1/2 and alpha < 0.75):
            # if (((delta + kappa) >= 1/2 and delta < (3 + 4*(alpha+kappa) - 4*(alpha+kappa)**2)/(16 * (alpha+kappa))) or ((delta + kappa) <= 1/2 and delta <= 1/3 + 1/3 * alpha - 1/3*sqrt(4*alpha** 2 + 2*alpha - 2) - epsilon )):
            if ((delta + kappa) < 1/2 and delta <= 1/3 + 1/3 * alpha - 1/3*sqrt(4*alpha** 2 + 2*alpha - 2) - epsilon ):
                # When d is full size and e is not.
                logging.info("Using Ernst (Appendix C_2)...")
                print("Using Ernst (Appendix C_2)...")
                gamma = max(alpha + delta - 1, alpha - 1 / 2)
                t = int(((1 / 2 - delta - gamma) / (2 * gamma)) * m) if t is None else t
                result = _ernst_Appendix_C_2(N, e, d_bit_length, d0, d0_bit_length, d1, d1_bit_length, m, t)
                if result is None:
                    logging.error(f"Failed to execute _ernst_Appendix_C_2: {result}")
                    print("Failed to execute _ernst_Appendix_C_2")
                else:
                    return result

        logging.info("No attacks were found to fit the provided parameters (known lsbs and msbs).")
        print("No attacks were found to fit the provided parameters (known lsbs and msbs).")
        return None

    if d0_bit_length > 0:
        # Known lsbs.
        M = 2 ** d0_bit_length
        delta = (d_bit_length - d0_bit_length) / n

        if e < RR(N) ** (7 / 8) and RR(N) ** (1 / 6 + 1 / 3 * sqrt(1 + 6 * alpha)) <= M:
            logging.info("Using Blomer-May (Section 6)...")
            print("Using Blomer-May (Section 6)...")

            t = int((2 / 3 * (1 - delta - alpha) / (2 * alpha - 1))) if t is None else t
            
            return _bm_6(N, e, d_bit_length, d0, d0_bit_length, M, m, t)

        if delta <= 5 / 6 - 1 / 3 * sqrt(1 + 6 * beta):
            logging.info("Using Ernst (Section 4.3)...")
            print("Using Ernst (Section 4.3)...")
            
            t = int((1 / 2 - delta) * m) if t is None else t
            
            return _ernst_4_3(N, e, d_bit_length, d0, d0_bit_length, M, m, t)

        # Last resort method: enumerate possible k values (very slow if e is too large).
        if d0_bit_length >= n / 4:
            logging.info("Using Boneh-Durfee-Frankel (Section 3)...")
            print("Using Boneh-Durfee-Frankel (Section 3)...")
            assert t is not None, "t can not be None for Boneh-Durfee-Frankel small roots."
            return _bdf_3(N, e, d_bit_length, d0, d0_bit_length, M, m, t)

        logging.info("No attacks were found to fit the provided parameters (known lsbs).")
        print("No attacks were found to fit the provided parameters (known lsbs).")
        return None

    if d1_bit_length > 0:
        # Known msbs.
        delta = (d_bit_length - d1_bit_length) / n

        if 0 <= t_ <= n / 2 and d1_bit_length >= n - t_:
            logging.info("Using Boneh-Durfee-Frankel (Section 4.2)...")
            print("Using Boneh-Durfee-Frankel (Section 4.2)...")
            return _bdf_4_2(N, e, d_bit_length, d1, d1_bit_length)

        if t_ <= n / 2 and d1_bit_length >= (3 / 4)*n:
            # From B-D-F full version theorem 3.3
            print("Using Boneh-Durfee-Frankel Full version (Theorem 3.3)...")
            t = int(1/2 * m) if t is None else t
            return _bdf_theorem_3_3(N, e, d_bit_length, d1, d1_bit_length, m, t)

        if n / 4 <= t_ <= n / 2 and d1_bit_length >= t_ and (is_prime(e) or factor_e):
            logging.info("Using Boneh-Durfee-Frankel (Section 4.1)...")
            print("Using Boneh-Durfee-Frankel (Section 4.1)...")
            assert t is not None, "t can not be None for Boneh-Durfee-Frankel small roots."
            return _bdf_4_1(N, e, d_bit_length, d1, d1_bit_length, m, t)

         # Blomer-May Section 4 is superseded by Ernst Section 4.2.
        if 1 / 2 < alpha <= (sqrt(6) - 1) / 2 and delta <= 1 / 8 * (5 - 2 * alpha - sqrt(36 * alpha ** 2 + 12 * alpha - 15)):
            logging.info("Using Blomer-May (Section 4)...")
            print("Using Blomer-May (Section 4)...")
            t = int((2 / 3 * (1 - delta - alpha) / (2 * alpha - 1))) if t is None else t
            return _bm_4(N, e, d_bit_length, d1, d1_bit_length, m, t)

        margin4_1_1 = 5 / 6 - 1 / 3 * sqrt(1 + 6 * beta) - delta
        margin4_1_2 = (3 / 16 - delta) if beta <= 11 / 16 else (1 / 3 + 1 / 3 * beta - 1 / 3 * sqrt(4 * beta ** 2 + 2 * beta - 2) - delta)
        if margin4_1_1 > max(0, margin4_1_2):
            logging.info("Using Ernst (Section 4.1.1)...")
            print("Using Ernst (Section 4.1.1)...")
            t = int((1 / 2 - delta) * m) if t is None else t
            # return True
            return _ernst_4_1_1(N, e, d_bit_length, d1, d1_bit_length, m, t)

        if margin4_1_2 > max(0, margin4_1_1):
            logging.info("Using Ernst (Section 4.1.2)...")
            print("Using Ernst (Section 4.1.2)...")
            gamma = max(delta, beta - 1 / 2)
            t = int(((1 / 2 - delta - gamma) / (2 * gamma)) * m) if t is None else t
            return _ernst_4_1_2(N, e, d_bit_length, d1, d1_bit_length, m, t)

        if alpha > 1 / 2 and delta <= 1 / 3 + 1 / 3 * alpha - 1 / 3 * sqrt(4 * alpha ** 2 + 2 * alpha - 2):
            logging.info("Using Ernst (Section 4.2)...")
            print("Using Ernst (Section 4.2)...")
            gamma = max(alpha + delta - 1, alpha - 1 / 2)
            t = int(((1 / 2 - delta - gamma) / (2 * gamma)) * m) if t is None else t
            return _ernst_4_2(N, e, d_bit_length, d1, d1_bit_length, m, t)

        logging.info("No attacks were found to fit the provided parameters (known msbs).")
        print("No attacks were found to fit the provided parameters (known msbs).")
        return None

    logging.info("No attacks were found to fit the provided parameters.")
    return None

