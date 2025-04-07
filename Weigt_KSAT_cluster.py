__author__ = 'david'

import numpy as np
from scipy.integrate import ode
from scipy.special import binom
from scipy.stats import poisson
import sys


observables = []


def delta(a, b):
    return int(a == b)


def ener(p_u_full, mean_c):
    e = 0
    for c in range(len(p_u_full)):
        for u in range(len(p_u_full[c])):
            e += u * p_u_full[c][u] * poisson.pmf(c, mean_c) / k
    return e


def rates_fms_time_ind(u0, s0, eta=1, n=1):
    delta_e = s0 - u0
    if delta_e > 0:
        return u0 / k / n * np.power(eta, delta_e)
    else:
        return u0 / k / n


def rates_mc(u0, s0, temp=1):
    delta_e = s0 - u0
    if delta_e > 0:
        return np.exp(-delta_e / temp)
    else:
        return 1


def flatten_p_u_full(p_u_full):
    x = []
    for c in range(len(p_u_full)):
        for u in range(c + 1):
            x.append(p_u_full[c][u])
    return x


def unflatten_p_u_full(x):
    p_u_full = []
    counter = 0
    for c in range(c_cutoff + 1):
        inner1 = []
        for u in range(c + 1):
            inner1.append(x[counter])
            counter += 1
        p_u_full.append(inner1)
    return p_u_full


def time_independent_mean(c, u, ratefunc=rates_fms_time_ind, ratefunc_kwargs=()):
    cumul = 0.0
    for n in range(c - u + 1):
        cumul += binom(c - u, n) * ((2 ** k - 2) ** (-n)) * ratefunc(u, n, **dict(ratefunc_kwargs))
    return cumul * ((1 - 1.0 / (2 ** k - 1)) ** (c - u))


def store_independent_means(ratefunc=rates_fms_time_ind, ratefunc_kwargs=()):
    ind_means = []
    for c in range(c_cutoff + 1):
        inner1 = []
        for u in range(c + 1):
            inner1.append(time_independent_mean(c, u, ratefunc=ratefunc, ratefunc_kwargs=ratefunc_kwargs))
        ind_means.append(inner1)
    return ind_means


def compute_all_bin_factors(ratefunc=rates_fms_time_ind, ratefunc_kwargs=()):
    all_rates = []
    for c in range(c_cutoff + 1):
        inner1 = []
        for u in range(c + 1):
            inner2 = []
            for n in range(c - u + 1):
                inner2.append((2 ** k - 2) ** (-u) * binom(c - n, u) * ((1 - 1.0 / (2 ** k - 1)) ** (c - n)) *
                              ratefunc(n, u, **dict(ratefunc_kwargs)))
            inner1.append(inner2)
        all_rates.append(inner1)
    return all_rates


def mean_u(p_u_full, c):
    cumul = 0.0
    for u in range(c + 1):
        cumul += u * p_u_full[c][u]
    return cumul


def mean_s(p_u_full, c):
    cumul = 0.0
    for u in range(c + 1):
        cumul += (c - u) * p_u_full[c][u]
    return cumul


def mean_rate_u(ind_means, p_u_full, c):
    cumul = 0.0
    for u in range(c + 1):
        cumul += ind_means[c][u] * u * p_u_full[c][u]
    return cumul


def mean_rate_s(ind_means, p_u_full, c):
    cumul = 0.0
    for u in range(c + 1):
        cumul += ind_means[c][u] * (c - u) * p_u_full[c][u]
    return cumul


def compute_needed_means(p_u_full, ind_means, mean_c):
    needed_means_rate_u = 0.0
    needed_means_u = 0.0
    for c in range(len(p_u_full)):
        needed_means_rate_u += mean_rate_u(ind_means, p_u_full, c) * poisson.pmf(c, mean_c)
        needed_means_u += mean_u(p_u_full, c) * poisson.pmf(c, mean_c)

    needed_means_rate_s = 0.0
    needed_means_s = 0.0
    for c in range(len(p_u_full)):
        needed_means_rate_s += mean_rate_s(ind_means, p_u_full, c) * poisson.pmf(c, mean_c)
        needed_means_s += mean_s(p_u_full, c) * poisson.pmf(c, mean_c)
    return needed_means_rate_u / needed_means_u, needed_means_rate_s / needed_means_s


def mean_over_n(c, u, p, bin_factors):
    cumul = 0.0
    for n in range(c - u + 1):
        cumul += bin_factors[c][u][n] * p[n]
    return cumul


def p_u_full_der(p_u_full, c, ind_means, bin_factors, needed_means_rates, e=1):
    u = 0
    der_u = [-ind_means[c][u] * p_u_full[c][u] / e + mean_over_n(c, u, p_u_full[c], bin_factors) / e -
             (k - 1) * needed_means_rates[0] * (-(u + 1) * p_u_full[c][u + 1]) / e -
             (k - 1) / (2 ** k - 1) * needed_means_rates[1] * ((c - u) * p_u_full[c][u]) / e]

    u += 1
    while u < c:

        der_u.append(-ind_means[c][u] * p_u_full[c][u] / e + mean_over_n(c, u, p_u_full[c], bin_factors) / e -
                     (k - 1) * needed_means_rates[0] *
                     (u * p_u_full[c][u] - (u + 1) * p_u_full[c][u + 1]) / e -
                     (k - 1) / (2 ** k - 1) * needed_means_rates[1] *
                     ((c - u) * p_u_full[c][u] - (c - u + 1) * p_u_full[c][u - 1]) / e)
        u += 1

    der_u.append(-ind_means[c][u] * p_u_full[c][u] / e + mean_over_n(c, u, p_u_full[c], bin_factors) / e -
                 (k - 1) * needed_means_rates[0] * (u * p_u_full[c][u]) / e -
                 (k - 1) / (2 ** k - 1) * needed_means_rates[1] *
                 (-(c - u + 1) * p_u_full[c][u - 1]) / e)

    return der_u


def p_u_full_der_c0(p_u_full, ind_means, bin_factors, e=0):
    c = 0
    u = 0
    return [-ind_means[c][u] * p_u_full[c][u] / e + mean_over_n(c, u, p_u_full[c], bin_factors) / e]


def all_p_u_full_der(p_u_full, ind_means, bin_factors, needed_means_rates, e=1):
    all_ders = [p_u_full_der_c0(p_u_full, ind_means, bin_factors, e=e)]
    for c in range(1, c_cutoff + 1):
        all_ders.append(p_u_full_der(p_u_full, c, ind_means, bin_factors, needed_means_rates, e=e))
    return all_ders


def build_ratefunc_kwargs(param, algorithm='FMS'):
    if algorithm == "FMS":
        return {"eta": param}
    elif algorithm == "MC":
        return {"temp": param}
    else:
        print("You must choose one of the available options for the algorithm")


def build_ders_kwargs(e, algorithm='FMS'):
    if algorithm == "FMS":
        return {"e": e}
    elif algorithm == "MC":
        return {"e": 1}
    else:
        print("You must choose one of the available options for the algorithm")


def derivative_full(t, x, mean_c, filew, ind_means, bin_factors, algorithm):
    p_u_full = unflatten_p_u_full(x)

    e = ener(p_u_full, alpha * k)
    print(t, e)

    if e < threshold:
        print("Threshold reached")
        return np.zeros(len(x))
    else:
        observables.append([t, e])
        filew.write(str(t) + "\t" + str(e) + "\n")
        ders_kwargs = build_ders_kwargs(e, algorithm)
        needed_means_rates = compute_needed_means(p_u_full, ind_means, mean_c)
        ders_p = all_p_u_full_der(p_u_full, ind_means, bin_factors, needed_means_rates, **ders_kwargs)

        return flatten_p_u_full(ders_p)


def solout(t, x):
    p_u_full = unflatten_p_u_full(x)
    e = ener(p_u_full, alpha * k)
    if e < threshold:
        print(t, e)
        observables.append([t, e])
        print("Threshold reached")
        return -1
    else:
        print(t, e)
        observables.append([t, e])
        return 0


def solution(t_limit, x, mean_c, method, param, filename, rtol, atol, ratefunc=rates_fms_time_ind, algorithm="FMS"):
    filew = open(filename, "w")
    ratefunc_kwargs = build_ratefunc_kwargs(param, algorithm=algorithm)
    ind_means = store_independent_means(ratefunc=ratefunc, ratefunc_kwargs=ratefunc_kwargs)
    bin_factors = compute_all_bin_factors(ratefunc=ratefunc, ratefunc_kwargs=ratefunc_kwargs)
    solver = ode(derivative_full).set_integrator(method, nsteps=10 ** 10, rtol=rtol, atol=atol)
    solver.set_f_params(mean_c, filew, ind_means, bin_factors, algorithm)
    # solver.set_solout(solout)
    solver.set_initial_value(x)

    global observables
    observables = []

    solver.integrate(t_limit)

    filew.close()


def print_solution(filename):
    w = open(filename, "w")
    for i in range(len(observables)):
        w.write(str(observables[i][0]))
        for j in range(1, len(observables[i])):
            w.write("\t" + str(observables[i][j]))
        w.write("\n")
    w.close()


def transition(x, method, t0, t_limit, eps, rtol, atol, param=0.0, ratefunc=rates_fms_time_ind, paramname="eta",
               algorithm="FMS"):
    mean_c = alpha * k
    print("\n", paramname, "=", param, "\t", "alpha = ", alpha, "\t", "cutoff = ", c_cutoff)
    fileout = "Weigt_KSAT_" + algorithm + "_alpha_" + str(alpha) + "_" + paramname + "_" + str(param) + \
              "_thr_" + str(threshold) + "_time_" + str(t_limit-t0) + "_eps_" + str(eps) + "_rtol_" + str(rtol) + \
              "_atol_" + str(atol) + ".txt"
    solution(t_limit, x, mean_c, method, param, fileout, rtol, atol, ratefunc=ratefunc, algorithm=algorithm)


def set_initial_p_u_full_rand():
    p_u_full = []
    lenght = 0
    for c in range(c_cutoff + 1):
        inner1 = []
        for u in range(c + 1):
            inner1.append(binom(c, u) * (2 ** (-k * u)) * ((1 - 2 ** (-k)) ** (c - u)))
            lenght += 1
        p_u_full.append(inner1)
    return p_u_full, lenght


def get_cutoff(mean_c, eps):
    cutoff = 0
    while poisson.cdf(cutoff, mean_c) < 1 - eps:
        cutoff += 1
    return cutoff - 1


def main():
    eta = float(sys.argv[1])
    # temp = 1.0
    global alpha
    alpha = float(sys.argv[2])
    global k
    k = 3
    eps = np.power(10.0, -6)
    global c_cutoff
    c_cutoff = get_cutoff(alpha * k, eps)
    t0 = 0.0
    t_limit = float(sys.argv[3])
    method = 'vode'
    global threshold
    threshold = 10 ** (-6)
    rtol = 1e-06
    atol = 1e-12

    p_u_full_0, len_p = set_initial_p_u_full_rand()
    x = flatten_p_u_full(p_u_full_0)
    transition(x, method, t0, t_limit, eps, rtol, atol, param=eta, ratefunc=rates_fms_time_ind, paramname="eta",
               algorithm="FMS")
    # transition(x, method, t0, t_limit, eps, rtol, atol, param=temp, ratefunc=rates_cav_mc, paramname="T",
    #            algorithm="MC")

    return 0


if __name__ == '__main__':
    main()
