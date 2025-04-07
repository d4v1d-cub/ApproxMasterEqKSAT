__author__ = 'david'

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import poisson, multinomial
import sys


def fms_rate_time_ind(u, s0, k, eta):
    if s0 > u:
        return u / k * eta ** (s0 - u)
    else:
        return u / k


def comp_all_rates(max_c, k, rate_func, *args):
    all_rates = []
    for u in range(max_c + 1):
        inner_u = np.zeros(max_c - u + 1)
        for s0 in range(max_c - u + 1):
            inner_u[s0] = rate_func(u, s0, k, *args)
        all_rates.append(inner_u)
    return all_rates


def all_poisson_weights(alpha, k, eps):
    poisson_w = []
    max_c = 0
    while 1 - poisson.cdf(max_c, alpha * k) > eps:
        poisson_w.append(poisson.pmf(max_c, alpha * k))
        max_c += 1
    poisson_w.append(poisson.pmf(max_c, alpha * k))
    poisson_w = np.array(poisson_w)
    return max_c, poisson_w


def init_probs(max_c, k, poisson_w):
    probs = [poisson_w[0]]
    for c in range(1, max_c + 1):
        inner_c = []
        for u in range(c + 1):
            inner_u = np.zeros(c - u + 2)
            for s0 in range(c - u + 1):
                inner_u[s0] = multinomial.pmf([u, s0, c - u - s0], c, [2 ** (-k), 2 ** (-k), 1 - 2 * 2 ** (-k)]) * \
                              poisson_w[c]
            inner_c.append(inner_u)
        inner_c.append(np.zeros(c + 1))
        probs.append(inner_c)
    return probs


def norm(probs, max_c):
    cumul = probs[0]
    for c in range(1, max_c + 1):
        for u in range(c + 1):
            cumul += np.sum(probs[c][u][:-1])
    return cumul


def flatten(probs):
    x = [probs[0]]
    for i in range(1, len(probs)):
        for j in range(len(probs[i])):
            for k in range(len(probs[i][j])):
                x.append(probs[i][j][k])
    return np.array(x)


def unflatten(x, max_c):
    probs = [x[0]]
    counter = 1
    for c in range(1, max_c + 1):
        inner_c = []
        for u in range(c + 1):
            inner_u = np.zeros(c - u + 2)
            for s0 in range(c - u + 2):
                inner_u[s0] = x[counter]
                counter += 1
            inner_c.append(inner_u)
        inner_c.append(np.zeros(c + 1))
        counter += c + 1
        probs.append(inner_c)
    return probs


def comp_av_u(probs, all_rates, max_c, poisson_w):
    cumul = 0
    for gamma in range(max_c):
        cumul_num = 0
        cumul_den = 0
        for u in range(gamma + 1):
            for s0 in range(gamma - u + 1):
                val = (u + 1) * probs[gamma + 1][u + 1][s0]
                cumul_num += all_rates[u + 1][s0] * val
                cumul_den += val
        cumul += poisson_w[gamma] * cumul_num / cumul_den
    return cumul


def comp_av_s0(probs, all_rates, max_c, poisson_w):
    cumul = 0
    for gamma in range(max_c):
        cumul_num = 0
        cumul_den = 0
        for u in range(gamma + 1):
            for s0 in range(gamma - u + 1):
                val = (s0 + 1) * probs[gamma + 1][u][s0 + 1]
                cumul_num += all_rates[u][s0 + 1] * val
                cumul_den += val
        cumul += poisson_w[gamma] * cumul_num / cumul_den
    return cumul


def comp_av_s(probs, all_rates, max_c, poisson_w):
    cumul = 0
    for gamma in range(max_c):
        cumul_num = 0
        cumul_den = 0
        for u in range(gamma + 1):
            for s0 in range(gamma - u + 1):
                val = (gamma - u - s0 + 1) * probs[gamma + 1][u][s0]
                cumul_num += all_rates[u][s0] * val
                cumul_den += val
        cumul += poisson_w[gamma] * cumul_num / cumul_den
    return cumul


def comp_all_avs(probs, all_rates, max_c, poisson_w):
    av_u = comp_av_u(probs, all_rates, max_c, poisson_w)
    av_s0 = comp_av_s0(probs, all_rates, max_c, poisson_w)
    av_s = comp_av_s(probs, all_rates, max_c, poisson_w)
    return av_u, av_s0, av_s


def derivative(probs, all_rates, av_u, av_s0, av_s, k, c, u, s0):
    der = -all_rates[u][s0] * probs[c][u][s0] + all_rates[s0][u] * probs[c][s0][u] - \
          (k - 1) * av_u * (u * probs[c][u][s0] - (u + 1) * probs[c][u + 1][s0]) - \
          (k - 1) * av_s0 / (2 ** k - 2) * ((c - u - s0) * probs[c][u][s0] -
                                            (c - u - s0 + 1) * probs[c][u - 1][s0]) - \
          (k - 1) * av_s * (s0 * probs[c][u][s0] - (s0 + 1) * probs[c][u][s0 + 1] +
                            ((c - u - s0) * probs[c][u][s0] -
                             (c - u - s0 + 1) * probs[c][u][s0 - 1]) / (2 ** k - 2))
    return der


def comp_ders(probs, all_rates, av_u, av_s0, av_s, k, max_c):
    all_ders = [0]
    for c in range(1, max_c + 1):
        inner_c = []
        for u in range(c + 1):
            inner_u = np.zeros(c - u + 2)
            for s0 in range(c - u + 1):
                inner_u[s0] = derivative(probs, all_rates, av_u, av_s0, av_s, k, c, u, s0)
            inner_c.append(inner_u)
        inner_c.append(np.zeros(c + 1))
        all_ders.append(inner_c)
    return all_ders


def energy(probs, k):
    e = 0
    for c in range(1, len(probs)):
        for u in range(c + 1):
            p_u = np.sum(probs[c][u])
            e += p_u * u
    return e / k


def der_full(t, x, all_rates, k, max_c, poisson_w, thr, thr_norm):
    probs = unflatten(x, max_c)
    av_u, av_s0, av_s = comp_all_avs(probs, all_rates, max_c, poisson_w)
    all_ders = comp_ders(probs, all_rates, av_u, av_s0, av_s, k, max_c)
    e = energy(probs, k)
    print(t, norm(probs, max_c), e)
    return flatten(all_ders) / e


def event_stop(t, x, all_rates, k, max_c, poisson_w, thr_e, thr_norm):
    probs = unflatten(x, max_c)
    e = energy(probs, k)
    n = norm(probs, max_c)
    if e < thr_e:
        return 0
    elif n < 1 - thr_norm:
        return 0
    else:
        return 1


event_stop.terminal = True


def solution(tl, x0, all_rates, k, max_c, poisson_w, thr_e, thr_norm, file_probs, file_energy, method, rtol, atol):
    sol = solve_ivp(der_full, [0, tl], x0, method=method, args=(all_rates, k, max_c, poisson_w, thr_e, thr_norm),
                    rtol=rtol, atol=atol, dense_output=True, events=event_stop)
    times = sol.t
    y_vals = sol.y
    fp = open(file_probs, "w")
    fe = open(file_energy, "w")
    eners = np.zeros((len(times)))

    probs = unflatten(x0, max_c)
    eners[0] = energy(probs, k)
    fp.write(str(times[0]))
    for val in x0:
        fp.write("\t" + str(val))
    fp.write("\n")
    fe.write(str(times[0]) + "\t" + str(eners[0]) + "\n")
    for ii in range(1, len(times)):
        x = y_vals[:, ii]
        probs = unflatten(x, max_c)
        eners[ii] = energy(probs, k)
        fp.write(str(times[ii]))
        for val in x:
            fp.write("\t" + str(val))
        fp.write("\n")
        fe.write(str(times[ii]) + "\t" + str(eners[ii]) + "\n")
    fp.close()
    fe.close()
    return times, eners


def several_alphas(tl, alpha, k, eta, eps, thr_e, thr_norm, file_probs, file_energy,
                   method="RK45", rtol=1e-6, atol=1e-12):
    max_c, poisson_w = all_poisson_weights(alpha, k, eps)
    all_rates = comp_all_rates(max_c, k, fms_rate_time_ind, eta)
    probs0 = init_probs(max_c, k, poisson_w)
    x0 = flatten(probs0)
    solution(tl, x0, all_rates, k, max_c, poisson_w, thr_e, thr_norm, file_probs, file_energy, method, rtol, atol)


def main():
    alpha = float(sys.argv[1])
    k = int(sys.argv[2])
    eta = float(sys.argv[3])
    eps = float(sys.argv[4])
    tl = float(sys.argv[5])

    thr_e = float(sys.argv[6])
    thr_norm = float(sys.argv[7])

    head_str = "CDA2_KSAT" + "_K_" + str(k)
    common_str = "alpha_" + str(alpha) + "_eta_" + str(eta) + "_eps_" + str(eps) + \
                 "_tl_" + str(tl) + "_thr_e_" + str(thr_e) + "_thr_n_" + str(thr_norm)

    file_probs = head_str + "_probs_" + common_str + ".txt"
    file_energy = head_str + "_eners_" + common_str + ".txt"

    several_alphas(tl, alpha, k, eta, eps, thr_e, thr_norm, file_probs, file_energy)

    return 0


if __name__ == '__main__':
    main()
