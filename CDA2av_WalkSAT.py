__author__ = 'david'

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import poisson, multinomial
import sys


def get_all_poisson_sums(max_c, pu_av, mean_c):
    poisson_probs = np.zeros(max_c + 1)
    poisson_sums = np.zeros(max_c + 1)
    for sj in range(max_c + 1):
        poisson_probs[sj] = poisson.pmf(sj, (1 - pu_av) * mean_c)
    
    for si in range(max_c + 1):
        poisson_sums[si] = 1 - poisson.cdf(si, (1 - pu_av) * mean_c)
    
    return poisson_probs, poisson_sums


def walksat_rate(u, K, q, nch_exc, poisson_prob_S, poisson_sum_S):
    if u > 0:
        cumul = 0
        pneigh = [poisson_prob_S, poisson_sum_S]
        for ch in range(nch_exc):
            prod = 1
            cumul_bits = 0
            for w in range(K - 1):
                bit = ((ch >> w) & 1)
                prod *= pneigh[bit]  
                cumul_bits += 1 - bit
            cumul += prod / (cumul_bits + 1);
        return u * (q / K + (1 - q) * cumul); 
    else:
        return 0


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


def comp_av_u(probs, max_c, poisson_w, K, q, nch_exc, poisson_probs, poisson_sums):
    cumul = 0
    for gamma in range(max_c):
        cumul_num = 0
        cumul_den = 0
        for u in range(gamma + 1):
            for s0 in range(gamma - u + 1):
                val = (u + 1) * probs[gamma + 1][u + 1][s0]
                cumul_num += walksat_rate(u + 1, K, q, nch_exc, poisson_probs[gamma - u], poisson_sums[gamma - u]) * val
                cumul_den += val
        cumul += poisson_w[gamma] * cumul_num / cumul_den
    return cumul


def comp_av_s0(probs, max_c, poisson_w, K, q, nch_exc, poisson_probs, poisson_sums):
    cumul = 0
    for gamma in range(max_c):
        cumul_num = 0
        cumul_den = 0
        for u in range(gamma + 1):
            for s0 in range(gamma - u + 1):
                val = (s0 + 1) * probs[gamma + 1][u][s0 + 1]
                cumul_num += walksat_rate(u, K, q, nch_exc, poisson_probs[gamma + 1 - u], poisson_sums[gamma + 1 - u]) * val
                cumul_den += val
        cumul += poisson_w[gamma] * cumul_num / cumul_den
    return cumul


def comp_av_s(probs, max_c, poisson_w, K, q, nch_exc, poisson_probs, poisson_sums):
    cumul = 0
    for gamma in range(max_c):
        cumul_num = 0
        cumul_den = 0
        for u in range(gamma + 1):
            for s0 in range(gamma - u + 1):
                val = (gamma - u - s0 + 1) * probs[gamma + 1][u][s0]
                cumul_num += walksat_rate(u, K, q, nch_exc, poisson_probs[gamma + 1 - u], poisson_sums[gamma + 1 - u]) * val
                cumul_den += val
        cumul += poisson_w[gamma] * cumul_num / cumul_den
    return cumul


def comp_all_avs(probs, max_c, poisson_w, K, q, nch_exc, poisson_probs, poisson_sums):
    av_u = comp_av_u(probs, max_c, poisson_w, K, q, nch_exc, poisson_probs, poisson_sums)
    av_s0 = comp_av_s0(probs, max_c, poisson_w, K, q, nch_exc, poisson_probs, poisson_sums)
    av_s = comp_av_s(probs, max_c, poisson_w, K, q, nch_exc, poisson_probs, poisson_sums)
    return av_u, av_s0, av_s


def derivative(probs, av_u, av_s0, av_s, K, c, u, s0, q, nch_exc, poisson_probs, poisson_sums):
    der = -walksat_rate(u, K, q, nch_exc, poisson_probs[c - u], poisson_sums[c - u]) * probs[c][u][s0] + \
          walksat_rate(s0, K, q, nch_exc, poisson_probs[c - s0], poisson_sums[c - s0]) * probs[c][s0][u] - \
          (K - 1) * av_u * (u * probs[c][u][s0] - (u + 1) * probs[c][u + 1][s0]) - \
          (K - 1) * av_s0 / (2 ** K - 2) * ((c - u - s0) * probs[c][u][s0] -
                                            (c - u - s0 + 1) * probs[c][u - 1][s0]) - \
          (K - 1) * av_s * (s0 * probs[c][u][s0] - (s0 + 1) * probs[c][u][s0 + 1] +
                            ((c - u - s0) * probs[c][u][s0] -
                             (c - u - s0 + 1) * probs[c][u][s0 - 1]) / (2 ** K - 2))
    return der


def comp_ders(probs, av_u, av_s0, av_s, K, max_c, q, nch_exc, poisson_probs, poisson_sums):
    all_ders = [0]
    for c in range(1, max_c + 1):
        inner_c = []
        for u in range(c + 1):
            inner_u = np.zeros(c - u + 2)
            for s0 in range(c - u + 1):
                inner_u[s0] = derivative(probs, av_u, av_s0, av_s, K, c, u, s0, q, nch_exc, poisson_probs, poisson_sums)
            inner_c.append(inner_u)
        inner_c.append(np.zeros(c + 1))
        all_ders.append(inner_c)
    return all_ders


def energy(probs, K):
    e = 0
    for c in range(1, len(probs)):
        for u in range(c + 1):
            p_u = np.sum(probs[c][u])
            e += p_u * u
    return e / K



def der_full(t, x, K, max_c, poisson_w, nch_exc, q, alpha, thr_e):
    probs = unflatten(x, max_c)
    e = energy(probs, K)
    pu_av = e / alpha
    poisson_probs, poisson_sums = get_all_poisson_sums(max_c, pu_av, alpha * K)

    av_u, av_s0, av_s = comp_all_avs(probs, max_c, poisson_w, K, q, nch_exc, poisson_probs, poisson_sums)
    all_ders = comp_ders(probs, av_u, av_s0, av_s, K, max_c, q, nch_exc, poisson_probs, poisson_sums)
    
    print(t, e)
    sys.stdout.flush()
    return flatten(all_ders) / e


def event_stop(t, x, K, max_c, poisson_w, nch_exc, q, alpha, thr_e):
    probs = unflatten(x, max_c)
    e = energy(probs, K)
    if e < thr_e:
        return 0
    else:
        return 1


event_stop.terminal = True


def solution(tl, x0, K, max_c, poisson_w, q, alpha, thr_e, file_probs, file_energy, method, rtol, atol):
    nch_exc = 2 ** (K - 1)
    sol = solve_ivp(der_full, [0, tl], x0, method=method, args=(K, max_c, poisson_w, nch_exc, q, alpha, thr_e),
                    rtol=rtol, atol=atol, dense_output=True, events=event_stop)
    times = sol.t
    y_vals = sol.y
    fe = open(file_energy, "w")

    probs = unflatten(x0, max_c)
    e = energy(probs, K)
    
    fe.write(str(times[0]) + "\t" + str(e) + "\n")
    for ii in range(1, len(times)):
        x = y_vals[:, ii]
        probs = unflatten(x, max_c)
        e = energy(probs, K)
        fe.write(str(times[ii]) + "\t" + str(e) + "\n")
    fe.close()

    fp = open(file_probs, "w")
    for val in x:
        fp.write("\t" + str(val))
    fp.write("\n")
    fp.close()


def several_alphas(tl, alpha, K, q, eps, thr_e, file_probs, file_energy,
                   method="RK23", rtol=1e-6, atol=1e-12):
    max_c, poisson_w = all_poisson_weights(alpha, K, eps)
    probs0 = init_probs(max_c, K, poisson_w)
    x0 = flatten(probs0)
    solution(tl, x0, K, max_c, poisson_w, q, alpha, thr_e, file_probs, file_energy, method, rtol, atol)


def main():
    alpha = float(sys.argv[1])
    K = int(sys.argv[2])
    q = float(sys.argv[3])
    eps = float(sys.argv[4])
    tl = float(sys.argv[5])

    thr_e = float(sys.argv[6])

    head_str = "CDA2av_WalkSAT" + "_K_" + str(K)
    common_str = "alpha_" + str(alpha) + "_q_" + str(q) + "_eps_" + str(eps) + \
                 "_tl_" + str(tl) + "_thr_e_" + str(thr_e)

    file_probs = head_str + "_probs_" + common_str + ".txt"
    file_energy = head_str + "_eners_" + common_str + ".txt"

    several_alphas(tl, alpha, K, q, eps, thr_e, file_probs, file_energy)

    return 0


if __name__ == '__main__':
    main()
