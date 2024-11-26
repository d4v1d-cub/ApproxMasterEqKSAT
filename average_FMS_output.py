__author__ = 'david'

import numpy as np
import sys

def average(K, N, eta, alpha, tl, dtN, nsamples):
    av_final = np.zeros(tl * N // dtN + 1)
    av_lowest = np.zeros(tl * N // dtN + 1)
    final_times = []
    final_eners = []
    true_nsamples = 0
    solved = 0
    for s in range(1, nsamples + 1):
        fname = f'out_FMS_K_{K}_N_{N}_eta_{eta}_alpha_{alpha}_idumgraph_{s}_tl_{tl}_dtN_{dtN}.txt'
        try:
            fin = open(fname, "r")
            for _ in range(15):
                fin.readline()
            counter = 0
            line = fin.readline()
            while line != '\n':
                line = line.split()
                if int(line[0]) > 0:
                    av_lowest[counter] += float(line[0])
                    av_final[counter] += float(line[1])
                    counter += 1
                if len(line) > 3:
                    t_final = float(line[5])
                    ef = int(line[1])
                    if ef == 0:
                        solved += 1
                line = fin.readline()
            true_nsamples += 1
            final_times.append(t_final)
            final_eners.append(ef)
        except(OSError, IOError):
            print(fname + "  not read")
            continue

    if true_nsamples > 0:
        av_lowest /= true_nsamples
        av_final /= true_nsamples
        fout_name = f'FMS_dynamics_K_{K}_N_{N}_eta_{eta}_alpha_{alpha}_tl_{tl}_dtN_{dtN}_nsamples_{true_nsamples}.txt'
        fout = open(fout_name, "w")
        fout.write("# time  av(final each step)  av(lowest each step)  \n")
        for i in range(len(av_final)):
            fout.write(str(float(i * dtN) / N) + "\t" + str(av_final[i]) + "\t" + str(av_lowest[i]) + "\n")
        fout.close()

        ffinal = f'FMS_final_stats_K_{K}_N_{N}_eta_{eta}_alpha_{alpha}_tl_{tl}_dtN_{dtN}_nsamples_{true_nsamples}.txt'
        ff = open(ffinal, "w")
        av_t = np.mean(final_times)
        std_t = np.std(final_times)
        av_e = np.mean(final_eners)
        std_e = np.std(final_eners)

        ff.write(f'# Psol={float(solved) / true_nsamples}\n')
        ff.write(f'# nsamples={true_nsamples}\n')
        ff.write(f'# av(tf)={av_t}\n')
        ff.write(f'# std(tf)={std_t}\n')
        ff.write(f'# av(ef)={av_e}\n')
        ff.write(f'# std(ef)={std_e}\n')

        for i in range(len(final_times)):
            ff.write(f'{final_eners[i]}\t{final_times[i]}\n')
        ff.close()
    else:
        print(f'No data for K={K}  N={N}  eta={eta}  alpha={alpha}  tl={tl}   dtN={dtN}')


def main():
    N = int(sys.argv[1])
    K = int(sys.argv[2])
    eta = sys.argv[3]
    alpha = sys.argv[4]
    tl = sys.argv[5]
    dtN = sys.argv[6]
    nsamples = sys.argv[7]

    average(K, N, eta, alpha, tl, dtN, nsamples)

    return 0


if __name__ == '__main__':
    main()
