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
        fname = 'out_FMS_K_'+ str(K) + '_N_' + str(N) + '_eta_' + str(eta) + '_alpha_' + str(alpha) + '_idumgraph_' + \
                str(s) + '_tl_' + str(tl) + '_dtN_' + str(dtN) + '.txt'
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
        fout_name = 'FMS_dynamics_K_'+ str(K) + '_N_' + str(N) + '_eta_' + str(eta) + '_alpha_' + str(alpha) + \
                    '_tl_' + str(tl) + '_dtN_' + str(dtN) + '_nsamples_' + str(true_nsamples) + '.txt'
        fout = open(fout_name, "w")
        fout.write("# time  av(final each step)  av(lowest each step)  \n")
        for i in range(len(av_final)):
            fout.write(str(float(i * dtN) / N) + "\t" + str(av_final[i]) + "\t" + str(av_lowest[i]) + "\n")
        fout.close()

        ffinal = 'FMS_final_stats_K_'+ str(K) + '_N_' + str(N) + '_eta_' + str(eta) + '_alpha_' + str(alpha) + \
                 '_tl_' + str(tl) + '_dtN_' + str(dtN) + '_nsamples_' + str(true_nsamples) + '.txt'
        ff = open(ffinal, "w")
        av_t = np.mean(final_times)
        std_t = np.std(final_times)
        av_e = np.mean(final_eners)
        std_e = np.std(final_eners)

        ff.write('# Psol=' + str(float(solved) / true_nsamples) + '\n')
        ff.write('# nsamples=' + str(true_nsamples) + '\n')
        ff.write('# av(tf)=' + str(av_t) + '\n')
        ff.write('# std(tf)=' + str(std_t) + '\n')
        ff.write('# av(ef)=' + str(av_e) + '\n')
        ff.write('# std(ef)=' + str(std_e) + '\n')

        for i in range(len(final_times)):
            ff.write(str(final_eners[i]) + '\t' + str(final_times[i]) + '\n')
        ff.close()
    else:
        print('No data for K=' + str(K) + '  N=' + str(N) + '    eta=' + str(eta) + '    alpha=' + str(alpha) + \
              '    tl=' + str(tl) + '     dtN=' + str(dtN))


def main():
    N = int(sys.argv[1])
    K = int(sys.argv[2])
    eta = sys.argv[3]
    alpha = sys.argv[4]
    tl = int(sys.argv[5])
    dtN = int(sys.argv[6])
    nsamples = int(sys.argv[7])

    average(K, N, eta, alpha, tl, dtN, nsamples)

    return 0


if __name__ == '__main__':
    main()
