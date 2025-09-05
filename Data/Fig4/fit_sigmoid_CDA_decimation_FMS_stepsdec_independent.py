__author__ = 'david'

import numpy as np
from scipy.optimize import curve_fit

def read_all(K, stepsdec, path):
    fin_name = f'{path}/Results/CDA_decimation_FMS_all_K_{K}_stepsdec_{stepsdec}.txt'
    fin = open(fin_name, "r")
    fin.readline()
    n_list = []
    data = {}
    min_alpha = 1000
    max_alpha = 0
    while True:
        j = fin.readline()
        if not j:
            break
        line = j.split()
        alpha = float(line[1])
        if alpha < min_alpha:
            min_alpha = alpha
        if alpha > max_alpha:
            max_alpha = alpha
        if int(line[2]) == stepsdec:
            n = int(line[0])
            if n not in n_list:
                n_list.append(n)
                data[n] = []
                data[n].append([alpha, float(line[5])])
            else:
                data[n].append([alpha, float(line[5])])

    fin.close()
    for n in n_list:
        data[n] = np.array(data[n])
    return n_list, data, min_alpha, max_alpha


def sigmoid(x, a, b):
    y = 1.0 / (1 + np.exp(a * x + b))
    return y


def make_fit(data, n, func_to_fit):
    answ = curve_fit(func_to_fit, data[n][:, 0], data[n][:, 1], 
                     bounds=((0, -np.inf), (np.inf, np.inf)), ftol=1e-10, xtol=1e-10, maxfev=100000)
    params = answ[0]
    error_params = np.sqrt(np.diag(answ[1]))    
    return params, error_params


def all_fits(n_list, data, func_to_fit):
    fits = {}
    for n in n_list:
        fits[n] = make_fit(data, n, func_to_fit)
    return fits


def print_params(fits, fileout, path):
    fout = open(f'{path}/{fileout}', "w")
    fout.write("#  N  a  err(a)  b  err(b)\n")
    for n in fits.keys():
        fout.write(str(n) + "\t")
        for i in range(len(fits[n][0])):
            fout.write(str(fits[n][0][i]) + "\t" + str(fits[n][1][i]) + "\n")
    fout.close()


def gen_data_for_plot(stepsdec, fits, n_list, func_to_fit, min_alpha, max_alpha, 
                      margin, npoints, path):
    alpha_list = np.linspace(min_alpha - margin, max_alpha + margin, npoints)
    for n in n_list:
        foutname = f'{path}/Results/sigmoid_data_N_{n}_stepsdec_{stepsdec}.txt'
        fout = open(foutname, "w")
        for alpha in alpha_list:
            estimate = func_to_fit(alpha, fits[n][0][0], fits[n][0][1])
            fout.write(f'{alpha}\t{estimate}\n')
        fout.close()

def main():
    K = 3
    stepsdec = 10
    margin = 0.5
    npoints = 1000

    path = "/media/david/Data/UH/Grupo_de_investigacion/Estudiantes/Jonathan/Decimation/CDA"

    n_list, data, min_alpha, max_alpha = read_all(K, stepsdec, path)
    fits = all_fits(n_list, data, sigmoid)
    print_params(fits, f'./Results/sigmoid_params_stepsdec_{stepsdec}.txt', path)
    gen_data_for_plot(stepsdec, fits, n_list, sigmoid, min_alpha, max_alpha, margin, npoints, path)

    return 0


if __name__ == '__main__':
    main()
