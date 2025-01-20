#include <iostream>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <vector>

using namespace std;

double binomial_coef(int n, int k){
    return gsl_ran_binomial_pdf(k, 1, n);
}


int get_max_gamma(double alpha, int K, double thr){
    int max_gamma = 0;
    double cdf_Q = gsl_cdf_poisson_Q(max_gamma, alpha * K);
    while (cdf_Q > thr){
        max_gamma++;
        cdf_Q = gsl_cdf_poisson_Q(max_gamma, alpha * K);
    }
    return max_gamma;
}


int init_indexes(int **&ind_2_vals, int **&vals_2_ind, int max_gamma){
    int nindex_single = (max_gamma + 1) * (max_gamma + 2) / 2;
    int ind = 0;
    ind_2_vals = new int *[nindex_single];
    vals_2_ind = new int *[max_gamma + 1];
    for (int gamma = 0; gamma < max_gamma + 1; gamma++){
        vals_2_ind[gamma] = new int [gamma + 1];
        for (int lp = 0; lp < gamma + 1; lp++){
            ind_2_vals[ind] = new int [2];
            ind_2_vals[ind][0] = gamma;
            ind_2_vals[ind][1] = lp;

            vals_2_ind[gamma][lp] = ind;

            ind++;
        }
    }

    return nindex_single;
}


// initializes all the joint and conditional probabilities
// prob_joint[index][ch], with ondex that goes over all the possible
// combinations of connectivities and links. If the maximum allowed value for gamma is 
// gamma_max, and there are K variables inside a clause, the number of combinations is 
// (gamma_max + 1)^{K} (gamma_max + 2)^{K} / 2^{K}. One could think
// index as a number with K digits in the base b=(gamma_max + 1) (gamma_max + 2) / 2.
// ch=0,...,2^{K}-1 is the combination of the s_a. ch=2^{K}-1 is the unsat combination 
// pu_cond[gamma][lp][sj], gamma = 0, ..., max_gamma; lp=0,...,gamma; sj=0,1 the state in the conditional
// pi_gamma_lp[sj] only needs two values, and is re-used for every gamma and lp when computing pu_cond
// pjoint_gamma_ln[sj] saves the joint probability of having one variable sj and the rest
// of the clause in the unsat combinatio. Is re-used for every gamma and lp when computing pu_cond  
void init_probs(double **&prob_joint, double ***&pu_cond, double *&pi_gamma_lp, 
                double *&pjoint_gamma_lp, double **&me_sum, int K, int nch_fn, double p0, 
                int max_gamma, long nindex_total, int nindex_single, double alpha, 
                int **ind_2_vals){
    double prod;
    int bit;
    prob_joint = new double *[nindex_total];
    me_sum = new double *[nindex_total];

    int ind_single, gamma_in, lp_in;

    for (long ind = 0; ind < nindex_total; ind++){
        prob_joint[ind] = new double [nch_fn];
        me_sum[ind] = new double [nch_fn];
        for (int ch = 0; ch < nch_fn; ch++){
            prod = 1;
            for (int w = 0; w < K; w++){
                ind_single = (ind / (long) pow(nindex_single, w)) % nindex_single;
                gamma_in = ind_2_vals[ind_single][0];
                lp_in = ind_2_vals[ind_single][1];
                
                bit = ((ch >> w) & 1);
                prod *= (bit + (1 - 2 * bit) * p0) * 
                         gsl_ran_poisson_pdf(gamma_in, alpha * K) * 
                         binomial_coef(gamma_in, lp_in) / pow(2, gamma_in);
            }
            prob_joint[ind][ch] = prod;
        }
    }
    
    pu_cond = new double **[max_gamma + 1];

    for (int gamma = 0; gamma < max_gamma + 1; gamma++){
        pu_cond[gamma] = new double *[gamma + 1];
        for (int lp = 0; lp < gamma + 1; lp++){
            pu_cond[gamma][lp] = new double [2];
        }
    }

    pi_gamma_lp = new double[2];
    pjoint_gamma_lp = new double [2];
}


// initializes all the joint and conditional probabilities
// prob_joint[index][ch], with ondex that goes over all the possible
// combinations of connectivities and links. If the maximum allowed value for gamma is 
// gamma_max, and there are K variables inside a clause, the number of combinations is 
// (gamma_max + 1)^{K} (gamma_max + 2)^{K} / 2^{K}. One could think
// index as a number with K digits in the base b=(gamma_max + 1) (gamma_max + 2) / 2.
// ch=0,...,2^{K}-1 is the combination of the s_a. ch=2^{K}-1 is the unsat combination 
// pu_cond[gamma][lp][sj], gamma = 0, ..., max_gamma; lp=0,...,gamma; sj=0,1 the state in the conditional
// pi_gamma_lp[sj] only needs two values, and is re-used for every gamma and lp when computing pu_cond
// pjoint_gamma_ln[sj] saves the joint probability of having one variable sj and the rest
// of the clause in the unsat combinatio. Is re-used for every gamma and lp when computing pu_cond  
void init_probs(vector < vector <double> > &prob_joint, double ***&pu_cond, double *&pi_gamma_lp, 
                double *&pjoint_gamma_lp, int K, int nch_fn, double p0, 
                int max_gamma, long nindex_total, int nindex_single, double alpha, 
                int **ind_2_vals){
    double prod;
    int bit;

    int ind_single, gamma_in, lp_in;

    for (long ind = 0; ind < nindex_total; ind++){
        for (int ch = 0; ch < nch_fn; ch++){
            prod = 1;
            for (int w = 0; w < K; w++){
                ind_single = (ind / (long) pow(nindex_single, w)) % nindex_single;
                gamma_in = ind_2_vals[ind_single][0];
                lp_in = ind_2_vals[ind_single][1];
                
                bit = ((ch >> w) & 1);
                prod *= (bit + (1 - 2 * bit) * p0) * 
                         gsl_ran_poisson_pdf(gamma_in, alpha * K) * 
                         binomial_coef(gamma_in, lp_in) / pow(2, gamma_in);
            }
            prob_joint[ind][ch] = prod;
        }
    }
    
    pu_cond = new double **[max_gamma + 1];

    for (int gamma = 0; gamma < max_gamma + 1; gamma++){
        pu_cond[gamma] = new double *[gamma + 1];
        for (int lp = 0; lp < gamma + 1; lp++){
            pu_cond[gamma][lp] = new double [2];
        }
    }

    pi_gamma_lp = new double[2];
    pjoint_gamma_lp = new double [2];
}


// initializes the auxiliary arrays for the Runge-Kutta integration
void init_RK_arr(double **&k1, double **&k2, double **&prob_joint_1, 
                 int nch_fn, long nindex_total){
    k1 = new double *[nindex_total];
    k2 = new double *[nindex_total];
    prob_joint_1 = new double *[nindex_total];

    for (long ind = 0; ind < nindex_total; ind++){
        k1[ind] = new double [nch_fn];
        k2[ind] = new double [nch_fn];
        prob_joint_1[ind] = new double [nch_fn];
        for (int ch = 0; ch < nch_fn; ch++){
            k1[ind][ch] = 0;
            k2[ind][ch] = 0;
            prob_joint_1[ind][ch] = 0;
        }
    }
}


// rate of the Focused Metropolis Search algorithm.
double rate_fms(int E0, int E1, int K, double eta){
    double dE = E1 - E0;
    if (dE > 0){
        return double(E0) / K * pow(eta, dE);
    }else{
        return double(E0) / K;
    }
}


void table_all_rates(int max_gamma, int K, double eta, double **&rates){
    rates = new double *[max_gamma + 1];
    for (int gamma = 0; gamma < max_gamma + 1; gamma++){
        rates[gamma] = new double [gamma + 1];
        for (int E0 = 0; E0 < gamma + 1; E0++){
            rates[gamma][E0] = rate_fms(E0, gamma - E0, K, eta);
        }
    }
}


double rate_fms(int E0, int E1, double **rates, double e_av){
    return rates[E0 + E1][E0] / e_av;
}


// it computes the conditional probabilities of having a partially unsatisfied clause, given the 
// value of one variable in the clause
void comp_pcond(double **prob_joint, double ***pu_cond, double *pi_gamma_lp, 
                double *pjoint_gamma_lp, int K, int nch_fn, int max_gamma, 
                long nindex_total, int nindex_single, int **vals_2_ind){
    double pu;
    int bit, ch_uns_flip, ind_single;
    for (int gamma = 0; gamma < max_gamma; gamma++){
        for (int lp = 0; lp < max_gamma; lp++){
            for (int s = 0; s < 2; s++){
                pi_gamma_lp[s] = 0;
                pjoint_gamma_lp[s] = 0;
            }
            ind_single = vals_2_ind[gamma][lp];
            for (long ind = ind_single; ind < nindex_total; ind+=nindex_single){
                for (int s = 0; s < 2; s++){
                    pjoint_gamma_lp[s] += prob_joint[ind][(nch_fn - 1) ^ (1 - s)];  
                    // when s = 0, it inverts the last bit of the unsat combination (nch_fn - 1)
                    // when s = 1, it does not touch the last bit of the unsat combination
                }
                for (int ch = 0; ch < nch_fn; ch++){
                    bit = (ch & 1);
                    pi_gamma_lp[bit] += prob_joint[ind][ch];
                }
            }

            for (int s = 0; s < 2; s++){
                pu_cond[gamma][lp][s] = pjoint_gamma_lp[s] / pi_gamma_lp[s]; 
            }
        }
    }
}


// it computes the conditional probabilities of having a partially unsatisfied clause, given the 
// value of one variable in the clause
void comp_pcond(vector <vector <double> > prob_joint, double ***pu_cond, double *pi_gamma_lp, 
                double *pjoint_gamma_lp, int K, int nch_fn, int max_gamma, 
                long nindex_total, int nindex_single, int **vals_2_ind){
    double pu;
    int bit, ch_uns_flip, ind_single;
    for (int gamma = 0; gamma < max_gamma; gamma++){
        for (int lp = 0; lp < max_gamma; lp++){
            for (int s = 0; s < 2; s++){
                pi_gamma_lp[s] = 0;
                pjoint_gamma_lp[s] = 0;
            }
            ind_single = vals_2_ind[gamma][lp];
            for (long ind = ind_single; ind < nindex_total; ind+=nindex_single){
                for (int s = 0; s < 2; s++){
                    pjoint_gamma_lp[s] += prob_joint[ind][(nch_fn - 1) ^ (1 - s)];  
                    // when s = 0, it inverts the last bit of the unsat combination (nch_fn - 1)
                    // when s = 1, it does not touch the last bit of the unsat combination
                }
                for (int ch = 0; ch < nch_fn; ch++){
                    bit = (ch & 1);
                    pi_gamma_lp[bit] += prob_joint[ind][ch];
                }
            }

            for (int s = 0; s < 2; s++){
                pu_cond[gamma][lp][s] = pjoint_gamma_lp[s] / pi_gamma_lp[s]; 
            }
        }
    }
}


void sum_fms(int K, long ind, int plc_he, double *prob_joint, double ***pu_cond, 
             double **rates, int nch_fn, int max_gamma, double e_av, double *me_sum_src,
             int nindex_single, int **ind_2_vals){

    int bit, ch_flip, uns, uns_flip;

    int ind_single = (ind / (long) pow(nindex_single, plc_he)) % nindex_single;
    int gamma = ind_2_vals[ind_single][0];
    int lp = ind_2_vals[ind_single][1];
    int ln = gamma - lp;

    double sums[2][2]; 
    for (int s = 0; s < 2; s++){
        for (int sat = 0; sat < 2; sat++){
            sums[s][sat] = 0;
        }
    }

    for (int up = 0; up < lp + 1; up++){
        for (int un = 0; un < ln + 1; un++){
            sums[0][0] +=  binomial_coef(lp, up) * binomial_coef(ln, un) * rate_fms(un, up, rates, e_av) * 
                           pow(pu_cond[lp][ln][0], up) * pow(1 - pu_cond[lp][ln][0], lp - up) * 
                           pow(pu_cond[ln - 1][lp + 1][1], un) * pow(1 - pu_cond[ln - 1][lp + 1][1], ln - un);
            sums[1][0] +=  binomial_coef(lp, up) * binomial_coef(ln, un) * rate_fms(up, un, rates, e_av) * 
                           pow(pu_cond[lp][ln][1], up) * pow(1 - pu_cond[lp][ln][1], lp - up) * 
                           pow(pu_cond[ln - 1][lp + 1][0], un) * pow(1 - pu_cond[ln - 1][lp + 1][0], ln - un);               

            sums[0][1] +=  binomial_coef(lp, up) * binomial_coef(ln, un) * rate_fms(un, up + 1, rates, e_av) * 
                           pow(pu_cond[lp][ln][0], up) * pow(1 - pu_cond[lp][ln][0], lp - up) * 
                           pow(pu_cond[ln - 1][lp + 1][1], un) * pow(1 - pu_cond[ln - 1][lp + 1][1], ln - un);
            sums[1][1] +=  binomial_coef(lp, up) * binomial_coef(ln, un) * rate_fms(up + 1, un, rates, e_av) * 
                           pow(pu_cond[lp][ln][1], up) * pow(1 - pu_cond[lp][ln][1], lp - up) * 
                           pow(pu_cond[ln - 1][lp + 1][0], un) * pow(1 - pu_cond[ln - 1][lp + 1][0], ln - un);

        }
    }

    for (int ch_src = 0; ch_src < nch_fn; ch_src++){
        bit = ((ch_src >> plc_he) & 1);
        ch_flip = (ch_src ^ (1 << plc_he));
        uns = (ch_src == nch_fn - 1);
        uns_flip = (ch_flip == nch_fn - 1);
        me_sum_src[ch_src] += -sums[bit][uns || uns_flip] * prob_joint[ch_src] + 
                              sums[1 - bit][uns || uns_flip] * prob_joint[ch_flip];
        // if any of the two, uns and uns_flip, is one, then one has to use the value in
        // sums[2]. One of them represents the probability of a jump when ch_src in unsat,
        // and therefore it goes from E[bit unsat] + 1 ----> E[bit sat]. The other jump makes
        // E[bit sat] ----> E[bit unsat] + 1
    }
    
}


void sum_fms(int K, long ind, int plc_he, vector <double> prob_joint, double ***pu_cond, 
             double **rates, int nch_fn, int max_gamma, double e_av, vector <double> &me_sum_src,
             int nindex_single, int **ind_2_vals){

    int bit, ch_flip, uns, uns_flip;

    int ind_single = (ind / (long) pow(nindex_single, plc_he)) % nindex_single;
    int gamma = ind_2_vals[ind_single][0];
    int lp = ind_2_vals[ind_single][1];
    int ln = gamma - lp;

    double sums[2][2]; 
    for (int s = 0; s < 2; s++){
        for (int sat = 0; sat < 2; sat++){
            sums[s][sat] = 0;
        }
    }

    for (int up = 0; up < lp + 1; up++){
        for (int un = 0; un < ln + 1; un++){
            sums[0][0] +=  binomial_coef(lp, up) * binomial_coef(ln, un) * rate_fms(un, up, rates, e_av) * 
                           pow(pu_cond[lp][ln][0], up) * pow(1 - pu_cond[lp][ln][0], lp - up) * 
                           pow(pu_cond[ln - 1][lp + 1][1], un) * pow(1 - pu_cond[ln - 1][lp + 1][1], ln - un);
            sums[1][0] +=  binomial_coef(lp, up) * binomial_coef(ln, un) * rate_fms(up, un, rates, e_av) * 
                           pow(pu_cond[lp][ln][1], up) * pow(1 - pu_cond[lp][ln][1], lp - up) * 
                           pow(pu_cond[ln - 1][lp + 1][0], un) * pow(1 - pu_cond[ln - 1][lp + 1][0], ln - un);               

            sums[0][1] +=  binomial_coef(lp, up) * binomial_coef(ln, un) * rate_fms(un, up + 1, rates, e_av) * 
                           pow(pu_cond[lp][ln][0], up) * pow(1 - pu_cond[lp][ln][0], lp - up) * 
                           pow(pu_cond[ln - 1][lp + 1][1], un) * pow(1 - pu_cond[ln - 1][lp + 1][1], ln - un);
            sums[1][1] +=  binomial_coef(lp, up) * binomial_coef(ln, un) * rate_fms(up + 1, un, rates, e_av) * 
                           pow(pu_cond[lp][ln][1], up) * pow(1 - pu_cond[lp][ln][1], lp - up) * 
                           pow(pu_cond[ln - 1][lp + 1][0], un) * pow(1 - pu_cond[ln - 1][lp + 1][0], ln - un);

        }
    }

    for (int ch_src = 0; ch_src < nch_fn; ch_src++){
        bit = ((ch_src >> plc_he) & 1);
        ch_flip = (ch_src ^ (1 << plc_he));
        uns = (ch_src == nch_fn - 1);
        uns_flip = (ch_flip == nch_fn - 1);
        me_sum_src[ch_src] += -sums[bit][uns || uns_flip] * prob_joint[ch_src] + 
                              sums[1 - bit][uns || uns_flip] * prob_joint[ch_flip];
        // if any of the two, uns and uns_flip, is one, then one has to use the value in
        // sums[2]. One of them represents the probability of a jump when ch_src in unsat,
        // and therefore it goes from E[bit unsat] + 1 ----> E[bit sat]. The other jump makes
        // E[bit sat] ----> E[bit unsat] + 1
    }
    
}


// it computes all the derivatives of the joint probabilities
void der_fms(double **prob_joint, double ***pu_cond, double **rates, int K, int l_max, 
             int nch_fn, long nindex_total, int nindex_single, double e_av, double **me_sum,
             int **ind_2_vals){
    for (long ind = 0; ind < nindex_total; ind++){
        for (int ch = 0; ch < nch_fn; ch++){
            me_sum[ind][ch] = 0;
        }
    }

    // candidate to be a parallel for
    #pragma omp parallel for
    for (long ind = 0; ind < nindex_total; ind++){
        for (int w = 0; w < K; w++){
            sum_fms(K, ind, w, prob_joint[ind], pu_cond, rates, nch_fn, l_max, e_av, 
                    me_sum[ind], nindex_single, ind_2_vals);
        }
    }
}


// it computes all the derivatives of the joint probabilities
void der_fms(vector <vector <double> > prob_joint, double ***pu_cond, double **rates, int K, int l_max, 
             int nch_fn, long nindex_total, int nindex_single, double e_av, vector <vector <double> > &me_sum,
             int **ind_2_vals){
    for (long ind = 0; ind < nindex_total; ind++){
        for (int ch = 0; ch < nch_fn; ch++){
            me_sum[ind][ch] = 0;
        }
    }

    // candidate to be a parallel for
    #pragma omp parallel for
    for (long ind = 0; ind < nindex_total; ind++){
        for (int w = 0; w < K; w++){
            sum_fms(K, ind, w, prob_joint[ind], pu_cond, rates, nch_fn, l_max, e_av, 
                    me_sum[ind], nindex_single, ind_2_vals);
        }
    }
}


double energy(double **prob_joint, long nindex_total, int nch_fn){
    double e = 0;
    for (long ind = 0; ind < nindex_total; ind++){
        e += prob_joint[ind][nch_fn - 1];
    }
    return e;
}


double energy(vector <vector <double> > prob_joint, long nindex_total, int nch_fn){
    double e = 0;
    for (long ind = 0; ind < nindex_total; ind++){
        e += prob_joint[ind][nch_fn - 1];
    }
    return e;
}


// peforms the integration of the differential equations with the 2nd order Runge-Kutta
// the method is implemented with adaptive step size
void RK2_fms(double alpha, int K, int nch_fn, double eta, int max_gamma, 
             double p0, char *fileener, double tl, double tol = 1e-2, double t0 = 0, double dt0 = 0.01, 
             double ef = 1e-6, double dt_min = 1e-7){
    double **rates;
    double **prob_joint, ***pu_cond, **me_sum, *pi_gamma_lp, *pjoint_gamma_lp;
    double e, error;                 
    int **ind_2_vals, **vals_2_ind;

    table_all_rates(max_gamma, K, eta, rates);
    
    int nindex_single = init_indexes(ind_2_vals, vals_2_ind, max_gamma);
    long nindex_total = (long) pow(nindex_single, K);

    init_probs(prob_joint, pu_cond, pi_gamma_lp, pjoint_gamma_lp, me_sum, K, 
               nch_fn, p0, max_gamma, nindex_total, nindex_single, alpha, ind_2_vals);

    // initialize auxiliary arrays for the Runge-Kutta integration
    double **k1, **k2, **prob_joint_1;
    init_RK_arr(k1, k2, prob_joint_1, nch_fn, nindex_total);

    ofstream fe(fileener);
    
    e = energy(prob_joint, nindex_total, nch_fn) * alpha;
    fe << t0 << "\t" << e << endl;   // it prints the energy density

    double dt1 = dt0;
    double t = t0;

    bool valid;

    // the time scale is already given in Monte Carlo steps. Inside the rates I am using 
    // the energy density e_av
    while (t < tl){
        if (e < ef){
            //  cout << "Final energy reached" << endl;
            break;
        }

        auto t1 = std::chrono::high_resolution_clock::now();

        comp_pcond(prob_joint, pu_cond, pi_gamma_lp, pjoint_gamma_lp, K, nch_fn, max_gamma, 
                   nindex_total, nindex_single, vals_2_ind);

        der_fms(prob_joint, pu_cond, rates, K, max_gamma, nch_fn, nindex_total, nindex_single, 
                 e, me_sum, ind_2_vals);   // in the rates, I use the energy density

        valid = true;
        for (int ind = 0; ind < nindex_total; ind++){
            for (int ch = 0; ch < nch_fn; ch++){
                k1[ind][ch] = dt1 * me_sum[ind][ch];
                prob_joint_1[ind][ch] = prob_joint[ind][ch] + k1[ind][ch];
                if (prob_joint_1[ind][ch] < 0){
                    valid = false;
                }
            }
        }

        while (!valid){
            //  cout << "joint probabilities became negative in the auxiliary step of RK2" << endl;
            dt1 /= 2;
            //  cout << "step divided by half    dt=" << dt1 << endl;
            if (dt1 < dt_min){
                dt_min /= 2;
                //  cout << "dt_min also halfed" << endl;
            }

            valid = true;
            for (int ind = 0; ind < nindex_total; ind++){
                for (int ch = 0; ch < nch_fn; ch++){
                    k1[ind][ch] = dt1 * me_sum[ind][ch];
                    prob_joint_1[ind][ch] = prob_joint[ind][ch] + k1[ind][ch];
                    if (prob_joint_1[ind][ch] < 0){
                        valid = false;
                    }
                }
            }
        }
        
        e = energy(prob_joint_1, nindex_total, nch_fn) * alpha;
        comp_pcond(prob_joint_1, pu_cond, pi_gamma_lp, pjoint_gamma_lp, K, nch_fn, max_gamma, 
                   nindex_total, nindex_single, vals_2_ind);

        der_fms(prob_joint_1, pu_cond, rates, K, max_gamma, nch_fn, nindex_total, nindex_single, e, 
                me_sum, ind_2_vals);
            
        valid = true;
        for (int ind = 0; ind < nindex_total; ind++){
            for (int ch = 0; ch < nch_fn; ch++){
                k2[ind][ch] = dt1 * me_sum[ind][ch];
                if (prob_joint[ind][ch] + (k1[ind][ch] + k2[ind][ch]) / 2 < 0){
                    valid = false;
                }
            }
        }
        
        if (!valid){
            //  cout << "Some probabilities would be negative if dt=" << dt1 << " is taken" << endl;
            dt1 /= 2;
            //  cout << "step divided by half    dt=" << dt1  << endl;
            if (dt1 < dt_min){
                dt_min /= 2;
                //  cout << "dt_min also halfed" << endl;
            }
            e = energy(prob_joint, nindex_total, nch_fn) * alpha;
        }else{
            error = 0;
            for (int ind = 0; ind < nindex_total; ind++){
                for (int ch = 0; ch < nch_fn; ch++){
                    error += fabs(k1[ind][ch] - k2[ind][ch]);
                }
            }

            error /= nch_fn * nindex_total;

            if (error < 2 * tol){
                //  cout << "step dt=" << dt1 << "  accepted" << endl;
                //  cout << "error=" << error << endl;
                t += dt1;
                for (int ind = 0; ind < nindex_total; ind++){
                    for (int ch = 0; ch < nch_fn; ch++){
                        prob_joint[ind][ch] += (k1[ind][ch] + k2[ind][ch]) / 2;
                    }
                }
                e = energy(prob_joint, nindex_total, nch_fn) * alpha;
                fe << t << "\t" << e << endl;

            }else{
                e = energy(prob_joint, nindex_total, nch_fn) * alpha;
                //  cout << "step dt=" << dt1 << "  rejected  new step will be attempted" << endl;
                //  cout << "error=" <<  error << endl;
            }

            dt1 = 4 * dt1 * sqrt(2 * tol / error) / 5;
            if(dt1 < dt_min){
                dt1 = dt_min;
            }

            //  cout << "Recommended step is dt=" << dt1 << endl;
        }

        auto t2 = std::chrono::high_resolution_clock::now();

        auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

        cout << endl << "iteration time:   " << ms_int.count() << "ms" << endl; 
    }

    fe.close();

}


// peforms the integration of the differential equations with the 2nd order Runge-Kutta
// the method is implemented with adaptive step size
void RK2_fms_vectors(double alpha, int K, int nch_fn, double eta, int max_gamma, 
                     double p0, char *fileener, double tl, double tol = 1e-2, double t0 = 0, double dt0 = 0.01, 
                     double ef = 1e-6, double dt_min = 1e-7){
    double **rates;
    double ***pu_cond, *pi_gamma_lp, *pjoint_gamma_lp;
    double e, error;                 
    int **ind_2_vals, **vals_2_ind;

    table_all_rates(max_gamma, K, eta, rates);
    
    init_indexes(ind_2_vals, vals_2_ind, max_gamma);

    vector < vector <double> > prob_joint(nindex_total, vector <double> (nch_fn));
    vector < vector <double> > me_sum(nindex_total, vector <double> (nch_fn));

    init_probs(prob_joint, pu_cond, pi_gamma_lp, pjoint_gamma_lp, K, 
               nch_fn, p0, max_gamma, nindex_total, nindex_single, alpha, ind_2_vals);

    // initialize auxiliary arrays for the Runge-Kutta integration
    vector < vector <double> > k1(nindex_total, vector <double> (nch_fn));
    vector < vector <double> > k2(nindex_total, vector <double> (nch_fn));
    vector < vector <double> > prob_joint_1(nindex_total, vector <double> (nch_fn));

    
    ofstream fe(fileener);
    
    e = energy(prob_joint, nindex_total, nch_fn) * alpha;
    fe << t0 << "\t" << e << endl;   // it prints the energy density

    double dt1 = dt0;
    double t = t0;

    bool valid;

    // the time scale is already given in Monte Carlo steps. Inside the rates I am using 
    // the energy density e_av
    while (t < tl){
        if (e < ef){
            //  cout << "Final energy reached" << endl;
            break;
        }

        auto t1 = std::chrono::high_resolution_clock::now();

        comp_pcond(prob_joint, pu_cond, pi_gamma_lp, pjoint_gamma_lp, K, nch_fn, max_gamma, 
                   nindex_total, nindex_single, vals_2_ind);

        der_fms(prob_joint, pu_cond, rates, K, max_gamma, nch_fn, nindex_total, nindex_single, 
                 e, me_sum, ind_2_vals);   // in the rates, I use the energy density

        valid = true;
        for (int ind = 0; ind < nindex_total; ind++){
            for (int ch = 0; ch < nch_fn; ch++){
                k1[ind][ch] = dt1 * me_sum[ind][ch];
                prob_joint_1[ind][ch] = prob_joint[ind][ch] + k1[ind][ch];
                if (prob_joint_1[ind][ch] < 0){
                    valid = false;
                }
            }
        }

        while (!valid){
            //  cout << "joint probabilities became negative in the auxiliary step of RK2" << endl;
            dt1 /= 2;
            //  cout << "step divided by half    dt=" << dt1 << endl;
            if (dt1 < dt_min){
                dt_min /= 2;
                //  cout << "dt_min also halfed" << endl;
            }

            valid = true;
            for (int ind = 0; ind < nindex_total; ind++){
                for (int ch = 0; ch < nch_fn; ch++){
                    k1[ind][ch] = dt1 * me_sum[ind][ch];
                    prob_joint_1[ind][ch] = prob_joint[ind][ch] + k1[ind][ch];
                    if (prob_joint_1[ind][ch] < 0){
                        valid = false;
                    }
                }
            }
        }
        
        e = energy(prob_joint_1, nindex_total, nch_fn) * alpha;
        comp_pcond(prob_joint_1, pu_cond, pi_gamma_lp, pjoint_gamma_lp, K, nch_fn, max_gamma, 
                   nindex_total, nindex_single, vals_2_ind);

        der_fms(prob_joint_1, pu_cond, rates, K, max_gamma, nch_fn, nindex_total, nindex_single, e, 
                me_sum, ind_2_vals);
            
        valid = true;
        for (int ind = 0; ind < nindex_total; ind++){
            for (int ch = 0; ch < nch_fn; ch++){
                k2[ind][ch] = dt1 * me_sum[ind][ch];
                if (prob_joint[ind][ch] + (k1[ind][ch] + k2[ind][ch]) / 2 < 0){
                    valid = false;
                }
            }
        }
        
        if (!valid){
            //  cout << "Some probabilities would be negative if dt=" << dt1 << " is taken" << endl;
            dt1 /= 2;
            //  cout << "step divided by half    dt=" << dt1  << endl;
            if (dt1 < dt_min){
                dt_min /= 2;
                //  cout << "dt_min also halfed" << endl;
            }
            e = energy(prob_joint, nindex_total, nch_fn) * alpha;
        }else{
            error = 0;
            for (int ind = 0; ind < nindex_total; ind++){
                for (int ch = 0; ch < nch_fn; ch++){
                    error += fabs(k1[ind][ch] - k2[ind][ch]);
                }
            }

            error /= nch_fn * nindex_total;

            if (error < 2 * tol){
                //  cout << "step dt=" << dt1 << "  accepted" << endl;
                //  cout << "error=" << error << endl;
                t += dt1;
                for (int ind = 0; ind < nindex_total; ind++){
                    for (int ch = 0; ch < nch_fn; ch++){
                        prob_joint[ind][ch] += (k1[ind][ch] + k2[ind][ch]) / 2;
                    }
                }
                e = energy(prob_joint, nindex_total, nch_fn) * alpha;
                fe << t << "\t" << e << endl;

            }else{
                e = energy(prob_joint, nindex_total, nch_fn) * alpha;
                //  cout << "step dt=" << dt1 << "  rejected  new step will be attempted" << endl;
                //  cout << "error=" <<  error << endl;
            }

            dt1 = 4 * dt1 * sqrt(2 * tol / error) / 5;
            if(dt1 < dt_min){
                dt1 = dt_min;
            }

            //  cout << "Recommended step is dt=" << dt1 << endl;
        }

        auto t2 = std::chrono::high_resolution_clock::now();

        auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

        cout << endl << "iteration time:   " << ms_int.count() << "ms" << endl; 
    }

    fe.close();

}


int main(int argc, char *argv[]) {
    double alpha = atof(argv[1]);
    int K = atoi(argv[2]);
    double eta = atof(argv[3]);
    double tl = atof(argv[4]);
    double tol = atof(argv[5]);
    int nthr = atoi(argv[6]);
    double eps_c = atof(argv[7]);

    int nch_fn = (1 << K);
    double p0 = 0.5;

    omp_set_num_threads(nthr);

    char fileener[300]; 
    sprintf(fileener, "CDA1av_lpln_FMS_ener_K_%d_alpha_%.4lf_eta_%.4lf_tl_%.2lf_tol_%.1e_epsc_%.e.txt", 
            K, alpha, eta, tl, tol, eps_c);


    int max_gamma = get_max_gamma(alpha, K, eps_c);
    
    // RK2_fms(alpha, K, nch_fn, eta, max_gamma, p0, fileener, tl, tol);
    RK2_fms_vectors(alpha, K, nch_fn, eta, max_gamma, p0, fileener, tl, tol);

    return 0;
}