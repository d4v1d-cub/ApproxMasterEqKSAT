#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <cmath>
#include <omp.h>
#include <chrono>

using namespace std;



int get_max_l(double alpha, int K, double thr){
    int l_max = 0;
    double cdf_Q = gsl_cdf_poisson_Q(l_max, alpha * K / 2);
    while (cdf_Q > thr){
        l_max++;
        cdf_Q = gsl_cdf_poisson_Q(l_max, alpha * K / 2);
    }
    return l_max;
}


// initializes all the joint and conditional probabilities
// prob_joint[in_lp][in_ln][ch], with in_lp and index that goes over all the possible
// combinations of l+ links. If the maximum allowed value is l_max, and there are
// K variables inside a clause, the number of combinations is (l_max)^{K}. One could think
// in_lp as a number with K digits in the base l_max.
// the index in_ln is the analogous for the links l-
// ch=0,...,2^{K}-1 is the combination of the s_a. ch=2^{K}-1 is the unsat combination 
// pu_cond[lp][ln][sj], lp = 0, ..., l_max; ln=0,...,l_max; sj=0,1 the state in the conditional
// pi[sj] only needs two values, and is re-used for every lp and ln in pu_cond
void init_probs(double ***&prob_joint, double ***&pu_cond, double *&pi_lp_ln, 
                double *&pjoint_lp_ln, double ***&me_sum, int K, int nch_fn, double p0, 
                int l_max, long nindex_l, double alpha){
    double prod;
    int bit;
    prob_joint = new double **[nindex_l];
    me_sum = new double **[nindex_l];

    int lp, ln;
    
    for (long i_p = 0; i_p < nindex_l; i_p++){
        prob_joint[i_p] = new double *[nindex_l];
        me_sum[i_p] = new double *[nindex_l];
        for (long i_n = 0; i_n < nindex_l; i_n++){
            prob_joint[i_p][i_n] = new double [nch_fn];
            me_sum[i_p][i_n] = new double [nch_fn];
            
            for (int ch = 0; ch < nch_fn; ch++){
                prod = 1;
                for (int w = 0; w < K; w++){
                    lp = (i_p / (int) pow(l_max, w)) % l_max;
                    ln = (i_n / (int) pow(l_max, w)) % l_max;

                    bit = ((ch >> w) & 1);
                    prod *= (bit + (1 - 2 * bit) * p0) * 
                            gsl_ran_poisson_pdf(lp, alpha * K / 2) * gsl_ran_poisson_pdf(ln, alpha * K / 2);
                }
                prob_joint[i_p][i_n][ch] = prod;
            }
        }
    }
    
    pu_cond = new double **[l_max];

    for (int lp = 0; lp < l_max; lp++){
        pu_cond[lp] = new double *[l_max];
        for (int ln = 0; ln < l_max; ln++){
            pu_cond[lp][ln] = new double [2];
        }
    }

    pi_lp_ln = new double[2];
    pjoint_lp_ln = new double [2];
}


// initializes the auxiliary arrays for the Runge-Kutta integration
void init_RK_arr(double ***&k1, double ***&k2, double ***&prob_joint_1, 
                int nch_fn, long nindex_l){
    k1 = new double **[nindex_l];
    k2 = new double **[nindex_l];
    prob_joint_1 = new double **[nindex_l];

    for (long i_p = 0; i_p < nindex_l; i_p++){
        k1[i_p] = new double *[nindex_l];
        k2[i_p] = new double *[nindex_l];
        prob_joint_1[i_p] = new double *[nindex_l];
        for (long i_n = 0; i_n < nindex_l; i_n++){
            k1[i_p][i_n] = new double [nch_fn];
            k2[i_p][i_n] = new double [nch_fn];
            prob_joint_1[i_p][i_n] = new double [nch_fn];
            for (int ch = 0; ch < nch_fn; ch++){
                k1[i_p][i_n][ch] = 0;
                k2[i_p][i_n][ch] = 0;
                prob_joint_1[i_p][i_n][ch] = 0;
            }
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


void table_all_rates(int l_max, int K, double eta, double **&rates){
    rates = new double *[l_max + 1];
    for (int E0 = 0; E0 < l_max + 1; E0++){
        rates[E0] = new double [l_max + 1];
        for (int E1 = 0; E1 < l_max + 1; E1++){
            rates[E0][E1] = rate_fms(E0, E1, K, eta);
        }
    }
}


double rate_fms(int E0, int E1, double **rates, double e_av){
    return rates[E0][E1] / e_av;
}


// it computes the conditional probabilities of having a partially unsatisfied clause, given the 
// value of one variable in the clause
void comp_pcond(double ***prob_joint, double ***pu_cond, double *pi_lp_ln, 
                double *pjoint_lp_ln, int K, int nch_fn, int l_max, long nindex_l){
    double pu;
    int bit;
    int ch_uns_flip;
    for (int lp = 0; lp < l_max; lp++){
        for (int ln = 0; ln < l_max; ln++){
            for (int s = 0; s < 2; s++){
                pi_lp_ln[s] = 0;
                pjoint_lp_ln[s] = 0;
            }

            for (long i_p = lp; i_p < nindex_l; i_p+=l_max){
                for (long i_n = ln; i_n < nindex_l; i_n+=l_max){
                    for (int s = 0; s < 2; s++){
                        pjoint_lp_ln[s] += prob_joint[i_p][i_n][(nch_fn - 1) ^ (1 - s)];  
                        // when s = 0, it inverts the last bit of the unsat combination (nch_fn - 1)
                        // when s = 1, it does not touch the last bit of the unsat combination
                    }
                    for (int ch = 0; ch < nch_fn; ch++){
                        bit = (ch & 1);
                        pi_lp_ln[bit] += prob_joint[i_p][i_n][ch];
                    }
                }
            }

            for (int s = 0; s < 2; s++){
                pu_cond[lp][ln][s] = pjoint_lp_ln[s] / pi_lp_ln[s]; 
            }
        }
    }
}


double binomial_coef(int n, int k){
    return gsl_ran_binomial_pdf(k, 1, n);
}


void sum_fms(int K, long i_p, long i_n, int plc_he, double *prob_joint, double ***pu_cond, 
             double **rates, int nch_fn, int l_max, double e_av, double *me_sum_src){

    int bit, ch_flip, uns, uns_flip;

    int lp = (i_p / (int) pow(l_max, plc_he)) % l_max;
    int ln = (i_n / (int) pow(l_max, plc_he)) % l_max;

    double sums[2][2]; 
    for (int s = 0; s < 2; s++){
        for (int sat = 0; sat < 2; sat++){
            sums[s][sat] = 0;
        }
    }

    for (int up = 0; up < lp; up++){
        for (int un = 0; un < ln; un++){
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

    delete [] sums;
    
}


// it computes all the derivatives of the joint probabilities
void der_fms(double ***prob_joint, double ***pu_cond, double **rates, int K, int l_max, 
             int nch_fn, long nindex_l, double e_av, double ***me_sum){
    for (long ip = 0; ip < nindex_l; ip++){
        for (long in = 0; in < nindex_l; in++){
            for (int ch = 0; ch < nch_fn; ch++){
                me_sum[ip][in][ch] = 0;
            }
        }
    }

    // candidate to be a parallel for
    #pragma omp parallel for
    for (long ip = 0; ip < nindex_l; ip++){
        for (long in = 0; in < nindex_l; in++){
            for (int w = 0; w < K; w++){
                sum_fms(K, ip, in, w, prob_joint[ip][in], pu_cond, rates, nch_fn, l_max, e_av, 
                        me_sum[ip][in]);
            }
        }
    }
}


double energy(double ***prob_joint, long nindex_l, int nch_fn){
    double e = 0;
    for (long ip = 0; ip < nindex_l; ip++){
        for (long in = 0; in < nindex_l; in++){
            e += prob_joint[ip][in][nch_fn - 1];
        }
    }
    return e;
}


// peforms the integration of the differential equations with the 2nd order Runge-Kutta
// the method is implemented with adaptive step size
void RK2_fms(double alpha, int K, int nch_fn, double eta, int l_max, 
             double p0, char *fileener, double tl, double tol = 1e-2, double t0 = 0, double dt0 = 0.01, 
             double ef = 1e-6, double dt_min = 1e-7){
    double **rates;
    double ***prob_joint, ***pu_cond, ***me_sum, *pi_lp_ln, *pjoint_lp_ln;
    double e, error;                 
    
    long nindex_l = (long) pow(l_max, K);
    
    table_all_rates(l_max, K, eta, rates);
    
    init_probs(prob_joint, pu_cond, pi_lp_ln, pjoint_lp_ln, me_sum, K, 
               nch_fn, p0, l_max, nindex_l, alpha);

    // initialize auxiliary arrays for the Runge-Kutta integration
    double ***k1, ***k2, ***prob_joint_1;
    init_RK_arr(k1, k2, prob_joint_1, nch_fn, nindex_l);

    ofstream fe(fileener);
    
    e = energy(prob_joint, nindex_l, nch_fn) * alpha;
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

        comp_pcond(prob_joint, pu_cond, pi_lp_ln, pjoint_lp_ln, K, nch_fn, l_max, nindex_l);

        der_fms(prob_joint, pu_cond, rates, K, l_max, nch_fn, nindex_l, e, me_sum);   // in the rates, I use the energy density

        valid = true;
        for (long ip = 0; ip < nindex_l; ip++){
            for (long in = 0; in < nindex_l; in++){
                for (int ch = 0; ch < nch_fn; ch++){
                    k1[ip][in][ch] = dt1 * me_sum[ip][in][ch];
                    prob_joint_1[ip][in][ch] = prob_joint[ip][in][ch] + k1[ip][in][ch];
                    if (prob_joint_1[ip][in][ch] < 0){
                        valid = false;
                    }
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
            for (long ip = 0; ip < nindex_l; ip++){
                for (long in = 0; in < nindex_l; in++){
                    for (int ch = 0; ch < nch_fn; ch++){
                        k1[ip][in][ch] = dt1 * me_sum[ip][in][ch];
                        prob_joint_1[ip][in][ch] = prob_joint[ip][in][ch] + k1[ip][in][ch];
                        if (prob_joint_1[ip][in][ch] < 0){
                            valid = false;
                        }
                    }
                }
            }
        }
        
        e = energy(prob_joint_1, nindex_l, nch_fn) * alpha;
        comp_pcond(prob_joint_1, pu_cond, pi_lp_ln, pjoint_lp_ln, K, nch_fn, l_max, nindex_l);

        der_fms(prob_joint_1, pu_cond, rates, K, l_max, nch_fn, nindex_l, e, me_sum);
            
        valid = true;
        for (long ip = 0; ip < nindex_l; ip++){
            for (long in = 0; in < nindex_l; in++){
                for (int ch = 0; ch < nch_fn; ch++){
                    k2[ip][in][ch] = dt1 * me_sum[ip][in][ch];
                    if (prob_joint[ip][in][ch] + (k1[ip][in][ch] + k2[ip][in][ch]) / 2 < 0){
                        valid = false;
                    }
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
            e = energy(prob_joint, nindex_l, nch_fn) * alpha;
        }else{
            error = 0;
            for (int ch_u = 0; ch_u < nch_fn; ch_u++){
                for (int ch = 0; ch < nch_fn; ch++){
                    error += fabs(k1[ch_u][ch] - k2[ch_u][ch]);
                }
            }

            error /= nch_fn * nch_fn;

            if (error < 2 * tol){
                //  cout << "step dt=" << dt1 << "  accepted" << endl;
                //  cout << "error=" << error << endl;
                t += dt1;
                for (int ch_u = 0; ch_u < nch_fn; ch_u++){
                    for (int ch = 0; ch < nch_fn; ch++){
                        prob_joint[ch_u][ch] += (k1[ch_u][ch] + k2[ch_u][ch]) / 2;
                    }
                }
                e = energy(prob_joint, nch_fn) * alpha;
                fe << t << "\t" << e << endl;

            }else{
                e = energy(prob_joint, nch_fn) * alpha;
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
    long nsamples = atol(argv[1]);
    double alpha = atof(argv[2]);
    int K = atoi(argv[3]);
    unsigned long seed_r = atol(argv[4]);
    double eta = atof(argv[5]);
    double tl = atof(argv[6]);
    double tol = atof(argv[7]);
    int nthr = atoi(argv[8]);
    double eps_c = atof(argv[9]);

    int nch_fn = (1 << K);
    double p0 = 0.5;

    omp_set_num_threads(nthr);

    gsl_rng * r;
    init_ran(r, seed_r);

    char fileener[300]; 
    sprintf(fileener, "CDA1av_FMS_ener_K_%d_alpha_%.4lf_eta_%.4lf_tl_%.2lf_seed_%li_tol_%.1e_nsamples_%li_epsc_%.e.txt", 
            K, alpha, eta, tl, seed_r, tol, nsamples, eps_c);


    int max_gamma = get_max_gamma(alpha, K, eps_c);
    
    RK2_fms(alpha, K, nch_fn, eta, max_gamma, nsamples, p0, fileener, tl, r, tol);

    return 0;
}