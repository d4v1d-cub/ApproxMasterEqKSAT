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
                double *&pjoint_lp_ln, double ***&me_sum, int nch_fn, double p0, 
                int l_max, long nindex_l){
    double prod;
    int bit;
    prob_joint = new double **[nindex_l];
    me_sum = new double **[nindex_l];
    
    for (long i_p = 0; i_p < nindex_l; i_p++){
        prob_joint[i_p] = new double *[nindex_l];
        me_sum[i_p] = new double *[nindex_l];
        for (long i_n = 0; i_n < nindex_l; i_n++){
            prob_joint[i_p][i_n] = new double [nch_fn];
            me_sum[i_p][i_n] = new double [nch_fn];
            for (int ch = 0; ch < nch_fn; ch++){
                prod = 1;
                for (int w = 0; w < K; w++){
                    bit = ((ch >> w) & 1);
                    prod *= (bit + (1 - 2 * bit) * p0);
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


void sum_fms(int K, long i_p, long i_n, int w, double *prob_joint, double ***pu_cond, 
             double **rates, int nch_fn, int l_max, double e_av, double *me_sum_src){

    int lp = (i_p / (int) pow(l_max, w)) % l_max;
    int ln = (i_n / (int) pow(l_max, w)) % l_max;


    long he;
    int ch_flip;
    bool bit, uns, uns_flip;

    for (E[0] = 0; E[0] < ln + 1; E[0]++){
        for (E[1] = 0; E[1] < lp + 1; E[1]++){

            terms[0][0] = rate_fms(E[0], E[1], rates, e_av) * fE[0][0][E[0]] * fE[1][0][E[1]];
            terms[1][0] = rate_fms(E[1], E[0], rates, e_av) * fE[0][1][E[0]] * fE[1][1][E[1]];

            bit = ((ch_u >> plc_he) & 1);

            terms[bit][1] = rate_fms(E[bit] + 1, E[1 - bit], rates, e_av) * 
                            fE[bit][bit][E[bit]] * fE[1 - bit][bit][E[1 - bit]];
            terms[1 - bit][1] = rate_fms(E[1 - bit], E[bit] + 1, rates, e_av) * 
                                fE[1 - bit][1 - bit][E[1 - bit]] * fE[bit][1 - bit][E[bit]];
            
            for (int ch_src = 0; ch_src < nch_fn; ch_src++){
                bit = ((ch_src >> plc_he) & 1);
                ch_flip = (ch_src ^ (1 << plc_he));
                uns = (ch_src == ch_u);
                uns_flip = (ch_flip == ch_u);
                me_sum_src[ch_src] += -terms[bit][uns || uns_flip] * prob_joint[ch_src] + 
                                      terms[1 - bit][uns || uns_flip] * prob_joint[ch_flip];
                // if any of the two, uns and uns_flip, is one, then one has to use the terms
                // in terms[1]. One of them represents the probability of a jump when ch_src in unsat,
                // and therefore it goes from E[bit unsat] + 1 ----> E[bit sat]. The other jump makes
                // E[bit sat] ----> E[bit unsat] + 1
            }

        }
    }
    delete_aux_arr(pu_l, fE, fEnew);

}


void comp_sums(double alpha, int K, int max_gamma, long nsamples, int ch_u, int plc_he, 
               double *prob_joint, double ***pu_cond, double **rates, int nch_fn, double e_av, 
               double *me_sum_src, gsl_rng * r){
    int gamma;
    for (long i = 0; i < nsamples; i++){
        gamma = gsl_ran_poisson(r, alpha * K);
        if (gamma > max_gamma){
            gamma = max_gamma;
        }       
        sum_fms(K, gamma, ch_u, plc_he, prob_joint, pu_cond, rates, nch_fn, e_av, me_sum_src, r);
    }
}


void comp_sums(double alpha, int K, int max_gamma, int ch_u, int plc_he, 
               double *prob_joint, double ***pu_cond, double **rates, int nch_fn, double e_av, 
               double *me_sum_src, gsl_rng * r){
    double *me_sum_gamma;
    me_sum_gamma = new double[nch_fn];
    for (int gamma = 0; gamma < max_gamma + 1; gamma++){
        for (int ch = 0; ch < nch_fn; ch++){
            me_sum_gamma[ch] = 0;
        }
        sum_fms(K, gamma, ch_u, plc_he, prob_joint, pu_cond, rates, nch_fn, e_av, me_sum_gamma, r);

        for (int ch = 0; ch < nch_fn; ch++){
            me_sum_src[ch] += me_sum_gamma[ch] * gsl_ran_poisson_pdf(gamma, alpha * K);
        }
    }
}


// it computes all the derivatives of the joint probabilities
void der_fms(double **prob_joint, double ***pu_cond, double **rates, long nsamples, 
             double alpha, int K, int max_gamma, int nch_fn, double e_av, double **me_sum, 
             gsl_rng * r){
    for (int ch_u = 0; ch_u < nch_fn; ch_u++){
        for (int ch = 0; ch < nch_fn; ch++){
            me_sum[ch_u][ch] = 0;
        }
    }

    // candidate to be a parallel for
    #pragma omp parallel for
    for (int ch_u = 0; ch_u < nch_fn; ch_u++){
        for (int w = 0; w < K; w++){
            comp_sums(alpha, K, max_gamma, nsamples, ch_u, w, prob_joint[ch_u], pu_cond, 
                      rates, nch_fn, e_av, me_sum[ch_u], r);
        }
        for (int ch = 0; ch < nch_fn; ch++){
            me_sum[ch_u][ch] /= nsamples;
        }
    }
}


// it computes all the derivatives of the joint probabilities
void der_fms(double **prob_joint, double ***pu_cond, double **rates, 
             double alpha, int K, int max_gamma, int nch_fn, double e_av, double **me_sum, 
             gsl_rng * r){
    for (int ch_u = 0; ch_u < nch_fn; ch_u++){
        for (int ch = 0; ch < nch_fn; ch++){
            me_sum[ch_u][ch] = 0;
        }
    }

    // candidate to be a parallel for
    #pragma omp parallel for
    for (int ch_u = 0; ch_u < nch_fn; ch_u++){
        for (int w = 0; w < K; w++){
            comp_sums(alpha, K, max_gamma, ch_u, w, prob_joint[ch_u], pu_cond, 
                      rates, nch_fn, e_av, me_sum[ch_u], r);
        }
    }
}


double energy(double **prob_joint, int nch_fn){
    double e = 0;
    for (int ch_u = 0; ch_u < nch_fn; ch_u++){
        e += prob_joint[ch_u][ch_u];
    }
    return e / nch_fn;
}


// peforms the integration of the differential equations with the 2nd order Runge-Kutta
// the method is implemented with adaptive step size
void RK2_fms(double alpha, int K, int nch_fn, double eta, int max_gamma, long nsamples, 
             double p0, char *fileener, double tl, gsl_rng * r, double tol = 1e-2, double t0 = 0, double dt0 = 0.01, 
             double ef = 1e-6, double dt_min = 1e-7){
    double **rates;
    double **prob_joint, ***pu_cond, **me_sum, **pi;
    double e, error;                 
    
    
    table_all_rates(max_gamma, K, eta, rates);
    
    init_probs(prob_joint, pu_cond, pi, me_sum, K, nch_fn, p0);

    // initialize auxiliary arrays for the Runge-Kutta integration
    double **k1, **k2, **prob_joint_1;
    init_RK_arr(k1, k2, prob_joint_1, nch_fn);

    ofstream fe(fileener);
    
    e = energy(prob_joint, nch_fn) * alpha;
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

        comp_pcond(prob_joint, pu_cond, pi, K, nch_fn);

        der_fms(prob_joint, pu_cond, rates, nsamples, alpha, K, max_gamma, nch_fn, e, me_sum, r);   // in the rates, I use the energy density

        valid = true;
        for (int ch_u = 0; ch_u < nch_fn; ch_u++){
            for (int ch = 0; ch < nch_fn; ch++){
                k1[ch_u][ch] = dt1 * me_sum[ch_u][ch];
                prob_joint_1[ch_u][ch] = prob_joint[ch_u][ch] + k1[ch_u][ch];
                if (prob_joint_1[ch_u][ch] < 0){
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
            for (int ch_u = 0; ch_u < nch_fn; ch_u++){
                for (int ch = 0; ch < nch_fn; ch++){
                    k1[ch_u][ch] = dt1 * me_sum[ch_u][ch];
                    prob_joint_1[ch_u][ch] = prob_joint[ch_u][ch] + k1[ch_u][ch];
                    if (prob_joint_1[ch_u][ch] < 0){
                        valid = false;
                    }
                }
            }
        }
        
        e = energy(prob_joint_1, nch_fn) * alpha;
        comp_pcond(prob_joint_1, pu_cond, pi, K, nch_fn);

        der_fms(prob_joint_1, pu_cond, rates, nsamples, alpha, K, max_gamma, nch_fn, e, me_sum, r);
            
        valid = true;
        for (int ch_u = 0; ch_u < nch_fn; ch_u++){
            for (int ch = 0; ch < nch_fn; ch++){
                k2[ch_u][ch] = dt1 * me_sum[ch_u][ch];
                if (prob_joint[ch_u][ch] + (k1[ch_u][ch] + k2[ch_u][ch]) / 2 < 0){
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
            e = energy(prob_joint, nch_fn) * alpha;
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