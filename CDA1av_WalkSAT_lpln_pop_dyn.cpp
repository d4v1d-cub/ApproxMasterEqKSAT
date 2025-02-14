#include <iostream>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <vector>
#include <map>


using namespace std;

void init_ran(gsl_rng * &r, unsigned long s){
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, s);
}

double binomial_coef(int n, int k){
    return gsl_sf_fact(n) / gsl_sf_fact(n - k) / gsl_sf_fact(k);
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


vector <double> get_probs(int nch_fn, double p0, int K){
    vector <double> probs(nch_fn);
    double prod;
    int bit;
    for (int ch = 0; ch < nch_fn; ch++){
        double prod = 1;
        int bit;
        for (int w = 0; w < K; w++){
            bit = ((ch >> w) & 1);
            prod *= (bit + (1 - 2 * bit) * p0);
        }
        probs[ch] = prod;
    }
    return probs;
}


vector <double> get_probs(int nch_fn, double p0, int K, vector <int> &lp_in, 
                          vector <int> &lp_in_2){
    vector <double> probs(nch_fn);
    double prod;
    int bit;
    for (int ch = 0; ch < nch_fn; ch++){
        double prod = 1;
        int bit;
        for (int w = 0; w < K; w++){
            bit = ((ch >> w) & 1);
            prod *= (bit + (1 - 2 * bit) * p0) * (lp_in_2[w] + 1) / (lp_in[w] + 1);
        }
        probs[ch] = prod;
    }
    return probs;
}


void update_map(map < pair <int, int>, vector <pair <long, int> > > &vals_2_ind, int gamma, int lp, int pop_index, int w){
    if (vals_2_ind.find(make_pair(gamma, lp)) == vals_2_ind.end()){
        vals_2_ind[make_pair(gamma, lp)] = vector <pair <long, int> > ();
    }
    vals_2_ind[make_pair(gamma, lp)].push_back(make_pair(pop_index, w));
}


// initializes all the joint and conditional probabilities
// prob_joint[i][ch], with i that goes over all population. Each element of the population.
// has two associated vectors (gamma_1, ..., gamma_K) with all the connectivities and (lp_1, ..., lp_K)
// with the number of satisfied variables in the clause.
// The population is such that for
// every pair (gamma, lp), one finds also the pair (gamma, gamma - lp - 1)
// ch=0,...,2^{K}-1 is the combination of the s_a. ch=2^{K}-1 is the unsat combination 
// pu_cond[gamma][lp][sj], gamma = 0, ..., max_gamma; lp=0,...,gamma; sj=0,1 the state in the conditional
// pi_gamma_lp[sj] only needs two values, and is re-used for every gamma and lp when computing pu_cond
// pjoint_gamma_ln[sj] saves the joint probability of having one variable sj and the rest
// of the clause in the unsat combinatio. Is re-used for every gamma and lp when computing pu_cond  
void init_probs(vector < vector <double> > &prob_joint, map < pair<long, int>, vector <double>  > &pu_cond, 
                vector < vector <double> > &me_sum, int K, int nch_fn, double p0, int max_gamma, double alpha, 
                long pop_size, map < pair <int, int>, vector <pair <long, int> > > &vals_2_ind,
                vector < pair < vector <int>, vector <int> > > &gamma_lp, gsl_rng *r){
    double prod;
    int bit;
    prob_joint = vector < vector <double> > ();
    me_sum = vector < vector <double> > ();

    gamma_lp = vector < pair < vector <int>, vector <int> > > ();


    for (long i = 0; i < pop_size / 2; i++){
        // Push first element

        me_sum.push_back(vector <double> (nch_fn));

        vector <int> gamma_in(K);
        vector <int> lp_in(K);
        for (int w = 0; w < K; w++){
            gamma_in[w] = gsl_ran_poisson(r, alpha * K);
            if (gamma_in[w] > max_gamma){
                gamma_in[w] = max_gamma;
            }
            lp_in[w] = gsl_ran_binomial(r, 0.5, gamma_in[w]);
            update_map(vals_2_ind, gamma_in[w], lp_in[w], 2 * i, w);
        }

        gamma_lp.push_back(make_pair(gamma_in, lp_in));
        prob_joint.push_back(get_probs(nch_fn, p0, K));

        // Push second element, guaranteeing that for all (gamma, lp), one also has the pair (gamma, gamma - lp - 1)

        me_sum.push_back(vector <double> (nch_fn));

        vector <int> lp_in_2(K);

        for (int w = 0; w < K; w++){
            if (lp_in[w] < gamma_in[w]){
                lp_in_2[w] = gamma_in[w] - lp_in[w] - 1;
            }else{
                lp_in_2[w] = lp_in[w];
            }   

            update_map(vals_2_ind, gamma_in[w], lp_in_2[w], 2 * i + 1, w);
        }

        gamma_lp.push_back(make_pair(gamma_in, lp_in_2));
        prob_joint.push_back(get_probs(nch_fn, p0, K, lp_in, lp_in_2));

    }

    for (const auto &pair: vals_2_ind){
        pu_cond[make_pair(pair.first.first, pair.first.second)] = vector <double> (2);
    }
}


// initializes the auxiliary arrays for the Runge-Kutta integration
void init_RK_arr(vector < vector <double> > &k1, vector < vector <double> > &k2, vector < vector <double> > &prob_joint_1, 
                 long real_pop_size, int nch_fn){
    k1 = vector < vector <double> > (real_pop_size, vector <double> (nch_fn));
    k2 = vector < vector <double> > (real_pop_size, vector <double> (nch_fn));
    prob_joint_1 = vector < vector <double> > (real_pop_size, vector <double> (nch_fn));
}


// initializes the arrays where the binomial weights will be stored. Those arrays are used 
// to compute the probability of selecting the variable belonging to the smallest number 
// of satisfied clauses
void init_poisson_probs(double *&poisson_probs, double *&poisson_sums, int max_c){
    poisson_probs = new double [max_c];
    poisson_sums = new double [max_c];
}


// it computes the binomial weights
void get_all_poisson_sums(int max_c, double pu_av, double *poisson_probs, 
                          double *poisson_sums, double mean_c){
    for (int sj = 0; sj < max_c; sj++){
        poisson_probs[sj] = gsl_ran_poisson_pdf(sj, (1 - pu_av) * mean_c);
    }
    
    for (int si = 0; si < max_c; si++){
        poisson_sums[si] = gsl_cdf_poisson_Q(si, (1 - pu_av) * mean_c);
    }
}



// it fills the array of the probabilities to be used when computing the walksat rate
void fill_pneigh(int S, double *poisson_probs, double *poisson_sums, 
                 double *pneigh){
    pneigh[0] = poisson_probs[S];
    pneigh[1] = poisson_sums[S];
}





// rate of the walksat algorithm used by Barthel et al. in 2003
// cj is a list of all the connectivities of the neighbors of node i that are 
// in unsatisfied clauses
// nch_exc is the number of possible combinations of the K-1 other variables in the clause
// therefore, nch_exc = 2^(K-1)
double rate_walksat(int E0, int S, int K, double q, double e_av, 
                    double *poisson_probs, double *poisson_sums, 
                    double *pneigh, int nch_exc){
    if (E0 > 0){
        bool cond;
        double cumul, prod;
        int cumul_bits, bit, k;
        cumul = 0;
    
        fill_pneigh(S, poisson_probs, poisson_sums, pneigh);
        for (int ch = 0; ch < nch_exc; ch++){
            prod = 1;
            cumul_bits = 0;
            for (int w = 0; w < K - 1; w++){
                bit = ((ch >> w) & 1);
                prod *= pneigh[bit];  
                cumul_bits += 1 - bit;
            }
            cumul += prod / (cumul_bits + 1);
        }
        return E0 * (q / K + (1 - q) * cumul) / e_av; 
    }else{
        return 0;
    }
}


// it computes the conditional probabilities of having a partially unsatisfied clause, given the 
// value of one variable in the clause
void comp_pcond(vector < vector <double> > &prob_joint, map < pair <long, int>, vector <double> > &pu_cond, 
                int K, int nch_fn, map < pair <int, int>, vector <pair <long, int> > > &vals_2_ind){
    int bit, plc, gamma, lp;
    long pop_index;
    long ind_pu_cond = 0;

    for (const auto& pair: vals_2_ind){
        double num[2] = {0, 0};
        double den[2] = {0, 0};

        for (long i = 0; i < pair.second.size(); i++){
            pop_index = pair.second[i].first;
            plc = pair.second[i].second;
            for (int s = 0; s < 2; s++){
                num[s] += prob_joint[pop_index][(nch_fn - 1) ^ ((1 - s) << plc)];  
                // when s = 0, it inverts the bit of the unsat combination (nch_fn - 1) at the position plc
                // when s = 1, it does not touch the unsat combination
            }

            for (int ch = 0; ch < nch_fn; ch++){
                bit = ((ch >> plc) & 1);
                den[bit] += prob_joint[pop_index][ch];  
            }
        }

        for (int s = 0; s < 2; s++){
            pu_cond[pair.first][s] = num[s] / den[s]; 
        }
    }
}


void sum_walksat(int K, int gamma, int lp, int plc_he, vector <double> &prob_joint, 
             map < pair <long, int>, vector <double> > &pu_cond, 
             double *poisson_probs, double *poisson_sums, 
             double q, int nch_fn, double e_av, vector <double> &me_sum_src){

    int bit, ch_flip, uns, uns_flip;

    int ln = gamma - lp;

    double sums[2][2]; 
    sums[0][0] = 0;
    sums[1][0] = 0;
    sums[1][1] = 0;

    double p_neigh[2];

    pair <long, int> pair_p = make_pair(gamma, lp);  

    if (ln > 0){
        pair <long, int> pair_n = make_pair(gamma, ln - 1);
        for (int un = 0; un < ln + 1; un++){
            sums[0][0] +=  binomial_coef(ln, un) * 
                        rate_walksat(un, gamma + 1 - un, K, q, e_av, poisson_probs, poisson_sums, p_neigh, nch_fn / 2) * 
                        pow(pu_cond[pair_n][1], un) * pow(1 - pu_cond[pair_n][1], ln - un);
        }
        for (int up = 0; up < lp + 1; up++){
            sums[1][0] += binomial_coef(lp, up) * 
                        rate_walksat(up, gamma + 1 - up, K, q, e_av, poisson_probs, poisson_sums, p_neigh, nch_fn / 2) * 
                        pow(pu_cond[pair_p][1], up) * pow(1 - pu_cond[pair_p][1], lp - up); 
                        
            sums[1][1] +=  binomial_coef(lp, up) * 
                        rate_walksat(up + 1, gamma - up, K, q, e_av, poisson_probs, poisson_sums, p_neigh, nch_fn / 2) * 
                        pow(pu_cond[pair_p][1], up) * pow(1 - pu_cond[pair_p][1], lp - up);
        }

    }else{
        sums[0][0] = rate_walksat(0, gamma + 1, K, q, e_av, poisson_probs, poisson_sums, p_neigh, nch_fn / 2);
        for (int up = 0; up < lp + 1; up++){
            sums[1][0] +=  binomial_coef(lp, up) * 
                        rate_walksat(up, gamma + 1 - up, K, q, e_av, poisson_probs, poisson_sums, p_neigh, nch_fn / 2) * 
                        pow(pu_cond[pair_p][1], up) * pow(1 - pu_cond[pair_p][1], lp - up);               
            sums[1][1] +=  binomial_coef(lp, up) * 
                        rate_walksat(up + 1, gamma - up, K, q, e_av, poisson_probs, poisson_sums, p_neigh, nch_fn / 2) * 
                        pow(pu_cond[pair_p][1], up) * pow(1 - pu_cond[pair_p][1], lp - up);
        }
    }


    for (int ch_src = 0; ch_src < nch_fn; ch_src++){
        bit = ((ch_src >> plc_he) & 1);
        ch_flip = (ch_src ^ (1 << plc_he));
        uns = (ch_src == nch_fn - 1);
        uns_flip = (ch_flip == nch_fn - 1);
        me_sum_src[ch_src] += -sums[bit][uns] * prob_joint[ch_src] + 
                            sums[1 - bit][uns_flip] * prob_joint[ch_flip];
        // If bit = 0, we can only have uns = 0. So in the first line we put sums[0][0], which is OK.
        // In the second line, we can have uns_flip = 0 or 1. There are two options: sums[1][0] or sums[1][1].
        // If bit = 1, we can have uns = 0 or 1. In the first line, we can have sums[1][0] or sums[1][1].
        // In the second line, we can only have uns_flip = 0. So we put sums[0][0].
        // We never need to use sums[0][1]. 
    }
    
}


// it computes all the derivatives of the joint probabilities
void der_walksat(vector <vector <double> > &prob_joint, map < pair <long, int>, vector <double> > &pu_cond, 
             double *poisson_probs, double *poisson_sums, double q, int K, int nch_fn, double e_av, vector <vector <double> > &me_sum,
             vector < pair < vector <int>, vector <int> > > &gamma_lp, 
             map < pair <int, int>, vector <pair <long, int> > > &vals_2_ind){
    for (long pop_ind = 0; pop_ind < me_sum.size(); pop_ind++){
        for (int ch = 0; ch < nch_fn; ch++){
            me_sum[pop_ind][ch] = 0;
        }
    }

    // candidate to be a parallel for
    #pragma omp parallel for
    for (long pop_ind = 0; pop_ind < prob_joint.size(); pop_ind++){
        for (int w = 0; w < K; w++){
            sum_walksat(K, gamma_lp[pop_ind].first[w], gamma_lp[pop_ind].second[w], w, 
                    prob_joint[pop_ind], pu_cond, poisson_probs, poisson_sums, q, nch_fn, 
                    e_av, me_sum[pop_ind]);
        }
    }
}


double energy(vector <vector <double> > &prob_joint, int nch_fn){
    double e = 0;
    for (long pop_ind = 0; pop_ind < prob_joint.size(); pop_ind++){
        e += prob_joint[pop_ind][nch_fn - 1];
    }
    return e / prob_joint.size();
}


// peforms the integration of the differential equations with the 2nd order Runge-Kutta
// the method is implemented with adaptive step size
void RK2_walksat(double alpha, int K, int nch_fn, double q, int max_gamma, long pop_size, gsl_rng *r,  
             double p0, char *fileener, double tl, double tol = 1e-2, double t0 = 0, double dt0 = 0.01, 
             double ef = 1e-6, double dt_min = 1e-7){
    double *poisson_probs, *poisson_sums;
    double **rates;
    double e, error, pu_av;                 

    init_poisson_probs(poisson_probs, poisson_sums, max_gamma + 1);
    
    vector <vector <double> > prob_joint;
    map < pair <long, int>, vector <double> > pu_cond;
    vector < vector <double> > me_sum;
    vector < pair < vector <int>, vector <int> > > gamma_lp;
    map < pair <int, int>, vector <pair <long, int> > > vals_2_ind;
    
    init_probs(prob_joint, pu_cond, me_sum, K, nch_fn, p0, 
               max_gamma, alpha, pop_size, vals_2_ind, gamma_lp, r);

    // initialize auxiliary arrays for the Runge-Kutta integration
    vector <vector <double> > k1, k2, prob_joint_1;
    init_RK_arr(k1, k2, prob_joint_1, prob_joint.size(), nch_fn);

    ofstream fe(fileener);
    
    pu_av = energy(prob_joint, nch_fn);
    e = pu_av * alpha;
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

        comp_pcond(prob_joint, pu_cond, K, nch_fn, vals_2_ind);
        get_all_poisson_sums(max_gamma + 1, pu_av, poisson_probs, poisson_sums, alpha * K);

        der_walksat(prob_joint, pu_cond, poisson_probs, poisson_sums, q, K, nch_fn, e, me_sum, 
                    gamma_lp, vals_2_ind);   // in the rates, I use the energy density

        valid = true;
        for (long pop_ind = 0; pop_ind < prob_joint.size(); pop_ind++){
            for (int ch = 0; ch < nch_fn; ch++){
                k1[pop_ind][ch] = dt1 * me_sum[pop_ind][ch];
                prob_joint_1[pop_ind][ch] = prob_joint[pop_ind][ch] + k1[pop_ind][ch];
                if (prob_joint_1[pop_ind][ch] < 0){
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
            for (long pop_ind = 0; pop_ind < prob_joint.size(); pop_ind++){
                for (int ch = 0; ch < nch_fn; ch++){
                    k1[pop_ind][ch] = dt1 * me_sum[pop_ind][ch];
                    prob_joint_1[pop_ind][ch] = prob_joint[pop_ind][ch] + k1[pop_ind][ch];
                    if (prob_joint_1[pop_ind][ch] < 0){
                        valid = false;
                    }
                }
            }
        }
        
        pu_av = energy(prob_joint_1, nch_fn);
        e = pu_av * alpha;
        comp_pcond(prob_joint_1, pu_cond, K, nch_fn, vals_2_ind);
        get_all_poisson_sums(max_gamma + 1, pu_av, poisson_probs, poisson_sums, alpha * K);

        der_walksat(prob_joint_1, pu_cond, poisson_probs, poisson_sums, q, K, nch_fn, e, 
                    me_sum, gamma_lp, vals_2_ind);
            
        valid = true;
        for (long pop_ind = 0; pop_ind < prob_joint.size(); pop_ind++){
            for (int ch = 0; ch < nch_fn; ch++){
                k2[pop_ind][ch] = dt1 * me_sum[pop_ind][ch];
                if (prob_joint[pop_ind][ch] + (k1[pop_ind][ch] + k2[pop_ind][ch]) / 2 < 0){
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
            pu_av = energy(prob_joint, nch_fn);
            e = pu_av * alpha;
        }else{
            error = 0;
            for (long pop_ind = 0; pop_ind < prob_joint.size(); pop_ind++){
                for (int ch = 0; ch < nch_fn; ch++){
                    error += fabs(k1[pop_ind][ch] - k2[pop_ind][ch]);
                }
            }

            error /= nch_fn * prob_joint.size();

            if (error < 2 * tol){
                //  cout << "step dt=" << dt1 << "  accepted" << endl;
                //  cout << "error=" << error << endl;
                t += dt1;
                for (long pop_ind = 0; pop_ind < prob_joint.size(); pop_ind++){
                    for (int ch = 0; ch < nch_fn; ch++){
                        prob_joint[pop_ind][ch] += (k1[pop_ind][ch] + k2[pop_ind][ch]) / 2;
                    }
                }
                pu_av = energy(prob_joint, nch_fn); 
                e = pu_av * alpha;
                fe << t << "\t" << e << endl;

            }else{
                pu_av = energy(prob_joint, nch_fn); 
                e = pu_av * alpha;
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
    long pop_size = atol(argv[1]);
    double alpha = atof(argv[2]);
    int K = atoi(argv[3]);
    unsigned long seed_r = atol(argv[4]);
    double q = atof(argv[5]);
    double tl = atof(argv[6]);
    double tol = atof(argv[7]);
    int nthr = atoi(argv[8]);
    double eps_c = atof(argv[9]);

    int nch_fn = (1 << K);
    double p0 = 0.5;

    gsl_rng * r;
    init_ran(r, seed_r);

    omp_set_num_threads(nthr);

    char fileener[300]; 
    sprintf(fileener, "CDA1av_lpln_popdyn_WalkSAT_ener_K_%d_alpha_%.4lf_q_%.4lf_tl_%.2lf_tol_%.1e_epsc_%.e_popsize_%li_seed_%li.txt", 
            K, alpha, q, tl, tol, eps_c, pop_size, seed_r);


    int max_gamma = get_max_gamma(alpha, K, eps_c);
    
    RK2_walksat(alpha, K, nch_fn, q, max_gamma, pop_size, r, p0, fileener, tl, tol);
    

    

    return 0;
}