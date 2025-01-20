#include <iostream>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
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


void generate_index_gamma(vector <vector <int> > &indexes, vector <int> vals, int gamma, int K, 
                          int w){
    if (w == K){
        indexes.push_back(vals);
    }else{
        for (int g = 0; g < gamma + 1; g++){
            vals.push_back(g);
            generate_index_gamma(indexes, vals, g, K, w + 1);
            vals.pop_back();
        }
    }
}


void generate_index_lp(vector <vector <int> > &indexes_ind, vector <int> ind_gamma, vector <int> ind_l,
                       int K, int w){
    if (w == K){
        indexes_ind.push_back(ind_l);
    }else{
        for (int lp = 0; lp < ind_gamma[w] + 1; lp++){
            ind_l.push_back(lp);
            generate_index_lp(indexes_ind, ind_gamma, ind_l, K, w + 1);
            ind_l.pop_back();
        }
    }
}


vector < vector < vector <int> > > generate_all_index_lp(vector <vector <int> > indexes_gamma, 
                                                         int K){
    vector < vector < vector < int > > > indexes_lp = vector < vector < vector <int> > > ();
    vector < int > vals = vector < int > (); 
    for (long ind = 0; ind < indexes_gamma.size(); ind++){
        vector < vector <int> > indexes_ind = vector < vector <int> > ();
        generate_index_lp(indexes_ind, indexes_gamma[ind], vals, K, 0);
        indexes_lp.push_back(indexes_ind);
    }
    return indexes_lp;
}


double count_independent_repetitions(vector <int> vals_g, vector <int> vals_l){
    double rep = 1;
    int w, counter;
    while(vals_g.size() > 1){
        w = 1;
        counter = 0;
        while (w < vals_g.size()){
            if (vals_g[w] == vals_g[0] && vals_l[w] == vals_l[0]){
                vals_g.erase(vals_g.begin() + w);
                vals_l.erase(vals_l.begin() + w);
                counter++;
            }else{
                w++;
            }
        }
        vals_g.erase(vals_g.begin());
        vals_l.erase(vals_l.begin());
        rep *= gsl_sf_fact(counter + 1);
    }
    return rep;
}


void prepare_index_maps(vector <vector <vector < pair <long, long> > > > &vals_2_ind, vector < vector <int> > &ncombs_full,
                        vector <vector <vector <int > > > &ncombs_glp, int max_gamma, int K,
                        vector <vector <int> > indexes_gamma, vector < vector < vector < int > > > indexes_lp,
                        vector <vector <vector <int > > > &place_there){
    ncombs_full = vector <vector <int> > (indexes_gamma.size(), vector <int> ());
    double rep;
    for (long ind_g = 0; ind_g < indexes_gamma.size(); ind_g++){
        for (long ind_l = 0; ind_l < indexes_lp[ind_g].size(); ind_l++){
            rep = count_independent_repetitions(indexes_gamma[ind_g], indexes_lp[ind_g][ind_l]);
            ncombs_full[ind_g].push_back(gsl_sf_fact(K) / rep);
        }
    }

    vals_2_ind = vector <vector <vector <pair <long, long> > > > (max_gamma + 1, vector <vector <pair <long, long> > > ());
    ncombs_glp = vector <vector <vector <int> > > (max_gamma + 1, vector <vector <int> > ());
    place_there = vector <vector <vector <int> > > (max_gamma + 1, vector <vector <int> > ());
    int w;
    bool cond;
    for (int gamma = 0; gamma < max_gamma + 1; gamma++){
        for (int lp = 0; lp < gamma + 1; lp++){
            vector <pair <long, long> > ind_in = vector <pair <long, long> > ();
            vector <int> comb_in = vector <int> ();
            vector <int> plc_in = vector <int> ();
            for (long ind_g = 0; ind_g < indexes_gamma.size(); ind_g++){
                for (long ind_l = 0; ind_l < indexes_lp[ind_g].size(); ind_l++){
                    w = 0;
                    cond = false;
                    while (w < K && !cond){
                        if (indexes_gamma[ind_g][w] == gamma && indexes_lp[ind_g][ind_l][w] == lp){
                            cond = true;
                        }
                        w++;
                    }
                    if (cond){
                        w--;
                        ind_in.push_back(pair <long, long> (ind_g, ind_l));

                        vector <int> indexes_g_exc(indexes_gamma[ind_g]);
                        indexes_g_exc.erase(indexes_g_exc.begin() + w);
                        vector <int> indexes_l_exc(indexes_lp[ind_g][ind_l]);
                        indexes_l_exc.erase(indexes_l_exc.begin() + w);

                        rep = count_independent_repetitions(indexes_g_exc, indexes_l_exc);
                        comb_in.push_back(gsl_sf_fact(K - 1) / rep);
                        
                        plc_in.push_back(w);
                    }
                }
            }

            vals_2_ind[gamma].push_back(ind_in);
            ncombs_glp[gamma].push_back(comb_in);
            place_there[gamma].push_back(plc_in);
        }
    }
}


// initializes all the joint and conditional probabilities
// prob_joint[index_g][index_lp][ch], with index_g that goes over all the possible
// combinations of connectivities and index_lp that goes over all the possible combinations
// of positive links lp. Each combination appears only once.
// ch=0,...,2^{K}-1 is the combination of the s_a. ch=2^{K}-1 is the unsat combination 
// pu_cond[gamma][lp][sj], gamma = 0, ..., max_gamma; lp=0,...,gamma; sj=0,1 the state in the conditional
// pi_gamma_lp[sj] only needs two values, and is re-used for every gamma and lp when computing pu_cond
// pjoint_gamma_ln[sj] saves the joint probability of having one variable sj and the rest
// of the clause in the unsat combinatio. Is re-used for every gamma and lp when computing pu_cond  
void init_probs(double ***&prob_joint, double ***&pu_cond, double *&pi_gamma_lp, 
                double *&pjoint_gamma_lp, double ***&me_sum, int K, int nch_fn, double p0, 
                int max_gamma, double alpha, vector <vector <int> > indexes_gamma, 
                vector < vector < vector < int > > > indexes_lp){
    double prod;
    int bit;
    prob_joint = new double **[indexes_gamma.size()];
    me_sum = new double **[indexes_gamma.size()];

    int ind_single, gamma_in, lp_in;

    for (long ind_g = 0; ind_g < indexes_gamma.size(); ind_g++){
        prob_joint[ind_g] = new double *[indexes_lp[ind_g].size()];
        me_sum[ind_g] = new double *[indexes_lp[ind_g].size()];
        for (long ind_l = 0; ind_l < indexes_lp[ind_g].size(); ind_l++){
            prob_joint[ind_g][ind_l] = new double [nch_fn];
            me_sum[ind_g][ind_l] = new double [nch_fn];
            for (int ch = 0; ch < nch_fn; ch++){
                prod = 1;
                for (int w = 0; w < K; w++){
                    gamma_in = indexes_gamma[ind_g][w];
                    lp_in = indexes_lp[ind_g][ind_l][w];
                    
                    bit = ((ch >> w) & 1);
                    prod *= (bit + (1 - 2 * bit) * p0) * 
                            gsl_ran_poisson_pdf(gamma_in, alpha * K) * 
                            binomial_coef(gamma_in, lp_in) / pow(2, gamma_in);
                }
                prob_joint[ind_g][ind_l][ch] = prod;
            }
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
void init_RK_arr(double ***&k1, double ***&k2, double ***&prob_joint_1, 
                 int nch_fn, vector <vector <int> > indexes_gamma, 
                 vector < vector < vector < int > > > indexes_lp){
    k1 = new double **[indexes_gamma.size()];
    k2 = new double **[indexes_gamma.size()];
    prob_joint_1 = new double **[indexes_gamma.size()];

    for (long ind_g = 0; ind_g < indexes_gamma.size(); ind_g++){
        k1[ind_g] = new double *[indexes_lp[ind_g].size()];
        k2[ind_g] = new double *[indexes_lp[ind_g].size()];
        prob_joint_1[ind_g] = new double *[indexes_lp[ind_g].size()];
        for (long ind_l = 0; ind_l < indexes_lp[ind_g].size(); ind_l++){
            k1[ind_g][ind_l] = new double [nch_fn];
            k2[ind_g][ind_l] = new double [nch_fn];
            prob_joint_1[ind_g][ind_l] = new double [nch_fn];
            for (int ch = 0; ch < nch_fn; ch++){
                k1[ind_g][ind_l][ch] = 0;
                k2[ind_g][ind_l][ch] = 0;
                prob_joint_1[ind_g][ind_l][ch] = 0;
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
void comp_pcond(double ***prob_joint, double ***pu_cond, double *pi_gamma_lp, 
                double *pjoint_gamma_lp, int K, int nch_fn, int max_gamma, 
                vector <vector <int> > indexes_gamma, vector < vector < vector < int > > > indexes_lp, 
                vector <vector <vector < pair <long, long> > > > vals_2_ind,
                vector <vector <vector <int > > > ncombs_glp, 
                vector <vector <vector <int > > > place_there){
    double pu;
    int bit, plc;
    long ind_g, ind_l;
    for (int gamma = 0; gamma < max_gamma + 1; gamma++){
        for (int lp = 0; lp < gamma + 1; lp++){
            for (int s = 0; s < 2; s++){
                pi_gamma_lp[s] = 0;
                pjoint_gamma_lp[s] = 0;
            }
            for (long i = 0; i < vals_2_ind[gamma][lp].size(); i++){
                ind_g = vals_2_ind[gamma][lp][i].first;
                ind_l = vals_2_ind[gamma][lp][i].second;
                plc = place_there[gamma][lp][i];
                for (int s = 0; s < 2; s++){
                    pjoint_gamma_lp[s] += prob_joint[ind_g][ind_l][(nch_fn - 1) ^ ((1 - s) << plc)] * 
                                          ncombs_glp[gamma][lp][i];  
                    // when s = 0, it inverts the last bit of the unsat combination (nch_fn - 1)
                    // when s = 1, it does not touch the last bit of the unsat combination
                }
                for (int ch = 0; ch < nch_fn; ch++){
                    bit = ((ch >> plc) & 1);
                    pi_gamma_lp[bit] += prob_joint[ind_g][ind_l][ch] * 
                                        ncombs_glp[gamma][lp][i];  
                }
            }

            for (int s = 0; s < 2; s++){
                pu_cond[gamma][lp][s] = pjoint_gamma_lp[s] / pi_gamma_lp[s]; 
            }
        }
    }
}


void sum_fms(int K, int gamma, int lp, int plc_he, double *prob_joint, double ***pu_cond, 
             double **rates, int nch_fn, int max_gamma, double e_av, double *me_sum_src){

    int bit, ch_flip, uns, uns_flip;

    int ln = gamma - lp;

    double sums[2][2]; 
    for (int s = 0; s < 2; s++){
        for (int sat = 0; sat < 2; sat++){
            sums[s][sat] = 0;
        }
    }

    if (ln > 0){
        for (int up = 0; up < lp + 1; up++){
            for (int un = 0; un < ln + 1; un++){
                sums[0][0] +=  binomial_coef(lp, up) * binomial_coef(ln, un) * rate_fms(un, up, rates, e_av) * 
                            pow(pu_cond[gamma][lp][0], up) * pow(1 - pu_cond[gamma][lp][0], lp - up) * 
                            pow(pu_cond[gamma][ln - 1][1], un) * pow(1 - pu_cond[gamma][ln - 1][1], ln - un);
                sums[1][0] +=  binomial_coef(lp, up) * binomial_coef(ln, un) * rate_fms(up, un, rates, e_av) * 
                            pow(pu_cond[gamma][lp][1], up) * pow(1 - pu_cond[gamma][lp][1], lp - up) * 
                            pow(pu_cond[gamma][ln - 1][0], un) * pow(1 - pu_cond[gamma][ln - 1][0], ln - un);               

                sums[0][1] +=  binomial_coef(lp, up) * binomial_coef(ln, un) * rate_fms(un, up + 1, rates, e_av) * 
                            pow(pu_cond[gamma][lp][0], up) * pow(1 - pu_cond[gamma][lp][0], lp - up) * 
                            pow(pu_cond[gamma][ln - 1][1], un) * pow(1 - pu_cond[gamma][ln - 1][1], ln - un);
                sums[1][1] +=  binomial_coef(lp, up) * binomial_coef(ln, un) * rate_fms(up + 1, un, rates, e_av) * 
                            pow(pu_cond[gamma][lp][1], up) * pow(1 - pu_cond[gamma][lp][1], lp - up) * 
                            pow(pu_cond[gamma][ln - 1][0], un) * pow(1 - pu_cond[gamma][ln - 1][0], ln - un);

            }
        }
    }else{
        for (int up = 0; up < lp + 1; up++){
            sums[0][0] +=  binomial_coef(lp, up) * rate_fms(0, up, rates, e_av) * 
                        pow(pu_cond[gamma][lp][0], up) * pow(1 - pu_cond[gamma][lp][0], lp - up);
            sums[1][0] +=  binomial_coef(lp, up) * rate_fms(up, 0, rates, e_av) * 
                        pow(pu_cond[gamma][lp][1], up) * pow(1 - pu_cond[gamma][lp][1], lp - up);               

            sums[0][1] +=  binomial_coef(lp, up) * rate_fms(0, up + 1, rates, e_av) * 
                        pow(pu_cond[gamma][lp][0], up) * pow(1 - pu_cond[gamma][lp][0], lp - up);
            sums[1][1] +=  binomial_coef(lp, up) * rate_fms(up + 1, 0, rates, e_av) * 
                        pow(pu_cond[gamma][lp][1], up) * pow(1 - pu_cond[gamma][lp][1], lp - up);
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
void der_fms(double ***prob_joint, double ***pu_cond, double **rates, int K, int l_max, 
             int nch_fn, double e_av, double ***me_sum,
             vector <vector <int> > indexes_gamma, vector < vector < vector < int > > > indexes_lp){
    for (long ind_g = 0; ind_g < indexes_gamma.size(); ind_g++){
        for (long ind_l = 0; ind_l < indexes_lp[ind_g].size(); ind_l++){
            for (int ch = 0; ch < nch_fn; ch++){
                me_sum[ind_g][ind_l][ch] = 0;
            }
        }
    }

    // candidate to be a parallel for
    #pragma omp parallel for
    for (long ind_g = 0; ind_g < indexes_gamma.size(); ind_g++){
        for (long ind_l = 0; ind_l < indexes_lp[ind_g].size(); ind_l++){
            for (int w = 0; w < K; w++){
                sum_fms(K, indexes_gamma[ind_g][w], indexes_lp[ind_g][ind_l][w], w, prob_joint[ind_g][ind_l], pu_cond, 
                        rates, nch_fn, l_max, e_av, me_sum[ind_g][ind_l]);
            }
        }
    }
}


double energy(double ***prob_joint, int nch_fn, vector < vector <int> > ncombs_full,
              vector <vector <int> > indexes_gamma, vector < vector < vector < int > > > indexes_lp){
    double e = 0;
    for (long ind_g = 0; ind_g < indexes_gamma[ind_g].size(); ind_g++){
        for (long ind_l = 0; ind_l < indexes_lp[ind_g].size(); ind_l++){
            e += prob_joint[ind_g][ind_l][nch_fn - 1] * ncombs_full[ind_g][ind_l];
        }
    }
    return e;
}


// peforms the integration of the differential equations with the 2nd order Runge-Kutta
// the method is implemented with adaptive step size
void RK2_fms(double alpha, int K, int nch_fn, double eta, int max_gamma, 
             double p0, char *fileener, double tl, double tol = 1e-2, double t0 = 0, double dt0 = 0.01, 
             double ef = 1e-6, double dt_min = 1e-7){
    double **rates;
    double ***prob_joint, ***pu_cond, ***me_sum, *pi_gamma_lp, *pjoint_gamma_lp;
    double e, error;                 
    
    vector <vector <vector < pair <long, long> > > > vals_2_ind;
    vector < vector <int> > ncombs_full;
    vector <vector <vector <int > > > ncombs_glp; 
    vector <vector <vector <int > > > place_there; 
    vector <vector <int> > indexes_gamma; 
    vector < vector < vector < int > > > indexes_lp;

    vector <int> vals;

    generate_index_gamma(indexes_gamma, vals, max_gamma, K, 0);
    indexes_lp = generate_all_index_lp(indexes_gamma, K);
    prepare_index_maps(vals_2_ind, ncombs_full, ncombs_glp, max_gamma, K, indexes_gamma, 
                       indexes_lp, place_there);
    
    long nindex_total = 0;
    for (long ind_g = 0; ind_g < indexes_gamma.size(); ind_g++){
        for (long ind_l = 0; ind_l < indexes_lp[ind_g].size(); ind_l++){
            nindex_total++;
        }
    }

    table_all_rates(max_gamma, K, eta, rates);
    
    init_probs(prob_joint, pu_cond, pi_gamma_lp, pjoint_gamma_lp, me_sum, K, nch_fn, p0, 
               max_gamma, alpha, indexes_gamma, indexes_lp);

    // initialize auxiliary arrays for the Runge-Kutta integration
    double ***k1, ***k2, ***prob_joint_1;
    init_RK_arr(k1, k2, prob_joint_1, nch_fn, indexes_gamma, indexes_lp);

    ofstream fe(fileener);
    
    e = energy(prob_joint, nch_fn, ncombs_full, indexes_gamma, indexes_lp) * alpha;
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
                   indexes_gamma, indexes_lp, vals_2_ind, ncombs_glp, place_there);

        der_fms(prob_joint, pu_cond, rates, K, max_gamma, nch_fn, e, me_sum, 
                indexes_gamma, indexes_lp);   // in the rates, I use the energy density

        valid = true;
        for (long ind_g = 0; ind_g < indexes_gamma.size(); ind_g++){
            for (long ind_l = 0; ind_l < indexes_lp[ind_g].size(); ind_l++){
                for (int ch = 0; ch < nch_fn; ch++){
                    k1[ind_g][ind_l][ch] = dt1 * me_sum[ind_g][ind_l][ch];
                    prob_joint_1[ind_g][ind_l][ch] = prob_joint[ind_g][ind_l][ch] + k1[ind_g][ind_l][ch];
                    if (prob_joint_1[ind_g][ind_l][ch] < 0){
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
            for (long ind_g = 0; ind_g < indexes_gamma.size(); ind_g++){
                for (long ind_l = 0; ind_l < indexes_lp[ind_g].size(); ind_l++){
                    for (int ch = 0; ch < nch_fn; ch++){
                        k1[ind_g][ind_l][ch] = dt1 * me_sum[ind_g][ind_l][ch];
                        prob_joint_1[ind_g][ind_l][ch] = prob_joint[ind_g][ind_l][ch] + k1[ind_g][ind_l][ch];
                        if (prob_joint_1[ind_g][ind_l][ch] < 0){
                            valid = false;
                        }
                    }
                }
            }
        }
        
        e = energy(prob_joint_1, nch_fn, ncombs_full, indexes_gamma, indexes_lp) * alpha;
        comp_pcond(prob_joint_1, pu_cond, pi_gamma_lp, pjoint_gamma_lp, K, nch_fn, max_gamma, 
                   indexes_gamma, indexes_lp, vals_2_ind, ncombs_glp, place_there);

        der_fms(prob_joint_1, pu_cond, rates, K, max_gamma, nch_fn, e, me_sum, 
                indexes_gamma, indexes_lp);
            
        valid = true;
        for (long ind_g = 0; ind_g < indexes_gamma.size(); ind_g++){
            for (long ind_l = 0; ind_l < indexes_lp[ind_g].size(); ind_l++){
                for (int ch = 0; ch < nch_fn; ch++){
                    k2[ind_g][ind_l][ch] = dt1 * me_sum[ind_g][ind_l][ch];
                    if (prob_joint[ind_g][ind_l][ch] + (k1[ind_g][ind_l][ch] + k2[ind_g][ind_l][ch]) / 2 < 0){
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
            e = energy(prob_joint, nch_fn, ncombs_full, indexes_gamma, indexes_lp) * alpha;
        }else{
            error = 0;
            for (long ind_g = 0; ind_g < indexes_gamma.size(); ind_g++){
                for (long ind_l = 0; ind_l < indexes_lp[ind_g].size(); ind_l++){
                    for (int ch = 0; ch < nch_fn; ch++){
                        error += fabs(k1[ind_g][ind_l][ch] - k2[ind_g][ind_l][ch]);
                    }
                }
            }

            error /= nch_fn * nindex_total;

            if (error < 2 * tol){
                //  cout << "step dt=" << dt1 << "  accepted" << endl;
                //  cout << "error=" << error << endl;
                t += dt1;
                for (long ind_g = 0; ind_g < indexes_gamma.size(); ind_g++){
                    for (long ind_l = 0; ind_l < indexes_lp[ind_g].size(); ind_l++){
                        for (int ch = 0; ch < nch_fn; ch++){
                            prob_joint[ind_g][ind_l][ch] += (k1[ind_g][ind_l][ch] + k2[ind_g][ind_l][ch]) / 2;
                        }
                    }
                }
                e = energy(prob_joint, nch_fn, ncombs_full, indexes_gamma, indexes_lp) * alpha;
                fe << t << "\t" << e << endl;

            }else{
                e = energy(prob_joint, nch_fn, ncombs_full, indexes_gamma, indexes_lp) * alpha;
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
    
    RK2_fms(alpha, K, nch_fn, eta, max_gamma, p0, fileener, tl, tol);
    

    

    return 0;
}