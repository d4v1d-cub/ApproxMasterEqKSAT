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
    


void init_ran(gsl_rng * &r, unsigned long s){
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, s);
}


typedef struct{
    int nfacn;      // number of factor nodes
    vector <long> fn_in; // list of factor nodes to which the nodes belongs
    vector <int> pos_fn;  // position in each factor node's list of nodes.
    long nch;   // the number of possible combinations of the states of 
                    // the neighboring clauses.
    long **fn_exc;  // remaining factor nodes after one removes a specific factor node
                    // first index: removed fn, second index: list of remaining fn.
    int **pos_fn_exc;   // position in the remaining fn in the same order as in fn_exc.
}Tnode;


typedef struct{
    int ch_unsat;  // the combination of the nodes that makes the clause unsatisfied.
    vector <long> nodes_in;  // nodes inside the factor node
    vector <int> links; // links to those nodes
    vector <int> pos_n;   // position in each node's list of factor nodes.
    long **nodes_exc;   // remaining nodes after one removes a specific node from the factor node.
}Thedge;


void init_graph(Tnode *&nodes, Thedge *&hedges, long N, long M){
    nodes = new Tnode[N];
    hedges = new Thedge[M];

    for (long i = 0; i < N; i++){
        nodes[i].nfacn = 0;
        nodes[i].fn_in = vector <long> ();
        nodes[i].pos_fn = vector <int> ();
    }

    for (long he = 0; he < M; he++){
        hedges[he].links = vector <int> ();
        hedges[he].nodes_in = vector <long> ();
        hedges[he].pos_n = vector <int> ();
    }
}


// This function reads all the information about the graph from a file.
void read_graph(char *filegraph, long N, long M, int K, Tnode *&nodes, Thedge *&hedges){
    
    init_graph(nodes, hedges, N, M);
    string trash_str;
    double trash_double;

    vector <long> nodes_in;
    for (int j = 0; j < K; j++){
        nodes_in.push_back(0);
    }

    ifstream fg(filegraph);

    long fn_count = 0;
    int he;
    bool new_fn;

    for (long i = 0; i < N; i++){
        fg >> trash_double;
        fg >> nodes[i].nfacn;
        nodes[i].nch = (1 >> nodes[i].nfacn); 
        getline(fg, trash_str);
        getline(fg, trash_str);
        nodes_in[0] = i;

        for (int k = 0; k < nodes[i].nfacn; k++){
            new_fn = true;
            for (int j = 0; j < K - 1; j++){
                fg >> trash_double;
                fg >> trash_double;
                fg >> nodes_in[j + 1];
                if (nodes_in[j + 1] < nodes_in[0]){
                    new_fn = false;
                }
                fg >> trash_double;
                fg >> trash_double;
            }
            // sort(nodes_in.begin(), nodes_in.end());
            if (new_fn){
                for (int w = 0; w < K; w++){
                    nodes[nodes_in[w]].fn_in.push_back(fn_count);
                    nodes[nodes_in[w]].pos_fn.push_back(w);
                    hedges[fn_count].nodes_in.push_back(nodes_in[w]);
                    hedges[fn_count].pos_n.push_back(nodes[nodes_in[w]].fn_in.size() - 1);
                }
                fn_count++;
            }
        }
    }

    for (long i = 0; i < N; i++){
        if (nodes[i].nfacn != nodes[i].fn_in.size()){
            cout << "Problem with node " << i << endl;
        }
    }

    fg.close();
}

// This function read the links from a file
void read_links(char *filelinks, long N, long M, int K, Tnode *nodes, Thedge *hedges){
    ifstream fl(filelinks);
    int trash_int, link;

    for (long he = 0; he < M; he++){
        hedges[he].ch_unsat = 0;
        for (int w = 0; w < K; w++){
           hedges[he].links.push_back(0);
        }
    }

    for (long i = 0; i < N; i++){
        fl >> trash_int;
        for (int hind = 0; hind < nodes[i].nfacn; hind++){
            fl >> link;
            hedges[nodes[i].fn_in[hind]].links[nodes[i].pos_fn[hind]] = link;
            hedges[nodes[i].fn_in[hind]].ch_unsat += (((1 + link) / 2) << nodes[i].pos_fn[hind]);
        }
    }

    fl.close();
}

double av(Thedge *hedges, long M){
    double answ = 0;
    for (long i = 0; i < M; i++){
        answ += hedges[i].nodes_in.size();
    }
    return answ / M;
}

// This function creates at run time.
void create_graph(long N, long M, int K, Tnode *&nodes, Thedge *&hedges, gsl_rng * r){
    init_graph(nodes, hedges, N, M);
    int w, h;
    long var;
    bool cond;
    for (long he = 0; he < M; he++){
        hedges[he].ch_unsat = 0;
        w = 0;
        while (w < K){
            var = gsl_rng_uniform_int(r, N);
            cond = true;
            h = 0;
            while (h < hedges[he].nodes_in.size() && cond){
                if (hedges[he].nodes_in[h] == var){
                    cond = false;
                }
                h++;
            }

            if (cond){
                hedges[he].nodes_in.push_back(var);
                if (gsl_rng_uniform_pos(r) < 0.5){
                    hedges[he].links.push_back(1);
                    hedges[he].ch_unsat += (1 << w);
                }else{
                    hedges[he].links.push_back(-1);
                }
                hedges[he].pos_n.push_back(nodes[var].nfacn);
                nodes[var].fn_in.push_back(he);
                nodes[var].pos_fn.push_back(w);
                nodes[var].nfacn++;
                w++;
            }
        }
    }

    for (long i = 0; i < N; i++){
            nodes[i].nch = (long) pow(2, nodes[i].nfacn);
    }
}


void get_info_exc(Tnode *nodes, Thedge *hedges, long N, long M, int K){
    int w, count;
    for (long he = 0; he < M; he++){
        hedges[he].nodes_exc = new long *[K];
        for (int j = 0; j < K; j++){
            hedges[he].nodes_exc[j] = new long [K - 1];
            count = 0;
            w = (j + 1) % K;
            while (w != j){
                hedges[he].nodes_exc[j][count] = hedges[he].nodes_in[w];
                w = (w + 1) % K;
                count++;
            }
        } 
    }

    int other;
    for (long i = 0; i < N; i++){
        nodes[i].fn_exc = new long *[nodes[i].nfacn];
        nodes[i].pos_fn_exc = new int *[nodes[i].nfacn];
        for (int hind = 0; hind < nodes[i].nfacn; hind++){
            nodes[i].fn_exc[hind] = new long [nodes[i].nfacn - 1];
            nodes[i].pos_fn_exc[hind] = new int [nodes[i].nfacn - 1];
            other = (hind + 1) % nodes[i].nfacn;
            count = 0;
            while (other != hind){
                nodes[i].fn_exc[hind][count] = nodes[i].fn_in[other];
                nodes[i].pos_fn_exc[hind][count] = nodes[i].pos_fn[other];
                count++;
                other = (other + 1) % nodes[i].nfacn;
            }
        }
    }
}


int get_max_c(Tnode *nodes, long N){
    int max_c = 0;
    for (long i = 0; i < N; i++){
        if (nodes[i].nfacn > max_c){
            max_c = nodes[i].nfacn;
        }
    }
    return max_c;
}


// initializes all the joint and conditional probabilities
void init_probs(double **&prob_joint, double ***&pu_cond, double **&pi, double **&me_sum, long M, int K, 
                int nch_fn, double p0){
    double prod;
    int bit;
    prob_joint = new double *[M];
    me_sum = new double *[M];
    pu_cond = new double **[M];
    for (long he = 0; he < M; he++){
        prob_joint[he] = new double [nch_fn];
        me_sum[he] = new double [nch_fn];
        for (int ch = 0; ch < nch_fn; ch++){
            prod = 1;
            for (int w = 0; w < K; w++){
                bit = ((ch >> w) & 1);
                prod *= (bit + (1 - 2 * bit) * p0);
            }
            prob_joint[he][ch] = prod;
        }

        pu_cond[he] = new double*[K];
        for (int w = 0; w < K; w++){
            pu_cond[he][w] = new double [2];
        }
    }

    pi = new double*[K];
    for (int w = 0; w < K; w++){
        pi[w] = new double[2];
    }
}


// initializes the auxiliary arrays for the Runge-Kutta integration
void init_RK_arr(double **&k1, double **&k2, double **&prob_joint_1, long M, 
                int nch_fn){
    k1 = new double *[M];
    k2 = new double *[M];
    prob_joint_1 = new double *[M];
    for (long he = 0; he < M; he++){
        k1[he] = new double [nch_fn];
        k2[he] = new double [nch_fn];
        prob_joint_1[he] = new double [nch_fn];
        for (int ch = 0; ch < nch_fn; ch++){
            k1[he][ch] = 0;
            k2[he][ch] = 0;
            prob_joint_1[he][ch] = 0;
        }
    }
}


// rate of the Focused Metropolis Search algorithm.
double rate_fms(int E0, int E1, int K, double eta, double e_av){
    double dE = E1 - E0;
    if (dE > 0){
        return E0 / K / e_av * pow(eta, -dE);
    }else{
        return E0 / K / e_av;
    }
}


// initializes the arrays where the binomial weights will be stored. Those arrays are used 
// to compute the probability of selecting the variable belonging to the smallest number 
// of satisfied clauses
void init_aux_arr(double **&binom_probs, double **&binom_sums, int max_c, int K){
    binom_probs = new double *[max_c];
    for (int gj = 0; gj < max_c; gj++){
        binom_probs[gj] = new double[gj + 1];
    }

    binom_sums = new double *[max_c];
    for (int s = 0; s < max_c; s++){
        binom_sums[s] = new double[max_c - s];
    }
}


// it computes the binomial weights
void get_all_binom_sums(int max_c, double pu_av, double **binom_probs, 
                        double **binom_sums){
    for (int gj = 0; gj < max_c; gj++){
        for (int sj = 0; sj < gj + 1; sj++){
            binom_probs[gj][sj] = gsl_ran_binomial_pdf(sj, 1 - pu_av, gj);
        }
    }
    
    for (int si = 0; si < max_c; si++){
        for (int gj = si; gj < max_c; gj++){
            binom_sums[si][gj - si] = gsl_cdf_binomial_Q(si, 1 - pu_av, gj);
        }
    }
}


// it fills the array of the probabilities to be used when computing the walksat rate
void fill_pneigh(int E0, int S, int *cj, double **binom_probs, double **binom_sums, 
                 double **pneigh, int K){
    for (int j = 0; j < K - 1; j++){
        pneigh[j][0] = binom_probs[cj[j] - 1][S];
        pneigh[j][1] = binom_sums[S][cj[j] - 1 - S];
    }
}


// rate of the walksat algorithm used by Barthel et al. in 2003
// cj is a list of all the connectivities of the neighbors of node i that are 
// in unsatisfied clauses
// nch_exc is the number of possible combinations of the K-1 other variables in the clause
// therefore, nch_exc = 2^(K-1)
double rate_walksat(int E0, int S, int K, double q, double e_av, int **cj, 
                    double **binom_probs, double **binom_sums, double **pneigh, 
                    int nch_exc){
    bool cond;
    double cumul, prod;
    int cumul_bits, bit, k;
    cumul = 0;
    for (int fn = 0; fn < E0; fn++){
        cond = true;
        k = 0;
        while (k < K - 1 && cond){
            if (cj[fn][k] - 1 < S){
                cond = false;    // if one of the other nodes in the clause has less 
                // than 'S' neighbors, it cannot belong to more than 'S' satisfied clauses
                // therefore, the variable for which one is computing the rate cannot be chosen
                // by walksat dynamics over that neighbor in that specific clause
            }
            k++;
        }
        if (cond){   
            fill_pneigh(E0, S, cj[fn], binom_probs, binom_sums, pneigh, K);
            for (int ch = 0; ch < nch_exc; ch++){
                prod = 1;
                cumul_bits = 0;
                for (int w = 0; w < K - 1; w++){
                    bit = ((ch >> w) & 1);
                    prod *= pneigh[w][bit];  // when bit=0, one uses the probability 
                    // of finding that the neighbor (fn, w) belongs to the exact same number of
                    // satisfied clauses. When bit=1, one uses the probability that the neighbor
                    // (fn, w) belongs to more than S satisfied clauses.
                    cumul_bits += 1 - bit;
                }
                cumul += prod / (cumul_bits + 1);  // When other neighbors have the same 
                // number of satisfied clauses S, the variable inside the rate will be picked 
                // uniformly at random among those variables with the same S.
            }
        }
    }
    return q * E0 / e_av / K + (1 - q) * cumul / e_av / K; 
}


// it computes the conditional probabilities of having a partially unsatisfied clause, given the 
// value of one variable in the clause
void comp_pcond(double **prob_joint, double ***pu_cond, double **pi, Thedge *hedges, long M, int K, 
                int nch_fn){
    double pu;
    int bit;
    for (long he = 0; he < M; he++){
        pu = prob_joint[he][hedges[he].ch_unsat];
        for (int w = 0; w < K; w++){
            for (int s = 0; s < 2; s++){
                pi[w][s] = 0;
            }
        }

        for (int ch = 0; ch < nch_fn; ch++){
            for (int w = 0; w < K; w++){
                bit = ((ch >> w) & 1);
                pi[w][bit] += prob_joint[he][ch];
            }
        }

        for (int w = 0; w < K; w++){
            for (int s = 0; s < 2; s++){
                pu_cond[he][w][s] = pu / pi[w][s];
            }
        }
    }
}


// It takes the product over all the conditional probabilities of the neighboring factor nodes,
// except for the origin fn_src
double prodcond(double ***pu_cond, int fn_src, Tnode node, int s, long ch){
    double prod = 1;
    int bit;
    for (int other = 0; other < node.nfacn - 1; other++){
        bit = ((ch >> other) & 1);
        prod *= 1 - bit - 
                (1 - 2 * bit) * pu_cond[node.fn_exc[fn_src][other]][node.pos_fn_exc[fn_src][other]][s];
    }
    return prod;
}


// it does the sum in the derivative of the CDA equations
// fn_src is the origin factor node where one is computing the derivative
// part_uns is 1 if the other variables in fn_src are partially
// unsatisfying their links, and is 0 otherwise. 
void sum_walksat(long node, int fn_src, Tnode *nodes, Thedge *hedges, 
                 double *prob_joint, double ***pu_cond, double **binom_probs, 
                 double **binom_sums, int K, int nch_fn, double q, 
                 double e_av, double *me_sum_src, int max_c){
    int he, plc_he;
    bool bit, uns, uns_flip;
    int ch_flip;
    double prod[2], r[2][2];
    int E[2];
    double **pneigh = (double **) malloc((K - 1) * sizeof(double *));
    for (int j = 0; j < K - 1; j++){
        pneigh[j] = (double*) malloc(2 * sizeof(double));
    }

    int ***cj = (int ***) malloc (2 * sizeof(int **));
    for (int s = 0; s < 2; s++){
        cj[s] = (int **) malloc (max_c * sizeof(int *));
        for (int h = 0; h < max_c; h++){
            cj[s][h] = (int *) malloc ((K - 1) * sizeof(int));
        }
    }

    

    for (long ch = 0; ch < nodes[node].nch / 2; ch++){
        prod[0] = prodcond(pu_cond, fn_src, nodes[node], 0, ch);
        prod[1] = prodcond(pu_cond, fn_src, nodes[node], 1, ch);
        E[0] = 0;
        E[1] = 0;
        for (int other = 0; other < nodes[node].nfacn - 1; other++){
            if ((ch >> other) & 1){
                he = nodes[node].fn_exc[fn_src][other];
                plc_he = nodes[node].pos_fn_exc[fn_src][other];
                bit = ((hedges[he].ch_unsat >> plc_he) & 1);  // if s=0 unsatisfies the clause, 
                // bit = 0 (false), otherwise bit = 1.
                for (int h = 0; h < K - 1; h++){
                    cj[bit][E[bit]][h] = 
                            nodes[hedges[he].nodes_exc[plc_he][h]].nfacn;
                    // It's the connectivity of one of the other nodes inside the factor node: 
                    // 'he = nodes[node].fn_in[hind]'
                }
                E[bit]++;
            }
        }

        r[0][0] = rate_walksat(E[0], nodes[node].nfacn - E[0], K, q, e_av, cj[0], binom_probs, 
                            binom_sums, pneigh, nch_fn / 2);
        r[1][0] = rate_walksat(E[1], nodes[node].nfacn - E[1], K, q, e_av, cj[1], binom_probs, 
                            binom_sums, pneigh, nch_fn / 2);

        he = nodes[node].fn_in[fn_src];
        plc_he = nodes[node].pos_fn[fn_src];
        bit = ((hedges[he].ch_unsat >> plc_he) & 1);  
        for (int h = 0; h < K - 1; h++){
            cj[bit][E[bit]][h] = 
                nodes[hedges[he].nodes_exc[plc_he][h]].nfacn;
                    // It's the connectivity of one of the other nodes inside the factor node: 
                    // 'he = nodes[node].fn_in[hind]'
        }

        r[bit][1] = rate_walksat(E[bit] + 1, nodes[node].nfacn - E[bit] - 1, K, q, e_av, cj[bit], binom_probs, 
                            binom_sums, pneigh, nch_fn);
        
        for (int ch_src = 0; ch_src < nch_fn; ch_src++){
            bit = ((ch_src >> plc_he) & 1);
            ch_flip = (ch_src ^ (1 << plc_he));
            uns = (ch_src == hedges[he].ch_unsat);
            uns_flip = (ch_flip == hedges[he].ch_unsat);
            me_sum_src[ch_src] += -r[bit][uns] * prod[bit] * prob_joint[ch_src] + 
                                  r[1 - bit][uns_flip] * prod[1 - bit] * prob_joint[ch_flip];
        }
    }

    // delete []pneigh;
    // delete []cj;
    free(pneigh);
    free(cj);
}


// it computes all the derivatives of the joint probabilities
void der_walksat(Tnode *nodes, Thedge *hedges, double **prob_joint, double ***pu_cond, 
                 double **binom_probs, double **binom_sums, long M, int K, int nch_fn, 
                 double q, double e_av, double **me_sum, int max_c){
    for (long he = 0; he < M; he++){
        for (int ch = 0; ch < nch_fn; ch++){
            me_sum[he][ch] = 0;
        }
    }

    // candidate to be a parallel for
    #pragma omp parallel for
    for (long he = 0; he < M; he++){
        for (int w = 0; w < K; w++){
            sum_walksat(hedges[he].nodes_in[w], hedges[he].pos_n[w], nodes, hedges, 
                        prob_joint[he], pu_cond, binom_probs, binom_sums, K, nch_fn,
                        q, e_av, me_sum[he], max_c);
        }
    }
}


double energy(double **prob_joint, Thedge *hedges, long M){
    double e = 0;
    for (long he = 0; he < M; he++){
        e += prob_joint[he][hedges[he].ch_unsat];
    }
    return e;
}


double norm(double *probs, int nelems){
    double n = 0;
    for (int i = 0; i < nelems; i++){
        n += probs[i];
    }
    return n;
}


// peforms the integration of the differential equations with the 2nd order Runge-Kutta
// the method is implemented with adaptive step size
void RK2_walksat(Tnode *nodes, Thedge *hedges, long N, long M, int K, int nch_fn, double q, int max_c, 
                 double p0, char *fileener, double tl, double tol = 1e-2, double t0 = 0, double dt0 = 0.01, 
                 double ef = 1e-6, double dt_min = 1e-7){
    double **binom_probs, **binom_sums;
    int ***cj;
    double **prob_joint, ***pu_cond, **me_sum, **pi;
    double e, pu_av, error, dif_norm;                 
                 
    // initalize all arrays that will be used inside the derivative
    init_aux_arr(binom_probs, binom_sums, max_c, K);
    init_probs(prob_joint, pu_cond, pi, me_sum, M, K, nch_fn, p0);



    // initialize auxiliary arrays for the Runge-Kutta integration
    double **k1, **k2, **prob_joint_1;
    init_RK_arr(k1, k2, prob_joint_1, M, nch_fn);

    ofstream fe(fileener);
    
    e = energy(prob_joint, hedges, M);
    pu_av = e / M;
    fe << t0 << "\t" << e / N << endl;   // it prints the energy density

    double dt1 = dt0;
    double t = t0;

    bool valid;

    // the time scale is already given in Monte Carlo steps. Inside the rates I am using 
    // the energy density e_av
    while (t < tl){
        if (e < ef){
            cout << "Final energy reached" << endl;
            break;
        }

        auto t1 = std::chrono::high_resolution_clock::now();

        comp_pcond(prob_joint, pu_cond, pi, hedges, M, K, nch_fn);
        get_all_binom_sums(max_c, pu_av, binom_probs, binom_sums);

        der_walksat(nodes, hedges, prob_joint, pu_cond, binom_probs, binom_sums, 
                    M, K, nch_fn, q, e / N, me_sum, max_c);   // in the rates, I use the energy density

        for (long he = 0; he < M; he++){
            for (int ch = 0; ch < nch_fn; ch++){
                k1[he][ch] = dt1 * me_sum[he][ch];
                prob_joint_1[he][ch] = prob_joint[he][ch] + k1[he][ch];
            }
        }

        e = energy(prob_joint_1, hedges, M);
        pu_av = e / M;
        comp_pcond(prob_joint_1, pu_cond, pi, hedges, M, K, nch_fn);
        get_all_binom_sums(max_c, pu_av, binom_probs, binom_sums);

        der_walksat(nodes, hedges, prob_joint_1, pu_cond, binom_probs, binom_sums, 
                    M, K, nch_fn, q, e / N, me_sum, max_c);
        
        valid = true;
        for (long he = 0; he < M; he++){
            for (int ch = 0; ch < nch_fn; ch++){
                k2[he][ch] = dt1 * me_sum[he][ch];
                if (prob_joint[he][ch] + (k1[he][ch] + k2[he][ch]) / 2 < 0){
                    valid = false;
                }
            }
        }

        if (!valid){
            cout << "Some probabilities would be negative if dt=" << dt1 << " is taken" << endl;
            dt1 /= 2;
            cout << "step divided by half" << endl;
            if (dt1 < dt_min){
                dt_min /= 2;
                cout << "dt_min also halfed" << endl;
            }
        }else{
            error = 0;
            for (long he = 0; he < M; he++){
                for (int ch = 0; ch < nch_fn; ch++){
                    error += fabs(k1[he][ch] - k2[he][ch]);
                }
            }

            error /= nch_fn * M;

            if (error < 2 * tol){
                cout << "step dt=" << dt1 << "  accepted" << endl;
                cout << "error=" << error << endl;
                t += dt1;
                for (long he = 0; he < M; he++){
                    for (int ch = 0; ch < nch_fn; ch++){
                        prob_joint[he][ch] += (k1[he][ch] + k2[he][ch]) / 2;
                    }
                }
                e = energy(prob_joint, hedges, M);
                pu_av = e / M;
                fe << t << "\t" << e / N << endl;

            }else{
                e = energy(prob_joint, hedges, M);
                pu_av = e / M;
                cout << "step dt=" << dt1 << "  rejected  new step will be attempted" << endl;
                cout << "error=" <<  error << endl;
            }

            dt1 = 4 * dt1 * sqrt(2 * tol / error) / 5;
            if (dt1 > M){
                dt1 = M;
            }else if(dt1 < dt_min){
                dt1 = dt_min;
            }

            cout << "Recommended step is dt=" << dt1 << endl;
        }

        auto t2 = std::chrono::high_resolution_clock::now();

        auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

        cout << endl << "iteration time:   " << ms_int.count() << "ms" << endl; 

    }

    fe.close();

}


int main(int argc, char *argv[]) {
    long N = atol(argv[1]);
    long M = atol(argv[2]);
    int K = atoi(argv[3]);
    unsigned long seed_r = atol(argv[4]);
    double q = atof(argv[5]);
    double tl = atof(argv[6]);
    double tol = atof(argv[7]);
    int nthr = atoi(argv[8]);

    int nch_fn = (1 << K);
    double p0 = 0.5;

    omp_set_num_threads(nthr);

    Tnode *nodes;
    Thedge *hedges;

    gsl_rng * r;
    init_ran(r, seed_r);

    // char filegraph[300];
    // char filelinks[300];
    // sprintf(filegraph, "KSATgraph_K_%d_N_%li_M_%li_simetric_1_model_1_idum1_%li_J_1_ordered.txt", K, N, M, seed_g);
    // sprintf(filelinks, "KSAT_K_%d_enlaces_N_%li_M_%li_idumenlaces_%li_idumgraph_%li_ordered.txt", K, N, M, seed_g, seed_g);

    char fileener[300]; 
    sprintf(fileener, "CDA_WalkSAT_ener_K_%d_N_%li_M_%li_q_%.3lf_tl_%.2lf_seed_%li_tol_%.1e.txt", 
            K, N, M, q, tl, seed_r, tol);

    create_graph(N, M, K, nodes, hedges, r);
    int max_c = get_max_c(nodes, N);
    get_info_exc(nodes, hedges, N, M, K);

    
    RK2_walksat(nodes, hedges, N, M, K, nch_fn, q, max_c, p0, fileener, tl, tol);

    return 0;
}