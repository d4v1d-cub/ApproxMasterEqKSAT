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
    int *ch_unsat_exc;   // an array containing the partially unsat configuration of the other 
                         // nodes inside the clause if one removes one node.
    int **ch_exc;        // an array containing the configuration of the other 
                         // nodes inside the clause if one removes one node.
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
void read_graph(char *filegraph, long N, long M, int K, 
                Tnode *&nodes, Thedge *&hedges){
    
    init_graph(nodes, hedges, N, M);
    string trash_str;
    double trash_double;

    vector <long> nodes_in;
    for (int j = 0; j < K; j++){
        nodes_in.push_back(0);
    }

    ifstream fg(filegraph);

    long fn_count = 0;
    bool new_fn;

    for (long i = 0; i < N; i++){
        fg >> trash_double;
        fg >> nodes[i].nfacn;
        nodes[i].nch = (long) pow(2, nodes[i].nfacn);
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


long find_fn(vector <long> nodes_in, Thedge *hedges, vector <long> fn_list){
    bool cond = true;
    int w;
    int i = 0;
    while (i < fn_list.size() && cond){
        w = 0;
        while (cond && w < nodes_in.size()){
            if (nodes_in[w] != hedges[fn_list[i]].nodes_in[w]){
                cond = false;
            }
            w++;
        }
        i++;
        cond = !cond;
    }
    return fn_list[i - 1];
}


int find_node(vector <long> nodes_in, long node){
    bool cond = true;
    int i = 0;
    while (i < nodes_in.size() && cond){
        if (nodes_in[i] == node){
            cond = false;
        }
        i++;
    }
    return i - 1;
}


void read_graph_old_order(char *filegraph, long N, long M, int K, 
                          Tnode *&nodes, Thedge *&hedges){
    
    init_graph(nodes, hedges, N, M);
    string trash_str;
    double trash_double;

    vector <long> nodes_in;
    for (int j = 0; j < K; j++){
        nodes_in.push_back(0);
    }

    ifstream fg(filegraph);

    long fn_count = 0;
    long he;
    bool new_fn;
    long neigh;

    for (long i = 0; i < N; i++){
        fg >> trash_double;
        fg >> nodes[i].nfacn;
        nodes[i].nch = (long) pow(2, nodes[i].nfacn);
        getline(fg, trash_str);
        getline(fg, trash_str);

        for (int k = 0; k < nodes[i].nfacn; k++){
            nodes_in[0] = i;
            new_fn = true;
            for (int j = 0; j < K - 1; j++){
                fg >> trash_double;
                fg >> trash_double;
                fg >> nodes_in[j + 1];
                if (nodes_in[j + 1] < nodes_in[0]){
                    new_fn = false;
                    neigh = nodes_in[j + 1];
                }
                fg >> trash_double;
                fg >> trash_double;
            }
            sort(nodes_in.begin(), nodes_in.end());
            if (new_fn){
                nodes[i].fn_in.push_back(fn_count);
                hedges[fn_count].pos_n.push_back(nodes[i].fn_in.size() - 1);
                nodes[i].pos_fn.push_back(find_node(nodes_in, i));    
                for (int w = 0; w < K; w++){
                    
                    hedges[fn_count].nodes_in.push_back(nodes_in[w]);
                }
                fn_count++;
            }else{
                he = find_fn(nodes_in, hedges, nodes[neigh].fn_in);
                nodes[i].fn_in.push_back(he);
                nodes[i].pos_fn.push_back(find_node(nodes_in, i));
                hedges[he].pos_n.push_back(nodes[i].fn_in.size() - 1);
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


void get_info_exc(Tnode *nodes, Thedge *hedges, long N, long M, int K, int nch_fn){
    int w, count, ch_exc;
    bool bit;
    for (long he = 0; he < M; he++){
        hedges[he].nodes_exc = new long *[K];
        hedges[he].ch_unsat_exc = new int [K];
        hedges[he].ch_exc = new int *[K];
        for (int j = 0; j < K; j++){
            ch_exc = 0;
            hedges[he].nodes_exc[j] = new long [K - 1];
            count = 0;
            w = (j + 1) % K;
            while (w != j){
                hedges[he].nodes_exc[j][count] = hedges[he].nodes_in[w];
                bit = ((hedges[he].ch_unsat >> w) & 1);
                ch_exc += (bit << count);
                w = (w + 1) % K;
                count++;
            }
            hedges[he].ch_unsat_exc[j] = ch_exc;

            hedges[he].ch_exc[j] = new int [nch_fn]; 
            for (int ch = 0; ch < nch_fn; ch++){ // translation of the whole clause chain 'ch'
                ch_exc = 0;                      // into the chain that sees one of the variables inside 
                count = 0;
                w = (j + 1) % K;
                while (w != j){
                    bit = ((ch >> w) & 1);
                    ch_exc += (bit << count);
                    w = (w + 1) % K;
                    count++;
                }
                hedges[he].ch_exc[j][ch] = ch_exc;
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
void init_probs(double ****&pcav, double ***&pu_cav, double *&pi, double ****&cme_sum, 
                double *&me_sum, double ***&sums_save, long N, long M, int K, int nch_exc, 
                double p0){
    double prod;
    int bit;
    me_sum = new double [N];
    pi = new double [N];
    sums_save = new double **[N];
    for (long i = 0; i < N; i++){
        pi[i] = p0;
        sums_save[i] = new double *[2];
        for (int s = 0; s < 2; s++){
            sums_save[i][s] = new double[2];
        }
    }

    pcav = new double ***[M];
    cme_sum = new double ***[M];
    pu_cav = new double **[M];
    for (long he = 0; he < M; he++){
        pcav[he] = new double **[K];
        cme_sum[he] = new double **[K];
        pu_cav[he] = new double *[K];
        for (int l = 0; l < K; l++){
            pcav[he][l] = new double *[2];
            cme_sum[he][l] = new double *[2];
            pu_cav[he][l] = new double [2];
            for (int s = 0; s < 2; s++){
                pcav[he][l][s] = new double [nch_exc];
                cme_sum[he][l][s] = new double [nch_exc];
                for (int ch = 0; ch < nch_exc; ch++){
                    prod = 1;
                    for (int j = 0; j < K - 1; j++){
                        bit = ((ch >> j) & 1);
                        prod *= (bit + (1 - 2 * bit) * p0);
                    }
                    pcav[he][l][s][ch] = prod;
                } 
            }
        }
    }
}


// initializes the auxiliary arrays for the Runge-Kutta integration
void init_RK_arr(double ****&k1c, double ****&k2c, double ****&pcav_1,
                 double *&k1, double *&k2, double *&pi_1, long N, long M, int K, 
                 int nch_exc){
    k1 = new double [N];
    k2 = new double [N];
    pi_1 = new double [N];
    for (long i = 0; i < N; i++){
        k1[i] = 0;
        k2[i] = 0;
        pi_1[i] = 0;
    }


    k1c = new double ***[M];
    k2c = new double ***[M];
    pcav_1 = new double ***[M];
    for (long he = 0; he < M; he++){
        k1c[he] = new double **[K];
        k2c[he] = new double **[K];
        pcav_1[he] = new double **[K];
        for (int l = 0; l < K; l++){
            k1c[he][l] = new double *[2];
            k2c[he][l] = new double *[2];
            pcav_1[he][l] = new double *[2];
            for (int s = 0; s < 2; s++){
                k1c[he][l][s] = new double[nch_exc];
                k2c[he][l][s] = new double[nch_exc];
                pcav_1[he][l][s] = new double[nch_exc];
                for (int ch = 0; ch < nch_exc; ch++){
                    k1c[he][l][s][ch] = 0;
                    k2c[he][l][s][ch] = 0;
                    pcav_1[he][l][s][ch] = 0;
                }
            }
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
void fill_pneigh(int S, int *cj, double **binom_probs, double **binom_sums, 
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
double rate_walksat(int E0, int S, int K, double q, double e_av, 
                    int **cj, double **binom_probs, double **binom_sums, 
                    double **pneigh, int nch_exc){
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
            fill_pneigh(S, cj[fn], binom_probs, binom_sums, pneigh, K);
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
    return q * E0 / e_av / K + (1 - q) * cumul / e_av; 
}


void get_pu_cav(double ****pcav, double ***pu_cav, Thedge *hedges, long M, int K){
    for (long he = 0; he < M; he++){
        for (int w = 0; w < K; w++){
            for (int s = 0; s < 2; s++){
                pu_cav[he][w][s] = pcav[he][w][s][hedges[he].ch_unsat_exc[w]];
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
                 double ***pcav, double ***pu_cav, double **binom_probs, 
                 double **binom_sums, int K, int nch_fn, double q, 
                 double e_av, double ***cme_sum_src, int max_c){
    int he, plc_he, plc_other;
    bool bit, uns, uns_flip, bit_other;
    int ch_flip, ch_exc, ch_exc_flip;
    
    double prod[2], r[2][2];
    int E[2];
    double **pneigh = new double *[K - 1];
    for (int j = 0; j < K - 1; j++){
        pneigh[j] = new double [2];
    }

    int ***cj = new int **[2];
    for (int s = 0; s < 2; s++){
        cj[s] = new int *[nodes[node].nfacn];
        for (int h = 0; h < nodes[node].nfacn; h++){
            cj[s][h] = new int [K - 1];
        }
    }
    

    for (long ch = 0; ch < nodes[node].nch / 2; ch++){
        prod[0] = prodcond(pu_cav, fn_src, nodes[node], 0, ch);
        prod[1] = prodcond(pu_cav, fn_src, nodes[node], 1, ch);
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
                            binom_sums, pneigh, nch_fn / 2);
        
        for (int ch_src = 0; ch_src < nch_fn; ch_src++){
            bit = ((ch_src >> plc_he) & 1);
            ch_flip = (ch_src ^ (1 << plc_he));
            uns = (ch_src == hedges[he].ch_unsat);
            uns_flip = (ch_flip == hedges[he].ch_unsat);
            for (int j = 0; j < K - 1; j++){
                plc_other = (plc_he + j + 1) % K;
                bit_other = ((ch_src >> plc_other) & 1); 
                ch_exc = hedges[he].ch_exc[plc_other][ch_src];
                ch_exc_flip = hedges[he].ch_exc[plc_other][ch_flip];
                cme_sum_src[plc_other][bit_other][ch_exc] += 
                    -r[bit][uns] * prod[bit] * pcav[plc_other][bit_other][ch_exc] + 
                    r[1 - bit][uns_flip] * prod[1 - bit] * pcav[plc_other][bit_other][ch_exc_flip];
            }
        }
    }

    for (int j = 0; j < K - 1; j++){
        delete []pneigh[j];
    }
    delete [] pneigh;
    pneigh = NULL;
    

    for (int s = 0; s < 2; s++){
        for (int h = 0; h < nodes[node].nfacn; h++){
            delete [] cj[s][h];
        }
        delete [] cj[s];
    }
    delete []cj;
    cj = NULL;

}



// it does the sum in the derivative of the CDA equations
// fn_src is the origin factor node where one is computing the derivative
// part_uns is 1 if the other variables in fn_src are partially
// unsatisfying their links, and is 0 otherwise. 
void sum_walksat(long node, int fn_src, Tnode *nodes, Thedge *hedges, 
                 double ***pcav, double ***pu_cav, double **binom_probs, 
                 double **binom_sums, int K, int nch_fn, double q, 
                 double e_av, double ***cme_sum_src, int max_c, double **sums_save){
    int he, plc_he, plc_other;
    bool bit, uns, uns_flip, bit_other;
    int ch_flip, ch_exc, ch_exc_flip;
    long neigh;
    
    double prod[2], r[2][2];
    int E[2];
    double **pneigh = new double *[K - 1];
    for (int j = 0; j < K - 1; j++){
        pneigh[j] = new double [2];
    }

    int ***cj = new int **[2];
    for (int s = 0; s < 2; s++){
        cj[s] = new int *[nodes[node].nfacn];
        for (int h = 0; h < nodes[node].nfacn; h++){
            cj[s][h] = new int [K - 1];
        }
    }
    
    sums_save[0][0] = 0;
    sums_save[1][0] = 0;
    sums_save[0][1] = 0;
    sums_save[1][1] = 0;

    for (long ch = 0; ch < nodes[node].nch / 2; ch++){
        prod[0] = prodcond(pu_cav, fn_src, nodes[node], 0, ch);
        prod[1] = prodcond(pu_cav, fn_src, nodes[node], 1, ch);
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

        sums_save[0][0] += r[0][0] * prod[0];
        sums_save[1][0] += r[1][0] * prod[1];

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
                            binom_sums, pneigh, nch_fn / 2);

        sums_save[bit][1] += r[bit][1] * prod[bit];
        
        for (int ch_src = 0; ch_src < nch_fn; ch_src++){
            bit = ((ch_src >> plc_he) & 1);
            ch_flip = (ch_src ^ (1 << plc_he));
            uns = (ch_src == hedges[he].ch_unsat);
            uns_flip = (ch_flip == hedges[he].ch_unsat);
            for (int j = 0; j < K - 1; j++){
                plc_other = (plc_he + j + 1) % K;
                bit_other = ((ch_src >> plc_other) & 1); 
                ch_exc = hedges[he].ch_exc[plc_other][ch_src];
                ch_exc_flip = hedges[he].ch_exc[plc_other][ch_flip];
                cme_sum_src[plc_other][bit_other][ch_exc] += 
                    -r[bit][uns] * prod[bit] * pcav[plc_other][bit_other][ch_exc] + 
                    r[1 - bit][uns_flip] * prod[1 - bit] * pcav[plc_other][bit_other][ch_exc_flip];
            }
        }
    }

    for (int j = 0; j < K - 1; j++){
        delete []pneigh[j];
    }
    delete [] pneigh;
    pneigh = NULL;
    

    for (int s = 0; s < 2; s++){
        for (int h = 0; h < nodes[node].nfacn; h++){
            delete [] cj[s][h];
        }
        delete [] cj[s];
    }
    delete []cj;
    cj = NULL;

}


// it takes the derivative for a single node with the auxiliary sums already computed
double der_single_node(double pi, double **sums_save, double *pu_cav, bool bit_uns){
    double der = 0;
    bool uns, uns_flip;
    for (int part_uns = 0; part_uns < 2; part_uns++){
        uns = part_uns * (1 - bit_uns);
        uns_flip = part_uns * bit_uns;
        der += -sums_save[0][uns] * (1 - part_uns - (1 - 2 * part_uns) * pu_cav[0]) * pi + 
               sums_save[1][uns_flip] * (1 - part_uns - (1 - 2 * part_uns) * pu_cav[1]) * (1 - pi);
    }
    return der;
}



// it computes all the derivatives of the joint probabilities
void der_walksat(Tnode *nodes, Thedge *hedges, double ****pcav, double ***pu_cav, double *pi, 
                 double **binom_probs, double **binom_sums, long N, long M, int K, int nch_fn, 
                 double q, double e_av, double ****cme_sum, double *me_sum, int max_c, 
                 double ***sums_save){
    for (long he = 0; he < M; he++){
        for (int w = 0; w < K; w++){
            for (int s = 0; s < 2; s++){
                for (int ch = 0; ch < nch_fn / 2; ch++){
                    cme_sum[he][w][s][ch] = 0;
                }
            }
        }     
    }

    // candidate to be a parallel for
    #pragma omp parallel for
    for (long he = 0; he < M; he++){
        for (int w = 0; w < K; w++){
            if (hedges[he].pos_n[w] == 0){
                sum_walksat(hedges[he].nodes_in[w], hedges[he].pos_n[w], nodes, hedges, 
                            pcav[he], pu_cav, binom_probs, binom_sums, K, nch_fn,
                            q, e_av, cme_sum[he], max_c, sums_save[hedges[he].nodes_in[w]]);
            }else{
                sum_walksat(hedges[he].nodes_in[w], hedges[he].pos_n[w], nodes, hedges, 
                            pcav[he], pu_cav, binom_probs, binom_sums, K, nch_fn,
                            q, e_av, cme_sum[he], max_c);
            }
            
        }
    }

    bool bit;
    for (long i = 0; i < N; i++){
        bit = ((hedges[nodes[i].fn_in[0]].ch_unsat >> nodes[i].pos_fn[0]) & 1);  
        me_sum[i] = der_single_node(pi[i], sums_save[i], 
                                    pu_cav[nodes[i].fn_in[0]][nodes[i].pos_fn[0]], bit);
    }

}


double energy(double ***pu_cav, double *pi, Thedge *hedges, long M){
    double e = 0;
    bool bit;
    for (long he = 0; he < M; he++){
        bit = (hedges[he].ch_unsat & 1);
        e += pu_cav[he][0][bit] * (bit + (1 - 2 * bit) * pi[hedges[he].nodes_in[0]]);
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
void RK2_walksat(Tnode *nodes, Thedge *hedges, long N, long M, int K, int nch_fn, double q, 
                 int max_c, double p0, char *fileener, double tl, double tol = 1e-2, 
                 double t0 = 0, double dt0 = 0.01, double ef = 1e-6, double dt_min = 1e-7){
    double **binom_probs, **binom_sums;
    double ****pcav, ***pu_cav, ****cme_sum, *pi, *me_sum, ***sums_save;
    double e, pu_av, error;                 
                 
    // initalize all arrays that will be used inside the derivative
    init_aux_arr(binom_probs, binom_sums, max_c, K);
    init_probs(pcav, pu_cav, pi, cme_sum, me_sum, sums_save, N, M, K, nch_fn / 2, p0);



    // initialize auxiliary arrays for the Runge-Kutta integration
    double ****k1c, ****k2c, ****pcav_1, *k1, *k2, *pi_1;
    init_RK_arr(k1c, k2c, pcav_1, k1, k2, pi_1, N, M, K, nch_fn / 2);

    ofstream fe(fileener);
    
    get_pu_cav(pcav, pu_cav, hedges, M, K);
    e = energy(pu_cav, pi, hedges, M);
    pu_av = e / M;
    fe << t0 << "\t" << e / N << endl;   // it prints the energy density

    double dt1 = dt0;
    double t = t0;

    bool valid;

    // the time scale is already given in Monte Carlo steps. Inside the rates I am using 
    // the energy density e_av
    while (t < tl){
        if (e / N < ef){
            cout << "Final energy reached" << endl;
            break;
        }

        auto t1 = std::chrono::high_resolution_clock::now();

        get_all_binom_sums(max_c, pu_av, binom_probs, binom_sums);

        der_walksat(nodes, hedges, pcav, pu_cav, pi, binom_probs, binom_sums, N, M, K, nch_fn, 
                    q, e / N, cme_sum, me_sum, max_c, sums_save);   // in the rates, I use the energy density

        valid = true;
        for (long he = 0; he < M; he++){
            for (int w = 0; w < K; w++){
                for (int s = 0; s < 2; s++){
                    for (int ch = 0; ch < nch_fn / 2; ch++){
                        k1c[he][w][s][ch] = dt1 * cme_sum[he][w][s][ch];
                        pcav_1[he][w][s][ch] = pcav[he][w][s][ch] + k1c[he][w][s][ch];
                        if (pcav_1[he][w][s][ch] < 0){
                            valid = false;
                        }
                    }
                }
            }
        }

        if (valid){
            for (long i = 0; i < N; i++){
                k1[i] = dt1 * me_sum[i];
                pi_1[i] = pi[i] + k1[i];
                if (pi_1[i] < 0){
                    valid = false;
                }
            }
        }

        while (!valid){
            cout << "some probabilities became negative in the auxiliary step of RK2" << endl;
            dt1 /= 2;
            cout << "step divided by half    dt=" << dt1 << endl;
            if (dt1 < dt_min){
                dt_min /= 2;
                cout << "dt_min also halfed" << endl;
            }

            valid = true;
            for (long he = 0; he < M; he++){
                for (int w = 0; w < K; w++){
                    for (int s = 0; s < 2; s++){
                        for (int ch = 0; ch < nch_fn / 2; ch++){
                            k1c[he][w][s][ch] = dt1 * cme_sum[he][w][s][ch];
                            pcav_1[he][w][s][ch] = pcav[he][w][s][ch] + k1c[he][w][s][ch];
                            if (pcav_1[he][w][s][ch] < 0){
                                valid = false;
                            }
                        }
                    }
                }
            }

            if (valid){
                for (long i = 0; i < N; i++){
                    k1[i] = dt1 * me_sum[i];
                    pi_1[i] = pi[i] + k1[i];
                    if (pi_1[i] < 0){
                        valid = false;
                    }
                }
            }
        }

        get_pu_cav(pcav_1, pu_cav, hedges, M, K);
        e = energy(pu_cav, pi_1, hedges, M);
        pu_av = e / M;
        get_all_binom_sums(max_c, pu_av, binom_probs, binom_sums);

        der_walksat(nodes, hedges, pcav_1, pu_cav, pi_1, binom_probs, binom_sums, N, M, K, nch_fn, 
                    q, e / N, cme_sum, me_sum, max_c, sums_save);

        valid = true;
        for (long he = 0; he < M; he++){
            for (int w = 0; w < K; w++){
                for (int s = 0; s < 2; s++){
                    for (int ch = 0; ch < nch_fn / 2; ch++){
                        k2c[he][w][s][ch] = dt1 * cme_sum[he][w][s][ch];
                        if (pcav[he][w][s][ch] + (k1c[he][w][s][ch] + k2c[he][w][s][ch]) / 2 < 0){
                            valid = false;
                        }
                    }
                }
            }
        }

        if (valid){
            for (long i = 0; i < N; i++){
                k2[i] = dt1 * me_sum[i];
                if (pi[i] + (k1[i] + k2[i]) / 2 < 0){
                    valid = false;
                }
            }
        }

        if (!valid){
            cout << "Some probabilities would be negative if dt=" << dt1 << " is taken" << endl;
            dt1 /= 2;
            cout << "step divided by half    dt=" << dt1  << endl;
            if (dt1 < dt_min){
                dt_min /= 2;
                cout << "dt_min also halfed" << endl;
            }
            get_pu_cav(pcav, pu_cav, hedges, M, K);
            e = energy(pu_cav, pi, hedges, M);
            pu_av = e / M;
        }else{
            error = 0;
            for (long he = 0; he < M; he++){
                for (int w = 0; w < K; w++){
                    for (int s = 0; s < 2; s++){
                        for (int ch = 0; ch < nch_fn / 2; ch++){
                            error += fabs(k2c[he][w][s][ch] - k1c[he][w][s][ch]);
                        }
                    }
                }
            }

            for (long i = 0; i < N; i++){
                error += fabs(k2[i] - k1[i]);
            }

            error /= (N + M * K * nch_fn);

            if (error < 2 * tol){
                cout << "step dt=" << dt1 << "  accepted" << endl;
                cout << "error=" << error << endl;
                t += dt1;

                for (long he = 0; he < M; he++){
                    for (int w = 0; w < K; w++){
                        for (int s = 0; s < 2; s++){
                            for (int ch = 0; ch < nch_fn / 2; ch++){
                                pcav[he][w][s][ch] += (k1c[he][w][s][ch] + k2c[he][w][s][ch]) / 2;
                            }
                        }
                    }
                }

                for (long i = 0; i < N; i++){
                    pi[i] += (k1[i] + k2[i]) / 2;
                }


                get_pu_cav(pcav, pu_cav, hedges, M, K);
                e = energy(pu_cav, pi, hedges, M);
                pu_av = e / M;
                fe << t << "\t" << e / N << endl;

            }else{
                get_pu_cav(pcav, pu_cav, hedges, M, K);
                e = energy(pu_cav, pi, hedges, M);
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
    // sprintf(filegraph, "KSATgraph_K_%d_N_%li_M_%li_simetric_1_model_1_idum1_-2_J_1_ordered.txt", 
    //                    K, N, M);
    // sprintf(filelinks, "KSAT_K_%d_enlaces_N_%li_M_%li_idumenlaces_-2_idumgraph_-2_ordered.txt", 
    //                    K, N, M);

    char fileener[300]; 
    sprintf(fileener, "CME_WalkSAT_ener_K_%d_N_%li_M_%li_q_%.4lf_tl_%.2lf_seed_%li_tol_%.1e.txt", 
            K, N, M, q, tl, seed_r, tol);

    create_graph(N, M, K, nodes, hedges, r);
    // read_graph_old_order(filegraph, N, M, K, nodes, hedges);
    // read_links(filelinks, N, M, K, nodes, hedges);
    int max_c = get_max_c(nodes, N);
    get_info_exc(nodes, hedges, N, M, K, nch_fn);

    
    RK2_walksat(nodes, hedges, N, M, K, nch_fn, q, max_c, p0, fileener, tl, tol);

    return 0;
}