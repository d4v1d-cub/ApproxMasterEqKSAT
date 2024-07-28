#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <cmath>

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
        nodes[i].nch = (1 << nodes[i].nfacn);
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


// rate of the Focused Metropolis Search algorithm.
double rate_fms(int E0, int E1, int K, double eta, double Eav){
    double dE = E1 - E0;
    if (dE > 0){
        return E0 / K / Eav * pow(eta, -dE);
    }else{
        return E0 / K / Eav;
    }
}


// initializes the arrays where the binomial weights will be stored. Those arrays are used 
// to compute the probability of selecting the variable belonging to the smallest number 
// of satisfied clauses
void init_aux_arr(double **&binom_probs, double **&binom_sums, double **&pneigh, int ***&cj, 
                  int max_c, int K){
    binom_probs = new double *[max_c];
    for (int gj = 0; gj < max_c; gj++){
        binom_probs[gj] = new double[gj + 1];
    }

    binom_sums = new double *[max_c + 1];
    for (int s = 0; s < max_c + 1; s++){
        binom_sums[s] = new double[max_c];
    }
    
    pneigh = new double *[K - 1];
    for (int j = 0; j < K - 1; j++){
        pneigh[j] = new double[2];
    }

    cj = new int **[2];
    for (int s = 0; s < 2; s++){
        cj[s] = new int *[max_c + 1];
        for (int h = 0; h < max_c + 1; h++){
            cj[s][h] = new int[K - 1];
        }
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
    
    for (int si = 0; si < max_c + 1; si++){
        for (int gj = 0; gj < max_c; gj++){
            binom_sums[si][gj] = gsl_cdf_binomial_Q(si, 1 - pu_av, gj);
        }
    }
}


// it fills the array of the probabilities to be used when computing the walksat rate
void fill_pneigh(int E0, int S, int *cj, double **binom_probs, double **binom_sums, 
                 double **pneigh, int K){
    for (int j = 0; j < K - 1; j++){
        pneigh[j][0] = binom_probs[cj[j] - 1][S];
        pneigh[j][1] = binom_sums[S][cj[j] - 1];
    }
}


// rate of the walksat algorithm used by Barthel et al. in 2003
// cj is a list of all the connectivities of the neighbors of node i that are 
// in unsatisfied clauses
double rate_walksat(int E0, int S, int K, double q, double Eav, int **cj, 
                    double **binom_probs, double **binom_sums, double **pneigh, 
                    int nchain){
    bool cond;
    double cumul, prod;
    int cumul_bits, bit, k;
    cumul = 0;
    for (int fn = 0; fn < E0; fn++){
        cond = true;
        k = 0;
        while (k < K - 1 && cond){
            if (cj[fn][k] - 1 < S){
                cond = false;
            }
            k++;
        }
        if (cond){   
            fill_pneigh(E0, S, cj[fn], binom_probs, binom_sums, pneigh, K);
            for (int ch = 0; ch < nchain; ch++){
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
    return q * E0 / Eav / K + (1 - q) * cumul / Eav / K; 
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
void sum_walksat(long node, int fn_src, int part_uns, Tnode *nodes, Thedge *hedges, 
                   double *prob_joint, double ***pu_cond, int ***cj, double **binom_probs, 
                   double **binom_sums, double **pneigh, int K, int nch_fn, double q, 
                   double Eav, double *me_sum_src){
    
    int E[2], he, plc_he;
    double r[2], prod[2];
    int plc_src = nodes[node].pos_fn[fn_src];
    bool bit;
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
                // cond = 0 (false), otherwise cond = 1.
                for (int h = 0; h < K - 1; h++){
                    cj[bit][E[bit]][h] = 
                            nodes[hedges[he].nodes_exc[plc_src][h]].nfacn;
                    // It's the connectivity of one of the other nodes inside the factor node: 
                    // 'he = nodes[node].fn_in[hind]'
                }
                E[bit]++;
            }
        }

        r[0] = rate_walksat(E[0], nodes[node].nfacn - E[0], K, q, Eav, cj[0], binom_probs, 
                            binom_sums, pneigh, nch_fn);
        r[1] = rate_walksat(E[1], nodes[node].nfacn - E[1], K, q, Eav, cj[1], binom_probs, 
                            binom_sums, pneigh, nch_fn);
        
        for (int ch_src = 0; ch_src < nch_fn; ch_src++){
            bit = ((ch_src >> plc_src) & 1);
            me_sum_src[ch_src] += -r[bit] * prod[bit] * prob_joint[ch_src] + 
                                  r[1 - bit] * prod[1 - bit] * prob_joint[ch_src ^ (1 << plc_src)];
        }
    }
}


int main(int argc, char *argv[]) {
    long N = atol(argv[1]);
    long M = atol(argv[2]);
    int K = atoi(argv[3]);
    long seed_g = atol(argv[4]);
    unsigned long seed_r = atol(argv[5]);
    double q = atof(argv[6]);

    int nch_fn = (1 << K);
    double p0 = 0.5;

    Tnode *nodes;
    Thedge *hedges;
    double **binom_probs, **binom_sums, **pneigh;
    int ***cj;
    double **prob_joint, ***pu_cond, **me_sum, **pi;

    gsl_rng * r;
    init_ran(r, seed_r);

    // char filegraph[300];
    // char filelinks[300];
    // sprintf(filegraph, "KSATgraph_K_%d_N_%li_M_%li_simetric_1_model_1_idum1_%li_J_1_ordered.txt", K, N, M, seed_g);
    // sprintf(filelinks, "KSAT_K_%d_enlaces_N_%li_M_%li_idumenlaces_%li_idumgraph_%li_ordered.txt", K, N, M, seed_g, seed_g);

    create_graph(N, M, K, nodes, hedges, r);
    int max_c = get_max_c(nodes, N);
    get_info_exc(nodes, hedges, N, M, K);

    init_aux_arr(binom_probs, binom_sums, pneigh, cj, max_c, K);
    double pu_av = 0.3;
    double Eav = pu_av * M;
    get_all_binom_sums(max_c, pu_av, binom_probs, binom_sums);

    int he_src = nodes[10].fn_in[0];
    init_probs(prob_joint, pu_cond, pi, me_sum, M, K, nch_fn, p0);
    comp_pcond(prob_joint, pu_cond, pi, hedges, M, K, nch_fn);

    sum_walksat(10, 0, 1, nodes, hedges, prob_joint[he_src], pu_cond, cj, binom_probs, 
                binom_sums, pneigh, K, nch_fn, q, Eav, me_sum[he_src]);

    return 0;
}