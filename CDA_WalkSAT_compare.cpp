#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <cmath>
#include <omp.h>
#include <chrono>
#include "CDA_WalkSAT_bib_ser_compare.h"

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


vector <long> check_fn(long M, Thedge *hedges, int K){
    vector <long> answ;
    for (long he = 0; he < M; he++){
        if(hedges[he].pos_n.size() < K){
            answ.push_back(he);
        }
    }
    return answ;
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



// it computes the conditional probabilities of having a partially unsatisfied clause, given the 
// value of one variable in the clause
void comp_pcond(double **prob_joint, double ***pu_cond, double **pi, Thedge *hedges, long M, int K, 
                int nch_fn){
    double pu;
    int bit;
    int ch_uns_flip;
    for (long he = 0; he < M; he++){
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
            bit = ((hedges[he].ch_unsat >> w) & 1); 
            ch_uns_flip = (hedges[he].ch_unsat ^ (1 << w));
            pu_cond[he][w][bit] = prob_joint[he][hedges[he].ch_unsat] / pi[w][bit];
            pu_cond[he][w][1 - bit] = prob_joint[he][ch_uns_flip] / pi[w][1 - bit];
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
                            binom_sums, pneigh, nch_fn / 2);
        
        for (int ch_src = 0; ch_src < nch_fn; ch_src++){
            bit = ((ch_src >> plc_he) & 1);
            ch_flip = (ch_src ^ (1 << plc_he));
            uns = (ch_src == hedges[he].ch_unsat);
            uns_flip = (ch_flip == hedges[he].ch_unsat);
            me_sum_src[ch_src] += -r[bit][uns] * prod[bit] * prob_joint[ch_src] + 
                                  r[1 - bit][uns_flip] * prod[1 - bit] * prob_joint[ch_flip];
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


double produccond(long long cadenaoffn, int spin, int length, int lugar, Tred * red, int factornode,
                  double ** probcond, int numcadenasxfacnode) {
    //produccond hace el producto de todas las probabilidades que condicionan en el vecino excepto la del factor node donde está "lugar"
    //la cadenaoffn viene en el mismo orden de los factor nodes en Tred * red, quitando factornode. Es importante aclarar
    //que las cifras menos significativas corresponden a los ultimos fator nodes, por eso cadenaoffn >> (length - w -1).
    //Aqui la cadena es de los factor nodes, y las probabilidades son las de SAT-UNSAT
    double prod = 1;
    int elemento;
    int h = 1;
    int w;
    int indice = (factornode + h) % red[lugar].numfactornodes;
    while (indice != factornode) {
        w = (factornode + h - 1) % length;
        elemento = int((cadenaoffn >> (length - w - 1)) & 1);
        prod *= (1 - elemento - (1 - 2* elemento) * probcond[indice][red[lugar].unsatisfied[indice] + spin * numcadenasxfacnode]);
        indice = (indice + 1) % red[lugar].numfactornodes;
        h++;
    }
    return prod;
}



double energy(int M, double ** probjoint, int *unsat){
    //Calcula E usando solo las probabilidades de los sitios unsatisfied

    double ener  = 0;

    for (int i = 0; i < M; i++){
        ener += probjoint[i][unsat[i]];
    }
    return ener;
}


void init_comb_pot_alphas(int maxconnect, int **&comb, double *&pot_alphas){
    comb=new int*[2];
    for(int k=0;k<2;k++)
    {
        comb[k]=new int[maxconnect+1];
        comb[k][0]=1;
    }
    pot_alphas=new double[maxconnect+1];
    pot_alphas[0]=1;
}

double rate_walksat_greedy_single_clause(vector <int> connectivities, double alpha_SAT, double alpha_UNSAT,
                                         int num_sat_clauses_i, int origin_clause_value, double Q, int **comb,
                                         double *pot_alphas, double *P_igual, double *P_mayor)
{
    if(!origin_clause_value)
        return 0;
    if(num_sat_clauses_i>=connectivities[0]||num_sat_clauses_i>=connectivities[1])
        return Q/3;

    for(int k=0;k<2;k++)
    {
        for(int s=1;s<=num_sat_clauses_i;s++)         // This seems to give comb[k][s] = (c_k - 2)! / (c_k - 2 - S)! / (S + 1)!
            comb[k][s]=comb[k][s-1]*(connectivities[k]-s)/s;   // comb[k][s] = binom(c_k - 1, S) * (c_k - 1 - S) / (S + 1) / (c_k - 1)
    }


    for(int s=1;s<=num_sat_clauses_i;s++)
        pot_alphas[s]=pot_alphas[s-1]*alpha_SAT/alpha_UNSAT;


    for(int k=0;k<2;k++)
    {
        P_igual[k]= (double) (comb[k][num_sat_clauses_i] * pot_alphas[num_sat_clauses_i] * pow(alpha_UNSAT, connectivities[k] - 1));
        P_mayor[k]=0;
        for(int s=0;s<=num_sat_clauses_i;s++)    // I think it has to be <=, not <
            P_mayor[k]+=comb[k][s]*pot_alphas[s];
        P_mayor[k]= (double) (1 - P_mayor[k] * pow(alpha_UNSAT, connectivities[k] - 1));
    }

//    return Q/3+(1-Q)*(P_mayor[1]*P_mayor[2]+(P_mayor[1]*P_igual[2]+P_mayor[2]*P_igual[1])/2+P_igual[1]*P_igual[2]/3);
    return Q/3+(1-Q)*(P_mayor[0]*P_mayor[1]+(P_mayor[0]*P_igual[1]+P_mayor[1]*P_igual[0])/2+P_igual[0]*P_igual[1]/3);
}


double total_rate_walksat_greedy(Tred * red, double alpha_SAT, double alpha_UNSAT, int num_sat_clauses_i, int lugar,
                                 long long cadenaoffn, int spin, int length, int factornode,
                                 int numvecxfacnode, int *** lugvecinoswfn, int elemfacnode, double Q, int **comb,
                                 double *pot_alphas, double *P_igual, double *P_mayor){
    double total_rate = 0;
    int origin_clause_value;
    vector <int> connectivities;
    for (int k = 0; k < numvecxfacnode; k++)
    {
        connectivities.push_back(red[lugvecinoswfn[lugar][factornode][k]].numfactornodes);
    }
    origin_clause_value = elemfacnode * (1 - red[lugar].valorenlaces[factornode] * (1 - 2 * spin)) / 2;
    total_rate += rate_walksat_greedy_single_clause(connectivities, alpha_SAT, alpha_UNSAT, num_sat_clauses_i,
                                                    origin_clause_value, Q, comb, pot_alphas, P_igual, P_mayor);

    int elemento;
    int h = 1;
    int w;
    int indice = (factornode + h) % red[lugar].numfactornodes;
    while (indice != factornode) {
        w = (factornode + h - 1) % length;
        elemento = int((cadenaoffn >> (length - w - 1)) & 1);
        for (int k = 0; k < numvecxfacnode; k++)
        {
            connectivities[k] = red[lugvecinoswfn[lugar][indice][k]].numfactornodes;
        }
        origin_clause_value = elemento * (1 - red[lugar].valorenlaces[indice] * (1 - 2 * spin)) / 2;
        total_rate += rate_walksat_greedy_single_clause(connectivities, alpha_SAT, alpha_UNSAT, num_sat_clauses_i,
                                                        origin_clause_value, Q, comb, pot_alphas, P_igual, P_mayor);
        indice = (indice + 1) % red[lugar].numfactornodes;
        h++;
    }
    return total_rate;
}


void ComputeProbcond(Tred * red, double ***probcond, double **probjoint, int N, int numcadenasxfactornode){
    double cumul_up, cumul_down;
    for (int j = 0; j < N; j++){
        for (int h = 0; h < red[j].numfactornodes; h++){
            cumul_up = 0, cumul_down = 0;
            for (int a = 0; a < numcadenasxfactornode; a++) {
                probcond[j][h][a] = probjoint[red[j].lugfacnodes[h]][red[j].cadfacnode[h][a]];
                cumul_up += probjoint[red[j].lugfacnodes[h]][red[j].cadfacnode[h][a]];
            }

            for (int a = numcadenasxfactornode; a < 2 * numcadenasxfactornode; a++){
                probcond[j][h][a] = probjoint[red[j].lugfacnodes[h]][red[j].cadfacnode[h][a]];
                cumul_down += probjoint[red[j].lugfacnodes[h]][red[j].cadfacnode[h][a]];
            }

            for (int a = 0; a < numcadenasxfactornode; a++) {
                probcond[j][h][a] /= cumul_up;
            }

            for (int a = numcadenasxfactornode; a < 2 * numcadenasxfactornode; a++){
                probcond[j][h][a] /= cumul_down;
            }
        }
    }
}


void UpdateCDAarrays(Tred * red, double *** probcond, double ** probjoint,
                     double ** mejoint, int **** tablaLocaleners, int *** lugvecinoswfn,
                     vector < vector <vector <int> > > &places,
                     int numcadenasxfactornode, int K, double E,
                     int q, int w, long long j, long long g, int M, double Q, int **comb,
                     double *pot_alphas, double *P_igual, double *P_mayor) {
//Esta función se encarga de añadir a los acumuladores "me" y "cme" las partes de las derivadas que
// corresponden a los nodos en places[q][w] y a la cadena formada por "j" y "g"

    int Elugar0 = 0, Elugar1 = 0, Efacnode0 = 0, Efacnode1 = 0, E0, E1,
            lugar, cad, cadfacnode, cadfacnode1;
    double rate0, rate1, prodcond0 = 0, prodcond1 = 0;
    long long cadena;

    cadena = g + (j <<
                  red[places[q][w][0]].numfacnodepartless[1]); // Se rearma la cadena. De izquierda a derecha se lee "jg"

    for (int k = 0; k < q; k++) {  // Va por todos los factor nodes
        Efacnode0 = (1 - red[places[q][w][0]].valorenlaces[k]) / 2; // Todas las energias se calculan para spin=0(1)
        Efacnode1 = (1 + red[places[q][w][0]].valorenlaces[k]) / 2; // Estas son las correspondientes al factor node q

        Elugar0 = tablaLocaleners[red[places[q][w][0]].numfacnodepartless[0]] // Este es el resto de la energía
                  [red[places[q][w][0]].numenlpospartless[k][0]][0][j] +
                  tablaLocaleners[red[places[q][w][0]].numfacnodepartless[1]]
                  [red[places[q][w][0]].numenlpospartless[k][1]][0][g];
        Elugar1 = tablaLocaleners[red[places[q][w][0]].numfacnodepartless[0]]
                  [red[places[q][w][0]].numenlpospartless[k][0]][1][j] +
                  tablaLocaleners[red[places[q][w][0]].numfacnodepartless[1]]
                  [red[places[q][w][0]].numenlpospartless[k][1]][1][g];


        for (int h = 0; h < numcadenasxfactornode; h++) {
            for (int e = 0; e < places[q][w].size(); e++) {
                lugar = places[q][w][e]; //Se va entonces por todos los lugares de la clase

                cad = (h + red[lugar].unsatisfied[k]) % numcadenasxfactornode; //se calcula la cadena
                E0 = Elugar0 + Efacnode0 * red[lugar].elemcad[k][cad];
                E1 = Elugar1 + Efacnode1 * red[lugar].elemcad[k][cad];

                rate0 = total_rate_walksat_greedy(red, 1 - E / M, E / M, red[lugar].numfactornodes - E0, lugar,
                                                  cadena, 0, red[lugar].numfactornodes - 1, k, K - 1, lugvecinoswfn,
                                                  red[lugar].elemcad[k][cad], Q, comb, pot_alphas,
                                                  P_igual, P_mayor) / E;

                rate1 = total_rate_walksat_greedy(red, 1 - E / M, E / M, red[lugar].numfactornodes - E1, lugar,
                                                  cadena, 1, red[lugar].numfactornodes - 1, k, K - 1, lugvecinoswfn,
                                                  red[lugar].elemcad[k][cad], Q, comb, pot_alphas,
                                                  P_igual, P_mayor) / E;

                prodcond0 = produccond(cadena, 0,  // Se calcula el producto de los factor noes vecinos excepto el k
                                       red[lugar].numfactornodes - 1,
                                       lugar, red, k, probcond[lugar], numcadenasxfactornode);
                prodcond1 = produccond(cadena, 1,
                                       red[lugar].numfactornodes - 1,
                                       lugar, red, k, probcond[lugar], numcadenasxfactornode);

                cadfacnode = red[lugar].cadfacnode[k][cad];
                cadfacnode1 = red[lugar].cadfacnode[k][(cad + numcadenasxfactornode) %
                                                       (2 * numcadenasxfactornode)];

                mejoint[red[lugar].lugfacnodes[k]][cadfacnode] +=
                            -rate0 * prodcond0 * probjoint[red[lugar].lugfacnodes[k]][cadfacnode];
                mejoint[red[lugar].lugfacnodes[k]][cadfacnode] +=
                            rate1 * prodcond1 * probjoint[red[lugar].lugfacnodes[k]][cadfacnode1];


                mejoint[red[lugar].lugfacnodes[k]][cadfacnode1] +=
                            -rate1 * prodcond1 * probjoint[red[lugar].lugfacnodes[k]][cadfacnode1];
                mejoint[red[lugar].lugfacnodes[k]][cadfacnode1] +=
                            rate0 * prodcond0 * probjoint[red[lugar].lugfacnodes[k]][cadfacnode];
            }
        }
    }
}



void CMEordering(Tred * red, double *** probcond, double ** probjoint,
                 double ** mejoint,
                 int **** tablaLocaleners, int *** lugvecinoswfn,
                 vector < vector <vector <int> > > &places,
                 int numcadenasxfactornode, int K, double E,
                 int M, double Q, int **comb,double *pot_alphas, double *P_igual, double *P_mayor) {
// Calcula los updates a probcond y prob
// Le llamé ordering porque usa un reordenamiento del grafo para almacenar las energias locales y calcular las sumas en la CME. La esencia
// está en agrupar a los nodos en clases que se identifican por cantidades de factor nodes y de enlaces positivos.


    for (int q = 1; q < places.size(); q++) { // Iteramos por todas las clases que contiene places
        for (int w = 0; w < places[q].size(); w++) {
            for (long long j = 0; j < red[places[q][w][0]].cantcadfacnodepartless[0]; j++) {
                for (long long g = 0; g < red[places[q][w][0]].cantcadfacnodepartless[1]; g++) {
                    UpdateCDAarrays(red, probcond, probjoint, mejoint, tablaLocaleners, lugvecinoswfn,
                                    places, numcadenasxfactornode,
                                    K, E, q, w, j, g, M, Q, comb, pot_alphas, P_igual, P_mayor);
                }
            }
        }
    }
}



void sacarcheck(double tiempo, char *filecheck2, double **probjoint, int M,
                int numcadenasxfactornode){
    // Esta función se encarga de imprimir las probabilidades cada cierto tiempo.
    // Estos ficheros pueden ser usados para iniciar otra integración que comience
    // donde terminó la última.

    ofstream fcheck2(filecheck2);
    fcheck2 << tiempo << endl;
    for (int i = 0; i < M; i++){
        for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
            fcheck2 << probjoint[i][a] << "\t";
        }
        fcheck2 << endl;
    }
    fcheck2.close();
}


double leercheck(char *filecheck2, double **probjoint, int M,
                 int numcadenasxfactornode){
    // Esta función se encarga de leer las probabilidades en los ficheros "filecheck1" y "filecheck2".
    // Estos ficheros pueden ser usados para iniciar otra integración que comience
    // donde terminó la última.
    ifstream fcheck2(filecheck2);
    double tiempo;
    fcheck2 >> tiempo;
    for (int i = 0; i < M; i++){
        for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
            fcheck2 >> probjoint[i][a];
        }
    }
    fcheck2.close();

    return tiempo;
}


void RK2(Tred * red, double *** probcond, int N, double step, double cutoff, double ** probjoint, char *fileenergy,
         int K , int numcadenasxfactornode, int **** tablaLocaleners, double tiempo, int *** lugvecinoswfn,
         vector < vector <vector <int> > > places, double ** mejoint,
         double tcheck, char *filecheck2, int M, int *unsat, double Q, int **comb,double *pot_alphas,
         double *P_igual, double *P_mayor, double tol, double step_min, double E_min) {
    double Nec = 2 * numcadenasxfactornode * M;

    ofstream fener(fileenergy);
    int t = 1;
    double **k1, **k2, **probjoint1; //Magnitudes auxiliares para el Runge Kutta de la Master Equation
    InitializeToZero(red, probjoint1, k1, k2, M, numcadenasxfactornode);

    double E = energy(M, probjoint, unsat);
//cout << tiempo << "\t" << E << endl;
    fener << tiempo / N << "\t" << E / N << endl;

    while (tiempo < cutoff && E > E_min) //Runge-Kutta
    {

        for (int j = 0; j < M; j++) {
            for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
                mejoint[j][a] = 0;
            }
        }

        time_t medtiempo1, medtiempo2;
        medtiempo1 = time(0);

        ComputeProbcond(red, probcond, probjoint, N, numcadenasxfactornode);
        CMEordering(red, probcond, probjoint, mejoint, tablaLocaleners, lugvecinoswfn, places,
                    numcadenasxfactornode, K, E, M, Q, comb, pot_alphas, P_igual, P_mayor);


        for (int j = 0; j < M; j++) {
            for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
                k1[j][a] = step * mejoint[j][a];
                mejoint[j][a] = 0;
                probjoint1[j][a] = probjoint[j][a] + k1[j][a];
            }
        }

        medtiempo2 = time(0);
        cout << endl;
        cout << "tiempo de un paso: " << double(medtiempo2 - medtiempo1) / 60 << endl;
        cout << endl;

        E = energy(M, probjoint, unsat); //probcond1 y prob1 fueron updateados en el primer paso de Runge Kutta

        ComputeProbcond(red, probcond, probjoint1, N, numcadenasxfactornode);
        CMEordering(red, probcond, probjoint1, mejoint, tablaLocaleners, lugvecinoswfn, places,
                    numcadenasxfactornode, K, E, M, Q, comb, pot_alphas, P_igual, P_mayor);

        for (int j = 0; j < M; j++) {
            for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
                k2[j][a] = step * mejoint[j][a];
            }
        }

        double error = 0;
        bool validez = true;

        for (int j = 0; j < M; j++) {
            for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
                error += fabs(k1[j][a] - k2[j][a]);
                if (probjoint[j][a] + (k1[j][a] + k2[j][a]) / 2 < 0) {
                    validez = false;
                    break;
                }
            }
        }

        if (!validez) { //Esta parte garantiza que ninguna probabilidad se volvio menor que cero
            tol = tol / 2;
            step_min = step_min / 2;
        }
        cout << "error=" << error << endl;
        cout << "Nec=" << Nec << endl;
        if (validez && (error < 2 * tol * Nec || step == step_min)) {
            for (int j = 0; j < M; j++) {
                for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
                    probjoint[j][a] += (k1[j][a] + k2[j][a]) / 2;
                }
            }


            tiempo += step;
            E = energy(M, probjoint, unsat);
            fener << tiempo / N << "\t" << E / N << endl;
            if (tiempo >= N * tcheck * t) {
                t++;
                sacarcheck(tiempo / N, filecheck2, probjoint, M, numcadenasxfactornode);
            }
        }
        else
            E = energy(M, probjoint, unsat);
        if (error <= 32 * step * step * Nec * tol / (25 * N * N))
            step = N;
        else
            step = 4 * step * sqrt(2 * Nec * tol / error) / 5;
        cout << "stepRecomendado=" << step << endl;
        if (step < step_min)
            step = step_min;
    }

    fener.close();
}


vector <double> sort_array(double *array, int len){
    vector <double> sorted = vector <double> (len, 0);
    for (int i = 0; i < len; i++){
        sorted[i] = array[i];
    }
    sort(sorted.begin(), sorted.end());
    return sorted;
}


vector <long> comparer(double **me_sum, double **mejoint, long N, long M, int nch_fn, double thr){
    vector <long> dif;
    vector <double> v1, v2;
    bool cond;
    for (long he = 0; he < M; he++){
        v1 = sort_array(me_sum[he], nch_fn);
        v2 = sort_array(mejoint[he], nch_fn);
        cond = true;
        for (int ch = 0; ch < nch_fn; ch++){
            if (fabs(v1[ch] - N * v2[ch]) > thr){
                dif.push_back(he);
            }
        }
    }
    return dif;
}


int main(int argc, char *argv[]) {
    long N = 100;
    long M = 240;
    int K = 3;
    unsigned long seed_r = 1;
    double q = 0.2;
    double tl = 20;
    double tol = 1e-3;
    int nthr = 1;

    int nch_fn = (1 << K);
    double p0 = 0.5;

    double dt1 = 0.01;
    double dt_min = 1e-7;

    omp_set_num_threads(nthr);

    Tnode *nodes;
    Thedge *hedges;

    gsl_rng * r;
    init_ran(r, seed_r);

    char filegraph[300];
    char filelinks[300];
    sprintf(filegraph, "KSATgraph_K_%d_N_%li_M_%li_simetric_1_model_1_idum1_-2_J_1_ordered.txt", 
                       K, N, M);
    sprintf(filelinks, "KSAT_K_%d_enlaces_N_%li_M_%li_idumenlaces_-2_idumgraph_-2_ordered.txt", 
                       K, N, M);

    char fileener[300]; 
    sprintf(fileener, "CDA_WalkSAT_ener_K_%d_N_%li_M_%li_q_%.4lf_tl_%.2lf_seed_%li_tol_%.1e.txt", 
            K, N, M, q, tl, seed_r, tol);

    // create_graph(N, M, K, nodes, hedges, r);
    read_graph_old_order(filegraph, N, M, K, nodes, hedges);
    read_links(filelinks, N, M, K, nodes, hedges);
    int max_c = get_max_c(nodes, N);
    get_info_exc(nodes, hedges, N, M, K);

    
    double **binom_probs, **binom_sums;
    double **prob_joint, ***pu_cond, **me_sum, **pi;
    double e, pu_av, error;                 
                 
    // initalize all arrays that will be used inside the derivative
    init_aux_arr(binom_probs, binom_sums, max_c, K);
    init_probs(prob_joint, pu_cond, pi, me_sum, M, K, nch_fn, p0);

    // initialize auxiliary arrays for the Runge-Kutta integration
    double **k1, **k2, **prob_joint_1;
    init_RK_arr(k1, k2, prob_joint_1, M, nch_fn);
    
    e = energy(prob_joint, hedges, M);
    pu_av = e / M;


 
    comp_pcond(prob_joint, pu_cond, pi, hedges, M, K, nch_fn);
    get_all_binom_sums(max_c, pu_av, binom_probs, binom_sums);

    der_walksat(nodes, hedges, prob_joint, pu_cond, binom_probs, binom_sums, 
                M, K, nch_fn, q, e / N, me_sum, max_c); 


    bool valid = true;
        for (long he = 0; he < M; he++){
            for (int ch = 0; ch < nch_fn; ch++){
                k1[he][ch] = dt1 * me_sum[he][ch];
                prob_joint_1[he][ch] = prob_joint[he][ch] + k1[he][ch];
                if (prob_joint_1[he][ch] < 0){
                    valid = false;
                }
            }
        }

        while (!valid){
            cout << "joint probabilities became negative in the auxiliary step of RK2" << endl;
            dt1 /= 2;
            cout << "step divided by half    dt=" << dt1 << endl;
            if (dt1 < dt_min){
                dt_min /= 2;
                cout << "dt_min also halfed" << endl;
            }

            valid = true;
            for (long he = 0; he < M; he++){
                for (int ch = 0; ch < nch_fn; ch++){
                    k1[he][ch] = dt1 * me_sum[he][ch];
                    prob_joint_1[he][ch] = prob_joint[he][ch] + k1[he][ch];
                    if (prob_joint_1[he][ch] < 0){
                        valid = false;
                    }
                }
            }
        }
        
        e = energy(prob_joint_1, hedges, M);
        pu_av = e / M;
        comp_pcond(prob_joint_1, pu_cond, pi, hedges, M, K, nch_fn);
        get_all_binom_sums(max_c, pu_av, binom_probs, binom_sums);

        der_walksat(nodes, hedges, prob_joint_1, pu_cond, binom_probs, binom_sums, 
                    M, K, nch_fn, q, e / N, me_sum, max_c);

    /////////////////////////////////////////////////////////////////////////////////////


    int idumgraph; //semilla del grafo
    int idumenlaces; //semilla de los enlaces
    long iduminicond; //semilla para la condición inicial de las probabilidades
    long tiempo; // tiempo de integración
    double Q; // Parámetro del ruido del algoritmo WalkSAT
    double step; // Paso de integración
    double tcheck; // tiempo en el que se imprimen los checks
    int readcheck; // Si es 0, no se leen los checks al inicio, si es uno, se leen. Ningún otro valor está permitido
    long tiempoprev; // tiempo de la integración previa. Si se van a leer los checks es necesario para leer los ficheros

    double step_min;
    double E_min;
    double t0;



    N = 100; // Número de nodos
    M = 240; // Número de factor nodes
    K = 3;  // Cantidad de nodos en un factor node
    idumgraph = 2; //semilla del grafo
    tiempo = 20; // tiempo de integración
    Q = 0.2; // Parámetro del ruido del algoritmo WalkSAT
    step = 0.01; // Paso de integración
    readcheck = 0; // Si es 0, no se leen los checks al inicio, si es uno, se leen. Ningún otro valor está permitido
    tol = 1e-3;



    idumenlaces = idumgraph; //semilla de los enlaces
    iduminicond = -2; //semilla para la condición inicial de las probabilidades
    tcheck = 0.25; // tiempo en el que se imprimen los checks
    tiempoprev = 1; // tiempo de la integración previa. Si se van a leer los checks es necesario para leer los ficheros
    step_min = 1e-7;
    E_min = 1e-6 * N;
    t0 = 0;

    int numcadenasxfactornode = ( int) pow (2, K - 1); //Número de combinaciones de los vecinos de un spin en cada factor node
    int numvecxfacnode = (K - 1);//Número de vecinos en cada factor node

    int maxconnect = -1;  // máxima conectividad

    Tred * red; // En esta estructura esta la mayor parte de la informacion del grafo. Ver la biblioteca.


    double *** probcond;
    double ** probjoint;
    double **mejoint;
    //meprev y cmeprev son para Adams Bashforth, probcond y prob son las probabilidades
    //tabladE es una tabla con la parte independiente del tiempo de los rates de FMS, para una cierta entrada de energía
    //me y cme se usan para acumular las sumas que updatean las probabilidades. Ver la funcion CMEordering

    int **** tablaLocaleners;
    int *** lugvecinoswfn;
    //tablaLocaleners es una tabla de las energias locales
    //en tablaLocaleners se divide la energia local en dos componentes, debidas a dos grupos que contienen aproximadamente la mitad
    //de los factor nodes involucrados. Eso permite ahorrar memoria. Ver FillRates.
    // lugvecinoswfn contiene los vecinos de cada lugar. Ver inicred() en la biblioteca "CavMasEq_KSAT_bib_par.h"
    //factornodesthere contiene los factorn nodes en que ve el vecino j al nodo i. Ver FillFactornodesThere.

    map <int, pair< int, int> > * diccionario;
    vector <vector <vector < int>  >  > places;
    vector <map <int, int> > dictfacnodes;
    //diccionario contiene, para un sitio y uno de sus vecinos dados, cuales son el factor node y el sitio dentro de ese factor node
    //en que el vecino contiene al sitio. Ver LlenarDiccionario() en la biblioteca "CavMasEq_KSAT_bib_par.h"
    // dictfacnodes contiene, dentro de cada factor node, las posiciones de los nodos que los componen
    //places contiene, separados en grupos, a todos los nodos. Todos aquellos con la misma cantidad de factor nodes, y la misma cantidad de
    //enlaces positivos, partenecen a la misma clase.

    int *unsat;

    double cutoff = (double) tiempo;

    // Files
    char filegrafo[200];
    char fileenergy[300];
    char filecheck2[300];
    char fileenlaces[200];

    int **comb;
    double *pot_alphas;
    double P_igual[2], P_mayor[2];

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////


    cout << "N=" << N << "\t" << "M=" << M << "\t" << "seed=" << idumgraph << "\t" << "q=" << Q << "\t"
    << "\t" << "readcheck=" << readcheck << "\t" << "tl=" << tiempo << "\t" << "tol="
    << tol << endl;



    // Nombres de los ficheros
    sprintf(filegrafo,"KSATgraph_K_%d_N_%li_M_%li_simetric_1_model_1_idum1_-%d_J_1_ordered.txt", K,N, M, idumgraph);

    sprintf(fileenlaces, "KSAT_K_%d_enlaces_N_%li_M_%li_idumenlaces_-%d_idumgraph_-%d_ordered.txt", K, N, M,
            idumenlaces, idumgraph);


    inicred(red, N, filegrafo, maxconnect, numvecxfacnode, lugvecinoswfn); // Lee el grafo
    LlenarDiccionario(red, N, numvecxfacnode, lugvecinoswfn, diccionario); // Llena el diccionario de la estructura red
    inienlaces(red, N, fileenlaces); // Lee los enlaces

    dictfacnodes =  InitFacNodeInfo(red, numvecxfacnode, lugvecinoswfn, N, M, K, numcadenasxfactornode,
                                    unsat, diccionario);

    FillSatisfiedUnsatisfied(N, red, numcadenasxfactornode, numvecxfacnode, lugvecinoswfn, diccionario); //Llena los arreglos por spin que dicen cuáles
    // son, para cada factor node, las combinaciones de sus vecinos que son unsatisfied y satisfied
    Initializecadvec(red, N, numvecxfacnode, numcadenasxfactornode, lugvecinoswfn, diccionario); //Llena cadvec, que dice cual es la cadena del vecino para
    // una cadena del lugar dada

    InitializeToZero(red, probcond, probjoint, mejoint, N, M, numcadenasxfactornode); // Toma estos arreglos y los inicializa a 0


    if (readcheck == 0){
        condinicrandom(probcond, probjoint, red, N, M, K, iduminicond, numcadenasxfactornode, numvecxfacnode,
                       lugvecinoswfn, dictfacnodes);
        //Inicializa las probablidades, con una condición inicial factorizada y random
    }else if(readcheck == 1)
    {
        t0 = leercheck(filecheck2, probjoint, M, numcadenasxfactornode);
        // Lee los checks y los pone de condición inicial
    }



    FillRates(N, red, tablaLocaleners, maxconnect); //Llena la tabla de las energias locales
    FillPlaces(N, red, maxconnect, places); //Llena places
    delete [] diccionario; // Borrar diccionario
    diccionario = NULL;


    init_comb_pot_alphas(maxconnect, comb, pot_alphas);


    double **k1ser, **k2ser, **probjoint1; //Magnitudes auxiliares para el Runge Kutta de la Master Equation
    InitializeToZero(red, probjoint1, k1ser, k2ser, M, numcadenasxfactornode);

    double E = energy(M, probjoint, unsat);



    for (int j = 0; j < M; j++) {
        for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
            mejoint[j][a] = 0;
        }
    }

    ComputeProbcond(red, probcond, probjoint, N, numcadenasxfactornode);
    CMEordering(red, probcond, probjoint, mejoint, tablaLocaleners, lugvecinoswfn, places,
                numcadenasxfactornode, K, E, M, Q, comb, pot_alphas, P_igual, P_mayor);


    for (int j = 0; j < M; j++) {
            for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
                k1ser[j][a] = step * mejoint[j][a];
                mejoint[j][a] = 0;
                probjoint1[j][a] = probjoint[j][a] + k1ser[j][a];
            }
        }

        E = energy(M, probjoint, unsat); //probcond1 y prob1 fueron updateados en el primer paso de Runge Kutta

        ComputeProbcond(red, probcond, probjoint1, N, numcadenasxfactornode);
        CMEordering(red, probcond, probjoint1, mejoint, tablaLocaleners, lugvecinoswfn, places,
                    numcadenasxfactornode, K, E, M, Q, comb, pot_alphas, P_igual, P_mayor);




    comparer(me_sum, mejoint, N, M, nch_fn, 1e-6);
    return 0;
}