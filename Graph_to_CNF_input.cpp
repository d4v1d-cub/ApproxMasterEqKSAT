#include <iostream>
#include <fstream>
#include <map>
#include <stdlib.h>
#include <vector>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

using namespace std;

void init_ran(gsl_rng * &r, unsigned long s){
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, s);
}


typedef struct{
    vector <long> nodes_in;  // nodes inside the factor node
    vector <int> links;
}Thedge;


void init_graph(Thedge *&hedges, long M){
    hedges = new Thedge[M];

    for (long he = 0; he < M; he++){
        hedges[he].links = vector <int> ();
        hedges[he].nodes_in = vector <long> ();
    }
}


long create_graph(long N, long M, int K, Thedge *&hedges, gsl_rng * r, long *&disc_before){
    disc_before = new long [N];
    int *numfn = (int *) calloc(N, sizeof(int));
    init_graph(hedges, M);
    int w, h;
    long var;
    bool cond;
    for (long he = 0; he < M; he++){
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
                }else{
                    hedges[he].links.push_back(-1);
                }
                numfn[var]++;
                w++;
            }
        }
    }

    long counter = 0;
    for (long i = 0; i < N; i++){
        disc_before[i] = counter;
        if (numfn[i] == 0){
            counter ++;
        }
    }
    return counter;
}


void PrintToInput(Thedge *hedges, long N, long M, long *disc_before, long ndisc){
    cout << "p" << "\t" << "cnf" << "\t" << N - ndisc << "\t" << M << endl;
    for (long he = 0; he < M; he++){
        for (int w = 0; w < hedges[he].nodes_in.size(); w++){
            cout << (hedges[he].nodes_in[w] - disc_before[hedges[he].nodes_in[w]]) * 
                    hedges[he].links[w] << "\t";
        }
        cout << 0 << endl;
    }
}



int main(int argc, char *argv[]) {
    long N = atol(argv[1]);
    long M = atol(argv[2]);
    int K = atoi(argv[3]);
    unsigned long seed = atol(argv[4]);

    Thedge *hedges;

    gsl_rng * r;
    init_ran(r, seed);

    long *disc_before;
    long ndisc = create_graph(N, M, K, hedges, r, disc_before);

    PrintToInput(hedges, N, M, disc_before, ndisc);

    return 0;
}