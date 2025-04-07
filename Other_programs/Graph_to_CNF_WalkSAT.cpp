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


long create_graph(long N, long M, int K, Thedge *&hedges, gsl_rng * r){
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
}


void PrintToInput(Thedge *hedges, long N, long M){
    cout << N << "\t" << M << endl;
    for (long he = 0; he < M; he++){
        cout << (hedges[he].nodes_in[0] + 1) * hedges[he].links[0];
        for (int w = 1; w < hedges[he].nodes_in.size(); w++){
            cout << "\t" << (hedges[he].nodes_in[w] + 1) * hedges[he].links[w];
        }
        cout << endl;
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

    create_graph(N, M, K, hedges, r);

    PrintToInput(hedges, N, M);

    return 0;
}