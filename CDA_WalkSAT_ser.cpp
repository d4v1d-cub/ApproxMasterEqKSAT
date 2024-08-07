#include <cmath>
#include <iostream>
#include "CDA_WalkSAT_bib_ser.h"

using namespace std;

double produccond(long long cadenaoffn, int spin, int length, int lugar, Tnode * red, int factornode,
                  double ** probcond, int numcadenasxfacnode) {
    //produccond hace el producto de todas las probabilidades que condicionan en el vecino excepto la del factor node donde está "lugar"
    //la cadenaoffn viene en el mismo orden de los factor nodes en Tnode * red, quitando factornode. Es importante aclarar
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


double total_rate_walksat_greedy(Tnode * red, double alpha_SAT, double alpha_UNSAT, int num_sat_clauses_i, int lugar,
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


void ComputeProbcond(Tnode * red, double ***probcond, double **probjoint, int N, int numcadenasxfactornode){
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


void UpdateCDAarrays(Tnode * red, double *** probcond, double ** probjoint,
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
                E0 = Elugar0 - Efacnode0 * red[lugar].elemcad[k][cad];
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



void CMEordering(Tnode * red, double *** probcond, double ** probjoint,
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


void RK2(Tnode * red, double *** probcond, int N, double step, double cutoff, double ** probjoint, char *fileenergy,
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


int main(int argc, char *argv[]) {

    int N; // Número de nodos
    int M; // Número de factor nodes
    int K;  // Cantidad de nodos en un factor node
    int idumgraph; //semilla del grafo
    int idumenlaces; //semilla de los enlaces
    long iduminicond; //semilla para la condición inicial de las probabilidades
    long tiempo; // tiempo de integración
    double Q; // Parámetro del ruido del algoritmo WalkSAT
    double step; // Paso de integración
    double tcheck; // tiempo en el que se imprimen los checks
    int readcheck; // Si es 0, no se leen los checks al inicio, si es uno, se leen. Ningún otro valor está permitido
    long tiempoprev; // tiempo de la integración previa. Si se van a leer los checks es necesario para leer los ficheros

    double tol;
    double step_min;
    double E_min;
    double t0;



    N = atoi(argv[1]); // Número de nodos
    M = atoi(argv[2]); // Número de factor nodes
    K = atoi(argv[3]);  // Cantidad de nodos en un factor node
    idumgraph = atoi(argv[4]); //semilla del grafo
    tiempo = atof(argv[5]); // tiempo de integración
    Q = atof(argv[6]); // Parámetro del ruido del algoritmo WalkSAT
    step = atof(argv[7]); // Paso de integración
    readcheck = atoi(argv[8]); // Si es 0, no se leen los checks al inicio, si es uno, se leen. Ningún otro valor está permitido
    tol = atof(argv[9]);


    /*N = 100;
    M = 240;
    K = 3;
    idumgraph = 2;
    tiempo = 20;
    Q = 1.0;
    step = 0.01;
    readcheck = 0;
    tol=1e-2;*/


    idumenlaces = idumgraph; //semilla de los enlaces
    iduminicond = -2; //semilla para la condición inicial de las probabilidades
    tcheck = 0.25; // tiempo en el que se imprimen los checks
    tiempoprev = 1; // tiempo de la integración previa. Si se van a leer los checks es necesario para leer los ficheros
    step_min = 1e-7;
    E_min = max(1e-6 * N, 1.0);
    t0 = 0;

    int numcadenasxfactornode = ( int) pow (2, K - 1); //Número de combinaciones de los vecinos de un spin en cada factor node
    int numvecxfacnode = (K - 1);//Número de vecinos en cada factor node

    int maxconnect = -1;  // máxima conectividad

    Tnode * red; // En esta estructura esta la mayor parte de la informacion del grafo. Ver la biblioteca.


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
    sprintf(filegrafo,"KSATgraph_K_%d_N_%d_M_%d_simetric_1_model_1_idum1_-%d_J_1_ordered.txt", K,N, M, idumgraph);

    sprintf(fileenlaces, "KSAT_K_%d_enlaces_N_%d_M_%d_idumenlaces_-%d_idumgraph_-%d_ordered.txt", K, N, M,
            idumenlaces, idumgraph);

    sprintf(filecheck2, "CDA_KSAT_WalkSAT_quenched_K_%d_ABo4_check_probjoint_N_%d_M_%d_q_%.4lf_step_%.3lf_time_%li_"
            "idumcondinic_%li_idumgraph_-%d_tol_%.2e.txt", K, N, M, Q, step, tiempoprev, iduminicond, idumgraph, tol);

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

    // Nombre del fichero de las energías
    sprintf(fileenergy, "CDA_KSAT_WalkSAT_quenched_K_%d_ABo4_ener_global_N_%d_M_%d_q_%.4lf_step_%.3lf_time_%li_"
                    "idumcondinic_%li_idumgraph_-%d_t0_%.3lf_tol_%.2e.txt", K, N, M, Q, step, tiempo, iduminicond,
            idumgraph, t0, tol);


    init_comb_pot_alphas(maxconnect, comb, pot_alphas);

    //Esta es la tipa
    RK2(red, probcond, N, step, cutoff * N, probjoint, fileenergy, K, numcadenasxfactornode,
        tablaLocaleners, t0 * N, lugvecinoswfn, places, mejoint, tcheck, filecheck2, M, unsat, Q,
        comb, pot_alphas, P_igual, P_mayor, tol, step_min, E_min);

    return 0;
}
