//
// Created by Francisco on 2/16/2017.
//

#ifndef CME_PSPINS_CAVMASEQ_PSPINS_BIB_H
#define CME_PSPINS_CAVMASEQ_PSPINS_BIB_H

#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <atomic>

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

using namespace std;

double ran3(long *idum) //generador de numeros aleatorios
{
    static int inext,inextp;
    static long ma[56];
    static int iff=0;
    long mj,mk;
    int i,ii,k;
    if (*idum < 0 || iff == 0) {
        iff=1;
        mj=labs(MSEED-labs(*idum));
        mj %= MBIG;
        ma[55]=mj;
        mk=1;
        for (i=1;i<=54;i++) {
            ii=(21*i) % 55;
            ma[ii]=mk;
            mk=mj-mk;
            if (mk < MZ) mk += MBIG;
            mj=ma[ii];
        }
        for (k=1;k<=4;k++)
            for (i=1;i<=55;i++) {
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }
        inext=0;
        inextp=31;
        *idum=1;
    }
    if (++inext == 56) inext=1;
    if (++inextp == 56) inextp=1;
    mj=ma[inext]-ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext]=mj;
    return mj*FAC;
}


typedef struct{
    int numvecinos;
    int numfactornodes;

    vector <int> valorenlaces;
    int numenlacespos;

    vector < vector <vector <int> > >  cadvec;
    vector <vector <int> > elemcad; //dice si la cadena es SAT-UNSAT

    vector <int> lugfacnodes; // lugfacnodes contiene los lugares de los factor nodes donde está el nodo
    vector <vector <int> > cadfacnode; // conversión de la cadena que ve el nodo a la que ve el factor node

    vector <int> lugaresupdate; //Aqui estan los vecinos j del nodo i que cumplen j>i
    int numlugupdate;

    vector <int> unsatisfied; //Dice cual es la cadena UNSAT
    long long cantcadfacnodes;

    vector <int> numfacnodepartless; //Ver LlenarLocaleners
    vector <vector <int> > numenlpospartless;
    vector <long long> cantcadfacnodepartless;

}Tnode;


int GetSmaller(int node, int *lugvecfnode, int numvecxfacnode){
    int smaller = node;
    for (int i = 0; i < numvecxfacnode; i++){
        if (lugvecfnode[i] < smaller){
            smaller = lugvecfnode[i];
        }
    }
    return smaller;
}


int WhichFacNode(vector <map <int, int> > dictfacnodes, int node, int *lugvecfnode, int numvecxfacnode, Tnode *red){
    int smaller = GetSmaller(node, lugvecfnode, numvecxfacnode);
    int h = 0;
    bool cond = false;
    while (h < red[smaller].numfactornodes && !cond){
        cond = dictfacnodes[red[smaller].lugfacnodes[h]].count(node) == 1;
        for (int i = 0; i < numvecxfacnode; i++){
            cond = cond && (dictfacnodes[red[smaller].lugfacnodes[h]].count(lugvecfnode[i]));
        }
        h++;
    }
    if (!cond){
        cout << "Node " << node << " is not in the neighborhood of node " << smaller << " in WhichFacNode" << endl;
        return -1;
    }else{
        return red[smaller].lugfacnodes[h - 1];
    }
}


bool CheckFacNodeNew(int *lugvecfnode, int node, int numvecxfacnode){
    bool answ = true;
    for (int i = 0; i < numvecxfacnode; i++){
        answ = (answ & (node < lugvecfnode[i]));
    }
    return answ;
}


map <int, int> PlacesInFacNode(int *lugvecfnode, int node, int numvecxfacnode, vector <int> &valor_enlaces,
                               map <int, pair< int, int> > * &diccionario, Tnode *red, int facnode){
    map <int, int> dict;
    dict[node] = 0;
    valor_enlaces.push_back(red[node].valorenlaces[facnode]);
    for (int i = 0; i < numvecxfacnode; i++){
        dict[lugvecfnode[i]] = i + 1;
        valor_enlaces.push_back(red[lugvecfnode[i]].valorenlaces[
                                        diccionario[lugvecfnode[i]][node * 1000 + facnode].first]);
    }
    return dict;
}

void InitCadFacNode(Tnode *red, int N, int K, int numvecxfacnode, int numcadenasxfactornode,
                    vector <map <int, int> > dictfacnodes, int ***lugvecinoswfn){
    // cadvec dice qué cadena ve el vecino cuando lugar ve la cadena "a"
    for (int j = 0; j < N; j++)
    {
        for (int h = 0; h < red[j].numfactornodes; h++) {
            red[j].cadfacnode.push_back(vector<int>());
            for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
                red[j].cadfacnode[h].push_back((a / numcadenasxfactornode) <<
                                               (K - dictfacnodes[red[j].lugfacnodes[h]][j] - 1));
                for (int k = 0; k < numvecxfacnode; k++){
                    red[j].cadfacnode[h][a] += ((a >> (numvecxfacnode - k - 1)) & 1) <<
                                               (K - dictfacnodes[red[j].lugfacnodes[h]][lugvecinoswfn[j][h][k]] - 1);
                }
            }
        }
    }
}

void FillUnSAT(int M, int numcadenasxfactornode, vector <int> * valor_enlaces, int *unsat, int K){
    //Llena los arreglos por nodo que dicen cuáles son, para cada factor node,
    // las combinaciones de sus vecinos que son unsatisfied y satisfied
    int prod, a;
    bool cond;

    for (int i = 0; i < M; i++){
        a = 0;
        cond = true;
        while (a < 2 * numcadenasxfactornode && cond){
            prod = 1;
            for(int k = 0; k < K; k++){
                prod *= 1 - (1 - 2 * ((a >> (K - k - 1)) & 1)) * valor_enlaces[i][k];
            }

            if (prod > 0){
                cond = false;
            }

            a++;
        }
        unsat[i] = a - 1;
    }
}

vector < map <int, int> > InitFacNodeInfo(Tnode *red, int numvecxfacnode, int *** lugvecinoswfn, int N, int M, int K,
                                          int numcadenasxfactornode, int *&unsat, map <int, pair< int, int> > * diccionario){
    vector < map <int, int> > dictfacnodes;
    vector <int> *valor_enlaces;
    valor_enlaces = new vector <int> [M];
    unsat = new int [M];
    int counter = 0;
    for (int i = 0; i < N; i++){
        for (int h = 0; h < red[i].numfactornodes; h++){
            if (CheckFacNodeNew(lugvecinoswfn[i][h], i, numvecxfacnode)){
                dictfacnodes.push_back(PlacesInFacNode(lugvecinoswfn[i][h], i, numvecxfacnode, valor_enlaces[counter],
                                                       diccionario, red, h));
                red[i].lugfacnodes.push_back(counter);
                counter++;
            } else{
                red[i].lugfacnodes.push_back(WhichFacNode(dictfacnodes, i, lugvecinoswfn[i][h], numvecxfacnode, red));
            }
        }
    }
    InitCadFacNode(red, N, K, numvecxfacnode, numcadenasxfactornode, dictfacnodes, lugvecinoswfn);
    FillUnSAT(M, numcadenasxfactornode, valor_enlaces, unsat, K);
    return dictfacnodes;
}



void inienlaces(Tnode *red, int N, char filename[])
{
    ifstream f(filename);
    int trash;
    for (int j = 0; j < N; j++)
    {
        red[j].numenlacespos = 0;
        f >> trash;
        for (int h = 0; h < red[j].numfactornodes; h++)
        {
            red[j].valorenlaces.push_back(0);
            f >> red[j].valorenlaces[h];
            if (red[j].valorenlaces[h] == 1)
            {
                red[j].numenlacespos++; //Aqui se guardan la cantidad de enlaces positivos
            }
        }
    }
    f.close();
}


void inicred(Tnode *&red,  int N, char filename[], int &maxconect, int numvecxfacnode, int *** &lugvecinoswfn) //Lee del grafo
{
    red = new Tnode[N];
    lugvecinoswfn = new int**[N];
    maxconect = 0;
    ifstream f(filename);
    for ( int i=0;i<N;i++)
    {
        int trash;
        f>>red[i].numvecinos;

        f>>red[i].numfactornodes;
        if (red[i].numfactornodes > maxconect) //Esto es para saber cuál es el máximo número de factor nodes al que pertenece un nodo
        {
            maxconect = red[i].numfactornodes;
        }


        red[i].cantcadfacnodes = (long long) pow(2, red[i].numfactornodes); // Es el número de combinaciones que pueden tener los factor
        // nodes si los clasificamos en unsatisfied y satisfied


        for ( int j=0;j<red[i].numfactornodes;j++)
        {
            f >> trash;
        }

        lugvecinoswfn[i] = new int *[red[i].numfactornodes];
        for ( int j=0;j<red[i].numfactornodes;j++)
        {
            lugvecinoswfn[i][j] = new int [numvecxfacnode];
            bool cond = true;
            for ( int h=0;h<numvecxfacnode;h++) {
                f >> trash;
                f >> trash;
                f >> lugvecinoswfn[i][j][h];

                if (lugvecinoswfn[i][j][h] < i) //Si pasé por el vecino ya, no lo updateo. Esto es para cuando se vaya a calcular la energía
                {
                    cond = false;
                }

                f >> trash;
                f >> trash;
            }
            if (cond){
                red[i].lugaresupdate.push_back(j);
            }
        }
        red[i].numlugupdate = (int) red[i].lugaresupdate.size();
    }
    f.close();
}


bool Estoy( int * lista, int length,  int elemento)
{
    bool answ = false;
    for (int i = 0; i < length; i++)
    {
        if (lista[i] == elemento)
        {
            answ = true;
        }
    }
    return answ;
}

bool Estoy( vector <int> lista,  int elemento)
{
    bool answ = false;
    for (int i = 0; i < lista.size(); i++)
    {
        if (lista[i] == elemento)
        {
            answ = true;
        }
    }
    return answ;
}


int FindPos( vector <int> lista,  int elemento)
{
    int pos = -1;
    int i = 0;
    bool cond = true;
    while (i < lista.size() && cond)
    {
        if (lista[i] == elemento)
        {
            pos = i;
            cond = false;
        }
        i++;
    }
    return pos;
}



bool EstoyFactorNode(int i,  int h,  int k,  int p,  int numvecxfacnode, int *** lugvecinoswfn)
{
    bool answ = true;
    for (int z = 0; z < numvecxfacnode; z++)
    {
        answ = answ && (Estoy(lugvecinoswfn[i][h], numvecxfacnode, lugvecinoswfn[lugvecinoswfn[i][h][k]][p][z])
                        || lugvecinoswfn[lugvecinoswfn[i][h][k]][p][z] == i);
    }
    return answ;
}


void LlenarDiccionario(Tnode *red, int N, int numvecxfacnode, int *** lugvecinoswfn,
                       map <int, pair< int, int> > * &diccionario)
{
    diccionario = new map <int, pair< int, int> > [N];
    for (int i = 0; i < N; i++)
    {
        pair < vector <int> , vector <int> > yapase; //Esto discrimina los vecinos por los que ya pase. En yapase.first estan los lugares,
        // y en yapase.second los correspondientes factor nodes en los que esos vecinos ven al
        //nodo i
        vector <int> indices; // Este garantiza que no repito una direccion
        for (int h = 0; h < red[i].numfactornodes; h++)
        {
            for ( int k = 0; k < numvecxfacnode; k++ )
            {
                if(!Estoy(yapase.first, lugvecinoswfn[i][h][k])) {
                    yapase.first.push_back(lugvecinoswfn[i][h][k]);
                    for (int p = 0; p < red[lugvecinoswfn[i][h][k]].numfactornodes; p++)
                    {
                        if (EstoyFactorNode(i, h, k, p, numvecxfacnode, lugvecinoswfn)){ //Busco en que factor node me ve el vecino
                            diccionario[i][lugvecinoswfn[i][h][k] * 1000 + p] = pair <int, int> (h, k);
                            indices.push_back(lugvecinoswfn[i][h][k] * 1000 + p); //Guardo la direccion para no repetirla
                            yapase.second.push_back(p);
                            break;
                        }
                    }
                }else{ //Si se repite
                    int position = FindPos(yapase.first, lugvecinoswfn[i][h][k]); //Busco donde estaba
                    int p = (yapase.second[position] + 1) % red[lugvecinoswfn[i][h][k]].numfactornodes;
                    while (p != yapase.second[position])
                    {
                        if (EstoyFactorNode(i, h, k, p, numvecxfacnode, lugvecinoswfn) && !Estoy(indices, lugvecinoswfn[i][h][k] * 1000 + p)){
                            //Si estoy en ese factor node, y no he guardado esa direccion, entonces:
                            diccionario[i][lugvecinoswfn[i][h][k] * 1000 + p] = pair <int, int> (h, k);
                            indices.push_back(lugvecinoswfn[i][h][k] * 1000 + p);
                            yapase.second[position] = p;
                            break;
                        }
                        p = (p + 1) % red[lugvecinoswfn[i][h][k]].numfactornodes;
                    }
                }
            }
        }
    }
}


void InitializeToZero(Tnode * red, double ***&probcond, double ** &probjoint, double ** &mejoint,
                      int N, int M, int numcadenasxfactornode) {

    probjoint = new double *[M];
    mejoint = new double *[M];
    for (int j = 0; j < M; j++) {
        probjoint[j] = new double[2 * numcadenasxfactornode];
        mejoint[j] = new double[2 * numcadenasxfactornode];
        for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
            probjoint[j][a] = 0;
            mejoint[j][a] = 0;
        }
    }

    probcond = new double **[N];
    for (int j = 0; j < N; j++) {
        probcond[j] = new double *[red[j].numfactornodes];
        for (int h = 0; h < red[j].numfactornodes; h++) {
            probcond[j][h] = new double[2 * numcadenasxfactornode];
            for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
                probcond[j][h][a] = 0; //probcond[lugar][factor node][cadenafn + spin*numcadenasxfactornode]
            }
        }
    }
}


void InitializeToZero(Tnode * red, double ** &probjoint, double ** &k1, double **&k2, int M,
                      int numcadenasxfactornode) {

    probjoint = new double *[M];
    k1 = new double *[M];
    k2 = new double *[M];
    for (int j = 0; j < M; j++) {
        probjoint[j] = new double[2 * numcadenasxfactornode];
        k1[j] = new double[2 * numcadenasxfactornode];
        k2[j] = new double[2 * numcadenasxfactornode];
        for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
            probjoint[j][a] = 0;
            k1[j][a] = 0;
            k2[j][a] = 0;
        }
    }
}


double Fprob( int dE, double temp)
{
    return (double) exp(-double(dE) / temp);
}

void TablaDeltaE (double * &tabladE, double T, int maxconnect) //Es la función que inicializa la tabla de los rats de Metrópolis
{
    tabladE = new double [maxconnect];
    for ( int i = 0; i < maxconnect; i++)
    {
        tabladE[i] = Fprob(i, T);
    }
}


void FillSatisfiedUnsatisfied(int N, Tnode * red,  int numcadenasxfactornode,
                              int numvecxfacnode, int *** lugvecinoswfn,
                              map <int, pair< int, int> > * diccionario){
    //Llena los arreglos por nodo que dicen cuáles son, para cada factor node,
    // las combinaciones de sus vecinos que son unsatisfied y satisfied
    int prod;
    for (int j =0; j < N; j++){
        for (int h = 0; h < red[j].numfactornodes; h++) {
            red[j].elemcad.push_back(vector <int> ());
            for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
                red[j].elemcad[h].push_back(-25);
            }
            for (int a = 0; a < numcadenasxfactornode; a++) {
                prod = 1;
                for (int k = 0; k < numvecxfacnode; k++) {
                    prod *= (1 - (1 - 2 * ((a >> (numvecxfacnode - k - 1)) & 1)) *
                                 red[lugvecinoswfn[j][h][k]].valorenlaces[
                                         diccionario[lugvecinoswfn[j][h][k]][j * 1000 + h].first]);
                }
                if (prod > 0) {
                    red[j].unsatisfied.push_back(a);
                    red[j].elemcad[h][a] = 1; //elemcad es 0 si la cadena es SAT y 0 si es UNSAT
                    red[j].elemcad[h][a + numcadenasxfactornode] = 1;
                }else{
                    red[j].elemcad[h][a] = 0;
                    red[j].elemcad[h][a + numcadenasxfactornode] = 0;
                }
            }
        }
    }
}

void Initializecadvec(Tnode *red, int N, int numvecxfacnode, int numcadenasxfactornode, int *** lugvecinoswfn,
                      map <int, pair< int, int> > * diccionario){
    // cadvec dice qué cadena ve el vecino cuando lugar ve la cadena "a"
    for (int j = 0; j < N; j++)
    {
        for (int h = 0; h < red[j].numfactornodes; h++)
        {
            red[j].cadvec.push_back(vector <vector <int> >  ());
            for (int k = 0; k < numvecxfacnode; k++)
            {
                red[j].cadvec[h].push_back(vector <int> ());
                for (int a = 0; a < 2 * numcadenasxfactornode; a++)
                {
                    red[j].cadvec[h][k].push_back(((a / numcadenasxfactornode) << (numvecxfacnode -
                                                                                   diccionario[lugvecinoswfn[j][h][k]][j * 1000 + h].second - 1)));
                    int r = (k + 1) % numvecxfacnode; //Necesitamos reordenar la cadena
                    while (r != k)
                    {
                        int factornodethere = diccionario[lugvecinoswfn[j][h][r]][j * 1000 + h].first;
                        red[j].cadvec[h][k][a] += (((a >> (numvecxfacnode - r - 1)) & 1) << (numvecxfacnode -
                                                                                             diccionario[lugvecinoswfn[j][h][k]][lugvecinoswfn[j][h][r] * 1000 + factornodethere].second - 1));
                        r = (r+1) %numvecxfacnode;
                    }

                    red[j].cadvec[h][k][a] +=  ((a >> (numvecxfacnode - k - 1)) & 1) * numcadenasxfactornode;
                }
            }
        }
    }

}


void condinicrandom(double *** probcond, double ** probjoint, Tnode * red, int N, int M, int K,
                    long idum,  int numcadenasxfactornode,  int numvecxfacnode, int *** lugvecinoswfn,
                    vector <map <int, int> > dictfacnodes){
    // Inicializa las probablidades, con una condición inicial factorizada y random

    double *prob;
    prob = new double [N];
    for (int j = 0; j < N; j++)
    {
        // prob[j] = ran3(&idum);
        prob[j] = 0.5;
    }

    double prod;
    for (int j = 0; j < N; j++)
    {
        for (int h = 0; h < red[j].numfactornodes; h++)
        {
            for (int k = 0; k < 2 * numcadenasxfactornode; k++)
            {
                prod = 1;
                for (int p = 0; p < numvecxfacnode; p++) { //El número binario de la cadena se lee de izquierda a derecha
                    prod *= ((k >> (numvecxfacnode - p - 1)) & 1) + //Hay que tener cuidado y multiplicar (1-Pi) cuando Si es -1
                            (1 - 2*((k >> (numvecxfacnode - p - 1)) & 1)) * prob[lugvecinoswfn[j][h][p]];
                }
                probcond[j][h][k] = prod;
            }
        }
    }

    int bits;
    for (int i = 0; i < M; i++){
        for (int a = 0; a < 2 * numcadenasxfactornode; a++) {
            prod = 1;
            for (map<int, int>::iterator it = dictfacnodes[i].begin(); it != dictfacnodes[i].end(); ++it) {
                bits = (a >> (K - it->second - 1)) & 1;
                prod *= bits - (2 * bits - 1) * prob[it->first];
            }
            probjoint[i][a] = prod;
        }
    }
}

void FillRates(int N, Tnode *red, int **** &tablaLocaleners, int maxconnect){
    // Llena la tabla de las energías locales tablaLocaleners[a][b][s][w]
    // Primer índice: cantidad de factor nodes
    // Segundo índice: cantidad de enlaces positivos
    // Tercer índice: valor del nodo
    // Cuarto índice: cadena de los factor nodes vecinos exceptuando uno
    tablaLocaleners = new int ***[maxconnect/ 2 + 2];
    tablaLocaleners[0] = new int **[1];
    tablaLocaleners[0][0] = new int *[2];
    for (int s = 0; s < 2; s++) {
        tablaLocaleners[0][0][s] = new int [1];
        tablaLocaleners[0][0][s][0] = 0; //Si no tienes vecinos, tu energia es por definicion 0
    }

    vector <vector <bool> > yaesta; //para saber por cuáles ya pasé y no repetir
    yaesta.push_back(vector <bool> ());
    yaesta[0].push_back(false);

    for (int j = 1; j < maxconnect/ 2 + 2; j++){
        tablaLocaleners[j] = new int **[j + 1];
        yaesta.push_back(vector <bool> ());
        for (int h = 0; h < j + 1; h++){
            tablaLocaleners[j][h] = new int *[2];
            yaesta[j].push_back(false);
        }
    }

    int Ej;
    int contador = 0;
    for (int j = 0; j < N; j++)
    {
        if (red[j].numfactornodes > 0) {
            int pos, fac = red[j].numfactornodes - 1;
            red[j].numfacnodepartless.push_back((fac) / 2); //quitando un factor node, cuántos quedan en la primera mitad
            red[j].numfacnodepartless.push_back(fac - red[j].numfacnodepartless[0]); //en la segunda mitad. Dividir en dos ahorra espacio
            red[j].cantcadfacnodepartless.push_back(pow(2, red[j].numfacnodepartless[0])); //cuántas cadenas representa esto
            red[j].cantcadfacnodepartless.push_back(pow(2, red[j].numfacnodepartless[1]));
            for (int q = 0; q < red[j].numfactornodes; q++) {  // Este ciclo se encarga de llenar los parámetros
                // correspondientes al caso en que falta un factor node (el q)
                red[j].numenlpospartless.push_back(vector<int>());

                if (q < red[j].numenlacespos) { // cuántos enlaces positivos hay en la primera mitad luego de quitar q
                    pos = red[j].numenlacespos - 1; // si q < numenlacespos entonces hay que quitar uno
                } else {
                    pos = red[j].numenlacespos;  // si no, q tiene un enlace negativo y no hay que quitar nada
                }


                if (pos >= red[j].numfacnodepartless[0]) {
                    // Si la cantidad de enlaces positivos es mayor que la cantidad de factor nodes
                    // en la primera mitad, entonces todos en la primera mitad son positivos, y en la
                    // segunda habrán pos - numfacnodepartless[0] positivos
                    red[j].numenlpospartless[q].push_back(red[j].numfacnodepartless[0]);
                    red[j].numenlpospartless[q].push_back(pos - red[j].numfacnodepartless[0]);
                } else {
                    // Si es menor, en la primera hay "pos" enlaces positivos y 0 negativos
                    red[j].numenlpospartless[q].push_back(pos);
                    red[j].numenlpospartless[q].push_back(0);
                }

                int r = 0;
                for (int p = 0; p < 2; p++) { // Este ciclo se encarga de calcular las energías
                    // p = 0, primera mitad, p = 1, segunda mitad
                    if (!yaesta[red[j].numfacnodepartless[p]][red[j].numenlpospartless[q][p]]) {
                        // yaesta se encarga de que no se calcule lo mismo dos veces
                        contador++; // cuenta la cantidad de clases distintas que hay
                        yaesta[red[j].numfacnodepartless[p]][red[j].numenlpospartless[q][p]] = true;
                        for (int s = 0; s < 2; s++) {
                            tablaLocaleners[red[j].numfacnodepartless[p]][red[j].numenlpospartless[q][p]][s] = new int [red[j].cantcadfacnodepartless[p]];

                            for (long long w = 0; w < red[j].cantcadfacnodepartless[p]; w++) {
                                Ej = 0;
                                for (int h = 0; h < red[j].numfacnodepartless[p]; h++) {
                                    if ((h + p * red[j].numfacnodepartless[0]) >= q) { // r se encarga de saltarse el factor node que estamos
                                        r = 1;                                         // quitando, que es el q
                                    }else{
                                        r = 0;
                                    }
                                    Ej += ((w >> (red[j].numfacnodepartless[p] - h - 1)) & 1) *
                                          (1 - (1 - 2 * s) *
                                               red[j].valorenlaces[h + p * red[j].numfacnodepartless[0] + r]) / 2;
                                }
                                tablaLocaleners[red[j].numfacnodepartless[p]][red[j].numenlpospartless[q][p]][s][w] = Ej;
                            }

                        }
                    }
                }
            }
        }
    }
    cout << contador << endl;
}

void FillPlaces(int N, Tnode * red, int maxconnect, vector < vector <vector <int> > > &places){
    // Llena places, que agrupa los nodos en clases. Cada clase se identifica por la cantidad de factor nodes
    // a la que pertenece y el numero de enlaces positivos que tiene. Eso se usa luego al calcular las sumas
    // de la CME, pues se pueden ahorrar iteraciones al incluir, en un mismo ciclo, todos los lugares que
    // pertenecen a una de estas clases

    places.push_back(vector <vector <int> > ());

    vector <vector <int> >  indices;
    indices.push_back(vector <int> ());
    indices[0].push_back(-1);

    for (int j = 1; j < maxconnect + 1; j++){
        places.push_back(vector <vector <int> > ());
        indices.push_back(vector <int>  ());
        for (int h = 0; h < j + 1; h++){
            indices[j].push_back(-1);
        }
    }

    for (int i = 0; i < N; i++)
    {
        if (indices[red[i].numfactornodes][red[i].numenlacespos] == -1)
        {
            places[red[i].numfactornodes].push_back(vector <int> ());
            indices[red[i].numfactornodes][red[i].numenlacespos] = (int) (places[red[i].numfactornodes].size() - 1);
        }
        places[red[i].numfactornodes][indices[red[i].numfactornodes][red[i].numenlacespos]].push_back(i);
    }
}


#endif //CME_PSPINS_CAVMASEQ_PSPINS_BIB_H
