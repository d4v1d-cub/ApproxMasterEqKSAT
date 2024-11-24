#include <iostream>
#include <fstream>
#include <map>
#include <stdlib.h>
#include <vector>

using namespace std;

typedef struct{
    int numvecinos;
    int numfactornodes;
    map < int, pair< int, int> > diccionario;
    int *valorenlaces;

    long **lugvecinoswfn;
    int **lugaresnoupdate;
    int numlugnoupdate;
    int *lugaresupdate;
    int numlugupdate;

    int *numvecxfacnode;



}TgrafoKSAT;

bool Estoy( long *lista, int length,  long elemento)
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

bool Estoy( vector <long> lista,  long elemento)
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

bool Estoy( vector <int> lista,  long elemento)
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


int FindPos( vector <long> lista,  long elemento)
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



bool EstoyFactorNode(TgrafoKSAT *red, long i,  int h,  int k,  int p,  int numvecxfacnode)
{
    bool answ = true;
    for (int z = 0; z < numvecxfacnode; z++)
    {
        answ = answ && (Estoy(red[i].lugvecinoswfn[h], numvecxfacnode, red[red[i].lugvecinoswfn[h][k]].lugvecinoswfn[p][z])
                        || red[red[i].lugvecinoswfn[h][k]].lugvecinoswfn[p][z] == i);
    }
    return answ;
}


void LlenarDiccionario(TgrafoKSAT *red, long N)
{
    for (long i = 0; i < N; i++)
    {
        pair < vector <long> , vector <int> > yapase;
        for (int h = 0; h < red[i].numfactornodes; h++)
        {
            for ( int k = 0; k < red[i].numvecxfacnode[h]; k++ )
            {
                if(!Estoy(yapase.first, red[i].lugvecinoswfn[h][k])) {
                    yapase.first.push_back(red[i].lugvecinoswfn[h][k]);
                    for (int p = 0; p < red[red[i].lugvecinoswfn[h][k]].numfactornodes; p++)
                    {
                        if (EstoyFactorNode(red, i, h, k, p, red[i].numvecxfacnode[h])){
                            red[i].diccionario[red[i].lugvecinoswfn[h][k] * 1000 + p] = pair <int, int> (h, k);
                            yapase.second.push_back(p);
                            break;
                        }
                    }
                }else{
                    int position = FindPos(yapase.first, red[i].lugvecinoswfn[h][k]);
                    for (int p = yapase.second[position] + 1; p < red[red[i].lugvecinoswfn[h][k]].numfactornodes; p++)
                    {
                        if (EstoyFactorNode(red, i, h, k, p, red[i].numvecxfacnode[h])){
                            red[i].diccionario[red[i].lugvecinoswfn[h][k] * 1000 + p] = pair <int, int> (h, k);
                            yapase.second[position] = p;
                            break;
                        }
                    }
                }
            }
        }
    }
}

void inienlaces(TgrafoKSAT *red, long N, char filename[])
{
    ifstream f(filename);
    int trash;
    for (int j = 0; j < N; j++)
    {
        red[j].valorenlaces = new int[red[j].numfactornodes];
        f >> trash;
        for (int h = 0; h < red[j].numfactornodes; h++)
        {
            f >> red[j].valorenlaces[h];
        }
    }
    f.close();
}



void inicred(TgrafoKSAT *red,  long N, char filename[], int &maxconect) //Lee del grafo
{
    maxconect = 0;
    ifstream f(filename);
    for ( int i=0;i<N;i++)
    {
        int trash=0;
        f>>red[i].numvecinos;

        f>>red[i].numfactornodes;
        if (red[i].numfactornodes > maxconect) //Esto es para saber cuál es el máximo número de factor nodes al que pertenece un nodo
        {
            maxconect = red[i].numfactornodes;
        }

        red[i].lugaresnoupdate =  (int **) malloc(red[i].numfactornodes* sizeof(int*)); //Son los lugares que no hay que updatear cuando se recalcula la energía local del sitio i, los j<i

        red[i].lugaresupdate = (int *) malloc(red[i].numfactornodes* sizeof(int)); //Son los que hay que updatear, los j>i

        red[i].lugvecinoswfn = (long int **) malloc(
                red[i].numfactornodes * sizeof(long int *)); //Aquí se guardan las posiciones de los vecinos

        red[i].numvecxfacnode = (int*) malloc(red[i].numfactornodes* sizeof(int));


        int cont1 = 0; //Para ir contando la cantidad de lugares a updatear
        int cont2 = 0; //y a no updatear

        for ( int j=0;j<red[i].numfactornodes;j++)
        {
            f >> red[i].numvecxfacnode[j];
        }

        for ( int j=0;j<red[i].numfactornodes;j++)
        {
            red[i].lugvecinoswfn[j] = (long int *) malloc(red[i].numvecxfacnode[j] * sizeof(long int));
            red[i].lugaresnoupdate[j] = (int*) malloc(2* sizeof(int));
            bool cond = true;
            int menor = N + 1; //para ubicar en qué nodo buscar el término de energía correspondiente al factor node j
            for ( int h=0;h<red[i].numvecxfacnode[j];h++) {
                f >> trash;
                f >> trash;
                f >> red[i].lugvecinoswfn[j][h];

                if (red[i].lugvecinoswfn[j][h] < i) //Si pasé por el vecino ya, no lo updateo. Esto es para cuando se vaya a calcular la energía
                {
                    if (menor > red[i].lugvecinoswfn[j][h]) //Me quedo con el menor
                    {
                        menor = red[i].lugvecinoswfn[j][h];
                    }
                    cond = false;
                }

                f >> trash;
                f >> trash;
            }
            if (cond){
                red[i].lugaresupdate[cont2] = j;
                cont2 ++;
            }else{
                red[i].lugaresnoupdate[cont1][0] = j;
                red[i].lugaresnoupdate[cont1][1] = menor;
                cont1 ++;
            }
        }
        red[i].numlugnoupdate = cont1;
        red[i].numlugupdate = cont2;

    }
    f.close();
}


vector <long> NumNotParticipating(TgrafoKSAT *red,  long N)
{
    vector <long> noconectados;
    for (long j = 0; j < N; j++)
    {
        if (red[j].numfactornodes == 0)
        {
            noconectados.push_back(j);
        }
    }
    return noconectados;
}

int restar(vector <long> noconectados, long lugar)
{
    int contador = 0;
    for (int i = 0; i < noconectados.size(); i++)
    {
        if (noconectados[i] < lugar)
        {
            contador++;
        }
    }
    return contador;
}

void TranslateAndPrint(TgrafoKSAT *red,  long N, long cantidadvariables, long M, vector <long> noconectados)
{

    cout << "p" << "\t" << "cnf" << "\t" << cantidadvariables << "\t" << M << endl;

    for (long j = 0; j < N; j++)
    {
        for (int h = 0; h < red[j].numlugupdate; h++)
        {
            cout << red[j].valorenlaces[red[j].lugaresupdate[h]] * (j+1 - restar(noconectados, j)) << "\t";
            for (int k = 0; k < red[j].numvecxfacnode[red[j].lugaresupdate[h]]; k++)
            {
                cout << red[red[j].lugvecinoswfn[red[j].lugaresupdate[h]][k]].valorenlaces[
                                red[red[j].lugvecinoswfn[red[j].lugaresupdate[h]][k]].diccionario[j * 1000 + red[j].lugaresupdate[h]].first]*
                        (red[j].lugvecinoswfn[red[j].lugaresupdate[h]][k] + 1 - restar(noconectados, red[j].lugvecinoswfn[red[j].lugaresupdate[h]][k])) << "\t";
            }
            cout << 0 << endl;
        }
    }
}


int main(int argc, char *argv[]) {
    long N = atol(argv[1]);
    long M = atol(argv[2]);
    int idumgraph = atoi(argv[3]); //semilla del grafo
    int idumenlaces = atoi(argv[4]); //semilla de los enlaces li
    int K = atoi(argv[5]);

    TgrafoKSAT *red;
    red = new TgrafoKSAT[N];

    int maxconect = -1;

    char filegrafo[200];
    sprintf(filegrafo,"KSATgraph_K_%d_N_%li_M_%li_simetric_1_model_1_idum1_-%d_J_1.txt", K,N, M, idumgraph);

    char fileenlaces[200];
    sprintf(fileenlaces, "KSAT_K_%d_enlaces_N_%li_M_%li_idumenlaces_-%d_idumgraph_-%d.txt", K, N, M,
            idumenlaces, idumgraph);

    inicred(red, N, filegrafo, maxconect); //Lee del grafo
    inienlaces(red, N, fileenlaces); //lee los enlaces
    LlenarDiccionario(red, N);

    vector <long> noconectados =  NumNotParticipating(red, N);

    TranslateAndPrint(red, N, long(N - noconectados.size()), M, noconectados);

    return 0;
}