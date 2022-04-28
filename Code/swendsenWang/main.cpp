// Algoritmo de clusters Swendsen-Wang para el modelo 2D de Ising

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <list>                 // salvar valores para la autocorrelacion


using namespace std;


double J = +1;                  // Fuerza de interaccion spin spin
int Lx, Ly;                     // Numero de spines en x y
int N;                          // Numero total de spines
int **s;                        // Arreglo de spines
double T;                       // Temperatura
double H = 0;                   // Campo Externo
int steps = 0;                  // Pasos de Monte Carlo hasta el momento
int seed=-59;                   // Variable para el generador de numeros aleatorios

double ran2 (int& idum) {       // Generador de numeros aleatorios de Numerical Recipes
    const int IM1 = 2147483563, IM2 = 2147483399;
    const double AM=(1.0/IM1);
    const int IMM1 = IM1-1;
    const int IA1 = 40014, IA2 = 40692, IQ1 = 53668, IQ2 = 52774;
    const int IR1 = 12211, IR2 = 3791, NTAB = 32;
    const int NDIV = 1+IMM1/NTAB;
    const double EPS = 3.0e-16, RNMX = 1.0-EPS;
    int j, k;
    static int idum2=123456789, iy = 0;
    static int iv[NTAB];
    double temp;
    if (idum <= 0) {
        idum = (idum == 0 ? 1 : -idum);
        idum2=idum;
        for (j=NTAB+7;j>=0;j--) {
            k=idum/IQ1;
            idum=IA1*(idum-k*IQ1)-k*IR1;
            if (idum < 0) idum += IM1;
            if (j < NTAB) iv[j] = idum;
        }
        iy=iv[0];
    }
    k=idum/IQ1;
    idum=IA1*(idum-k*IQ1)-k*IR1;
    if (idum < 0) idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}



void initialize ( ) {
    s = new int* [Lx];

    for (int i = 0; i < Lx; i++)
        s[i] = new int [Ly];
    for (int i = 0; i < Lx; i++)
        for (int j = 0; j < Ly; j++)
            s[i][j] = ran2(seed) < 0.5 ? +1 : -1;   // Distribucion inicial aleatoria
    steps = 0;
}

bool **iBondFrozen, **jBondFrozen;  // Red de enlaces> Dos enlaces por spin
double freezeProbability;           // 1 - e^(-2J/kT)
int **cluster;                      // Numeros de cluster para los spines
int *labelLabel;                    // Para determinar numeros de cluster  "propios"
bool *sNewChosen;                   // Se selecciono el nuevo numero?
int *sNew;                          // Nuevos numeros de spin aleatorios en cada cluster

void initializeClusterVariables() {

    // crear las redes de enlaces en ambas direcciones
    iBondFrozen = new bool* [Lx];
    jBondFrozen = new bool* [Lx];
    for (int i = 0; i < Lx; i++) {
        iBondFrozen[i] = new bool [Ly];
        jBondFrozen[i] = new bool [Ly];
    }

    // calcular la probabilidad decongelamiento de enlace
    freezeProbability = 1 - exp(-2*J/T);

    // crear un arreglo para los numeros de cluster
    cluster = new int* [Lx];
    for (int i = 0; i < Lx; i++)
        cluster[i] = new int [Ly];

    // crear arreglos de tamano N para:
    labelLabel = new int [N];        // punteros a numeros de cluster "propios"
    sNewChosen = new bool [N];       // crear nuevos valores de spin cluster
    sNew = new int [N];              // nuevos valores de numeros de cluster
}

// declarar funciones
void freezeOrMeltBonds();
int properLabel(int label);
void labelClusters();
void flipClusterSpins();

void oneMonteCarloStep() {

    // primero construye una red de enlaces con enlaces congelados
    freezeOrMeltBonds();

    // Usar el algoritmo de Hoshen-Kopelman para identificar los numeros de cluster
    labelClusters();

    // Voltea aleatoriamente los clusters de spin
    flipClusterSpins();

    ++steps;
}

void freezeOrMeltBonds() {

    // visita todos los spines de la red
    for (int i = 0; i < Lx; i++)
    for (int j = 0; j < Ly; j++) {

        // congela o derrite  los dos  enlaces conectados a este spin
        // usando un criterio que depende del factor de Boltzmann
        iBondFrozen[i][j] = jBondFrozen[i][j] = false;

        // enlace en la direccion i
        int iNext = i == Lx-1 ? 0 : i+1;
        if (s[i][j] == s[iNext][j] && ran2(seed) < freezeProbability)
            iBondFrozen[i][j] = true;

        // enlace en la  direccion j
        int jNext = j == Ly-1 ? 0 : j+1;
        if (s[i][j] == s[i][jNext] && ran2(seed) < freezeProbability)
            jBondFrozen[i][j] = true;
    }
}

int properLabel(int label) {
    while (labelLabel[label] != label)
        label = labelLabel[label];
    return label;
}

void labelClusters() {

    int label = 0;

    // visita todos los lugares de la red
    for (int i = 0; i < Lx; i++)
    for (int j = 0; j < Ly; j++) {

        // encuentra sitios previamentes visitados conectados a i,j por enlaces congelados
        int bonds = 0;
        int iBond[4], jBond[4];

        // revisar en enlace a i-1,j
        if (i > 0 && iBondFrozen[i - 1][j]) {
            iBond[bonds] = i - 1;
            jBond[bonds++] = j;
        }

        // aplicar la condicion de fronteras periodicas al borde:
        // si i,j es el ultimo sitio, revisa i+1,j
        if (i == Lx - 1 && iBondFrozen[i][j]) {
            iBond[bonds] = 0;
            jBond[bonds++] = j;
        }

        // revisa el enlace a i,j-1
        if (j > 0 && jBondFrozen[i][j - 1]) {
            iBond[bonds] = i;
            jBond[bonds++] = j - 1;
        }

        // revisar condiciones de frontera periodicas
        if (j == Ly - 1 && jBondFrozen[i][j]) {
            iBond[bonds] = i;
            jBond[bonds++] = 0;
        }

        // revisar el numero de enlaces hacia sitios previamente visitados
        if (bonds == 0) { // se necesita empezar un nuevo cluster
            cluster[i][j] = label;
            labelLabel[label] = label;
            ++label;
        } else {          // re-numera los spines enlazados con el menor de los numeros
            int minLabel = label;
            for (int b = 0; b < bonds; b++) {
                int pLabel = properLabel(cluster[iBond[b]][jBond[b]]);
                if (minLabel > pLabel)
                    minLabel = pLabel;
            }

            // establecer el numero del sitio al menor de los numeros
            cluster[i][j] = minLabel;

            // reestablecer los punteros a numeros propios
            for (int b = 0; b < bonds; b++) {
                int pLabel = cluster[iBond[b]][jBond[b]];
                labelLabel[pLabel] = minLabel;

                // reestablecer numeros de cluster en los sitios conectados
                cluster[iBond[b]][jBond[b]] = minLabel;
            }
        }
    }
}

void flipClusterSpins() {

    for (int i = 0; i < Lx; i++)
    for (int j = 0; j < Ly; j++) {

        // los valores de numeros de cluster aleatorios aun no se han definido
        int n = i * Lx + j;
        sNewChosen[n] = false;

        // reemplazar todos los numeros impropios por los propios
        cluster[i][j] = properLabel(cluster[i][j]);
    }

    int flips = 0;    // para contar el numero de spines volteados
    for (int i = 0; i < Lx; i++)
    for (int j = 0; j < Ly; j++) {

        // encontrar el numero propio del cluster
        int label = cluster[i][j];

        // seleccionar un nuevo valor para todos los spines del cluster aleatoriamente
        // solo si esto ya no se ha hecho
        if (!sNewChosen[label]) {
            sNew[label] = ran2(seed) < 0.5 ? +1 : -1;
            sNewChosen[label] = true;
        }

        // evaluar el  valor de spin y contar el numero de flips (volteos)
        if (s[i][j] != sNew[label]) {
            s[i][j] = sNew[label];
            ++flips;
        }
    }
}

double eSum;                // acumulador de la energia por spin
double eSqdSum;             // acumulador del cuadrado de la energia por spin
int nSum;                   // numer de terminos de la suma

void initializeObservables() {
    eSum = eSqdSum = 0;     // acumuladores inicializados
    nSum = 0;               // cero terminos sumados
}

void measureObservables() {
    int sSum = 0, ssSum = 0;
    for (int i = 0; i < Lx; i++)
    for (int j = 0; j < Ly; j++) {
        sSum += s[i][j];
        int iNext = i == Lx-1 ? 0 : i+1;
        int jNext = j == Ly-1 ? 0 : j+1;
        ssSum += s[i][j]*(s[iNext][j] + s[i][jNext]);
    }
    double e = -(J*ssSum + H*sSum)/N;
    eSum += e;
    eSqdSum += e * e;
    ++nSum;
}

double eAve;                // energia media por spin
double eError;              // incertidumbre por Monte Carlo

void computeAverages() {
    eAve = eSum / nSum;
    eError = eSqdSum / nSum;
    eError = sqrt(eError - eAve*eAve);
    eError /= sqrt(double(nSum));
}

int main() {

    cout << " Modelo Bidimensional de Ising - Alg. de Swendsen-Wang\n"
         << " -----------------------------------------------------\n"
         << " Introduzca numero de spines L en cada direccion: ";
    cin >> Lx;
    Ly = Lx;
    N = Lx * Ly;
    cout << " Introduzca temperatura T: ";
    cin >> T;
    cout << " Introduzca numero de Pasos de Monte Carlo a realizar: ";
    int MCSteps;
    cin >> MCSteps;

    initialize();
    initializeClusterVariables();

    int thermSteps = MCSteps / 5;
    cout << " Realizando " << thermSteps
         << " pasos de termalizacion ..." << flush;
    for (int i = 0; i < thermSteps; i++)
        oneMonteCarloStep();
    cout << " Listo!\n Realizando pasos de produccion ..." << flush;

    initializeObservables();
    for (int i = 0; i < MCSteps; i++) {
        oneMonteCarloStep();
        measureObservables();
    }
    cout << " Listo!" << endl;
    computeAverages();
    cout << " Energia por spin = " << eAve << " +- " << eError << endl;
    system("pause");
}

