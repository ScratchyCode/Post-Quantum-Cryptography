// Coded by Pietro Squilla
// Analisi e test computazionale dello Shortest Vector Problem
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

#define DIMENSIONI 200
#define MAX_DIMENSIONI 1000000
#define MAX_THREAD 10
#define INIZIO -100
#define FINE 100
#define TOTAL_POINTS ((FINE - INIZIO + 1) * (FINE - INIZIO + 1) - 1)

typedef struct{
    int *vettore;
    double distanza;
} Vettore;

typedef struct{
    double **base;
    int **punti;
    int num_punti;
    int dimensioni;
    Vettore *min_vettore;
} Dati;

double calcolaDistanzaEuclidea(double *vettore, int dimensioni){
    double somma = 0.0;
    
    for(int i = 0; i < dimensioni; i++){
        somma += vettore[i] * vettore[i];
    }
    
    return sqrt(somma);
}

double *trasformaInBase(double **base, int *vettore, int dimensioni){
    double *vettoreTrasformato = malloc(dimensioni * sizeof(double));
    
    for(int i = 0; i < dimensioni; i++){
        vettoreTrasformato[i] = 0;
        for(int j = 0; j < dimensioni; j++){
            vettoreTrasformato[i] += base[i][j] * vettore[j];
        }
    }
    
    return vettoreTrasformato;
}

void aggiornaVettoreMinimo(double **base, int *vettoreCorrente, Vettore *minVettore, int dimensioni){
    // verifica se il vettore corrente è l'origine
    int isOrigine = 1;
    
    for(int i = 0; i < dimensioni; i++){
        if(vettoreCorrente[i] != 0){
            isOrigine = 0;
            break;
        }
    }
    
    // se il vettore corrente non è l'origine, aggiorna il vettore minimo
    if(!isOrigine){
        double *vettoreTrasformato = trasformaInBase(base, vettoreCorrente, dimensioni);
        double distanza = calcolaDistanzaEuclidea(vettoreTrasformato, dimensioni);
        free(vettoreTrasformato);
        if(distanza < minVettore->distanza){
            minVettore->distanza = distanza;
            for(int i = 0; i < dimensioni; i++){
                minVettore->vettore[i] = vettoreCorrente[i];
            }
        }
    }
}

void *trovaMinDistanza(void *arg){
    Dati *dati = (Dati*)arg;
    double **base = dati->base;
    int dimensioni = dati->dimensioni;
    Vettore *minVettore = dati->min_vettore;
    
    for(int i = 0; i < dati->num_punti; i++){
        aggiornaVettoreMinimo(base, dati->punti[i], minVettore, dimensioni);
    }
    
    return NULL;
}

int LUPDecompose(double **A, int N, double Tol, int *P){
    int i, j, k, imax; 
    double maxA, *ptr, absA;
    
    for(i = 0; i <= N; i++)
        P[i] = i;
    
    for(i = 0; i < N; i++){
        maxA = 0.0;
        imax = i;
        
        for(k = i; k < N; k++)
            if((absA = fabs(A[k][i])) > maxA){
                maxA = absA;
                imax = k;
            }
        
        if (maxA < Tol) return 0; 
        
        if(imax != i){
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;
            
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;
            
            P[N]++;
        }
        
        for(j = i + 1; j < N; j++){
            A[j][i] /= A[i][i];
            
            for(k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }
    
    return 1;
}

int main() {
    srand(time(NULL));
    int dimensioni = DIMENSIONI;
    double **base = malloc(dimensioni * sizeof(double*));
    
    // genero una base casuale
    for(int i = 0; i < dimensioni; i++){
        base[i] = malloc(dimensioni * sizeof(double));
        for(int j = 0; j < dimensioni; j++){
            base[i][j] = (double)(rand()) / RAND_MAX;
        }
    }
    
    int *P = malloc((dimensioni + 1) * sizeof(int));
    if(!LUPDecompose(base, dimensioni, 1e-7, P)){
        printf("Errore nella decomposizione LUP\n");
        return -1;
    }
    
    // preparo i dati per i thread
    int **punti = malloc(TOTAL_POINTS * sizeof(int*));
    int pointCounter = 0;
    
    for(int x = INIZIO; x <= FINE; x++){
        for(int y = INIZIO; y <= FINE; y++){
            if(x == 0 && y == 0) continue; // esclude l'origine
            punti[pointCounter] = malloc(dimensioni * sizeof(int));
            punti[pointCounter][0] = x;
            punti[pointCounter][1] = y;
            pointCounter++;
        }
    }
    
    Dati *dati = malloc(MAX_THREAD * sizeof(Dati));
    pthread_t *threads = malloc(MAX_THREAD * sizeof(pthread_t));
    
    Vettore *minVettori = malloc(MAX_THREAD * sizeof(Vettore));
    for(int i = 0; i < MAX_THREAD; i++){
        minVettori[i].vettore = malloc(dimensioni * sizeof(int));
        minVettori[i].distanza = INFINITY;
    }
    
    for(int i = 0; i < MAX_THREAD; i++){
        dati[i].base = base;
        dati[i].punti = &punti[i * (TOTAL_POINTS / MAX_THREAD)];
        dati[i].num_punti = TOTAL_POINTS / MAX_THREAD;
        dati[i].dimensioni = dimensioni;
        dati[i].min_vettore = &minVettori[i];
        pthread_create(&threads[i], NULL, trovaMinDistanza, (void*)&dati[i]);
    }
    
    for(int i = 0; i < MAX_THREAD; i++){
        pthread_join(threads[i], NULL);
    }
    
    Vettore minVettore;
    minVettore.vettore = malloc(dimensioni * sizeof(int));
    minVettore.distanza = INFINITY;
    
    for(int i = 0; i < MAX_THREAD; i++){
        if(minVettori[i].distanza < minVettore.distanza){
            minVettore.distanza = minVettori[i].distanza;
            for(int j = 0; j < dimensioni; j++){
                minVettore.vettore[j] = minVettori[i].vettore[j];
            }
        }
    }
    
    printf("Il vettore con la minima distanza Euclidea è: ");
    printf("[");
    for(int i = 0; i < dimensioni; i++){
        printf("%d,", minVettore.vettore[i]);
    }
    printf("]");
    printf("\nLa sua distanza è: %lf\n", minVettore.distanza);
    
    // libero la memoria
    for(int i = 0; i < dimensioni; i++){
        free(base[i]);
    }
    free(base);
    free(P);
    
    for(int i = 0; i < TOTAL_POINTS; i++){
        free(punti[i]);
    }
    free(punti);
    
    free(dati);
    free(threads);

    for(int i = 0; i < MAX_THREAD; i++){
        free(minVettori[i].vettore);
    }
    free(minVettori);
    
    free(minVettore.vettore);
    
    return 0;
}

// provare a stampare le componenti del vettore moltiplicato per la base
