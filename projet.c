#include <stdio.h>
#include <string.h>
#include <stdlib.h>
// matrice au format Matrix Market 
typedef struct {
    double *values;  // Tableau des valeurs non nulles de la matrice creuse
    int *rows;       // Tableau des indices des lignes des valeurs non nulles
    int *cols;       // Tableau des indices des colonnes des valeurs non nulles
    int nRows;     // Nombre de lignes de la matrice 
    int nCols;     // Nombre de colonnes de la matrice
    int nnz;         // Nombre de valeurs non nulles
    //autre chose ? 
} Matrix;

// vecteur au format Matrix Market
typedef struct {
    double *values;  // Valeurs non nulles du vecteur
    int size;        // Taille du vecteur
    // autre chose ?
} Vector;

int main(int argc, char **argv){

    // Si on a pas fournit assez d'argument au lancement du programme
    if(argc != 4){
        fprintf(stderr, "Donnez le bon nombre d'arguments");
        return 1;
    }
    // Ici je stock les noms de fichiers envoyer en argument lors du lancement du programme dans des chaînes de caractères
    const char *sysMtx = argv[1];
    const char *bMtx = argv[2];
    double precision = atof(argv[3]);

    //On initialise les matrices nécessaires
    Matrix A;
    Vector b, x;

    //A partir d'ici, c'est le code pour créer la matrice creuse.
    FILE *mtxFile = fopen(sysMtx, "r");
    if(!mtxFile){
        fprintf(stderr, "Erreur lors de l'ouverture du fichier %s", sysMtx);
        return 1;
    }
    // Pour lire l'en-tête qui contient les dimensions de la matrice et le nombre d'élément non nuls (voir énoncé)
    fscanf(mtxFile, "%d %d %d", &A.nRows, &A.nCols, &A.nnz);

    // Pour allouer de la mémoire pour stocker les éléments non nuls
    A.values = malloc(A.nnz * sizeof(double));
    A.rows = malloc(A.nnz * sizeof(int));
    A.cols = malloc(A.nnz * sizeof(int));

    fprintf(stderr, "%d\n", A.nnz);
    //On stocke les valeurs non nuls dans notre structure.
    for(int i = 0; i < A.nnz; i++){
        fscanf(mtxFile, "%d %d %lf", &A.rows[i], &A.cols[i], &A.values[i]);
        fprintf(stderr, "Ligne : %d, colonne : %d, valeur : %lf\n", A.rows[i], A.cols[i], A.values[i]);
    }


    fclose(mtxFile);

    //reste à gérer le vecteur (et à tester le code) et je dois aussi trouver un moyen de créer le makefile. 
    // Si tu sais comment faire fait moi signe sinon je fais ça dans la journée.
    
}