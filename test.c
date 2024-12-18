#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

// Sructure des matrices creuses au format column compress
typedef struct {
    double *values;     // tableaux des valeurs non nulles
    int *rowIndexes;    // tableaux des indices des lignes des valeurs non nulles triés par colonne
    int *colPointers;   // Début de chaque colonne dans le tableau rowIndexes
    int nRows;          // Nombre de lignes dans la matrice
    int nCols;          // Nombre de colonnes dans la matrice
    int nnz;            // Nombre de valeurs non nulles
} SparseMatrix;

int extractTriangularMtx(SparseMatrix *M, SparseMatrix *U, bool Upper){
    
    U->nCols = M->nCols;
    U->nRows = M->nCols;
    U->nnz = 0;
    U->values = malloc(M->nnz * sizeof(double));
    if(U->values == NULL){
        printf("Erreur allocation mémoire\n");
        return 1;
    }
    U->rowIndexes = malloc(M->nnz * sizeof(int));
    if(U->rowIndexes == NULL){
        printf("Erreur allocation mémoire\n");
        return 1;
    }
    U->colPointers = malloc((M->nCols + 1) * sizeof(int));
    if(U->colPointers == NULL){
        printf("Erreur allocation mémoire\n");
        return 1;
    }
    //On met uniquement le premier pointeur à 0 parce que la liste des valeurs de la colonne 0 commence forcément à la cellule 0 de rowIndexes.
    U->colPointers[0] = 0;
    //Je parcours chaque pointeur de colonne 
    for(int col = 0; col < U->nCols; col++){
        //Je boucle sur les valeurs de cette colonne. 
        //Dès qu'on dépasse l'indice de cellule à partir duquel on change de colonne, on sort de la boucle pour passer à la colonne suivante.
        for (int i = M->colPointers[col]; i < M->colPointers[col+1]; i++){
            //La valeur contenu dans la cellule i de rowIndexes, c'est la ligne où se situe l'élément courant de la colonne qu'on est entrain de parcourir.
            int row = M->rowIndexes[i];
            double value = M->values[i];
            //On ne veut garder dans notre matrice que ce qui se trouve sur la diagonale et en dessous donc row >= col
            //printf(" i %d, row %d, col %d, value %lf\n", i,row,col,value);
            if(!Upper){
                if(row >= col){
                U->values[U->nnz] = value;
                U->rowIndexes[U->nnz] = row;
                U->nnz++;
            }
            }else{
                if(row < col){
                    U->values[U->nnz] = value;
                    U->rowIndexes[U->nnz] = row;
                    U->nnz++;
                }
            }
        }
        //A chaque fin d'itération de valeurs d'une colonne, on met l'indice à partir duquel on compte les éléments dans la colonne suivante dans le pointeur suivant.
        U->colPointers[col + 1] = U->nnz;
    }
    // Comme on a pas le même nombre de valeur non nulles que dans la matrice d'origine (vu qu'on garde que les valeurs sous la diag et la diag) on doit réallouer le nombre correcte de cellule.
    U->values = realloc(U->values, U->nnz * sizeof(double));
    U->rowIndexes = realloc(U->rowIndexes, U->nnz * sizeof(int));

    return 0;

}

int loadSparseMatrix(SparseMatrix *matrix, const char *sysMtx){

    //Initialise matrice
    matrix->values = NULL;
    matrix->rowIndexes = NULL;
    matrix->colPointers = NULL;
    matrix->nRows = 0;
    matrix->nCols = 0;
    matrix->nnz = 0;

    // Ouverture du fichier qui renvoie 1 si erreur (voir énoncé)
    FILE *mtxFile = fopen(sysMtx, "r");
    if (!mtxFile) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier %s\n", sysMtx);
        return 1;
    }

    // récupérer les dimensions de la matrice et le nombre de valeur non nulles
    fscanf(mtxFile, "%d %d %d", &matrix->nRows, &matrix->nCols, &matrix->nnz);

    // Allocation de mémoire pour les variables qui implémentent la structure column compress
    matrix->values = malloc(matrix->nnz * sizeof(double));
    matrix->rowIndexes = malloc(matrix->nnz * sizeof(int));
    matrix->colPointers = malloc((matrix->nCols + 1) * sizeof(int));

    // Tableau temporaire pour compter les éléments dans chaque colonne et suivre le remplissage du tableau rowIndices
    int *colCounts = calloc(matrix->nCols, sizeof(int));

    // on lit les triplets (ligne, colonne, valeur) et on compte le nombre d'élément dans chaque colonne
    // de la matrice parce qu'il faut pouvoir les placer à la bonne cellule dans la tableau des indices trié par colonne (je sais c'est pas évident pose moi des questions)
    int *rowTmp = malloc(matrix->nnz * sizeof(int));
    int *colTmp = malloc(matrix->nnz * sizeof(int));
    int i = 0;
    for (i = 0; i < matrix->nnz; i++) {
        fscanf(mtxFile, "%d %d %lf", &rowTmp[i], &colTmp[i], &matrix->values[i]);
        rowTmp[i]--; // Passer à un indexage 0
        colTmp[i]--;
        colCounts[colTmp[i]]++;
    }
    fclose(mtxFile);

    // Ici on veut que le pointeur vers les colonnes indique ou commence chaque colonne dans le tableau rowIndices
    matrix->colPointers[0] = 0;
    for (int i = 0; i < matrix->nCols; i++) {
        matrix->colPointers[i + 1] = matrix->colPointers[i] + colCounts[i];
    }
    // matrix->colPointers[matrix->nCols]++;

    // Remplir rowIndices en triant les valeurs par colonne 
    int *colOffsets = calloc(matrix->nCols, sizeof(int));
    for (int i = 0; i < matrix->nnz; i++) {
        int col = colTmp[i];
        int dest = matrix->colPointers[col] + colOffsets[col];
        matrix->rowIndexes[dest] = rowTmp[i];
        colOffsets[col]++;
    }

    free(colCounts);
    free(colOffsets);
    free(rowTmp);
    free(colTmp);

    return 0;

}

void freeSparseMatrix(SparseMatrix *matrix){
    if(matrix->values)
        free(matrix->values);
    if(matrix->rowIndexes)
        free(matrix->rowIndexes);
    if(matrix->colPointers)
        free(matrix->colPointers);

}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <inputMatrix.mtx> \n", argv[0]);
        return 1;
    }

    const char *sysMtx = argv[1];

    // Initialisation de la matrice et des vecteurs
    SparseMatrix A;
    loadSparseMatrix(&A, sysMtx);
    

   /*  // Affichage pour vérifier si mon code est bon
    printf("Matrice au format CCS de taile %dx%d:\n", A.nRows, A.nCols);
    for (int i = 0; i < A.nnz; i++) {
        printf("Valeur : %lf, Ligne: %d et i vaut %d \n", A.values[i], A.rowIndexes[i], i);
    }
    for (int i = 0; i <= A.nCols; i++) {
        printf("ColPointer[%d]: %d\n", i, A.colPointers[i]);
    } */

    SparseMatrix L; 
    extractTriangularMtx(&A, &L, false);

    // Affichage pour lower matrix
    printf("LOWER MATRIX\n");
    printf("Matrice au format CCS de taile %dx%d:\n", L.nRows, L.nCols);
    for (int i = 0; i < L.nnz; i++) {
        printf("Valeur : %lf, Ligne: %d et i vaut %d \n", L.values[i], L.rowIndexes[i], i);
    }
    for (int i = 0; i <= L.nCols; i++) {
        printf("ColPointer[%d]: %d\n", i, L.colPointers[i]);
    }
    SparseMatrix U;
    extractTriangularMtx(&A, &U, true);

    // Affichage pour upper matrix
    printf("UPPER MATRIX\n");
    printf("Matrice au format CCS de taile %dx%d:\n", U.nRows, U.nCols);
    for (int i = 0; i < U.nnz; i++) {
        printf("Valeur : %lf, Ligne: %d et i vaut %d \n", U.values[i], U.rowIndexes[i], i);
    }
    for (int i = 0; i <= U.nCols; i++) {
        printf("ColPointer[%d]: %d\n", i, U.colPointers[i]);
    }
    

    freeSparseMatrix(&A);
    freeSparseMatrix(&L);
    freeSparseMatrix(&U);
    
}