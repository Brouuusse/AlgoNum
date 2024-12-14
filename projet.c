#include <stdio.h>
#include <stdlib.h>

// Sructure des matrices creuses au format column compress
typedef struct {
    double *values;     // tableaux des valeurs non nulles
    int *rowIndexes;    // tableaux des indices des lignes des valeurs non nulles triés par colonne
    int *colPointers;   // Début de chaque dans le tableau rowIndexes
    int nRows;          // Nombre de lignes dans la matrice
    int nCols;          // Nombre de colonnes dans la matrice
    int nnz;            // Nombre de valeurs non nulles
} SparseMatrix;

void solveUpperTriangularCCS(SparseMatrix *U, SparseMatrix *b, SparseMatrix *x) {
    // Initialisation du vecteur solution à 0
    for (int i = 0; i < U->nCols; i++) {
        x->values[i] = 0.0;
    }

    // Résolution colonne par colonne (triangulaire supérieur)
    for (int col = U->nCols - 1; col >= 0; col--) {
        // Chercher l'index où commence cette colonne
        int startIdx = U->colPointers[col];
        int endIdx = U->colPointers[col + 1];

        // Trouver la valeur diagonale (doit exister dans une matrice triangulaire supérieure)
        double diagValue = 0.0;
        for (int idx = startIdx; idx < endIdx; idx++) {
            if (U->rowIndexes[idx] == col) {
                diagValue = U->values[idx];
                break;
            }
        }

        if (diagValue == 0.0) {
            fprintf(stderr, "Erreur : la matrice est singulière (élément diagonal nul en colonne %d)\n", col);
            exit(1);
        }

        // Calcul de x[col]
        double sum = 0.0;
        for (int idx = startIdx; idx < endIdx; idx++) {
            int row = U->rowIndexes[idx];
            if (row > col) {
                sum += U->values[idx] * x->values[row];
            }
        }

        // Rechercher si b contient une valeur pour cette ligne
        double bValue = 0.0;
        for (int i = 0; i < b->nnz; i++) {
            if (b->rowIndexes[i] == col) {
                bValue = b->values[i];
                break;
            }
        }

        x->values[col] = (bValue - sum) / diagValue;
    }
}

void solveLowerTriangular(SparseMatrix A, SparseMatrix b, SparseMatrix *x) {
    // Initialisation de la solution
    double *sum = malloc(sizeof(double) * A.nCols);
    for (int i = 0; i < A.nRows; i++) {
        x->values[i] = 0.0;
        sum[i] = 0.0;
    }

    // Résolution par substitution avant
    for (int col = 0; col < A.nCols; col++) {

        // Parcourir les éléments de la colonne courante
        for (int i = A.colPointers[col]; i < A.colPointers[col + 1]; i++) {
            int row = A.rowIndexes[i];

            if (row > col) {
                // Accumuler les contributions pour les termes précédents
                sum[row] += A.values[i] * x->values[col];
                printf("La ligne actuelle est %d et la colonne actuelle est %d \n", row, col);
                printf("La somme de la ligne %d vaut : %lf et on est en train de lui ajouter la valeur %lf multipliee par x[%d] = %lf \n\n", row, sum[row], A.values[i], row, x->values[row]);
            } else if (row == col) {
                // Résoudre pour la diagonale
                printf("La valeur du tableau de valeur est %lf, la valeur de b est %lf et la somme vaut %lf\n", A.values[i], b.values[col], sum[col]);
                x->values[row] = (b.values[col] - sum[col]) / A.values[i];
                printf("x[%d] = %lf \n\n", row, x->values[row]);
            }
        }
    }
    printf("Les solutions sont:\n");
    for (int i = 0; i < A.nRows; i++){
        printf("- x[%d] = %lf \n", i, x->values[i]);
    }
}


int loadSparseMatrix(SparseMatrix *matrix, const char* sysMtx){

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
    for (int i = 0; i < matrix->nnz; i++) {
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

// int loadSparseVector(SparseVector *vector, const char* bMtx, SparseMatrix *matrix){


//     //Création du vecteur creux.
//     FILE *vctFile = fopen(bMtx, "r");
//     if (!vctFile) {
//         fprintf(stderr, "Erreur lors de l'ouverture du fichier %s\n", bMtx);
//         return 1;
//     }
//     //lecture de l'en-tête (j'ai ignoré la valeur du milieu parce que comme c'est un vecteur la colonne vaudra toujours 1 donc pas utile de la stocker)
//     fscanf(vctFile, "%d %*d %d", &vector->size, &vector->nnz);
//     vector->indexes = malloc(vector->nnz*sizeof(int));
//     vector->values = malloc(vector->nnz*sizeof(double));

//     // je garnis les valeurs et les lignes associés à ses valeurs dans le tableau d'indice et le tableau de valeur de la structure vecteur
//     for (int i = 0; i < vector->nnz; i++){
//         fscanf(vctFile, "%d %*d %lf", &vector->indexes[i], &vector->values[i]);
//     }
//     fclose(vctFile);

//     return 0;

// }

// void freeSparseVector(SparseVector *vector){
//     if(vector->indexes)
//         free(vector->indexes);
//     if(vector->values)
//         free(vector->values);

// }


int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <inputMatrix.mtx> <inputRhs.mtx> <precision>\n", argv[0]);
        return 1;
    }

    const char *sysMtx = argv[1];
    const char *bMtx = argv[2];
    double precision = atof(argv[3]);

    // Initialisation de la matrice et des vecteurs
    SparseMatrix A;
    loadSparseMatrix(&A, sysMtx);
    // Initialisation du vecteur
    SparseMatrix b;
    loadSparseMatrix(&b, bMtx);

    // Affichage pour vérifier si mon code est bon
    printf("Matrice au format CCS :\n");
    for (int i = 0; i < A.nnz; i++) {
        printf("Valeur : %lf, Ligne: %d\n", A.values[i], A.rowIndexes[i]);
    }
    for (int i = 0; i <= A.nCols; i++) {
        printf("ColPointer[%d]: %d\n", i, A.colPointers[i]);
    }
    printf("Vecteur\n");
    for(int i = 0; i < b.nnz; i++){
        printf("Valeur : %lf, Ligne : %d\n", b.values[i], b.rowIndexes[i]);
    }
    SparseMatrix x;
    x.values = malloc(sizeof(double)*A.nRows);
    solveLowerTriangular(A, b, &x);

    freeSparseMatrix(&A);
    freeSparseMatrix(&b);
    freeSparseMatrix(&x);

}
