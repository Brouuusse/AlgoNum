#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

// Sructure des matrices creuses au format column compress
typedef struct {
    double *values;     // tableaux des valeurs non nulles
    int *rowIndexes;    // tableaux des indices des lignes des valeurs non nulles triés par colonne
    int *colPointers;   // Début de chaque colonne dans le tableau rowIndexes
    int nRows;          // Nombre de lignes dans la matrice
    int nCols;          // Nombre de colonnes dans la matrice
    int nnz;            // Nombre de valeurs non nulles
} SparseMatrix;


void writeSolutionToFile(SparseMatrix *x) {
    FILE *file = fopen("solution.txt", "w");
    if (!file) {
        fprintf(stderr, "Erreur : impossible d'ouvrir le fichier solution.txt en écriture.\n");
        return;
    }
    fprintf(file, "%d %d %d\n", x->nRows, x->nCols, x->nnz);

    for (int i = 0; i < x->nnz; i++) {
        fprintf(file, "%d %d %.16g\n", x->rowIndexes[i] + 1, 1, x->values[i]);
    }

    fclose(file);

    printf("Solution écrite dans le fichier 'solution.txt'.\n");
}

void freeSparseMatrix(SparseMatrix *matrix){
    if(matrix->values)
        free(matrix->values);
    if(matrix->rowIndexes)
        free(matrix->rowIndexes);
    if(matrix->colPointers)
        free(matrix->colPointers);

}

bool findRowIndex(SparseMatrix m, int rowToSearch, int *matchedIndex) {
    // Effectue une recherche dans rowIndexes et retourne true si trouvé, avec l'indice correspondant
    for (int i = 0; i < m.nnz; i++) {
        if (m.rowIndexes[i] == rowToSearch) {
            *matchedIndex = i;
            return true;
        }
    }
    return false;
}

void solveLowerTriangular(SparseMatrix A, SparseMatrix b, SparseMatrix *x) {
    // Initialisation des solutions et des sommes
    double *sum = calloc(A.nRows, sizeof(double));
    if (!sum) {
        fprintf(stderr, "Erreur : allocation mémoire échouée.\n");
        return;
    }

    // Initialiser les valeurs de la solution
    for (int i = 0; i < A.nRows; i++) {
        x->values[i] = 0.0;
    }

    // Résolution par substitution avant
    for (int col = 0; col < A.nCols; col++) {
        for (int i = A.colPointers[col]; i < A.colPointers[col + 1]; i++) {
            int row = A.rowIndexes[i];

            if (row > col) {
                // Accumuler les contributions pour les termes précédents
                sum[row] += A.values[i] * x->values[col];
            } else if (row == col) {
                // Résoudre pour la diagonale
                int bRowIndex = -1;
                double bValue = 0.0;

                if (findRowIndex(b, row, &bRowIndex)) {
                    bValue = b.values[bRowIndex];
                }

                x->values[row] = (bValue - sum[col]) / A.values[i];
            }
        }
    }

    // Libérer la mémoire utilisée pour les sommes
    free(sum);
}

int extractTriangularMtx(SparseMatrix *M, SparseMatrix *U, bool Upper) {
    // Initialisation de la matrice U
    U->nCols = M->nCols;
    U->nRows = M->nRows;
    U->nnz = 0;

    // Allouer suffisamment de mémoire pour les valeurs et les index
    U->values = malloc(M->nnz * sizeof(double));
    U->rowIndexes = malloc(M->nnz * sizeof(int));
    U->colPointers = malloc((M->nCols + 1) * sizeof(int));

    // Vérifier si l'allocation a échoué
    if (!U->values || !U->rowIndexes || !U->colPointers) {
        fprintf(stderr, "Erreur d'allocation mémoire\n");
        free(U->values);
        free(U->rowIndexes);
        free(U->colPointers);
        return 1;
    }

    // Initialiser le premier pointeur de colonne
    U->colPointers[0] = 0;

    // Parcourir chaque colonne de la matrice M
    for (int col = 0; col < M->nCols; col++) {
        for (int i = M->colPointers[col]; i < M->colPointers[col + 1]; i++) {
            int row = M->rowIndexes[i];
            double value = M->values[i];

            // Filtrer les éléments selon la matrice triangulaire (Upper ou Lower)
            if ((Upper && row < col) || (!Upper && row >= col)) {
                U->values[U->nnz] = value;
                U->rowIndexes[U->nnz] = row;
                U->nnz++;
            }
        }
        // Mettre à jour les colPointers pour la colonne courante
        U->colPointers[col + 1] = U->nnz;
    }

    // Réduire la taille des allocations aux dimensions exactes
    U->values = realloc(U->values, U->nnz * sizeof(double));
    U->rowIndexes = realloc(U->rowIndexes, U->nnz * sizeof(int));

    // Vérifier si la réallocation a échoué
    if (!U->values || !U->rowIndexes) {
        fprintf(stderr, "Erreur de réallocation mémoire\n");
        free(U->values);
        free(U->rowIndexes);
        free(U->colPointers);
        return 1;
    }

    return 0;
}

void multiplyUx(SparseMatrix *U, double *x, double *Ux) {
    // Réinitialiser Ux à 0, sans boucle inutile
    memset(Ux, 0, U->nRows * sizeof(double));

    // Effectuer la multiplication
    for (int col = 0; col < U->nCols; col++) {
        double x_col = x[col];
        for (int idx = U->colPointers[col]; idx < U->colPointers[col + 1]; idx++) {
            int row = U->rowIndexes[idx];
            Ux[row] += U->values[idx] * x_col;
        }
    }
}

void substractbUx(double *b, double *Ux, int vectorLength, double *b_Ux){
    for(int i = 0; i < vectorLength; i++){
        b_Ux[i] = b[i] - Ux[i];
    }
}

void fromSparsetoDouble(SparseMatrix *matrix, double *vector) {
    // Réinitialiser le vecteur à 0 en utilisant memset
    memset(vector, 0, matrix->nRows * sizeof(double));

    // Remplir le vecteur avec les valeurs non nulles de la matrice
    for (int i = 0; i < matrix->nnz; i++) {
        vector[matrix->rowIndexes[i]] = matrix->values[i];
    }
}

void fromDoubletoSparse(double *vector, SparseMatrix *matrix, int vectorLength) {
    // Initialiser les dimensions de la matrice creuse
    matrix->nCols = 1;
    matrix->nRows = vectorLength;

    // Étape 1 : Compter les éléments non nuls (nnz) en une seule boucle
    int nnz = 0;
    for (int i = 0; i < vectorLength; i++) {
        if (vector[i] != 0.0) {
            nnz++;
        }
    }
    matrix->nnz = nnz;

    // Allocation mémoire pour les valeurs, indices de lignes et pointeurs de colonnes
    matrix->values = malloc(nnz * sizeof(double));
    matrix->rowIndexes = malloc(nnz * sizeof(int));
    matrix->colPointers = malloc(2 * sizeof(int));

    // Vérification d'allocation
    if (!matrix->values || !matrix->rowIndexes || !matrix->colPointers) {
        printf("Erreur allocation mémoire\n");
        free(matrix->values);
        free(matrix->rowIndexes);
        free(matrix->colPointers);
        return;
    }
    // Étape 2 : Remplir `values` et `rowIndexes` en une seule boucle
    int idx = 0;
    for (int i = 0; i < vectorLength; i++) {
        if (vector[i] != 0.0) {
            matrix->values[idx] = vector[i];
            matrix->rowIndexes[idx] = i;
            idx++;
        }
    }

    // Étape 3 : Remplir les pointeurs de colonnes
    matrix->colPointers[0] = 0;
    matrix->colPointers[1] = nnz;
}

bool converge(double *x_curr, double *x_next, double precision, int vectorLength) {
    double squaredNorm = 0.0;

    // Calcul de la norme au carré pour éviter l'utilisation de sqrt si ce n'est pas nécessaire
    for (int i = 0; i < vectorLength; i++) {
        double diff = x_next[i] - x_curr[i];
        squaredNorm += diff * diff;
        // Arrêt anticipé si la norme dépasse le seuil de non-convergence (200^2)
        if (squaredNorm > 40000.0) { 
            return true;
        }
    }

    // Comparer la norme au carré à la précision pour éviter une racine carrée inutile
    return (squaredNorm < precision * precision);
}

// Fonction pour libérer toutes les ressources
void freeAll(SparseMatrix *L, SparseMatrix *U, double *x_curr, double *b_vector, 
             double *Ux, double *b_Ux, SparseMatrix *b_prime, SparseMatrix *x_next_tmp) {
    if (L) freeSparseMatrix(L);
    if (U) freeSparseMatrix(U);
    if (x_curr) free(x_curr);
    if (b_vector) free(b_vector);
    if (Ux) free(Ux);
    if (b_Ux) free(b_Ux);
    if (b_prime) freeSparseMatrix(b_prime);
    if (x_next_tmp) freeSparseMatrix(x_next_tmp);
}

int resolutionGS(SparseMatrix *A, SparseMatrix *b, double precision) {
    bool hasConverged = false;

    // Extraction des matrices triangulaires L et U
    SparseMatrix L, U;
    if (extractTriangularMtx(A, &L, false)) {
        fprintf(stderr, "Erreur allocation mémoire pour L.\n");
        return 1;
    }
    if (extractTriangularMtx(A, &U, true)) {
        fprintf(stderr, "Erreur allocation mémoire pour U.\n");
        freeSparseMatrix(&L);
        return 1;
    }

    // Allocation des vecteurs nécessaires
    double *x_curr = malloc(A->nRows * sizeof(double));
    if (!x_curr) {
        fprintf(stderr, "Erreur allocation mémoire pour x_curr.\n");
        freeSparseMatrix(&L);
        freeSparseMatrix(&U);
        return 1;
    }
    memset(x_curr, 0, A->nRows * sizeof(double)); 

    double *b_vector = malloc(b->nRows * sizeof(double));
    if (!b_vector) {
        fprintf(stderr, "Erreur allocation mémoire pour b_vector.\n");
        freeSparseMatrix(&L);
        freeSparseMatrix(&U);
        free(x_curr);
        return 1;
    }
    fromSparsetoDouble(b, b_vector);

    double *Ux = malloc(U.nRows * sizeof(double));
    if (!Ux) {
        fprintf(stderr, "Erreur allocation mémoire pour Ux.\n");
        freeSparseMatrix(&L);
        freeSparseMatrix(&U);
        free(x_curr);
        free(b_vector);
        return 1;
    }
    memset(Ux, 0, U.nRows * sizeof(double));

    double *b_Ux = malloc(U.nRows * sizeof(double));
    if (!b_Ux) {
        fprintf(stderr, "Erreur allocation mémoire pour b_Ux.\n");
        freeSparseMatrix(&L);
        freeSparseMatrix(&U);
        free(x_curr);
        free(b_vector);
        free(Ux);
        return 1;
    }
    memset(b_Ux, 0, U.nRows * sizeof(double));

    // Préallocation des matrices temporaires
    SparseMatrix b_prime = {0};
    SparseMatrix x_next_tmp = {0};
    x_next_tmp.nRows = A->nRows;
    x_next_tmp.nCols = 1;
    x_next_tmp.values = malloc(A->nRows * sizeof(double));
    if (!x_next_tmp.values) {
        fprintf(stderr, "Erreur allocation mémoire pour x_next_tmp.values.\n");
        freeAll(&L, &U, x_curr, b_vector, Ux, b_Ux, &b_prime, &x_next_tmp);
        return 1;
    }
    x_next_tmp.rowIndexes = malloc(A->nRows * sizeof(int));
    if (!x_next_tmp.rowIndexes) {
        fprintf(stderr, "Erreur allocation mémoire pour x_next_tmp.rowIndexes.\n");
        freeAll(&L, &U, x_curr, b_vector, Ux, b_Ux, &b_prime, &x_next_tmp);
        return 1;
    }
    x_next_tmp.colPointers = malloc(2 * sizeof(int));
    if (!x_next_tmp.colPointers) {
        fprintf(stderr, "Erreur allocation mémoire pour x_next_tmp.colPointers.\n");
        freeAll(&L, &U, x_curr, b_vector, Ux, b_Ux, &b_prime, &x_next_tmp);
        return 1;
    }

    // Boucle principale de Gauss-Seidel
    while (!hasConverged) {
        // Calculer U * x^k
        multiplyUx(&U, x_curr, Ux);

        // Calculer b' = b - U * x^k
        substractbUx(b_vector, Ux, A->nRows, b_Ux);

        // Remettre "vecteur" b' en "matrice creuse"
        fromDoubletoSparse(b_Ux, &b_prime, A->nRows);

        // Résoudre système L * x^(k+1) = b'
        solveLowerTriangular(L, b_prime, &x_next_tmp);

        // Mise à jour de x_curr et vérification de la convergence
        double *x_next = x_next_tmp.values;
        hasConverged = converge(x_curr, x_next, precision, A->nRows);
        memcpy(x_curr, x_next, A->nRows * sizeof(double));
    }

    // Sauvegarde de la solution finale
    SparseMatrix x_final;
    fromDoubletoSparse(x_curr, &x_final, A->nRows);
    writeSolutionToFile(&x_final);

    // Libération des ressources
    freeAll(&L, &U, x_curr, b_vector, Ux, b_Ux, &b_prime, &x_next_tmp);
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
        rowTmp[i]--;
        colTmp[i]--;
        colCounts[colTmp[i]]++;
    }
    fclose(mtxFile);

    matrix->colPointers[0] = 0;
    for (int i = 0; i < matrix->nCols; i++) {
        matrix->colPointers[i + 1] = matrix->colPointers[i] + colCounts[i];
    }

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

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <inputMatrix.mtx> <inputRhs.mtx> <precision>\n", argv[0]);
        return 1;
    }
    int result;
    const char *sysMtx = argv[1];
    const char *bMtx = argv[2];
    double precision = atof(argv[3]);

    // Initialisation de la matrice et des vecteurs
    SparseMatrix A;
    loadSparseMatrix(&A, sysMtx);
    // Initialisation du vecteur
    SparseMatrix b;
    loadSparseMatrix(&b, bMtx);

    SparseMatrix x;
    x.values = malloc(sizeof(double)*A.nRows);
    x.colPointers = malloc(sizeof(int) * A.nCols);
    x.rowIndexes = malloc(sizeof(int) * A.nnz);

    result = resolutionGS(&A, &b, precision);
    freeSparseMatrix(&A);
    freeSparseMatrix(&b);
    freeSparseMatrix(&x);
    return result;
}