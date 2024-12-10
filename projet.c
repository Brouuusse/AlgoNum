#include <stdio.h>
#include <stdlib.h>

// Sructure des matrices creuses au format column compress
typedef struct {
    double *values;     // tableaux des valeurs non nulles
    int *rowIndices;    // tableaux des indices des lignes des valeurs non nulles triés par colonne
    int *colPointers;   // Début de chaque dans le tableau rowIndices
    int nRows;          // Nombre de lignes dans la matrice
    int nCols;          // Nombre de colonnes dans la matrice
    int nnz;            // Nombre de valeurs non nulles
} SparseMatrix;

typedef struct {
    double *values;    // Valeurs non nulles
    int *indices;      // Indices des valeurs non nulles dans le vecteur
    int nnz;           // Nombre de valeurs non nulles
    int size;          // Taille totale du vecteur (utile pour référence)
} SparseVector;

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
    SparseVector b, x;

    // Ouvertur du fichier qui renvoie 1 si erreur (voir énoncé)
    FILE *mtxFile = fopen(sysMtx, "r");
    if (!mtxFile) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier %s\n", sysMtx);
        return 1;
    }

    // récupérer les dimensions de la matrice et le nombre de valeur non nulles
    fscanf(mtxFile, "%d %d %d", &A.nRows, &A.nCols, &A.nnz);

    // Allocation de mémoire pour les variables qui implémentent la structure column compress
    A.values = malloc(A.nnz * sizeof(double));
    A.rowIndices = malloc(A.nnz * sizeof(int));
    A.colPointers = malloc((A.nCols + 1) * sizeof(int));

    // Tableau temporaire pour compter les éléments dans chaque colonne et suivre le remplissage du tableau rowIndices
    int *colCounts = calloc(A.nCols, sizeof(int));

    // on lit les triplets (ligne, colonne, valeur) et on compte le nombre d'élément dans chaque colonne
    // de la matrice parce qu'il faut pouvoir les placer à la bonne cellule dans la tableau des indices trié par colonne (je sais c'est pas évident pose moi des questions)
    int *rowTmp = malloc(A.nnz * sizeof(int));
    int *colTmp = malloc(A.nnz * sizeof(int));
    for (int i = 0; i < A.nnz; i++) {
        fscanf(mtxFile, "%d %d %lf", &rowTmp[i], &colTmp[i], &A.values[i]);
        rowTmp[i]--; // Passer à un indexage 0
        colTmp[i]--;
        colCounts[colTmp[i]]++;
    }
    fclose(mtxFile);

    // Ici on veut que le pointeur vers les colonnes indique ou commence chaque colonne dans le tableau rowIndices
    A.colPointers[0] = 0;
    for (int i = 0; i < A.nCols; i++) {
        A.colPointers[i + 1] = A.colPointers[i] + colCounts[i];
    }

    // Remplir rowIndices en triant les valeurs par colonne 
    int *colOffsets = calloc(A.nCols, sizeof(int));
    for (int i = 0; i < A.nnz; i++) {
        int col = colTmp[i];
        int dest = A.colPointers[col] + colOffsets[col];
        A.rowIndices[dest] = rowTmp[i];
        colOffsets[col]++;
    }

    free(colCounts);
    free(colOffsets);
    free(rowTmp);
    free(colTmp);

    //Création du vecteur creux.
    FILE *vctFile = fopen(bMtx, "r");
    if (!vctFile) {
        fprintf(stderr, "Erreur lors de l'ouverture du fichier %s\n", bMtx);
        return 1;
    }
    //lecture de l'en-tête (j'ai ignoré la valeur du milieu parce que comme c'est un vecteur la colonne vaudra toujours 1 donc pas utile de la stocker)
    fscanf(vctFile, "%d %*d %d", &b.size, &b.nnz);
    b.indices = malloc(b.nnz*sizeof(int));
    b.values = malloc(b.nnz*sizeof(double));

    // je garnis les valeurs et les lignes associés à ses valeurs dans le tableau d'indice et le tableau de valeur de la structure vecteur
    for (int i = 0; i < b.nnz; i++){
        printf("Ah\n");
        fscanf(vctFile, "%d %*d %lf", &b.indices[i], &b.values[i]);
    }
    // Affichage pour vérifier si mon code est bon
    fprintf(stderr, "Matrice au format CCS :\n");
    for (int i = 0; i < A.nnz; i++) {
        fprintf(stderr, "Valeur: %lf, Ligne: %d\n", A.values[i], A.rowIndices[i]);
    }
    for (int i = 0; i <= A.nCols; i++) {
        fprintf(stderr, "ColPointer[%d]: %d\n", i, A.colPointers[i]);
    }
    fprintf(stderr, "Vecteur\n");
    for(int i = 0; i < b.nnz; i++){
        fprintf(stderr, "Valeur : %lf, Ligne : %d\n", b.values[i], b.indices[i]);
    }

    free(A.values);
    free(A.rowIndices);
    free(A.colPointers);
    free(b.indices);
    free(b.values);
    return 0;
}
