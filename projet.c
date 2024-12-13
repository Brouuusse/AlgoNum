#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

// Sructure des matrices creuses au format column compress
typedef struct {
    double *values;     // tableaux des valeurs non nulles
    int *rowIndexes;    // tableaux des indices des lignes des valeurs non nulles triés par colonne
    int *colPointers;   // Début de chaque colonne dans le tableau rowIndexes
    int nRows;          // Nombre de lignes dans la matrice
    int nCols;          // Nombre de colonnes dans la matrice
    int nnz;            // Nombre de valeurs non nulles
} SparseMatrix;

// void solveUpperTriangularCCS(SparseMatrix *U, SparseMatrix *b, SparseMatrix *x) {
//     // Initialisation du vecteur solution à 0
//     for (int i = 0; i < U->nCols; i++) {
//         x->values[i] = 0.0;
//     }

//     // Résolution colonne par colonne (triangulaire supérieur)
//     for (int col = U->nCols - 1; col >= 0; col--) {
//         // Chercher l'index où commence cette colonne
//         int startIdx = U->colPointers[col];
//         int endIdx = U->colPointers[col + 1];

//         // Trouver la valeur diagonale (doit exister dans une matrice triangulaire supérieure)
//         double diagValue = 0.0;
//         for (int idx = startIdx; idx < endIdx; idx++) {
//             if (U->rowIndexes[idx] == col) {
//                 diagValue = U->values[idx];
//                 break;
//             }
//         }

//         if (diagValue == 0.0) {
//             fprintf(stderr, "Erreur : la matrice est singulière (élément diagonal nul en colonne %d)\n", col);
//             exit(1);
//         }

//         // Calcul de x[col]
//         double sum = 0.0;
//         for (int idx = startIdx; idx < endIdx; idx++) {
//             int row = U->rowIndexes[idx];
//             if (row > col) {
//                 sum += U->values[idx] * x->values[row];
//             }
//         }

//         // Rechercher si b contient une valeur pour cette ligne
//         double bValue = 0.0;
//         for (int i = 0; i < b->nnz; i++) {
//             if (b->rowIndexes[i] == col) {
//                 bValue = b->values[i];
//                 break;
//             }
//         }

//         x->values[col] = (bValue - sum) / diagValue;
//     }
// }
void freeSparseMatrix(SparseMatrix *matrix){
    if(matrix->values)
        free(matrix->values);
    if(matrix->rowIndexes)
        free(matrix->rowIndexes);
    if(matrix->colPointers)
        free(matrix->colPointers);

}

bool isInRowIndexes(SparseMatrix m, int rowToSearch){
    bool isInside = false;
    for (int i = 0; i < m.nnz; i++){
        if(m.rowIndexes[i]==rowToSearch){
            isInside = true;
        }
    }
    return isInside;
}

int matchedRow(SparseMatrix m, int rowToSearch){
    int i =0;
    while(i < m.nnz && m.rowIndexes[i]!=rowToSearch){
        i++;
    }
    return i;
}

void solveLowerTriangular(SparseMatrix A, SparseMatrix b, SparseMatrix *x) {
    // Initialisation de la solution
    double *sum = malloc(sizeof(double) * A.nCols);
    for (int i = 0; i < A.nRows; i++) {
        x->values[i] = 0.0;
        sum[i] = 0.0;
    }
    int bRow;
    // Le problème actuel c'est qu'on ne parcours pas b correctement. 
    // Résolution par substitution avant
    for (int col = 0; col < A.nCols; col++) {
        // Parcourir les éléments de la colonne courante
        for (int i = A.colPointers[col]; i < A.colPointers[col + 1]; i++) {
            //printf("Colonne %d \n", col);
            int row = A.rowIndexes[i];
            //printf("La ligne actuelle est %d et la colonne actuelle est %d et i vaut %d\n", row, col, i);
            if (row > col) {
                // Accumuler les contributions pour les termes précédents
                sum[row] += A.values[i] * x->values[col];
                //printf("La ligne actuelle est %d et la colonne actuelle est %d \n", row, col);
                //printf("La somme de la ligne %d vaut : %lf et on est en train de lui ajouter la valeur %lf multipliee par x[%d] = %lf \n\n", row, sum[row], A.values[i], col, x->values[col]);
            } else if (row == col) {
                // Résoudre pour la diagonale
                if(isInRowIndexes(b, row)){
                    bRow = matchedRow(b, row);
                    //printf("b vaut = %lf à la cellule %d\n", b.values[bRow], bRow);
                    //printf("La ligne actuelle est %d et la colonne actuelle est %d. La valeur du tableau de valeur est %lf, la valeur de b est %lf et la somme vaut %lf\n",row, col, A.values[i], b.values[bRow], sum[col]);
                    x->values[row] = (b.values[bRow] - sum[col]) / A.values[i];
                    //printf("x[%d] = %lf \n\n", row, x->values[row]);
                }else{
                    //printf("La valeur du tableau de valeur est %lf, b vaut 0 car l'indice de ligne de A (%d) et b (%d) ne match pas et la somme vaut %lf \n", A.values[i], A.rowIndexes[i], b.rowIndexes[i], sum[col]);
                    x->values[row] = (0 - sum[col]) / A.values[i];
                    //printf("x[%d] = %lf \n\n", row, x->values[row]);
                }
            }
        }
    }
    printf("Les solutions sont:\n");
    for (int i = 0; i < A.nRows; i++){
        printf("- x[%d] = %.17lf \n", i, x->values[i]);
    }
    free(sum);
}

int extractTriangularMtx(SparseMatrix *M, SparseMatrix *U, bool Upper){
    //Upper = true si on extrait la mat triang sup 
    //      = false si on extraut la mat triang inf (avec diag)
    //printf("est rentré ds extr\n");
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
    //printf("avant boucle extract\n");
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

void multiplyUx(SparseMatrix *U, double *x, double *Ux){
    
    for(int i = 0; i < U->nRows; i++){
        Ux[i] = 0.0;
    }

    for (int col = 0; col < U->nCols; col++){
        for(int idx = U->colPointers[col]; idx < U->colPointers[col + 1]; idx++){
            int row = U->rowIndexes[idx];
            Ux[row] += U->values[idx] * x[col]; 
        }
    }
}

void substractbUx(double *b, double *Ux, int vectorLength, double *b_Ux){
    for(int i = 0; i < vectorLength; i++){
        b_Ux[i] = b[i] - Ux[i];
    }
}

void fromSparsetoDouble(SparseMatrix *matrix, double *vector){
    for(int i = 0; i < matrix->nRows; i++){
        vector[i] = 0.0;
    }

    for(int i = 0; i < matrix->nnz; i++){
        vector[matrix->rowIndexes[i]] = matrix->values[i];
    }
}

void fromDoubletoSparse(double *vector, SparseMatrix *matrix, int vectorLength) {

    int count = 0;
    for(int i = 0; i < vectorLength; i++){
        if(vector[i] != 0.0){
            count++;
        }
    }
    matrix->nCols = 1;
    matrix->nRows = vectorLength;
    matrix->nnz = count;

    matrix->values = malloc(matrix->nnz * sizeof(double));
    if (matrix->values == NULL) {
        printf("Erreur allocation mémoire pour les valeurs\n");
        return;
    }
    matrix->rowIndexes = malloc(matrix->nnz * sizeof(int));
    if (matrix->rowIndexes == NULL) {
        printf("Erreur allocation mémoire pour les indices de ligne\n");
        free(matrix->values);
        return;
    }
    matrix->colPointers = malloc((matrix->nCols + 1) * sizeof(int));
    if (matrix->colPointers == NULL) {
        printf("Erreur allocation mémoire pour les pointeurs de colonne\n");
        free(matrix->values);
        free(matrix->rowIndexes);
        return;
    }
    // on parcours le vecteur en incrémentant au fur et à mesure le nombre de valeurs non nulles et en garnissant le tableau valeur et rowIndexes
    int idx = 0;
    for (int i = 0; i < vectorLength; i++) {
        if (vector[i] != 0.0) {
            matrix->values[idx] = vector[i];
            matrix->rowIndexes[idx] = i;
            idx++;
        }
    }

    // La première colonne commence à l'indice 0
    matrix->colPointers[0] = 0;
    // Une seule colonne contient tous les éléments non nuls
    matrix->colPointers[1] = matrix->nnz; 
}

bool converge(double *x_curr, double *x_next, double precision, int vectorLength){
    double norm = 0.0;
    for (int i = 0; i < vectorLength; i++){
        double diff = x_next[i] - x_curr[i];
        norm += diff * diff;
    }

    norm = sqrt(norm);

    if(norm > 200){
        //stopper si converge pas
        //ne pas print les sol dans un fichier si on passe par ici parce converge pas
        return true;
    }

    return (norm<precision);
}

int resolutionGS(SparseMatrix *A, SparseMatrix *b, double precision){
    //printf("est rentré dans res\n");
    bool hasConverged = false;
    //2) on crée nos matrices (L + D) et U et un vecteur x et b
    SparseMatrix L;
    if(extractTriangularMtx(A, &L, false)){
        printf("Erreur allocation mémoire\n");
        return 1;
    }
    SparseMatrix U;
    if(extractTriangularMtx(A, &U, true)){
        printf("Erreur allocation mémoire\n");
        freeSparseMatrix(&L);
        return 1;
    }
    //printf("a fini extrait\n");
    double *x_curr = malloc(A->nRows * sizeof(double));
    if(x_curr == NULL){
        freeSparseMatrix(&L);
        freeSparseMatrix(&U);
        printf("Erreur allocation mémoire\n");
        return 1;
    }
    for(int i = 0; i < A->nRows; i++){
        x_curr[i] = 0.0; // on inititalise à zéro, ceci est x^0
    }

    double *b_vector = malloc(b->nRows * sizeof(double));
    if(b_vector == NULL){
        printf("Erreur allocation mémoire\n");
        freeSparseMatrix(&L);
        freeSparseMatrix(&U);
        free(x_curr);
        return 1;
    }
    fromSparsetoDouble(b, b_vector);

    double *Ux = malloc(U.nRows * sizeof(double));
    if(Ux == NULL){
        printf("Erreur allocation mémoire\n");
        freeSparseMatrix(&L);
        freeSparseMatrix(&U);
        free(x_curr);
        free(b_vector);
        return 1;
    }
    for(int i = 0; i < U.nRows; i++){
        Ux[i] = 0.0; // on inititalise à zéro
    }

    double *b_Ux = malloc(U.nRows * sizeof(double));
    if(b_Ux == NULL){
        printf("Erreur allocation mémoire\n");
        freeSparseMatrix(&L);
        freeSparseMatrix(&U);
        free(x_curr);
        free(b_vector);
        free(Ux);
        return 1;
    }
    for(int i = 0; i < U.nRows; i++){
        b_Ux[i] = 0.0; // on inititalise à zéro
    }
    //boucle tant que ça converge pas
    while (!hasConverged) {
        //printf("est rentré boucle conv\n");
        // Calculer U * x^k
        multiplyUx(&U, x_curr, Ux);

        // Calculer b' = b - U * x^k
        substractbUx(b_vector, Ux, A->nRows, b_Ux);

        // Remettre "vecteur" b' en "matrice creuse"
        SparseMatrix b_prime;
        fromDoubletoSparse(b_Ux, &b_prime, A->nRows);
        // Initialiser x_next_tmp
        SparseMatrix x_next_tmp;
        x_next_tmp.nRows = A->nRows;
        x_next_tmp.nCols = 1;
        x_next_tmp.nnz = 0; //
        x_next_tmp.values = malloc(A->nRows * sizeof(double));
        x_next_tmp.rowIndexes = malloc(A->nRows * sizeof(int));
        x_next_tmp.colPointers = malloc((x_next_tmp.nCols + 1) * sizeof(int));

        // Résoudre système L * x^(k+1) = b'
        solveLowerTriangular(L, b_prime, &x_next_tmp);

        //actualisez donnée de matrice x_next_tmp
        x_next_tmp.nnz = A->nRows;
        for(int i = 0; i < A->nRows; i++){
            x_next_tmp.rowIndexes[i] = i;
        }
        x_next_tmp.colPointers[0] = 0;
        x_next_tmp.colPointers[1] = A->nRows;

        double *x_next = malloc(A->nRows * sizeof(double));
        if(x_next == NULL){
            printf("Erreur allocation mémoire\n");
            freeSparseMatrix(&L);
            freeSparseMatrix(&U);
            freeSparseMatrix(&b_prime);
            freeSparseMatrix(&x_next_tmp);
            free(x_curr);
            free(b_vector);
            free(Ux);
            return 1;
        } 
        for(int i = 0; i < A->nRows; i++){
            x_next[i] = 0.0;
        }
        fromSparsetoDouble(&x_next_tmp, x_next);

        // Vérifier la convergence
        hasConverged = converge(x_curr, x_next, precision, A->nRows);

        // Si convergence pas atteinte, mettre à jour x_curr avec x_next_tmp
        if (!hasConverged) {
           for(int i = 0; i < A->nRows; i++){
            x_curr[i] = x_next[i];
           }
        }

        // Libérer b_prime pour éviter les fuites mémoire
        freeSparseMatrix(&b_prime);
        freeSparseMatrix(&x_next_tmp);
        free(x_next);
    }


    freeSparseMatrix(&L);
    freeSparseMatrix(&U);
    free(x_curr);
    free(b_vector);
    free(Ux);
    free(b_Ux);

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

    /* // Affichage pour vérifier si mon code est bon
    printf("Matrice au format CCS de taile %dx%d:\n", A.nRows, A.nCols);
    for (int i = 0; i < A.nnz; i++) {
        printf("Valeur : %lf, Ligne: %d et i vaut %d \n", A.values[i], A.rowIndexes[i], i);
    }
    for (int i = 0; i <= A.nCols; i++) {
        printf("ColPointer[%d]: %d\n", i, A.colPointers[i]);
    }
    printf("Vecteur de taille %dx%d\n", b.nRows, b.nCols);
    for(int i = 0; i < b.nnz; i++){
        printf("Valeur : %lf, Ligne : %d\n", b.values[i], b.rowIndexes[i]);
    } */
    SparseMatrix x;
    x.values = malloc(sizeof(double)*A.nRows);
    x.colPointers = malloc(sizeof(int) * A.nCols);
    x.rowIndexes = malloc(sizeof(int) * A.nnz);

    resolutionGS(&A, &b, precision);
    

    freeSparseMatrix(&A);
    freeSparseMatrix(&b);
    freeSparseMatrix(&x);
}
