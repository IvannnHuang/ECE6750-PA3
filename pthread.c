#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <pthread.h>

#define BILLION 1000000000L

/*
Implements a CSR sparse matrix format.
*/
typedef struct {
    int m; // The number of rows
    int n; // The number of columns
    int nnz; // The number of nonzero values
    int *rowOffsets; // The row pointer array
    int *cindexes; // The column index array
    int *x; // Both this and the line below are coordinates used when reading from the .mtx file, and should not be used in computation.
    int *y;
    double *values; // The values array
}CSRdata; 

/*
Implements a COO sparse matrix format.
*/
typedef struct {
    int m; // The number of rows
    int n; // The number of columns
    int nnz; // The number of nonzero values
    int *row; // The row index array
    int *column; // The column index array
    double *values; // The values array
    int pointer; // A pointer used to indicate where the next value should be placed
}COOdata; 

/*
Worker thread all inputs
*/
typedef struct {
    CSRdata *a;
    CSRdata *b;
    COOdata *c;
    int *row_ptr;   // prefix sum of nnz_per_row
    int tid;        // thread id
    int p;          // total number of threads
} ThreadArg;

/*
Counts the number of nonzero values in each result row. Output:
1. nnz_per_row[i] = number of valid outputs in row i
2. row_ptr[i]     = starting offset of row i in COO output
   row_ptr[m]     = total nnz of result matrix
*/
void resNNZCount(CSRdata *a, CSRdata *b, int* nnz_per_row, int* row_ptr) {
    int n = b -> n;
    bool * marker = (bool *) calloc(n, sizeof(bool));

    // porefix sum starts at 0
    row_ptr[0] = 0;

    for (int i = 0; i < a -> m; i++) {
        int nnzCount = 0;
        for (int j = a -> rowOffsets[i]; j < a -> rowOffsets[i + 1]; j++) {
            int aCol = a -> cindexes[j];
            for (int k = b -> rowOffsets[aCol]; k < b -> rowOffsets[aCol + 1]; k++) {
                int bCol = b -> cindexes[k];
                if (!marker[bCol]) {
                    marker[bCol] = true;
                    nnzCount++;
                }
            }
        }

        nnz_per_row[i] = nnzCount;
        row_ptr[i + 1] = row_ptr[i] + nnzCount;
        
        memset(marker, 0, n * sizeof(bool));
    }
    free (marker);
}

/*
Fills in the CSR arrays in the CSRData struct. Used when converting the .mtx files read in to their 
CSR format for Gustavson's algorithm.
*/
void convert(CSRdata *a) {
    // initialize rowOffsets
    memset(a->rowOffsets, 0, (a->m + 1) * sizeof(int));

    for (int i = 0; i < a -> nnz; i++) {
        a->rowOffsets[a -> x[i] + 1] += 1;
    }
    int sum = 0;
    for (int i = 0; i <= a -> m; i++) {
        sum += a-> rowOffsets[i];
        a->rowOffsets[i] = sum;
    }
    int *temp = (int *)malloc((a->m + 1) * sizeof(int));
    memcpy(temp, a->rowOffsets, (a->m + 1) * sizeof(int));
    double *values_temp = (double *)malloc(a->nnz * sizeof(double));

    for (int i = 0; i < a->nnz; i++) {
        int row = a->x[i];
        int dest = temp[row]++;
        a->cindexes[dest] = a->y[i];
        values_temp[dest] = a->values[i];
    }
    free(temp);
    free(a->values);
    a->values = values_temp;
}

/*
Thread wroker for pthread:
- alternating rows
- each row writes to preassigned COO segment
*/
void *worker(void *arg) {
    ThreadArg * targ = (ThreadArg *) arg;
    CSRdata *a = targ -> a;
    CSRdata *b = targ -> b;
    COOdata *c = targ -> c;
    int * row_ptr = targ -> row_ptr; // prefix sum of nnz_per_row
    int tid = targ -> tid;  // thread id
    int p = targ -> p;      // total number of threads
    
    int n = b -> n;

    for (int i = tid; i < a->m; i += p) {
        // A temporary row to store the values collected in this iteration
        double *tempRow = (double *)calloc(n, sizeof(double));

        // Marks which elements of the row are valid
        int *marker = (int *)calloc(n, sizeof(int));

        // Gustavson accumulation for row i
        for (int j = a->rowOffsets[i]; j < a->rowOffsets[i+1]; j++) {
            int y = a->cindexes[j];     // k idx
            double val = a->values[j];

            for (int j2 = b->rowOffsets[y]; j2 < b->rowOffsets[y+1]; j2++) {
                int k = b->cindexes[j2];
                double val2 = b->values[j2];

                if (!marker[k]) {
                    marker[k] = 1;
                    tempRow[k] = val*val2;
                }
                else {
                    tempRow[k] += val * val2;
                }
            }
        }

        // Write row i into reserved COO segment
        int pos = row_ptr[i];
        for (int col = 0; col < n; col++) {
            if (marker[col]) {
                c -> row[pos] = i;
                c -> column[pos] = col;
                c -> values[pos] = tempRow[col];
                pos++;
            }
        }

        free(tempRow);
        free(marker);        
    }

    return NULL;

}

/*
Implements a pthread version of Gustavson's algorithm.
*/
void sparseMM_pthread(CSRdata *a, CSRdata *b, COOdata *c, int p) {
    int m = a -> m;
    int n = b -> n;
 
    // count output nnz and prefix sum per row
    int *nnz_per_row = (int *) malloc(m * sizeof(int));
    int *row_ptr = (int *) malloc((m + 1) * sizeof(int));

    resNNZCount(a, b, nnz_per_row, row_ptr);    
    int total_nnz = row_ptr[m];

    // allocate COO result
    c->m = m;
    c->n = n;
    c->pointer = 0;
    c->nnz = total_nnz;
    c->row = (int *)malloc(total_nnz * sizeof(int));
    c->column = (int *)malloc(total_nnz * sizeof(int));
    c->values = (double *)malloc(total_nnz * sizeof(double));

    // create threads
    pthread_t *threads = (pthread_t *) malloc(p * sizeof(pthread_t));
    ThreadArg * args = (ThreadArg *) malloc(p * sizeof(ThreadArg));

    for (int t = 0; t < p; t++) {
        args[t].a = a;
        args[t].b = b;
        args[t].c = c;
        args[t].row_ptr = row_ptr;
        args[t].tid = t;
        args[t].p = p;
        
        pthread_create(&threads[t], NULL, worker, &args[t]);
    }

    // join threads
    for (int t = 0; t < p; t++) {
        pthread_join(threads[t], NULL);
    }

    free(nnz_per_row);
    free(row_ptr);
    free(threads);
    free(args);
}

int main(int argc, char **argv) {
    FILE *fptr;

    if (argc < 3) {
        printf("Usage: %s matrix_file num_threads\n", argv[0]);
        return 1;
    }

    // Determines which matrix will be multiplied by itself.
    fptr = fopen(argv[1], "r");
    if (!fptr) {
        perror("fopen");
        return 1;
    }

    int p = atoi(argv[2]);

    // Initialize structs
    CSRdata *a = (CSRdata *)malloc(sizeof(CSRdata));
    CSRdata *b = (CSRdata *)malloc(sizeof(CSRdata));
    COOdata *c = (COOdata *)malloc(sizeof(COOdata));

    char buffer[256];
    int counter = 0;
    int xp = 0;
    int yp = 0;

    // Reads in a .mtx file to the x and y arrays of the `a` CSRData struct
    while (fgets(buffer, sizeof(buffer), fptr) != NULL) {
        if (buffer[0] == '%') {
            continue;
        }
        if (counter == 0) {
            int rowValue, colValue, numnnz;
            sscanf(buffer, "%d %d %d", &rowValue, &colValue, &numnnz);
            a->m = rowValue;
            a->n = colValue;
            a->nnz = numnnz;

            a->rowOffsets = (int *)calloc((rowValue + 1), sizeof(int));
            a->cindexes = (int *)malloc(numnnz * sizeof(int));
            a->values = (double *)malloc(numnnz * sizeof(double));

            a -> x = (int *)malloc(numnnz * sizeof(int));
            a -> y = (int *)malloc(numnnz * sizeof(int));
        }

        else {
           int i, j;
           double val;
            sscanf(buffer, "%d %d %lf", &i, &j, &val);
            
            a -> x[xp] = i - 1;
            a -> y[yp] = j - 1;
            a -> values[yp] = val;
            xp += 1;
            yp += 1;
        }
        counter += 1;
    }
    fclose(fptr);

    a->nnz = xp;
    convert(a);

    // Copies the `a` matrix to the `b` matrix
    b->m = a->m;
    b->n = a->n;
    b->nnz = a->nnz;
    int rowCount = a->m + 1;
    b->rowOffsets = (int *)malloc(rowCount * sizeof(int));
    memcpy(b->rowOffsets, a->rowOffsets, rowCount * sizeof(int));
    b->cindexes = (int *)malloc(a->nnz * sizeof(int));
    memcpy(b->cindexes, a->cindexes, a->nnz * sizeof(int));
    b->values = (double *)malloc(a->nnz * sizeof(double));
    memcpy(b->values, a->values, a->nnz * sizeof(double));

    // Times the matrix multiplication method
    struct timespec start, end;
    double time;
    clock_gettime(CLOCK_MONOTONIC, &start);
    sparseMM_pthread(a, b, c, p);
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    time =
    BILLION *(end.tv_sec - start.tv_sec) +(end.tv_nsec - start.tv_nsec);
    time = time / BILLION;

    // Prints elapsed time
    printf("Elapsed: %lf seconds\n", time);

    // Prints result in COO format for debugging; should be commented out on large matrices.
    for (int i = 0; i < c->nnz; i++) {
        printf("%d\t%d\t%.2f\n", c -> row[i], c -> column[i], c -> values[i]);
    }
    return 0;
}


