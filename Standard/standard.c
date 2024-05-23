#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "string.h"
#include "math.h"
#include "time.h"
#include "sys/stat.h"
#include "stdbool.h"
#include "../Huffman.c"
#include "../SFDC-delta.c"
#include "../SFDC-gamma.c"

// verifica l'esistenza del file all'interno del filesystem
bool file_exists (char *filename) { 
    struct stat   buffer;
    return (stat (filename, &buffer) == 0);
}

// scrive il contenuto del file in text
void load(char *filename, int n, int *text) {
    FILE *fp = fopen(filename, "r");
    char c;
    int i = 0;
    while (i<n && (c=getc(fp))) {
        text[i++] = c;
    }
    text[i] = '\0';
    fclose(fp);
}

// restituisce la dimensione in byte del file
int get_file_size(char *filename) {
    FILE *fp = fopen(filename, "r");
    fseek(fp, 0L, SEEK_END);
    int sz = ftell(fp);
    fclose(fp);
    return sz;
}

// legge una sequenza di interi da un file, li memorizza in text e
// restituisce il numero totale di interi letti.
int load_int(char *filename, int n, int *text) {
    FILE *fp = fopen(filename, "r");
    int i = 0;
    while (i<n && fread(&text[i], 4, 1, fp)==1) i++;
    fclose(fp);
    return i;
}

// returns the number of integers read
int load_int_block(char *filename, int *text, int block_size, int block_index)
{
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("Failed to open file: %s\n", filename);
        return 0;
    }
    /*
    printf("\033[1;31m"); 
    printf("Loaded block %d of file %s\n", block_index, filename); fflush(stdout);
    printf("\033[0m"); 
    */
    fseek(fp, block_index * block_size * sizeof(unsigned int), SEEK_SET);
    int i = 0;
    while (i < block_size && !feof(fp))
    {
        int read_count = fread(&text[i], sizeof(int), 1, fp);
        if (read_count < 1)
            break;
        i++;
    }

    fclose(fp);
    return i;
}

// verifica se x e y sono uguali nella finestra [i...j]
int check_decoding(int* x, int*y, int i, int j) {
    for(int k=i; k<=j; k++)
        if(x[k]!=y[k]) return k;
    return -1;
}

int main(int argc, char **argv) {
    
    printf("\n");
    if(argc!=4) {
        printf("Error in input arguments!\n\n");
        printf("Usage: %s <delta|gamma> <file> <lambda>\n", argv[0]);
        return 0;
    }

    char *variant = argv[1];
    char *file = argv[2];
    int nl = atoi(argv[3]);
    nl--;

    if(strcmp(variant,"delta") && strcmp(variant,"gamma")) {
        printf("Error: first argument must be 'delta' or 'gamma'!\n\n");
        return 0;
    }
    if(!file_exists(file)) {
        printf("Error: file does not exist!\n\n");
        return 0;
    }

    // LOAD FILE AND COMPUTE SOME STATS ********************************
    int n = 104857600;
    int sz = get_file_size(file);
    printf("File size : %d Bytes\n", sz);fflush(stdout); // dimensione del file in byte
    if(n<=0 || n>sz) n=sz;
    int *text = (int*) malloc (sizeof(int)*(n+1));
    printf("Loading sequence\n");fflush(stdout);
    n = load_int(file, n, text);
    printf("Length of the sequence : %d\n",n);fflush(stdout); // numero di interi nella sequenza letta
    int maxval = 0; double avgval = 0.0;
    for(int i=0; i<n; i++) {
        if(text[i]>maxval) maxval = text[i];
        avgval += text[i];
    }
    maxval++;
    printf("Maximum value in the sequence: %d\n",maxval-1);fflush(stdout);
    printf("Average vale in the sequence: %.2f\n",avgval/(double)n);fflush(stdout);
    printf("Compute frquencies\n");fflush(stdout);
    int *freq = (int*) malloc (sizeof(int)*maxval);;
    int nchars = compute_freq(text, n, freq, maxval);
    printf("Size of the alphabet: %d\n",nchars);fflush(stdout);
    // *****************************************************************
    
    
    // RUN HUFFMAN ALGORITHM  ******************************************
    HNODE **hset = (HNODE**) malloc (sizeof(HNODE*)*maxval);
    uint64_t *map = (uint64_t*) malloc (sizeof(uint64_t)*maxval);
    int *codelen = (int*) malloc (sizeof(int)*maxval);
    build_huffman_codes(text, n, maxval, hset, map, codelen, freq);
    int maxcodelen = 0; double avgcodelen = 0.0;
    for(int i=0; i<n; i++) {
        if(codelen[text[i]]>maxcodelen) maxcodelen = codelen[text[i]];
        avgcodelen += codelen[text[i]];
    }
    printf("Maximum code length: %d\n",maxcodelen);fflush(stdout);
    printf("Average code length: %.2f\n",avgcodelen/(double)n);fflush(stdout);
    HNODE *root = hset[0];
    for(int i=0; i<maxval; i++) codelen[i] = map[i] = 0;
    build_codeset(map, codelen, root, 0, 0);
    //*****************************************************************

    STACKNODE* stack = (STACKNODE*) malloc(sizeof(STACKNODE)*n);
    double time, avgdelay, avgbits;
    clock_t init, end;

    int max_char, idx=0;
    printf("\nCONDING WITH %d LAYERS ... \n", nl+1);fflush(stdout);
    uint32_t **sfdc;
    if(strcmp(variant,"delta")==0)
        sfdc = delta_encode(text, n, map, codelen, nl, &avgdelay, &avgbits, &max_char, idx);
    else 
        sfdc = gamma_encode(text, n, map, codelen, nl, &avgdelay, &avgbits, &max_char, idx);

    printf("\033[1;32m");
    printf("Average delay: %.2f\n",avgdelay);fflush(stdout);
    printf("\033[0m");
    printf("Average number of bits/e: %.2f\n",avgbits);fflush(stdout);

    printf("Making fast map\n");fflush(stdout);
    HNODE **fmap = (HNODE**) malloc (sizeof(HNODE*)*maxval);
    for(int i=0; i<maxval; i++) fmap[i] = NULL;
    make_fast_map(root, nl, fmap, 0, 0);

    printf("Allocating blind text\n");fflush(stdout);
    int *dt = (int*) malloc (sizeof(int)*(n+1));
    for(int i=0; i<=n; i++) dt[i] = '\0';
    int *shtab = (int*) malloc (sizeof(int)*maxval);
    for(int i=0; i<maxval; i++) shtab[i]=0;
    int c;
    int totsh=0;
    int maxsh = 0;

    double compressed = get_compressed_size(root, map, codelen);
    printf("Compressed text : %.2f bytes\n", compressed);fflush(stdout);
    printf("Compression ratio : %.2f\n", compressed/(double)n);fflush(stdout);

    /*
    // Print Huffman tree
    print_huffman_tree(root, map, codelen);
    */
    
    printf("Start decoding\n");fflush(stdout);
    init = clock();
    if(strcmp(variant,"delta")==0)
        delta_decode(sfdc, 0, n-1, nl, root, n, dt, stack);//, fmap);
    else 
        gamma_decode(sfdc, 0, n-1, nl, root, n, dt, stack);
    end = clock();
    time = (end-init)*1000.0 / (double)CLOCKS_PER_SEC;
    printf("Full decoding time : %.2f\n", time);fflush(stdout);
    int k = check_decoding(text,dt,0,n-1);
    if(k>-1)
        printf("ERROR IN DECODING (%d)",k);

    if(avgdelay<10000) {
        time = 0.0;
        avgdelay = 0.0;
        int j;
        for(int k=0; k<10000; k++) {
            int i = rand()%n;
            init = clock();
            if(strcmp(variant,"delta")==0)
                j = delta_decode(sfdc, i, i, nl, root, n, dt, stack);//, fmap);
            else 
                j = gamma_decode(sfdc, i, i, nl, root, n, dt, stack);
            end = clock();
            avgdelay += j-i;
            time += (end-init)*1000000.0 / (double)CLOCKS_PER_SEC; //microsecondi
            if(text[i]!=dt[i]) printf("ERROR IN DECODING (%d)",k);
            for(int h=i; h<=j; h++) dt[h]='\0';
        }
        time /= 10000.0;
        avgdelay /= 10000.0;
        printf("Average access time : %.2f microsec.\n", time);
    }
    else
        printf("Average access time : --- microsec.\n");

    for(int i=0; i<=nl; i++)  free(sfdc[i]);
    free(sfdc);
    free(fmap);
    free(shtab);
    printf("\n");
    free(dt);
    free(stack);
    return 1;
}
