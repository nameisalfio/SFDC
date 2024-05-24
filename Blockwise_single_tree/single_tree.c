#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include "string.h"
#include "math.h"
#include "time.h"
#include "sys/stat.h"
#include "stdbool.h"
#include "../Huffman.c"
#include "../SFDC-delta.c"
#include "../SFDC-gamma.c"

#define BLOCK_SIZE 1000 * 1

int get_print_interval(int total_blocks)
{
    if (total_blocks < 100)
        return 10;
    else if (total_blocks < 1000)
        return 100;
    else if (total_blocks < 10000)
        return 1000;
    else if (total_blocks < 100000)
        return 10000;
    else
        return 100000;
}

// verifies the existence of the file in the filesystem
bool file_exists(char *filename)
{
    struct stat buffer;
    return (stat(filename, &buffer) == 0);
}

// write the content of the file in text
void load(char *filename, int n, int *text)
{
    FILE *fp = fopen(filename, "r");
    char c;
    int i = 0;
    while (i < n && (c = getc(fp)))
        text[i++] = c;
    text[i] = '\0';
    fclose(fp);
}

// get size of the file in bytes
int get_file_size(char *filename)
{
    FILE *fp = fopen(filename, "r");
    fseek(fp, 0L, SEEK_END);
    int sz = ftell(fp);
    fclose(fp);
    return sz;
}

// verifies if x and y are equal in the window [i...j]
int check_decoding(int *x, int *y, int i, int j)
{
    for (int k = i; k <= j; k++)
    {
        if (x[k] != y[k])
        {
            printf("Error: text[%d]: %d, decoded[%d]: %d\n", k, x[k], k, y[k]);
            return k;
        }
    }
    return -1;
}

// reads a sequence of integers from a file, stores them in text
int load_int(char *filename, int n, unsigned int *text)
{
    FILE *fp = fopen(filename, "r");
    int i = 0;
    while (i < n && fread(&text[i], sizeof(unsigned int), 1, fp) == 1)
        i++;
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

int main(int argc, char **argv)
{
    printf("\n");
    if (argc != 4)
    {
        printf("Error in input arguments!\n\n");
        printf("Usage: %s <delta|gamma> <file> <lambda>\n", argv[0]);
        return 0;
    }

    char *variant = argv[1];
    char *file = argv[2];
    int nl = atoi(argv[3]);
    nl--;

    if (strcmp(variant, "delta") && strcmp(variant, "gamma"))
    {
        printf("Error: first argument must be 'delta' or 'gamma'!\n\n");
        return 0;
    }
    if (!file_exists(file))
    {
        printf("Error: file does not exist!\n\n");
        return 0;
    }

    // LOAD FILE AND COMPUTE SOME STATS ********************************
    int n = INT_MAX;
    int sz = get_file_size(file);
    printf("File size : %d Bytes\n", sz);
    fflush(stdout); // size of the file in bytes
    printf("Block size : %d Bytes\n\n", BLOCK_SIZE);
    fflush(stdout); // size of the block in bytes
    if (n <= 0 || n > sz)
        n = sz;
    int *text = (int *)malloc(sizeof(int) * (n + 1));
    printf("Loading sequence...\n");
    fflush(stdout);
    int N = n = load_int(file, n, text);
    printf("Length of the sequence : %d\n", n);
    fflush(stdout); // number of integers in the read sequence
    double total_avgdelay = 0.0;
    double total_avgbits = 0.0;
    double maxdelay = 0.0;
    int maxblock = 0;
    int max_char;
    double time, avgdelay, avgbits;
    int total_blocks = (sz + BLOCK_SIZE * sizeof(int) - 1) / (BLOCK_SIZE * sizeof(int)); // total number of blocks
    int print_interval = get_print_interval(total_blocks);

    int maxval = 0;
    double avgval = 0.0;
    for (int i = 0; i < n; i++)
    {
        if (text[i] > maxval)
            maxval = text[i];
        avgval += text[i];
    }
    maxval++;
    printf("Maximum value in the sequence: %d\n", maxval - 1);
    fflush(stdout);
    printf("Average value in the sequence: %.2f\n", avgval / (double)n);
    fflush(stdout);
    int *freq = (int *)malloc(sizeof(int) * maxval);
    int *prev_freq = (int *)malloc(sizeof(int) * maxval);
    for (int i = 0; i < maxval; i++)
    {
        freq[i] = 0;
        prev_freq[i] = 0;
    }
    int nchars = compute_freq(text, n, freq, maxval);
    printf("Size of the alphabet: %d\n", nchars);
    fflush(stdout);
    //*****************************************************************

    // RUN HUFFMAN ALGORITHM  ******************************************
    HNODE **hset = (HNODE **)malloc(sizeof(HNODE *) * maxval);
    uint64_t *map = (uint64_t *)malloc(sizeof(uint64_t) * maxval);
    int *codelen = (int *)malloc(sizeof(int) * maxval);
    build_huffman_codes(text, n, maxval, hset, map, codelen, freq);
    int maxcodelen = 0;
    double avgcodelen = 0.0;
    for (int i = 0; i < n; i++)
    {
        if (codelen[text[i]] > maxcodelen)
            maxcodelen = codelen[text[i]];
        avgcodelen += codelen[text[i]];
    }
    printf("Maximum code length: %d\n", maxcodelen);
    fflush(stdout);
    printf("Average code length: %.2f\n", avgcodelen / (double)n);
    fflush(stdout);
    //*****************************************************************

    printf("\nCODING WITH %d LAYERS ... \n", nl + 1);
    fflush(stdout);
    STACKNODE *stack = (STACKNODE *)malloc(sizeof(STACKNODE) * n);
    HNODE **huffman_trees = (HNODE **)malloc(sizeof(HNODE *) * total_blocks);
    uint32_t ***encoded_blocks = (uint32_t ***)malloc(sizeof(uint32_t **) * total_blocks);
    uint32_t **original_blocks = (uint32_t **)malloc(sizeof(uint32_t *) * total_blocks);
    HNODE *root = NULL;

    // build Huffman tree
    build_huffman_codes(text, n, maxval, hset, map, codelen, freq);
    root = hset[0];
    build_codeset(map, codelen, root, 0, 0);
    int tree_size = get_tree_size(root, 0);
    long unsigned total_compressed_size = 0;
    double total_weighted_avgdelay = 0.0;
    double total_weighted_avgbits = 0.0;

    // Encoding
    for (int block_index = 0; block_index < total_blocks; block_index++)
    {
        // load current block
        int compressed_size = 0;
        int block_length = (block_index == total_blocks - 1 && sz % BLOCK_SIZE != 0) ? sz % BLOCK_SIZE : BLOCK_SIZE;
        n = load_int_block(file, text, block_length, block_index);
        unsigned int *block_data = (int *)malloc(sizeof(int) * block_length);

        // compute maxval, avgval and frequencies
        maxval = 0;
        avgval = 0.0;
        for (int i = 0; i < n; i++)
        {
            if (text[i] > maxval)
                maxval = text[i];
            avgval += text[i];
        }
        maxval++;
        nchars = compute_freq(text, n, freq, maxval);

        // compute maxcodelen and avgcodelen
        maxcodelen = 0;
        avgcodelen = 0.0;
        for (int i = 0; i < n; i++)
        {
            if (codelen[text[i]] > maxcodelen)
                maxcodelen = codelen[text[i]];
            avgcodelen += codelen[text[i]];
        }

        // compute the SFDC encoding
        if (strcmp(variant, "delta") == 0)
            encoded_blocks[block_index] = delta_encode(text, block_length, map, codelen, nl, &avgdelay, &avgbits, &max_char, block_index);
        else
            encoded_blocks[block_index] = gamma_encode(text, block_length, map, codelen, nl, &avgdelay, &avgbits, &max_char, block_index);

        huffman_trees[block_index] = root;
        compressed_size = get_compressed_size(root, map, codelen);

        // print block's stats
        // print block's stats
        if (block_index % print_interval == 0)
        {
            printf("\n_________________________________  [Block %d]  _________________________________\n\n", block_index);
            fflush(stdout);
            printf("Maximum value: %d\n", maxval - 1);
            fflush(stdout);
            printf("Average value: %.3f\n", avgval / (double)n);
            fflush(stdout);
            printf("Size of the alphabet: %d\n", nchars);
            fflush(stdout);
            printf("Number of symbols: %d\n", n);
            fflush(stdout);
            printf("Maximum code length: %d\n", maxcodelen);
            fflush(stdout);
            printf("Average code length: %.2f\n", avgcodelen / (double)n);
            fflush(stdout);
            printf("Compressed size: %d\n", compressed_size);
            fflush(stdout);
            printf("Tree size: %d\n", tree_size);
            fflush(stdout);
            printf("Avg-delay: %f\n", avgdelay);
            fflush(stdout);
            printf("Avg-bit: %f\n", avgbits);
            fflush(stdout);
            // print_huffman_tree(root, map, codelen);
        }

        // store the original block, the encoded block, the huffman tree
        memcpy(block_data, text, sizeof(int) * block_length);
        original_blocks[block_index] = block_data;

        // update stats
        total_weighted_avgdelay += avgdelay * block_length;
        total_weighted_avgbits += avgbits * block_length;   
        total_avgdelay += avgdelay;
        total_avgbits += avgbits;
        total_compressed_size += compressed_size;

        if (avgdelay > maxdelay)
        {
            maxdelay = avgdelay;
            maxblock = block_index;
        }

        // update previous frequencies
        memcpy(prev_freq, freq, sizeof(int) * maxval);
    }

    double total_delay = total_weighted_avgdelay / N;
    double total_bits = total_weighted_avgbits / N;
    double tree_size_percentage = ((double)tree_size / (double)total_compressed_size) * 100;

    printf("_____________________________________________________________________________________\n\n");
    fflush(stdout);
    printf("\033[1;31m");
    printf("Maximum delay over all blocks: %.2f in block %d caused by char %d\n", maxdelay, maxblock, max_char);
    fflush(stdout);
    printf("\033[1;32m");
    printf("Average delay over all blocks: %.2f\n", total_delay);
    fflush(stdout);
    printf("Average number of bits/e for blocks: %.2f\n", total_bits);
    fflush(stdout);
    printf("\033[0m");
    printf("Total compressed size: %lu bytes\n", total_compressed_size);
    fflush(stdout);
    printf("Size of the Huffman tree: %d bytes\n", tree_size);
    fflush(stdout);
    printf("\033[1;32m");
    printf("Huffman trees size as a percentage of the compressed text: %.2f%%\n", tree_size_percentage);
    fflush(stdout);
    printf("\033[0m");
    printf("Number of blocks: %d\n", total_blocks);
    fflush(stdout);

    // Decoding
    double total_decode_time = 0.0;
    double total_access_time = 0.0;
    int access_count = 0;

    printf("\n\nStart Decoding...\n\n");
    fflush(stdout);
    for (int block_index = 0; block_index < total_blocks; block_index++)
    {
        int block_length = (block_index == total_blocks - 1 && sz % BLOCK_SIZE != 0) ? sz % BLOCK_SIZE : BLOCK_SIZE;

        // load the encoded block
        unsigned int *block_decoded_text = (unsigned int *)malloc(sizeof(unsigned int) * block_length);
        memset(block_decoded_text, 0, sizeof(unsigned int) * block_length);
        clock_t decode_start = clock();

        // decode the block
        if (strcmp(variant, "delta") == 0)
            delta_decode(encoded_blocks[block_index], 0, block_length - 1, nl, huffman_trees[block_index], block_length, block_decoded_text, stack);
        else
            gamma_decode(encoded_blocks[block_index], 0, block_length - 1, nl, huffman_trees[block_index], block_length, block_decoded_text, stack);

        clock_t decode_end = clock();

        double decode_time = (decode_end - decode_start) * 1000.0 / (double)CLOCKS_PER_SEC; // millisecondsS
        total_decode_time += decode_time;

        // compute average acces time
        for (int k = 0; k < 10000 && k < block_length; k++)
        {
            int i = rand() % block_length; // random access
            clock_t init = clock();

            if (strcmp(variant, "delta") == 0)
                delta_decode(encoded_blocks[block_index], i, i, nl, huffman_trees[block_index], block_length, block_decoded_text, stack);
            else
                gamma_decode(encoded_blocks[block_index], i, i, nl, huffman_trees[block_index], block_length, block_decoded_text, stack);

            clock_t end = clock();
            double access_time = (end - init) * 1000000.0 / (double)CLOCKS_PER_SEC; // microseconds
            total_access_time += access_time;
            access_count++;
        }

        // check decoding
        int block_size = min(BLOCK_SIZE, n - block_index * BLOCK_SIZE);
        int error_position = check_decoding(original_blocks[block_index], block_decoded_text, 0, block_size - 1);
        if (error_position != -1)
        {
            printf("Decoding error at position %d in block %d\n", error_position, block_index);
            fflush(stdout);
            exit(EXIT_FAILURE);
        }
        else
        {
            if (block_index % print_interval == 0)
            {
                printf("_________________________________   [Block %d]   ________________________________\n", block_index);
                fflush(stdout);
                printf("\nDecoding of block %d successful\n", block_index);
                fflush(stdout);
                printf("Decoding time: %f\n", decode_time);
                fflush(stdout);
            }
        }
    }

    double avg_decode_time = total_decode_time / total_blocks;
    double avg_access_time = total_access_time / access_count;

    printf("____________________________________________________________________________________\n\n");
    fflush(stdout);
    printf("\033[1;31m");
    printf("\nFull decoding time: %.2f millisec.\n", total_decode_time);
    fflush(stdout);
    printf("\033[1;32m");
    printf("Average decoding time over all blocks: %.2f millisec.\n", avg_decode_time);
    fflush(stdout);
    printf("Average access time: %.2f microsec.\n", avg_access_time);
    fflush(stdout);
    printf("\033[0m");
    printf("\n____________________________________________________________________________________\n\n\n");
    fflush(stdout);

    // Free memory
    for (int i = 0; i < total_blocks; i++)
    {
        free(original_blocks[i]);
        free(encoded_blocks[i]);
    }
    free(huffman_trees);
    free(original_blocks);
    free(encoded_blocks);
    free(text);
    free(hset);
    free(map);
    free(codelen);
    free(stack);
    free(prev_freq);
    free(freq);

    return EXIT_SUCCESS;
}