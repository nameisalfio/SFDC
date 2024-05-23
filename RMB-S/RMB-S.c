#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <stdint.h>
#include "string.h"
#include "math.h"
#include "time.h"
#include "sys/stat.h"
#include "stdbool.h"
#include "../Huffman.c"
#include "../SFDC-delta.c"
#include "../SFDC-gamma.c"
#include "list.c"

#define RANK 2 // Number of rare markers
#define K 20   // Common chars sequence length

// int.dblp.txt
// ------------------------------------------------
// #define THRESHOLD 0.045 // 2
// #define THRESHOLD 0.053 // 4
// #define THRESHOLD 0.073 // 6
// #define THRESHOLD 0.050 // 8
// #define THRESHOLD 0.043 // 10

// int.english.txt
// ------------------------------------------------
// #define THRESHOLD 0.162 // 2
// #define THRESHOLD 0.125 // 4
// #define THRESHOLD 0.113 // 6
// #define THRESHOLD 0.118 // 8
// #define THRESHOLD 0.126 // 10

// int.protein.txt
// ------------------------------------------------
#define THRESHOLD 0.058 // 2
// #define THRESHOLD 0.124 // 4
// #define THRESHOLD 0.037 // 6
// #define THRESHOLD 0.037 // 8
//#define THRESHOLD 0.125 // 10

typedef struct
{
    int character;
    int frequency;
} CharFreq;

int compare(const void *a, const void *b)
{
    CharFreq *cf1 = (CharFreq *)a;
    CharFreq *cf2 = (CharFreq *)b;
    return cf1->frequency - cf2->frequency;
}

CharFreq *find_rare_markers(int n, int *freq, int maxval)
{
    CharFreq *min_freqs = malloc(sizeof(CharFreq) * maxval);
    int count = 0;

    for (int i = 0; i < maxval; i++)
    {
        if (freq[i] != 0)
        {
            min_freqs[count].character = i;
            min_freqs[count].frequency = freq[i];
            count++;
        }
    }

    qsort(min_freqs, count, sizeof(CharFreq), compare);

    CharFreq *result = malloc(sizeof(CharFreq) * RANK);
    memcpy(result, min_freqs, sizeof(CharFreq) * RANK);

    free(min_freqs);
    return result;
}

bool is_rare_marker(CharFreq *infrequents, int char_to_check)
{
    for (int i = 0; i < RANK; i++)
        if (infrequents[i].character == char_to_check)
            return true;
    return false;
}

void print_block(List *original_blocks, int block_index, int block_length)
{
    unsigned int *block_data = (unsigned int *)get_by_idx(original_blocks, block_index);
    printf("[Block %d] length: %d --> [%d", block_index, block_length, block_data[0]);
    fflush(stdout);
    for (int i = 1; i < block_length; i++)
        printf(" - %d", block_data[i]);
    fflush(stdout);
    printf("]\n\n");
}

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

int check_splitting_consistency(List *original_blocks, List *block_sizes, int total_len, int *Text)
{
    int *concatenated_blocks = malloc(total_len * sizeof(int));
    if (concatenated_blocks == NULL)
    {
        printf("Failed to allocate memory for concatenated blocks\n");
        return 0;
    }

    // Concatenate the blocks
    int offset = 0;
    for (int i = 0; i < list_length(block_sizes); i++)
    {
        int block_len = (int)(intptr_t)get_by_idx(block_sizes, i);
        int *block_data = (int *)get_by_idx(original_blocks, i);
        memcpy(concatenated_blocks + offset, block_data, block_len * sizeof(int));
        offset += block_len;
    }

    // Compare the concatenated blocks with the original text
    int is_consistent = memcmp(Text, concatenated_blocks, total_len * sizeof(int)) == 0;

    free(concatenated_blocks);
    return is_consistent;
}

bool check_block_consistency(int block_length, int block_index, List *block_sizes, CharFreq *infrequents, int *text, int next_block_start)
{
    if (block_length == 0)
    {
        printf("Error: Block %d has length 0\n", block_index);
        return false;
    }
    if (block_index < list_length(block_sizes) - 1 && !is_rare_marker(infrequents, text[next_block_start - 1]))
    {
        printf("Block %d has a different character at the end: %d\n", block_index, text[next_block_start - 1]);
        return false;
    }
    return true;
}

// cosine distance between two frequency vectors
double cosine_distance(int *freq1, int *freq2, int size)
{
    double dot_product = 0.0;
    double norm1 = 0.0;
    double norm2 = 0.0;

    for (int i = 0; i < size; i++)
    {
        dot_product += (double)freq1[i] * freq2[i];
        norm1 += (double)freq1[i] * freq1[i];
        norm2 += (double)freq2[i] * freq2[i];
    }

    if (norm1 == 0.0 || norm2 == 0.0)
        return 1.0;

    double cosine_similarity = dot_product / (sqrt(norm1) * sqrt(norm2));

    double cosine_distance = 1.0 - cosine_similarity;
    return cosine_distance;
}

bool has_significant_change(int *freq1, int *freq2, int size, double *dist_cum)
{
    double cos_dist = cosine_distance(freq1, freq2, size);
    *dist_cum += cos_dist;
    return cos_dist > THRESHOLD;
}

bool char_in_tree(int c, HNODE *node)
{
    if (node == NULL)
        return false;
    if (node->is_char)
        return node->c == c;

    return char_in_tree(c, node->left) || char_in_tree(c, node->right);
}

bool all_chars_in_prev_tree(int *block, HNODE *prev_tree, int size)
{
    for (int i = 0; i < size; i++)
        if (!char_in_tree(block[i], prev_tree))
            return false;
    return true;
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
            printf("\nError: text[%d]: %d, decoded[%d]: %d\n", k, x[k], k, y[k]);
            return k;
        }
    }
    return -1;
}

// reads a sequence of integers from a file and stores them in text
int load_int(char *filename, int n, unsigned int *text)
{
    FILE *fp = fopen(filename, "r");
    int i = 0;
    while (i < n && fread(&text[i], sizeof(unsigned int), 1, fp) == 1)
        i++;
    fclose(fp);
    return i;
}

// reads an indexed block of size block_lenght
int load_block(char *filename, unsigned int *text, int n, int block_start)
{
    FILE *fp = fopen(filename, "r");

    // Seek to the offset
    if (fseek(fp, block_start * sizeof(int), SEEK_SET) != 0)
    {
        printf("Failed to seek to position %d in file: %s\n", block_start, filename);
        fclose(fp);
        return 0;
    }

    int i = 0;
    while (i < n && fread(&text[i], sizeof(unsigned int), 1, fp) == 1)
        i++;
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
    printf("File size : %d Bytes\n\n", sz);
    fflush(stdout); // size of the file in bytes
    if (n <= 0 || n > sz)
        n = sz;
    int *Text = (int *)malloc(sizeof(int) * (n + 1));
    int *text = (int *)malloc(sizeof(int) * (n + 1));
    printf("LOADING SEQUENCE ...\n\n");
    fflush(stdout);
    int N = n = load_int(file, n, text);
    load_int(file, n, Text);
    printf("Length of the sequence : %d\n", n);
    fflush(stdout); // number of integers in the read sequence
    double total_avgdelay = 0.0;
    double total_avgbits = 0.0;
    double maxdelay = 0.0;
    double dist_cum = 0.0;
    double avgval = 0.0;
    double time, avgdelay, avgbits;
    int maxblock = 0;
    int max_char;
    int maxval = 0;
    int ntrees = 0;

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

    // Find the first k characters with the least frequency
    CharFreq *infrequents = find_rare_markers(n, freq, maxval);
    printf("\nThe %s with the least frequency %s:\n\n", RANK > 1 ? "characters" : "character",
           RANK > 1 ? "are" : "is");
    fflush(stdout);
    for (int i = 0; i < RANK; i++)
        printf("Char: %d \tfreq: %d\n", infrequents[i].character, infrequents[i].frequency);
    printf("\n");
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

    List *block_sizes = new_list();
    List *encoded_blocks = new_list();
    List *original_blocks = new_list();
    STACKNODE *stack = (STACKNODE *)malloc(sizeof(STACKNODE) * n);
    HNODE *root = NULL;
    unsigned int *block_data = NULL;

    int block_start = 0;
    int block_length = 0;
    int block_index = 0;
    int next_block_start = 0;
    int total_blocks = 0;
    int min_block_size = 0;
    int max_block_size = 0;

    // RARE MARKER BLOCK-SEGMENTATION ******************************************
    printf("\nRARE MARKER BLOCK-SEGMENTATION (RANK %d)...\n", RANK); // migiorare la suddivisione con una soglia
    fflush(stdout);
    while (block_start < N)
    {
        next_block_start = block_start;
        while (next_block_start < N)
        {
            if (is_rare_marker(infrequents, text[next_block_start]))
            {
                do
                {
                    next_block_start++;
                } while (next_block_start < N && is_rare_marker(infrequents, text[next_block_start]));

                // include rare markers within K characters distance
                int lookahead = next_block_start;
                while (lookahead < N && lookahead - next_block_start <= K)
                {
                    if (is_rare_marker(infrequents, text[lookahead]))
                    {
                        next_block_start = lookahead;
                        do
                        {
                            next_block_start++;
                        } while (next_block_start < N && is_rare_marker(infrequents, text[next_block_start]));
                    }
                    lookahead++;
                }

                break;
            }
            next_block_start++;
        }
        block_length = next_block_start - block_start;

        list_append(block_sizes, (void *)(intptr_t)block_length);

        // update min_block_size and max_block_size
        if (block_length < min_block_size || min_block_size == 0)
            min_block_size = block_length;

        if (block_length > max_block_size)
            max_block_size = block_length;

        // load the block data and check its consistency
        n = load_block(file, text, block_length, block_start);
        if (n != block_length)
        {
            printf("Error: block %d -> read %d instead of %d bytes\n", block_index, n, block_length);
            exit(EXIT_FAILURE);
        }

        block_data = (int *)malloc(sizeof(int) * block_length);
        memcpy(block_data, text, sizeof(unsigned int) * block_length);
        list_append(original_blocks, block_data);

        // check consistency of the block
        if (!check_block_consistency(block_length, block_index, block_sizes, infrequents, text, next_block_start))
            exit(EXIT_FAILURE);

        block_index++;
        block_start = next_block_start;
    }
    total_blocks = list_length(block_sizes);
    int print_interval = get_print_interval(total_blocks);
    //*****************************************************************

    // check consistency of text splitting
    if (!check_splitting_consistency(original_blocks, block_sizes, total_block_size(block_sizes), Text))
        exit(EXIT_FAILURE);

    HNODE **huffman_trees = (HNODE **)malloc(sizeof(HNODE *) * total_blocks);
    double total_weighted_avgdelay = 0.0;
    double total_weighted_avgbits = 0.0;
    int total_compressed_size = 0;
    int total_tree_size = 0;
    int tree_size = 0;

    // ENCODING PHASE *************************************************
    printf("\n\nCODING WITH %d LAYERS ... \n", nl + 1);
    fflush(stdout);
    for (block_index = 0; block_index < total_blocks; block_index++)
    {
        // load current block
        int compressed_size = 0;
        n = block_length = (int)(intptr_t)get_by_idx(block_sizes, block_index);
        text = get_by_idx(original_blocks, block_index);

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

        // build new huffman tree if there is a significant change in the frequencies
        if (block_index == 0 || has_significant_change(prev_freq, freq, maxval, &dist_cum) || !all_chars_in_prev_tree(text, root, n))
        {
            build_huffman_codes(text, n, maxval, hset, map, codelen, freq);
            root = hset[0];
            build_codeset(map, codelen, root, 0, 0);

            tree_size = get_tree_size(root);
            total_tree_size += tree_size;
            ntrees++;
        }
        else
            root = huffman_trees[block_index - 1];

        huffman_trees[block_index] = root;
        compressed_size = get_compressed_size(root, map, codelen);

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
            list_append(encoded_blocks, delta_encode(text, block_length, map, codelen, nl, &avgdelay, &avgbits, &max_char, block_index));
        else
            list_append(encoded_blocks, gamma_encode(text, block_length, map, codelen, nl, &avgdelay, &avgbits, &max_char, block_index));

        // update total_avgdelay and total_avgbits
        total_weighted_avgdelay += avgdelay * block_length;
        total_weighted_avgbits += avgbits * block_length;
        total_compressed_size += compressed_size;
        total_avgdelay += avgdelay;
        total_avgbits += avgbits;

        if (avgdelay > maxdelay)
        {
            maxdelay = avgdelay;
            maxblock = block_index;
        }

        // update previous frequencies
        memset(prev_freq, 0, sizeof(int) * maxval);
        memcpy(prev_freq, freq, sizeof(int) * maxval);

        // print some block stats
        if (block_index % print_interval == 0)
        {
            printf("\n_________________________________  [Block %d]  _________________________________\n\n", block_index);
            fflush(stdout);
            printf("Block length: %d\n", block_length);
            fflush(stdout);
            printf("Residual bytes yet to encoded: %d\n", N - cum_block_size(block_sizes, block_index));
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
            printf("Cosine distance from previous block: %f\n", dist_cum);
            fflush(stdout);
            printf("Avg-delay: %f\n", avgdelay);
            fflush(stdout);
            printf("Avg-bit: %f\n", avgbits);
            fflush(stdout);
            // print_block(original_blocks, block_index, block_length);
            if (block_index > 0)
                printf(are_equal(huffman_trees[block_index], huffman_trees[block_index - 1]) ? "Huffman tree was reused\n"
                                                                                             : "Huffman tree was recalculated\n");
            // print_huffman_tree(huffman_trees[block_index], map, codelen);
        }
    }
    //*****************************************************************

    double total_delay = total_weighted_avgdelay / N;
    double total_bits = total_weighted_avgbits / N;
    double avg_distance = dist_cum / total_blocks;
    double tree_size_percentage = ((double)total_tree_size / (double)total_compressed_size) * 100;

    printf("__________________________________________________________________________________\n\n");
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
    printf("Total compressed size: %d bytes\n", total_compressed_size);
    fflush(stdout);
    printf("Compression ratio : %.2f\n", total_compressed_size/(double)N);fflush(stdout);
    printf("Size of the Huffman tree: %d bytes\n", total_tree_size);
    fflush(stdout);
    printf("\033[1;32m");
    printf("Huffman trees size as a percentage of the compressed text: %.2f%%\n", tree_size_percentage);
    fflush(stdout);
    printf("\033[0m");
    printf("Average cosine distance: %.3f\n", avg_distance);
    fflush(stdout);
    printf("Minimum block size: %d\n", min_block_size);
    fflush(stdout);
    printf("Maximum block size: %d\n", max_block_size);
    fflush(stdout);
    printf("Average block size: %.2f\n", average_block_size(block_sizes));
    fflush(stdout);
    printf("Number of blocks: %d\n", total_blocks);
    printf("Number of trees created: %d\n", ntrees);
    fflush(stdout);

    int access_count = 0;
    block_index = 0;
    double total_decode_time = 0.0;
    double total_access_time = 0.0;

    // DECODING PHASE *************************************************
    printf("\n\nDECODING ...\n\n");
    fflush(stdout);
    for (block_index = 0; block_index < total_blocks; block_index++)
    {
        unsigned int *block_data = (unsigned int *)get_by_idx(original_blocks, block_index);
        int block_length = (int)(intptr_t)get_by_idx(block_sizes, block_index);

        unsigned int *block_decoded_text = (unsigned int *)malloc(sizeof(unsigned int) * block_length);
        memset(block_decoded_text, 0, sizeof(unsigned int) * block_length);

        // decode the block
        clock_t decode_start = clock();

        if (strcmp(variant, "delta") == 0)
            delta_decode(get_by_idx(encoded_blocks, block_index), 0, block_length - 1, nl, huffman_trees[block_index], block_length, block_decoded_text, stack);
        else
            gamma_decode(get_by_idx(encoded_blocks, block_index), 0, block_length - 1, nl, huffman_trees[block_index], block_length, block_decoded_text, stack);

        clock_t decode_end = clock();

        double decode_time = (decode_end - decode_start) * 1000.0 / (double)CLOCKS_PER_SEC; // milliseconds
        total_decode_time += decode_time;

        // check decoding
        int error_position = check_decoding(block_data, block_decoded_text, 0, block_length - 1);
        if (error_position == -1)
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
        else
        {
            printf("Decoding error at position %d in block %d\n", error_position, block_index);
            fflush(stdout);
            printf("Size of block %d: %d\n", block_index, block_length);
            fflush(stdout);
            exit(EXIT_FAILURE);
        }

        // compute average acces time
        for (int k = 0; k < 10000 && k < block_length; k++)
        {
            int i = rand() % block_length; // random access
            clock_t init = clock();

            if (strcmp(variant, "delta") == 0)
                delta_decode(get_by_idx(encoded_blocks, block_index), i, i, nl, huffman_trees[block_index], block_length, block_decoded_text, stack);
            else
                gamma_decode(get_by_idx(encoded_blocks, block_index), i, i, nl, huffman_trees[block_index], block_length, block_decoded_text, stack);

            clock_t end = clock();
            double access_time = (end - init) * 1000000.0 / (double)CLOCKS_PER_SEC;

            if (block_data[i] != block_decoded_text[i])
            {
                printf("Error: block %d --> text[%d]: %d, decoded[%d]: %d\n", block_index, i, block_data[i], i, block_decoded_text[i]);
                exit(EXIT_FAILURE);
            }
            total_access_time += access_time;
            access_count++;
        }
    }
    //*****************************************************************

    double avg_decode_time = total_decode_time / total_blocks;
    double avg_access_time = total_access_time / access_count;

    printf("__________________________________________________________________________________\n\n");
    fflush(stdout);
    printf("\033[1;31m");
    printf("Full decoding time: %.2f millisec.\n", total_decode_time);
    fflush(stdout);
    printf("\033[1;32m");
    printf("Average decoding time over all blocks: %.2f millisec.\n", avg_decode_time);
    fflush(stdout);
    printf("Average access time: %.2f microsec.\n", avg_access_time);
    fflush(stdout);
    printf("\033[0m");
    printf("\n____________________________________________________________________________________\n\n\n");
    fflush(stdout);

    // free memory
    for (int i = 0; i < total_blocks; i++)
        huffman_trees[i];
    free(huffman_trees);
    free(text);
    free(hset);
    free(map);
    free(codelen);
    free(stack);
    free(prev_freq);
    free(freq);
    free(infrequents);

    // destroy lists
    list_destroy(block_sizes);
    list_destroy(encoded_blocks);
    list_destroy(original_blocks);

    return EXIT_SUCCESS;
}