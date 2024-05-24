#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stddef.h> 
#include "string.h"
#include "math.h"
#include "time.h"

uint32_t I32 = 1;
uint64_t I64 = 1;

// estraggono il bit h da code
#define getbit32(code, h)  (((code) >> (31-(h))) & I32)
#define getbit64(code, h)  (((code) >> (63-(h))) & I64)

typedef struct node {
    int c;
    int freq;
    struct node *left, *right;
    int is_char;
    int pos;
} HNODE;

typedef struct stack_node {
    HNODE *node;
    int pos;
} STACKNODE;

typedef struct stackcode {
    int i;
    int j;
} STACKCODE;

void print_queue(STACKCODE* stack, int top, int* y) {
    printf("HEAD-");
    for(int k=0; k<=top; k++) {
        printf("[%d]-", y[stack[k].i]);
    }
    printf("NIL\n");fflush(stdout);
}

// calcola la frequenza di ogni carattere all'interno del testo
// ritorna il numero di caratteri distinti
int compute_freq(int* text, int n, int *freq, int maxval) {
    for(int i=0; i<maxval; i++) freq[i]=0;
    int nchars = 0;
    for(int i=0; i<n; i++) {
        // Verifica che text[i] sia un indice valido per freq
        if (text[i] < 0 || text[i] >= maxval) {
            printf("Invalid value in text: %d\n", text[i]);
            exit(1);
        }
        if(freq[text[i]]==0) nchars++;
        freq[text[i]]++;
    }
    return nchars;
}

// crea un insieme di nodi per la codifica di Huffman a partire dalle frequenze degli elementi
void make_charset(int *freq, HNODE **hset, int nchars, int maxval) {
    for(int i=0; i<maxval; i++) hset[i] = NULL;
    int pos = 0;
    for(int i=0; i<maxval; i++) {
        if(freq[i]>0) {
            HNODE *node = (HNODE *) malloc(sizeof(HNODE));
            node->c = i;
            node->freq = freq[i];
            node->left = node->right = NULL;
            node->is_char = 1;
            hset[pos++] = node;
        }
    }
    // se il numero di caratteri non coincide con il numero di nodi creati
    if(nchars!=pos) printf("\nERROR CHARSET (%d-%d)\n\n",pos, nchars);
}

// stampa il set di caratteri
void print_charset(HNODE **hset, int nchars) {
    for(int i=0; i<nchars; i++) {
        printf("[%c] : %d\n", hset[i]->c, hset[i]->freq);
    }
}

// usa insertion-sort per ordinare i nodi
void order_chars(HNODE **hset, int nchars) {
    for(int i=1; i<nchars; i++) {
        int j=i;
        while(j>0 && hset[j]->freq > hset[j-1]->freq) {
            HNODE *tmp = hset[j];
            hset[j] = hset[j-1];
            hset[j-1] = tmp;
            j--;
        }
    }
}

// esplora l'albero costruendo i codici binari per i caratteri 
void build_codeset(uint64_t *map, int *codelen, HNODE *node, uint64_t code, int depth) {
    if(node == NULL) {
        printf("Error: node is NULL\n");
        //exit(EXIT_FAILURE);
    }
    if(node->is_char) {
        map[node->c] = code; // codifica del carattere 
        codelen[node->c] = depth; // lunghezza della codifica
    }
    else {
        code = code << 1;
        build_codeset(map, codelen, node->left, code, depth+1);
        code = code | I64;
        build_codeset(map, codelen, node->right, code, depth+1);
    }
}

//pone a 1 l'iesimo bit più significativo di map (parte da 0)
void setbit32(uint32_t *map, int i) {
    (*map) |= (1<<(31-i));
}

//restituisce il valore del bit all'indice i in map
int getbit(int map, int i) {
    map >>= i;
    return (map & 1);
}


//stampa il codice binario rappresentato da code con la lunghezza specificata
void print_code(uint64_t code, int codelen) {
    for(int i=64-codelen; i<=63; i++) {
        int bit = getbit64(code, i);
        printf("%d",bit);
    }
}

//costruisce una mappa 'veloce' per l'accesso ai nodi dell'albero di Huffman
void make_fast_map(HNODE* node, int nl, HNODE **fmap, char bitstream, int depth) {
    if(node->is_char || depth==nl) {
        fmap[bitstream] = node;
        return;
    }
    make_fast_map(node->left, nl, fmap, bitstream, depth+1);
    bitstream = bitstream | (1<<depth);
    make_fast_map(node->right, nl, fmap, bitstream, depth+1);
}

// costruisce i codici di Huffman per i caratteri in text
int build_huffman_codes(int *text,
                        int n,
                        int maxval,
                        HNODE **hset,
                        uint64_t *map,
                        int *codelen,
                        int *freq
                        ) {

    int nchars = compute_freq(text, n, freq, maxval);
    make_charset(freq, hset, nchars, maxval);
    order_chars(hset, nchars);
    
    int nnodes = nchars;
    while(nnodes>1) {
        HNODE *a = hset[nnodes-1];
        HNODE *b = hset[nnodes-2];
        HNODE *node = (HNODE *) malloc(sizeof(HNODE));
        node->c = 'X';
        node->freq = a->freq + b->freq;
        node->left = a;
        node->right = b;
        node->is_char = 0;
        hset[nnodes-2] = node;
        hset[nnodes-1] = NULL;
        nnodes--;
        int j = nnodes-1;
        while(j>0 && hset[j]->freq > hset[j-1]->freq) {
            HNODE *tmp = hset[j];
            hset[j] = hset[j-1];
            hset[j-1] = tmp;
            j--;
        }
    }

    HNODE *root = hset[0];
    for(int i=0; i<maxval; i++) codelen[i] = map[i] = 0;
    build_codeset(map, codelen, root, 0, 0);
    return 1;
}

void collect_nodes(HNODE *root, HNODE **nodes, int *index) {
    if (root == NULL) return;

    if (root->is_char) {
        nodes[*index] = root;
        (*index)++;
    }

    collect_nodes(root->left, nodes, index);
    collect_nodes(root->right, nodes, index);
}

int compare_nodes(const void *a, const void *b) {
    HNODE *nodeA = *(HNODE**)a;
    HNODE *nodeB = *(HNODE**)b;
    return nodeA->c - nodeB->c;
}

void print_huffman_tree(HNODE *root, uint64_t *map, int* codelen) {
    HNODE *nodes[256];
    int index = 0;

    collect_nodes(root, nodes, &index);

    qsort(nodes, index, sizeof(HNODE*), compare_nodes);

    printf("\nHuffman Tree:\n\n");fflush(stdout);
    for (int i = 0; i < index; i++) {
        printf("Char: %-10d Freq: %-10d Code-length: %-10d Code: ", nodes[i]->c, nodes[i]->freq, codelen[nodes[i]->c]);
        for(int j=64-codelen[nodes[i]->c]; j<=63; j++) {
            int bit = getbit64(map[nodes[i]->c], j);
            printf("%d",bit);
        }
        printf("\n");
    }
    printf("\n");fflush(stdout);
}

int get_compressed_size(HNODE *root, uint64_t *map, int* codelen) {
    HNODE *nodes[256];
    int index = 0;
    int total_bits = 0;

    collect_nodes(root, nodes, &index);

    for (int i = 0; i < index; i++) 
        total_bits += codelen[nodes[i]->c] * nodes[i]->freq;

    return (total_bits + 7) / 8; // round up to nearest byte
}

size_t get_tree_size(HNODE *root, size_t depth) {
    if (root == NULL) return 0;
    if (root->left == NULL && root->right == NULL) return depth; // leaf
    return depth + get_tree_size(root->left, depth + 1) + get_tree_size(root->right, depth + 1); 
}

bool are_equal(HNODE *root1, HNODE *root2) {

    if (root1 == NULL && root2 == NULL) return 1;

    if (root1 == NULL || root2 == NULL) return 0;

    return (root1->c == root2->c &&
            root1->freq == root2->freq &&
            root1->is_char == root2->is_char &&
            are_equal(root1->left, root2->left) &&
            are_equal(root1->right, root2->right));
}