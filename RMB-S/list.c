#include <stdlib.h>

typedef struct Node {
    void *data;
    struct Node *next;
} Node;

typedef struct List {
    Node *head;
    int length;
} List;

List *new_list() {
    List *list = malloc(sizeof(List));
    list->head = NULL;
    list->length = 0;
    return list;
}

void list_append(List *list, void *data) {
    Node *new_node = malloc(sizeof(Node));
    new_node->data = data;
    new_node->next = NULL;

    if (list->head == NULL) {
        list->head = new_node;
    } else {
        Node *current = list->head;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = new_node;
    }

    list->length++;
}

int list_length(List *list) {
    return list->length;
}

void *get_last_element(List *list) {
    if (list == NULL || list->head == NULL) {
        printf("Error: List is empty\n");
        exit(EXIT_FAILURE);
    }

    Node *current = list->head;
    while (current->next != NULL) {
        current = current->next;
    }
    return current->data;
}

void *get_by_idx(List *list, int idx) {
    if (list == NULL || idx < 0 || idx >= list_length(list)) {
        printf("Error: Invalid index\n");
        exit(EXIT_FAILURE);
    }

    Node *current = list->head;
    for (int i = 0; i < idx && current != NULL; i++) {
        current = current->next;
    }
    return current != NULL ? current->data : NULL;
}

int cum_block_size(List *block_sizes, int block_index) {
    if (block_sizes == NULL || block_sizes->head == NULL) {
        printf("Error: List is empty\n");
        exit(EXIT_FAILURE);
    }

    Node *current = block_sizes->head;
    int sum = 0;
    int count = 0;
    while (current != NULL && count <= block_index) {
        sum += (int)(intptr_t)current->data;
        current = current->next;
        count++;
    }
    return sum;
}

int total_block_size(List *block_sizes) {
    return cum_block_size(block_sizes, list_length(block_sizes) - 1);
}

double average_block_size(List *block_sizes) {
    return total_block_size(block_sizes) / list_length(block_sizes);
}

void list_destroy(List *list) {
    Node *current = list->head;
    Node *next;

    while (current != NULL) {
        next = current->next;
        free(current);
        current = next;
    }

    free(list);
}