// SFDC code - gamma variant
// pending bits are maintained on all avalilable positions
// tuned for working of file of 4bytes values

#define min(a,b) a<b?a:b
#define max(a,b) a>b?a:b

uint32_t **gamma_encode(int *y, int n, uint64_t *map, int *codelen, int nlayers, double *avgdelay, double *avgbits, int *max_char, int idx) {
    int idlebits = 0; // number of idle bits
    uint32_t **s;
    int STACKSIZE = 100*n;
    STACKCODE* stack = (STACKCODE*) malloc(sizeof(STACKCODE)*n);
    int top = -1;
    int pos;
    *avgdelay = 0.0;
    s = (uint32_t**) malloc (sizeof(uint32_t*)*(nlayers+1));
    int nsample = (n/32)+1;
    for(int i=0; i<=nlayers; i++) {
        s[i] = (uint32_t*) malloc (sizeof(uint32_t)*10*nsample);
        for(int j=0; j<10*nsample; j++) s[i][j] = 0;
    }

    int *delay = (int *)malloc(sizeof(int) * n);
    for (int i = 0; i < n; i++) delay[i] = 0.0;
    int max_delay = 0.0;
    int max_delay_index = 0;

    for(int i=0; i<n; i++) {
        int j = 0;
        while(j<=nlayers && j<codelen[y[i]]) {
            if(getbit64(map[y[i]],64-codelen[y[i]]+j))
                setbit32(&s[j][i/32],i%32);
            j++;
        }
        if(j<codelen[y[i]]) {
            top++;
            stack[top].i = i;
            stack[top].j = j;
        }
        while(j<=nlayers && top>=0) {
            int h = stack[top].j;
            pos = stack[top].i;
            if(getbit64(map[y[pos]],64-codelen[y[pos]]+h))
                setbit32(&s[j][i/32],i%32);
            stack[top].j++;
            if(stack[top].j == codelen[y[pos]]) {
                delay[pos] += i-pos;
                if(delay[pos] > max_delay) {
                    max_delay = delay[pos];
                    max_delay_index = pos;
                }
                (*avgdelay) += i-pos;
                top--;
            }
            j++;
        }
    }
    double totalbits = n*(nlayers+1);
    int i=n;
    int j=0;
    while(top>=0) {
        totalbits++;
        int h = stack[top].j;
        pos = stack[top].i;
        if(getbit64(map[y[pos]],64-codelen[y[pos]]+h))
            setbit32(&s[j][i/32],i%32);
        stack[top].j++;
        if(stack[top].j == codelen[y[pos]]) {
            delay[pos] += min(i,n-1)-pos;
            if(delay[pos] > max_delay) {
                max_delay = delay[pos];
                max_delay_index = pos;
            }
            (*avgdelay) += min(i,n-1)-pos;
            top--;
        }
        j++;
        if(j>nlayers) {
            i++;
            j=0;
        }
    }
    *avgbits = totalbits/(double)n;
    *avgdelay = (*avgdelay)/(double)n;
    *max_char = y[max_delay_index];

/*
    // print stats for block
    for (int i = 0; i < n; i++)
        printf("Text[%d] = %d \tdelay %d\n", i, y[i], delay[i]);
    // printf("\nText[%d] = %d has highest delay of %d on a text of length %d\n", max_delay_index, y[max_delay_index], max_delay, n);
*/

    free(stack);
    return s;
}


int gamma_decode(
                uint32_t **blce, // the BLCE representation of the text
                int i, // the starting position of window we want to decode
                int j, // the starting position of window we want to decode (for a single char i=j)
                int nlayers, // the number of fixed layers in the BLCE
                HNODE *root, // the root of the Huffman tree
                int n, // the length of the text (position of the last char)
                int* dt, //decoded text
                STACKNODE* stack//,
                //HNODE **fmap
                )
{
    int pos = i;
    int top = -1;
    do {
        if(i<n) {
            ++top;
            stack[top].node = root;
            stack[top].pos = i;
        }
        int b=0;
        while(b<=nlayers && top>=0) {
            if(getbit32(blce[b][i/32], (i%32)))
                stack[top].node = stack[top].node->right;
            else
                stack[top].node = stack[top].node->left;
            if(stack[top].node->is_char) {
                dt[stack[top].pos] = stack[top].node->c;
                top--;
                if(top<0 && dt[j]!='\0') {
                    return (i<n?i:(n-1));
                }
            }
            b++;
        }
        i++;
    } while(1);
}


uint32_t **gamma_encode2(int *y, int n, uint64_t *map, int *codelen, int nlayers, double *avgdelay, double *avgbits) {
    int idlebits = 0; // number of idle bits
    uint32_t **s;
    int STACKSIZE = 100*n;
    STACKCODE* stack = (STACKCODE*) malloc(sizeof(STACKCODE)*n);
    int top = -1;
    int pos;
    *avgdelay = 0.0;
    s = (uint32_t**) malloc (sizeof(uint32_t*)*(nlayers+1));
    int nsample = (n/32)+1;
    for(int i=0; i<=nlayers; i++) {
        s[i] = (uint32_t*) malloc (sizeof(uint32_t)*10*nsample);
        for(int j=0; j<10*nsample; j++) s[i][j] = 0;
    }

    for(int i=0; i<n; i++) {
        int j = 0;
        while(j<=nlayers && j<codelen[y[i]]) {
            if(getbit64(map[y[i]],64-codelen[y[i]]+j))
                setbit32(&s[j][i/32],i%32);
            j++;
        }
        if(j<codelen[y[i]]) {
            top++;
            stack[top].i = i;
            stack[top].j = j;
        }
        while(j<=nlayers && top>=0) {
            int h = stack[top].j;
            pos = stack[top].i;
            if(getbit64(map[y[pos]],64-codelen[y[pos]]+h))
                setbit32(&s[j][i/32],i%32);
            stack[top].j++;
            if(stack[top].j == codelen[y[pos]]) {
                (*avgdelay) += i-pos;
                top--;
            }
            j++;
        }
    }
    double totalbits = n*(nlayers+1);
    int i=n;
    int j=0;
    while(top>=0) {
        totalbits++;
        int h = stack[top].j;
        pos = stack[top].i;
        if(getbit64(map[y[pos]],64-codelen[y[pos]]+h))
            setbit32(&s[j][i/32],i%32);
        stack[top].j++;
        if(stack[top].j == codelen[y[pos]]) {
            (*avgdelay) += min(i,n-1)-pos;
            top--;
        }
        j++;
        if(j>nlayers) {
            i++;
            j=0;
        }
    }
    *avgbits = totalbits/(double)n;
    *avgdelay = (*avgdelay)/(double)n;
    free(stack);
    return s;
}
