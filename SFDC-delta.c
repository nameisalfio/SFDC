// SFDC code - delta variant
// all pending bits are maintained on the dynamic layer
// tuned for working of file of 4bytes values

uint32_t **delta_encode(int *y, int n, uint64_t *map, int *codelen, int nlayers, double *avgdelay, double *avgbits, int *max_char, int idx)
{
    int idlebits = 0; // idle bits
    uint32_t **s;
    STACKCODE *stack = (STACKCODE *)malloc(sizeof(STACKCODE) * n);
    int top = -1;
    int pos;
    *avgdelay = 0.0;
    s = (uint32_t **)malloc(sizeof(uint32_t *) * (nlayers + 1));
    int nsample = (n / 32) + 1;
    for (int i = 0; i < nlayers; i++)
    {
        s[i] = (uint32_t *)malloc(sizeof(uint32_t) * nsample);
        for (int j = 0; j < nsample; j++)
            s[i][j] = 0;
    }
    s[nlayers] = (uint32_t *)malloc(sizeof(uint32_t) * (20 * nsample));
    for (int j = 0; j < 20 * nsample; j++)
        s[nlayers][j] = 0;

    int *delay = (int *)malloc(sizeof(int) * n);
    for (int i = 0; i < n; i++) delay[i] = 0.0;
    int max_delay = 0.0;
    int max_delay_index = 0;

    for (int i = 0; i < n; i++)
    {
        int j = 0;
        while (j < nlayers && j < codelen[y[i]])
        {
            if (getbit64(map[y[i]], 64 - codelen[y[i]] + j))
                setbit32(&s[j][i / 32], i % 32);
            j++;
        }
        if (j < codelen[y[i]])
        {
            top++;
            stack[top].i = i;
            stack[top].j = j;
        }
        if (top >= 0)
        {
            j = stack[top].j;
            pos = stack[top].i;
            if (getbit64(map[y[pos]], 64 - codelen[y[pos]] + j))
                setbit32(&s[nlayers][i / 32], i % 32);
            stack[top].j++;
            if (stack[top].j == codelen[y[pos]])
            {
                delay[pos] += i - pos;
                if (delay[pos] > max_delay)
                {
                    max_delay = delay[pos];
                    max_delay_index = pos;
                }
                (*avgdelay) += i - pos;
                top--;
            }
        }
    }
    double totalbits = n * (nlayers + 1);
    int i = n;
    while (top >= 0)
    {
        totalbits++;
        int j = stack[top].j;
        pos = stack[top].i;
        if (getbit64(map[y[pos]], 64 - codelen[y[pos]] + j))
            setbit32(&s[nlayers][i / 32], i % 32);
        stack[top].j++;
        if (stack[top].j == codelen[y[pos]])
        {
            delay[pos] += n - 1 - pos;
            if (delay[pos] > max_delay)
            {
                max_delay = delay[pos];
                max_delay_index = pos;
            }
            (*avgdelay) += n - 1 - pos;
            top--;
        }
        i++;
    }
    *avgbits = totalbits / (double)n;
    *avgdelay = (*avgdelay) / (double)n;
    *max_char = y[max_delay_index];

/*
    // print stats for block
    for (int i = 0; i < n; i++)
        printf("Text[%d] = %d \tdelay %d\n", i, y[i], delay[i]);
*/
    //printf("\nText[%d] = %d has highest delay of %d on a text of length %d\n", max_delay_index, y[max_delay_index], max_delay, n);

    free(stack);
    return s;
}

int delta_decode(
    uint32_t **blce, // the BLCE representation of the text
    int i,           // the starting position of window we want to decode
    int j,           // the starting position of window we want to decode (for a single char i=j)
    int nlayers,     // the number of fixed layers in the BLCE
    HNODE *root,     // the root of the Huffman tree
    int n,           // the length of the text (position of the last char)
    int *dt,         // decoded text
    STACKNODE *stack //,
    // HNODE **fmap
)
{
    int pos = i;
    int top = -1;
    HNODE *node;
    int idiv32 = i / 32;
    int imod32 = i % 32;
    do
    {
        if (i < n)
        {
            int b = 0;
            node = root;
            while (b < nlayers && !node->is_char)
            {
                if (getbit32(blce[b][idiv32], imod32)) // gets the i-th bit
                    node = node->right;
                else
                    node = node->left;
                b++;
            }
            if (node->is_char)
            {
                dt[i] = node->c;
                if (top >= 0)
                {
                    if (getbit32(blce[nlayers][idiv32], imod32))
                        stack[top].node = stack[top].node->right;
                    else
                        stack[top].node = stack[top].node->left;
                    if (stack[top].node->is_char)
                    {
                        dt[stack[top].pos] = stack[top].node->c;
                        top--;
                    }
                }
            }
            else
            {
                if (getbit32(blce[nlayers][idiv32], imod32))
                    node = node->right;
                else
                    node = node->left;
                if (node->is_char)
                    dt[i] = node->c;
                else
                {
                    ++top;
                    stack[top].node = node;
                    stack[top].pos = i;
                }
            }
        }
        else
        {
            if (top >= 0)
            {
                if (getbit32(blce[nlayers][idiv32], imod32))
                    stack[top].node = stack[top].node->right;
                else
                    stack[top].node = stack[top].node->left;
                if (stack[top].node->is_char)
                {
                    dt[stack[top].pos] = stack[top].node->c;
                    top--;
                }
            }
        }
        i++;
        imod32++;
        if (imod32 == 32)
        {
            idiv32++;
            imod32 = 0;
        }
    } while (top >= 0 || dt[j] == '\0');
    i--;
    return (i < n ? i : (n - 1));
}
