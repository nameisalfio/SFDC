# SFDC Text Compression Algorithm

The SFDC algorithm is designed for efficient text compression, particularly optimized for files with 4-byte values. The algorithm maintains all pending bits on the dynamic layer, ensuring effective compression and decompression of text data.

## Features

- **Dynamic Layer Management**: All pending bits are maintained on the dynamic layer for optimized performance.
- **4-Byte Value Optimization**: Specifically tuned for files containing 4-byte values.
- **Efficient Encoding and Decoding**: Utilizes a stack-based approach to manage encoding and decoding processes.

## Functions

### `uint32_t **delta_encode(int *y, int n, uint64_t *map, int *codelen, int nlayers, double *avgdelay, double *avgbits, int *max_char, int idx)`

Encodes the input text using the SFDC delta variant algorithm.

#### Parameters

- `int *y`: Input array of integers representing the text to be encoded.
- `int n`: Length of the input text.
- `uint64_t *map`: Mapping of characters to their corresponding codes.
- `int *codelen`: Array containing the length of each code.
- `int nlayers`: Number of fixed layers in the compression scheme.
- `double *avgdelay`: Output parameter for the average delay in encoding.
- `double *avgbits`: Output parameter for the average bits per character.
- `int *max_char`: Output parameter for the character with the maximum delay.
- `int idx`: Index parameter for encoding context.

#### Returns

- `uint32_t **`: A pointer to the encoded representation of the text.

### `int delta_decode(uint32_t **blce, int i, int j, int nlayers, HNODE *root, int n, int *dt, STACKNODE *stack)`

Decodes the encoded text using the SFDC delta variant algorithm.

#### Parameters

- `uint32_t **blce`: The encoded representation of the text.
- `int i`: Starting position of the window to decode.
- `int j`: Ending position of the window to decode (for a single character, `i=j`).
- `int nlayers`: Number of fixed layers in the encoded text.
- `HNODE *root`: Root of the Huffman tree used for decoding.
- `int n`: Length of the text.
- `int *dt`: Array to store the decoded text.
- `STACKNODE *stack`: Stack used for managing the decoding process.

#### Returns

- `int`: The position of the last decoded character.


