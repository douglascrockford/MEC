/* MEC.h Matrix Exponentiation Cryptography
    Douglas Crockford
    Public Domain
    2017-04-14
*/

typedef unsigned char      mec8;
typedef unsigned long int  mec32;
typedef unsigned long long mec64;

extern mec32 mec_base[25];

extern void mec_generate(
    mec32 output[25],
    mec32 input[25],
    mec8 *private_key,
    mec32 private_key_length
);
