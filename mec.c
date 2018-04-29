/* MEC.c Matrix Exponentiation Cryptography
    Douglas Crockford
    2017-04-14
    Public Domain
*/

#include "MEC.h"

/*
    The modulus is the largest prime number that is less than 2**32.
*/

static const mec64 modulus = 4294967291;

/*
    The modulus correction is (2**64) % modulus. It is used to correct
    for mec64 additions that overflow.
*/

static const mec64 modulus_correction = 25;

/*
    The mec base matrix is next_prime(magic_square(Mars) * 160000000).
*/

mec32 mec_base[25] = {
    1760000027, 3840000041, 1120000031, 3200000087,  480000019,
     640000019, 1920000037, 4000000007, 1280000003, 2560000019,
    2720000011,  800000011, 2080000007, 3360000041, 1440000043,
    1600000009, 2880000007,  160000003, 2240000069, 3520000007,
    3680000003,  960000011, 3040000009,  320000077, 2400000011
};

static mec32 x[25];
static mec32 z[25];

static mec64 madd2(mec64 a, mec64 b) {

/*
    Modular addition of two terms. If the addition overflows (which is likely)
    then correct the sum.
*/

    mec64 sum = a + b;
    if (sum < a) {
        sum = sum + modulus_correction;
    }
    return sum;
}


static mec32 madd5(mec64 a, mec64 b, mec64 c, mec64 d, mec64 e) {

/*
    Modular addition of five terms.
*/

    return (mec32)((
        madd2(a, madd2(b, madd2(c, madd2(d, e))))
    ) % modulus);
}

static void mm(mec32 y[25]) {

/*
    Modular multiplication[5x5]. The multiplication is unrolled.
    z = x * y
    x = z
*/

    z[0] = madd5(
        (mec64)x[0] * (mec64)y[0], (mec64)x[1] * (mec64)y[5],
        (mec64)x[2] * (mec64)y[10], (mec64)x[3] * (mec64)y[15],
        (mec64)x[4] * (mec64)y[20]
    );
    z[1] = madd5(
        (mec64)x[0] * (mec64)y[1], (mec64)x[1] * (mec64)y[6],
        (mec64)x[2] * (mec64)y[11], (mec64)x[3] * (mec64)y[16],
        (mec64)x[4] * (mec64)y[21]
    );
    z[2] = madd5(
        (mec64)x[0] * (mec64)y[2], (mec64)x[1] * (mec64)y[7],
        (mec64)x[2] * (mec64)y[12], (mec64)x[3] * (mec64)y[17],
        (mec64)x[4] * (mec64)y[22]
    );
    z[3] = madd5(
        (mec64)x[0] * (mec64)y[3], (mec64)x[1] * (mec64)y[8],
        (mec64)x[2] * (mec64)y[13], (mec64)x[3] * (mec64)y[18],
        (mec64)x[4] * (mec64)y[23]
    );

    z[4] = madd5(
        (mec64)x[0] * (mec64)y[4], (mec64)x[1] * (mec64)y[9],
        (mec64)x[2] * (mec64)y[14], (mec64)x[3] * (mec64)y[19],
        (mec64)x[4] * (mec64)y[24]
    );
    z[5] = madd5(
        (mec64)x[5] * (mec64)y[0], (mec64)x[6] * (mec64)y[5],
        (mec64)x[7] * (mec64)y[10], (mec64)x[8] * (mec64)y[15],
        (mec64)x[9] * (mec64)y[20]
    );
    z[6] = madd5(
        (mec64)x[5] * (mec64)y[1], (mec64)x[6] * (mec64)y[6],
        (mec64)x[7] * (mec64)y[11], (mec64)x[8] * (mec64)y[16],
        (mec64)x[9] * (mec64)y[21]
    );
    z[7] = madd5(
        (mec64)x[5] * (mec64)y[2], (mec64)x[6] * (mec64)y[7],
        (mec64)x[7] * (mec64)y[12], (mec64)x[8] * (mec64)y[17],
        (mec64)x[9] * (mec64)y[22]
    );
    z[8] = madd5(
        (mec64)x[5] * (mec64)y[3], (mec64)x[6] * (mec64)y[8],
        (mec64)x[7] * (mec64)y[13], (mec64)x[8] * (mec64)y[18],
        (mec64)x[9] * (mec64)y[23]
    );
    z[9] = madd5(
        (mec64)x[5] * (mec64)y[4], (mec64)x[6] * (mec64)y[9],
        (mec64)x[7] * (mec64)y[14], (mec64)x[8] * (mec64)y[19],
        (mec64)x[9] * (mec64)y[24]
    );
    z[10] = madd5(
        (mec64)x[10] * (mec64)y[0], (mec64)x[11] * (mec64)y[5],
        (mec64)x[12] * (mec64)y[10], (mec64)x[13] * (mec64)y[15],
        (mec64)x[14] * (mec64)y[20]
    );
    z[11] = madd5(
        (mec64)x[10] * (mec64)y[1], (mec64)x[11] * (mec64)y[6],
        (mec64)x[12] * (mec64)y[11], (mec64)x[13] * (mec64)y[16],
        (mec64)x[14] * (mec64)y[21]
    );
    z[12] = madd5(
        (mec64)x[10] * (mec64)y[2], (mec64)x[11] * (mec64)y[7],
        (mec64)x[12] * (mec64)y[12], (mec64)x[13] * (mec64)y[17],
        (mec64)x[14] * (mec64)y[22]
    );
    z[13] = madd5(
        (mec64)x[10] * (mec64)y[3], (mec64)x[11] * (mec64)y[8],
        (mec64)x[12] * (mec64)y[13], (mec64)x[13] * (mec64)y[18],
        (mec64)x[14] * (mec64)y[23]
    );
    z[14] = madd5(
        (mec64)x[10] * (mec64)y[4], (mec64)x[11] * (mec64)y[9],
        (mec64)x[12] * (mec64)y[14], (mec64)x[13] * (mec64)y[19],
        (mec64)x[14] * (mec64)y[24]
    );
    z[15] = madd5(
        (mec64)x[15] * (mec64)y[0], (mec64)x[16] * (mec64)y[5],
        (mec64)x[17] * (mec64)y[10], (mec64)x[18] * (mec64)y[15],
        (mec64)x[19] * (mec64)y[20]
    );
    z[16] = madd5(
        (mec64)x[15] * (mec64)y[1], (mec64)x[16] * (mec64)y[6],
        (mec64)x[17] * (mec64)y[11], (mec64)x[18] * (mec64)y[16],
        (mec64)x[19] * (mec64)y[21]
    );
    z[17] = madd5(
        (mec64)x[15] * (mec64)y[2], (mec64)x[16] * (mec64)y[7],
        (mec64)x[17] * (mec64)y[12], (mec64)x[18] * (mec64)y[17],
        (mec64)x[19] * (mec64)y[22]
    );
    z[18] = madd5(
        (mec64)x[15] * (mec64)y[3], (mec64)x[16] * (mec64)y[8],
        (mec64)x[17] * (mec64)y[13], (mec64)x[18] * (mec64)y[18],
        (mec64)x[19] * (mec64)y[23]
    );
    z[19] = madd5(
        (mec64)x[15] * (mec64)y[4], (mec64)x[16] * (mec64)y[9],
        (mec64)x[17] * (mec64)y[14], (mec64)x[18] * (mec64)y[19],
        (mec64)x[19] * (mec64)y[24]
    );
    z[20] = madd5(
        (mec64)x[20] * (mec64)y[0], (mec64)x[21] * (mec64)y[5],
        (mec64)x[22] * (mec64)y[10], (mec64)x[23] * (mec64)y[15],
        (mec64)x[24] * (mec64)y[20]
    );
    z[21] = madd5(
        (mec64)x[20] * (mec64)y[1], (mec64)x[21] * (mec64)y[6],
        (mec64)x[22] * (mec64)y[11], (mec64)x[23] * (mec64)y[16],
        (mec64)x[24] * (mec64)y[21]
    );
    z[22] = madd5(
        (mec64)x[20] * (mec64)y[2], (mec64)x[21] * (mec64)y[7],
        (mec64)x[22] * (mec64)y[12], (mec64)x[23] * (mec64)y[17],
        (mec64)x[24] * (mec64)y[22]
    );
    z[23] = madd5(
        (mec64)x[20] * (mec64)y[3], (mec64)x[21] * (mec64)y[8],
        (mec64)x[22] * (mec64)y[13], (mec64)x[23] * (mec64)y[18],
        (mec64)x[24] * (mec64)y[23]
    );
    z[24] = madd5(
        (mec64)x[20] * (mec64)y[4], (mec64)x[21] * (mec64)y[9],
        (mec64)x[22] * (mec64)y[14], (mec64)x[23] * (mec64)y[19],
        (mec64)x[24] * (mec64)y[24]
    );

    x[0] = z[0];
    x[1] = z[1];
    x[2] = z[2];
    x[3] = z[3];
    x[4] = z[4];

    x[5] = z[5];
    x[6] = z[6];
    x[7] = z[7];
    x[8] = z[8];
    x[9] = z[9];

    x[10] = z[10];
    x[11] = z[11];
    x[12] = z[12];
    x[13] = z[13];
    x[14] = z[14];

    x[15] = z[15];
    x[16] = z[16];
    x[17] = z[17];
    x[18] = z[18];
    x[19] = z[19];

    x[20] = z[20];
    x[21] = z[21];
    x[22] = z[22];
    x[23] = z[23];
    x[24] = z[24];
}

void mec_generate(
    mec32 output[25],
    mec32 input[25],
    mec8 *private_key,
    mec32 private_key_length
) {

/*
    Crunch the input and the private key with Exponentiation by squaring and
    multiplying.

    The caller must provide an array of mec32[25] to receive the output.

    If the input is mec_base, then the output is a public key.

    If the input is a public key, then the output is a shared secret.

    The private key is an array of random bytes.

    The recommended private key length is at least 64.
*/

    mec8 bit;
    mec8 byte;
    mec32 byte_nr;

    x[0] = input[0];
    x[1] = input[1];
    x[2] = input[2];
    x[3] = input[3];
    x[4] = input[4];

    x[5] = input[5];
    x[6] = input[6];
    x[7] = input[7];
    x[8] = input[8];
    x[9] = input[9];

    x[10] = input[10];
    x[11] = input[11];
    x[12] = input[12];
    x[13] = input[13];
    x[14] = input[14];

    x[15] = input[15];
    x[16] = input[16];
    x[17] = input[17];
    x[18] = input[18];
    x[19] = input[19];

    x[20] = input[20];
    x[21] = input[21];
    x[22] = input[22];
    x[23] = input[23];
    x[24] = input[24];

    for (byte_nr = 0; byte_nr < private_key_length; byte_nr += 1) {
        byte = private_key[byte_nr];
        for (bit = 128; bit != 0; bit >>= 1) {
            mm(x);
            if ((bit & byte) != 0) {
                mm(input);
            }
        }
    }

    output[0] = x[0];
    output[1] = x[1];
    output[2] = x[2];
    output[3] = x[3];
    output[4] = x[4];

    output[5] = x[5];
    output[6] = x[6];
    output[7] = x[7];
    output[8] = x[8];
    output[9] = x[9];

    output[10] = x[10];
    output[11] = x[11];
    output[12] = x[12];
    output[13] = x[13];
    output[14] = x[14];

    output[15] = x[15];
    output[16] = x[16];
    output[17] = x[17];
    output[18] = x[18];
    output[19] = x[19];

    output[20] = x[20];
    output[21] = x[21];
    output[22] = x[22];
    output[23] = x[23];
    output[24] = x[24];

    x[0] = 0;
    x[1] = 0;
    x[2] = 0;
    x[3] = 0;

    x[4] = 0;
    x[5] = 0;
    x[6] = 0;
    x[7] = 0;

    x[8] = 0;
    x[9] = 0;
    x[10] = 0;
    x[11] = 0;

    x[12] = 0;
    x[13] = 0;
    x[14] = 0;
    x[15] = 0;

    z[0] = 0;
    z[1] = 0;
    z[2] = 0;
    z[3] = 0;

    z[4] = 0;
    z[5] = 0;
    z[6] = 0;
    z[7] = 0;

    z[8] = 0;
    z[9] = 0;
    z[10] = 0;
    z[11] = 0;

    z[12] = 0;
    z[13] = 0;
    z[14] = 0;
    z[15] = 0;

    byte = 0;
}
