title   mec.x64.asm: mec for x64.

; 2017-04-16
; Public Domain

; No warranty expressed or implied. Use at your own risk. You have been warned.

; MEC is a public key cryptographic system built on modular matrix
; exponentiation. MEC has not been adequately reviewed, so its security
; properties are currently unknown.

; This is an implementation in X64 assembly language.

UNIX    equ 0                   ; calling convention: 0 for Win64, 1 for Unix

; Unix passes parameters in
;   r7  r6  r2  r1
; Win64 passes parameters in
;   r1  r2  r8  r9
; What a world, what a world.

NO_LEAK equ 1                   ; 0 for speed, 1 for leak resistence

; This file can be configured to run fast, or to run in near constant time to
; minimize key leakage.

; This file publishes two public symbols:

public mec_base; array of 25 mec32

public mec_generate;(
;   output: array of 25 mec32,
;   input: array of 25 mec32,
;   public_key: array of mec8,
;   public_key_length: mec32
;)

MODULUS equ 4294967291          ; Largest prime less than 2^32

; Repair the register names. Over the long and twisted evolution of x86, the
; register names have picked up some weird, inconsistent conventions. We can
; simplify, naming them r0 thru r15. (We don't use rsp or rbp.)

r0  equ     rax
r1  equ     rcx
r2  equ     rdx
r3  equ     rbx
r6  equ     rsi
r7  equ     rdi

e0  equ     eax
e2  equ     edx

r3l equ     bl
r3h equ     bh

pad macro
    align   16
    endm

;  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

mec_data segment page read

mec_base:
            dword   1760000027, 3840000041, 1120000031, 3200000087,  480000019
            dword    640000019, 1920000037, 4000000007, 1280000003, 2560000019
            dword   2720000011,  800000011, 2080000007, 3360000041, 1440000043
            dword   1600000009, 2880000007,  160000003, 2240000069, 3520000007
            dword   3680000003,  960000011, 3040000009,  320000077, 2400000011

mec_data ends

;  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

mec_state segment page write

save:                       ; Preserved registers

save_r3     qword   0
save_r6     qword   0
save_r7     qword   0
save_r12    qword   0
save_r13    qword   0
save_r14    qword   0
save_r15    qword   0
            pad

; matrix contains the current value of the matrix.

matrix:
            qword   0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            qword   0, 0, 0
            pad

if NO_LEAK

; decoy will take an unneeded product.

decoy:
            qword   0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            qword   0, 0, 0
            pad

endif

; product is the result of matrix multiplication.

product:
            qword   0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            qword   0, 0, 0

mec_state ends

;  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

copy macro to,from

;; to and from identify two 5x5 mec32 matrixes. Copy matrix from to matrix to.
;; Move two cells at a time.

    mov     r0,[from]
    mov     [to],r0
    mov     r0,[from+1*8]
    mov     [to+1*8],r0
    mov     r0,[from+16]
    mov     [to+16],r0
    mov     r0,[from+24]
    mov     [to+24],r0
    mov     r0,[from+32]
    mov     [to+32],r0
    mov     r0,[from+40]
    mov     [to+40],r0
    mov     r0,[from+48]
    mov     [to+48],r0
    mov     r0,[from+56]
    mov     [to+56],r0
    mov     r0,[from+64]
    mov     [to+64],r0
    mov     r0,[from+72]
    mov     [to+72],r0
    mov     r0,[from+80]
    mov     [to+80],r0
    mov     r0,[from+88]
    mov     [to+88],r0
    mov     r0,[from+96]
    mov     [to+96],e0
    endm

clear macro it

;; Set 25 dwords to zero, two at a time.

    xor     r0,r0
    mov     [it],r0
    mov     [it+8],r0
    mov     [it+16],r0
    mov     [it+24],r0
    mov     [it+32],r0
    mov     [it+40],r0
    mov     [it+48],r0
    mov     [it+56],r0
    mov     [it+64],r0
    mov     [it+72],r0
    mov     [it+80],r0
    mov     [it+88],r0
    mov     [it+96],e0
    endm

mul___ macro  i, j

;; Multiply two mec32 cells, producing a 64 bit product. r1 and r13 contain
;; pointers to the two matrixes. The product of r1[i] * r13[j] goes in r15.

    mov     e0,[r1+i*4]
    mov     e2,[r13+j*4]
    mul     r2
    mov     r15,r0
    endm

muladd macro  i, j

;; Multiply two mec32 cells, producing a 64 bit product. r1 and r13 contain
;; pointers to the two matrixes. Add the product of r1[i] * r13[j] to r15,
;; taking care to avoid overflow errors.

    mov     e0,[r1+i*4]
    mov     e2,[r13+j*4]
    mul     r2

; Add the product to the accumulator. If it overflows, make a modulo
; correction. The bit that is lost is the 64th bit. That bit's contribution
; to the modulo total is 25, because 2^64 mod MODULUS is 25. So if the add
; produces a carry, then add 25.

    add     r15,r0
    mov     r0,r12
    cmovnc  r0,r2
    add     r15,r0
    endm

modstore macro i

;; Divide the sum by the modulus and store the remainder in the product matrix.
;; r2 had better be zero. Clear r15.

    mov     r0,r15
    div     r14
    mov     [r10+i*4],e2
    xor     r15,r15
    endm

matrixmultiply macro

;;  r1      left matrix
;;  r13     right matrix
;;  r10     product matrix

    mul___  0,0
    muladd  1,5
    muladd  2,10
    muladd  3,15
    muladd  4,20
    modstore 0
    mul___  0,1
    muladd  1,6
    muladd  2,11
    muladd  3,16
    muladd  4,21
    modstore 1
    mul___  0,2
    muladd  1,7
    muladd  2,12
    muladd  3,17
    muladd  4,22
    modstore 2
    mul___  0,3
    muladd  1,8
    muladd  2,13
    muladd  3,18
    muladd  4,23
    modstore 3
    mul___  0,4
    muladd  1,9
    muladd  2,14
    muladd  3,19
    muladd  4,24
    modstore 4
    mul___  5,0
    muladd  6,5
    muladd  7,10
    muladd  8,15
    muladd  9,20
    modstore 5
    mul___  5,1
    muladd  6,6
    muladd  7,11
    muladd  8,16
    muladd  9,21
    modstore 6
    mul___  5,2
    muladd  6,7
    muladd  7,12
    muladd  8,17
    muladd  9,22
    modstore 7
    mul___  5,3
    muladd  6,8
    muladd  7,13
    muladd  8,18
    muladd  9,23
    modstore 8
    mul___  5,4
    muladd  6,9
    muladd  7,14
    muladd  8,19
    muladd  9,24
    modstore 9
    mul___  10,0
    muladd  11,5
    muladd  12,10
    muladd  13,15
    muladd  14,20
    modstore 10
    mul___  10,1
    muladd  11,6
    muladd  12,11
    muladd  13,16
    muladd  14,21
    modstore 11
    mul___  10,2
    muladd  11,7
    muladd  12,12
    muladd  13,17
    muladd  14,22
    modstore 12
    mul___  10,3
    muladd  11,8
    muladd  12,13
    muladd  13,18
    muladd  14,23
    modstore 13
    mul___  10,4
    muladd  11,9
    muladd  12,14
    muladd  13,19
    muladd  14,24
    modstore 14
    mul___  15,0
    muladd  16,5
    muladd  17,10
    muladd  18,15
    muladd  19,20
    modstore 15
    mul___  15,1
    muladd  16,6
    muladd  17,11
    muladd  18,16
    muladd  19,21
    modstore 16
    mul___  15,2
    muladd  16,7
    muladd  17,12
    muladd  18,17
    muladd  19,22
    modstore 17
    mul___  15,3
    muladd  16,8
    muladd  17,13
    muladd  18,18
    muladd  19,23
    modstore 18
    mul___  15,4
    muladd  16,9
    muladd  17,14
    muladd  18,19
    muladd  19,24
    modstore 19
    mul___  20,0
    muladd  21,5
    muladd  22,10
    muladd  23,15
    muladd  24,20
    modstore 20
    mul___  20,1
    muladd  21,6
    muladd  22,11
    muladd  23,16
    muladd  24,21
    modstore 21
    mul___  20,2
    muladd  21,7
    muladd  22,12
    muladd  23,17
    muladd  24,22
    modstore 22
    mul___  20,3
    muladd  21,8
    muladd  22,13
    muladd  23,18
    muladd  24,23
    modstore 23
    mul___  20,4
    muladd  21,9
    muladd  22,14
    muladd  23,19
    muladd  24,24
    modstore 24

    endm

;  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

mec_code segment page execute

mec_generate:

; Registers r0, r1, r2, r8, r9, r10, and r11 are clobbered.
; Registers r6 and r7 are clobbered on Unix. The other registers
; are not disturbed.

;   r0      low product, quotient
;   r1      matrix
;   r2      high product, remainder
;   r3      bit, byte
;   r6      input
;   r7      output
;   r8      private key
;   r9      private key length
;   r10     product matrix
;   r11     decoy
;   r12     25
;   r13     right matrix
;   r14     modulus
;   r15     accumulator

    mov     save_r3,r3

if UNIX
    mov     r8,r2               ;; UNIX
    mov     r9,r1               ;; UNIX
else
    mov     save_r6,r6          ;; WIN64
    mov     save_r7,r7          ;; WIN64
    mov     r6,r2               ;; WIN64
    mov     r7,r1               ;; WIN64
endif

    mov     save_r12,r12
    mov     save_r13,r13
    mov     save_r14,r14
    mov     save_r15,r15

    mov     r1,matrix
    copy    r1,r6               ; matrix <- input
    mov     r10,product

if NO_LEAK
    mov     r11,decoy
endif

    mov     r12,25
    mov     r14,MODULUS
    xor     r15,r15
    pad

outer:

; Get the next bit of the private key.
;
    mov     r3l,[r8]            ; the next byte of the private key
    mov     r3h,128             ; a byte with the high bit set
    add     r8,1                ; r8 points to the next byte of the private key
    pad

inner:

    mov     r13,r1
    matrixmultiply              ; product <- matrix * matrix
    copy    r1,r10              ; matrix <- product

if NO_LEAK

    mov     r13,r6
    matrixmultiply              ; product <- matrix * input
    mov     r2,r1
    test    r3l,r3h             ; if next bit
    cmovz   r2,r11              ;      matrix <- product
    copy    r2,r10              ; else decoy <- product
                                ; fi
else

    test    r3l,r3h             ; if next bit
    jz      skip
    mov     r13,r6
    matrixmultiply              ;     product <- matrix * input
    copy    r1,r10              ;     matrix <- product
    pad

skip:                           ; fi

endif

    shr     r3h,1               ; is there another bit in this byte?
    jnz     inner

    sub     r9,1                ; is there another byte of private key?
    jnz     outer

    copy    r7,r1               ; output <- matrix

; Erasure. Leave nothing in mec_state or the clobbered registers.

    clear   r1                  ; clear matrix
    mov     r3,[save_r3]

if UNIX
    xor     r6,r6               ; UNIX
    xor     r7,r7               ; UNIX
else
    mov     r6,[save_r6]        ; WIN64
    mov     r7,[save_r7]        ; WIN64
endif

    xor     r8,r8
    xor     r9,r9
    clear   r10                 ; clear product

if NO_LEAK
    clear   r11                 ; clear decoy
endif

    mov     r12,[save_r12]
    mov     r13,[save_r13]
    mov     r14,[save_r14]
    mov     r15,[save_r15]
    mov     r2,save
    clear   r2                  ; clear save
    ret

mec_code ends
    end
