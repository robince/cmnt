/*
 * =============================================================
  Header File for MATLAB 64 bit Integer Modulo Arithmetic Implementation

  Robin Ince

 * =============================================================
 */

#include <stdint.h>
typedef uint64_t mod_t;
typedef uint32_t uint32;
typedef uint64_t uint64;
typedef int64_t int64;

// 2^61 - 1 (I think!) 
#define M61   0x1FFFFFFFFFFFFFFF
#define M32   0x00000000FFFFFFFF
#define SIGN_BIT  0x1000000000000000
