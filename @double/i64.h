
/*
 * =============================================================
  Header File for MATLAB 64 bit Integer Modulo Arithmetic Implementation

  Robin Ince

 * =============================================================
 */

typedef unsigned __int64 mod_t;
typedef unsigned __int32 uint32;
typedef unsigned __int64 uint64;
typedef __int64 int64;

// 2^61 - 1 (I think!) 
#define M61   0x1FFFFFFFFFFFFFFF
#define M32   0x00000000FFFFFFFF
#define SIGN_BIT  0x1000000000000000