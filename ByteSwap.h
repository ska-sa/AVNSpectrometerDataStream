//Manual implementation of byte swap functions for old versions of GCC which don't have them

#ifndef BYTE_SWAP_H
#define BYTE_SWAP_H

//System includes
#ifdef _WIN32
#include <stdint.h>
#else
#include <inttypes.h>
#endif

//Library includes

//Local includes

//For *nix flavours
#ifndef _WIN32

//Get GCC version
#define GCC_VERSION (__GNUC__ * 10000 \
    + __GNUC_MINOR__ * 100 \
    + __GNUC_PATCHLEVEL__)

//For x86 platforms (assuming nothing else weird is being used...)
#ifndef __amd64__

//For GCC <4.8
#if GCC_VERSION < 40800

#pragma message "Detected GCC < 4.8 on x86, Defining missing 'uint16_t __builtin_bswap16(uint16_t)' function"

//There is no __builtin_bswap16 for GCC <4.8 in x86
static inline uint16_t __builtin_bswap16(uint16_t u16Value)
{
    return ( u16Value << 8 )|( u16Value >> 8 );
}

#endif // GCC_VERSION < 40800

#endif // ! __amd64__

#endif // ! WIN32

#endif //BYTE_SWAP_H
