#ifndef RINTINTIN_SQRT_H
#define RINTINTIN_SQRT_H


// see this is the whole stupid thing
// i can't just link <math.h> becuase the linker will pull in acos
// which isn't inline ASM, rather its part of libm
// so then you're linked to all of libm, and OOPS if you have a different version
// than the platform you're shipping on, and now it can't start, 
// so you need to package the fuckin libm into a fuckin appimage and 
// what the fuck is going on? 

#if RINTINTIN_USE_MATH_H
#include <math.h>
#else

#define fabs rintintin_abs
#define fmin rintintin_min
#define fmax rintintin_max
#define sqrt rintintin_sqrt

#ifdef __GNUC__
    #define FORCE_INLINE __attribute__((always_inline)) inline
#elif defined(_MSC_VER)
    #define FORCE_INLINE __forceinline
#else
    #define FORCE_INLINE inline
#endif

FORCE_INLINE double rintintin_abs(double a) { return a < 0? -a : a; }
FORCE_INLINE double rintintin_min(double a, double b) { return a < b? a : b; }
FORCE_INLINE double rintintin_max(double a, double b) { return a > b? a : b; }


// Use inline ASM on x86/x64 because why the hell not
#if defined(__x86_64__) || defined(__i386__)
    #define HAS_X86_ASM 1
#elif defined(__aarch64__) || defined(__arm__)
    #define HAS_ARM_ASM 1
#endif

static FORCE_INLINE double rintintin_sqrt(double x) {
#ifdef HAS_X86_ASM
    // Use the beautiful fsqrt instruction because Intel isn't completely insane
    double result;
    __asm__ volatile (
        "fsqrt"
        : "=t" (result)
        : "0" (x)
    );
    return result;
    
#elif defined(HAS_ARM_ASM) && defined(__aarch64__)
    // ARM64 has a sqrt instruction too, bless their hearts
    double result;
    __asm__ volatile (
        "fsqrt %d0, %d1"
        : "=w" (result)
        : "w" (x)
    );
    return result;
    
#else
    // FINE. We'll do it the hard way with Newton-Raphson
    // because apparently we can't have nice things
    
    if (x < 0.0) {
        // NaN for negative inputs, just like the "real" sqrt
        return 0.0/0.0; // This generates NaN without libm
    }
    
    if (x == 0.0 || x == 1.0) {
        return x; // Easy cases
    }

    // Quick and dirty initial guess using bit manipulation
    // This is the "magic number" method, don't ask me to explain it again
    union {
        double d;
        long long i;
    } u;
    u.d = x;
    u.i = (1LL << 61) + (u.i >> 1) - (1LL << 51);
    
    double guess = u.d;
    
    // Newton-Raphson: x_n+1 = 0.5 * (x_n + (S / x_n))
    // Do it 4 times because that's usually enough for double precision
    // and I'm tired of this bullshit
    for (int i = 0; i < 4; i++) {
        guess = 0.5 * (guess + x / guess);
    }
    
    return guess;
#endif
}
#endif

#endif // RINTINTIN_SQRT_H
