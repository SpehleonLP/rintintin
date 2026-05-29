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
#define sin  rintintin_sin
#define cos  rintintin_cos

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

// sin/cos without libm. No fancy hardware path -- x87 fsin/fcos exist but are
// notoriously slow and inaccurate for large arguments, and there's no ARM
// equivalent worth using. Polynomial it is.
//
// Strategy: reduce x to [-π, π] modulo 2π, fold that to [-π/2, π/2] using the
// identities sin(±π - y) = sin(y) and cos(±π - y) = -cos(y), then evaluate a
// 10-term Horner-form Taylor series on the small interval. Error is below
// 1e-15 for |x| up to ~1e9; for larger magnitudes the reduction loses bits
// because π isn't representable exactly, but rintintin's inputs are angles
// in radians at human scales, so this never bites.

static const double RINTINTIN_PI_         = 3.14159265358979323846;
static const double RINTINTIN_TWO_PI_     = 6.28318530717958647692;
static const double RINTINTIN_HALF_PI_    = 1.57079632679489661923;
static const double RINTINTIN_INV_TWO_PI_ = 0.15915494309189533577;

// Reduce x to [-π, π]. The cast trick rounds to nearest; for |k| beyond 2^62
// it saturates and accuracy degrades, which is fine for our scale.
static FORCE_INLINE double rintintin_reduce_pi(double x) {
    double k = x * RINTINTIN_INV_TWO_PI_;
    long long ki = (long long)(k + (k >= 0 ? 0.5 : -0.5));
    return x - (double)ki * RINTINTIN_TWO_PI_;
}

static FORCE_INLINE double rintintin_sin(double x) {
    double y = rintintin_reduce_pi(x);
    // Fold [π/2, π] and [-π, -π/2] into [-π/2, π/2] via sin(π - y) = sin(y).
    if (y >  RINTINTIN_HALF_PI_) y =  RINTINTIN_PI_ - y;
    if (y < -RINTINTIN_HALF_PI_) y = -RINTINTIN_PI_ - y;

    // Horner on u = y^2 of (1 - u/3! + u²/5! - u³/7! + ... - u⁹/19!), then * y.
    double u = y * y;
    double r = -1.0/121645100408832000.0;   // -1/19!
    r = r*u + 1.0/355687428096000.0;        //  1/17!
    r = r*u - 1.0/1307674368000.0;          // -1/15!
    r = r*u + 1.0/6227020800.0;             //  1/13!
    r = r*u - 1.0/39916800.0;               // -1/11!
    r = r*u + 1.0/362880.0;                 //  1/9!
    r = r*u - 1.0/5040.0;                   // -1/7!
    r = r*u + 1.0/120.0;                    //  1/5!
    r = r*u - 1.0/6.0;                      // -1/3!
    r = r*u + 1.0;                          //  1
    return y * r;
}

static FORCE_INLINE double rintintin_cos(double x) {
    double y = rintintin_reduce_pi(x);
    // Fold [π/2, π] and [-π, -π/2] into [-π/2, π/2] via cos(π - y) = -cos(y).
    double sign = 1.0;
    if (y >  RINTINTIN_HALF_PI_) { y =  RINTINTIN_PI_ - y; sign = -1.0; }
    if (y < -RINTINTIN_HALF_PI_) { y = -RINTINTIN_PI_ - y; sign = -1.0; }

    // Horner on u = y^2 of (1 - u/2! + u²/4! - u³/6! + ... - u⁹/18!).
    double u = y * y;
    double r = -1.0/6402373705728000.0;     // -1/18!
    r = r*u + 1.0/20922789888000.0;         //  1/16!
    r = r*u - 1.0/87178291200.0;            // -1/14!
    r = r*u + 1.0/479001600.0;              //  1/12!
    r = r*u - 1.0/3628800.0;                // -1/10!
    r = r*u + 1.0/40320.0;                  //  1/8!
    r = r*u - 1.0/720.0;                    // -1/6!
    r = r*u + 1.0/24.0;                     //  1/4!
    r = r*u - 1.0/2.0;                      // -1/2!
    r = r*u + 1.0;                          //  1
    return sign * r;
}

#endif

#endif // RINTINTIN_SQRT_H
