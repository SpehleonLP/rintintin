#ifndef SIMD_DOUBLE_PTR_H
#define SIMD_DOUBLE_PTR_H

#if defined(__GNUC__) || defined(__clang__)
	#ifdef __has_builtin
		#if __has_builtin(__atomic_add_fetch)
			#define HAS_ATOMIC 1
			#define HAS_SYNC 0
		#else
			#define HAS_ATOMIC 0
			#define HAS_SYNC 1
		#endif
	#else
		// GCC without __has_builtin, use version check
		#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 8))
			#define HAS_ATOMIC 1
			#define HAS_SYNC 0
		#else
			#define HAS_ATOMIC 0
			#define HAS_SYNC 1
		#endif
	#endif
#else
#define HAS_ATOMIC 0
#define HAS_SYNC 0
#endif


#if __STDC_VERSION__ >= 199901L
	#define RESTRICT restrict
#else
	#define RESTRICT 
#endif

/// indicate a variable is aligned
#if __STDC_VERSION__ >= 201112L
	#define ALIGN_AS(x) _Alignas(x)
#elif defined(__GNUC__) || defined(__clang__)
    #define ALIGN_AS(x) __attribute__((aligned(x)))
#elif defined(_MSC_VER)
    #define ALIGN_AS(x) __declspec((align(x)))
#elif defined(__INTEL_COMPILER) || defined(__ICC)
    #define ALIGN_AS(x) __declspec(align(x))
#else
    #define ALIGN_AS(x) 
#endif

/// indicate what a pointer is pointing to is aligned
#if defined(__GNUC__) || defined(__clang__)
    #define ASSUME_ALIGNED_T(x) __attribute__((assume_aligned(x)))
#elif defined(_MSC_VER)
    #define ASSUME_ALIGNED_T(x)
#elif defined(__INTEL_COMPILER) || defined(__ICC)
    #define ASSUME_ALIGNED_T(x)
#else
    #define ASSUME_ALIGNED_T(x) 
#endif

/// runtime hints to indicate what a pointer is pointing to is aligned
#if defined(__GNUC__) || defined(__clang__)
    #define ASSUME_ALIGNED(data, SIMD_WIDTH) (data = __builtin_assume_aligned(data, SIMD_WIDTH))
#elif defined(_MSC_VER)
    #define ASSUME_ALIGNED(data, SIMD_WIDTH)  __assume((uintptr_t)data % SIMD_WIDTH == 0)
#elif defined(__INTEL_COMPILER) || defined(__ICC)
    #define ASSUME_ALIGNED(data, SIMD_WIDTH)  __assume_aligned(data, SIMD_WIDTH)
#else
    #define ASSUME_ALIGNED(data, SIMD_WIDTH) 
#endif


/* Hierarchical SIMD instruction set detection and alignment */
#if defined(__AVX512F__)
    /* AVX512: 512-bit vectors, 64-byte alignment optimal */
    #define SIMD_ALIGN ALIGN_AS(64)
    #define SIMD_WIDTH 64
    #define SIMD_DOUBLES_PER_VEC 8
    #define SIMD_LEVEL "AVX512"
    #define SIMD_TARGET_FUNC __attribute__((target("avx512f")))
    

#elif defined(__AVX2__)
    /* AVX2: 256-bit vectors, 32-byte alignment optimal */
    #define SIMD_ALIGN ALIGN_AS(32)
    #define SIMD_WIDTH 32
    #define SIMD_DOUBLES_PER_VEC 4
    #define SIMD_LEVEL "AVX2"
    #define SIMD_TARGET_FUNC __attribute__((target("avx2")))
    

#elif defined(__AVX__)
    /* AVX: 256-bit vectors, 32-byte alignment optimal */
    #define SIMD_ALIGN ALIGN_AS(32)
    #define SIMD_WIDTH 32
    #define SIMD_DOUBLES_PER_VEC 4
    #define SIMD_LEVEL "AVX"
    #define SIMD_TARGET_FUNC __attribute__((target("avx")))

#elif defined(__SSE4_2__)
    /* SSE4.2: 128-bit vectors, 16-byte alignment optimal */
    #define SIMD_ALIGN ALIGN_AS(16)
    #define SIMD_WIDTH 16
    #define SIMD_DOUBLES_PER_VEC 2
    #define SIMD_LEVEL "SSE4.2"
    #define SIMD_TARGET_FUNC __attribute__((target("sse4.2")))

#elif defined(__SSE4_1__)
    /* SSE4.1: 128-bit vectors, 16-byte alignment optimal */
    #define SIMD_ALIGN ALIGN_AS(16)
    #define SIMD_WIDTH 16
    #define SIMD_DOUBLES_PER_VEC 2
    #define SIMD_LEVEL "SSE4.1"
    #define SIMD_TARGET_FUNC __attribute__((target("sse4.1")))

#elif defined(__SSSE3__)
    /* SSSE3: 128-bit vectors, 16-byte alignment optimal */
    #define SIMD_ALIGN ALIGN_AS(16)
    #define SIMD_WIDTH 16
    #define SIMD_DOUBLES_PER_VEC 2
    #define SIMD_LEVEL "SSSE3"
    #define SIMD_TARGET_FUNC __attribute__((target("ssse3")))

#elif defined(__SSE3__)
    /* SSE3: 128-bit vectors, 16-byte alignment optimal */
    #define SIMD_ALIGN ALIGN_AS(16)
    #define SIMD_WIDTH 16
    #define SIMD_DOUBLES_PER_VEC 2
    #define SIMD_LEVEL "SSE3"
    #define SIMD_TARGET_FUNC __attribute__((target("sse3")))

#elif defined(__SSE2__) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2) || defined(_M_X64)
    /* SSE2: 128-bit vectors, 16-byte alignment optimal - baseline for x64 */
    #define SIMD_ALIGN ALIGN_AS(16)
    #define SIMD_WIDTH 16
    #define SIMD_DOUBLES_PER_VEC 2
    #define SIMD_LEVEL "SSE2"
    #define SIMD_TARGET_FUNC __attribute__((target("sse2")))
    
#else
    /* No SIMD support detected - use natural alignment */
    #define SIMD_ALIGN
    #define SIMD_WIDTH 8
    #define SIMD_DOUBLES_PER_VEC 1
    #define SIMD_LEVEL "SCALAR"
    #define SIMD_TARGET_FUNC
#endif

/* Fix for MSVC - it doesn't support target attributes */
#ifdef _MSC_VER
    #undef SIMD_TARGET_FUNC
    #define SIMD_TARGET_FUNC
#endif

#define TO_SIMD_PTR(ptr) ((simd_double_ptr_t)(ptr))

/* Compile-time information macros */
#define GET_SIMD_LEVEL() SIMD_LEVEL
#define GET_SIMD_WIDTH() SIMD_WIDTH
#define GET_SIMD_DOUBLES_PER_VEC() SIMD_DOUBLES_PER_VEC

#endif /* SIMD_DOUBLE_PTR_H */
