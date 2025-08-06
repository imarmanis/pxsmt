#ifndef __VERIFIER_H__
#define __VERIFIER_H__

#include <stdbool.h>
#include <stddef.h>
#include <stdatomic.h>
#include <immintrin.h>

bool __VERIFIER_assume(bool);
#define assume(b) __VERIFIER_assume(b)

bool __VERIFIER_assert(bool);
#define assert(b) __VERIFIER_assert(b)

#define faa(x, u) __atomic_fetch_add(x, u, 0)

#define cas(x, u, w) ({    \
    typeof(u) __dummy = u; \
    __atomic_compare_exchange_n(x, &__dummy, w, false, 2, 2); \
})
/* false : strong */

int __VERIFIER_nondet_int();

void __VERIFIER_parallel(void (*)(void), ...);
#define par(...) __VERIFIER_parallel(__VA_ARGS__)

#define flush(x) ({         \
            _mm_clflush(x); \
        })

#define flush_opt(x) ({     \
            _mm_clflushopt(x); \
        })

#define mfence() _mm_mfence();

#define sfence() _mm_sfence();


#define rec __VERIFIER_recovery

#define pthread_t int
#define pthread_mutex_t int
#define pthread_mutex_init(l) { *l = 0; } ((void) 0)
#define pthread_mutex_lock(l) { assume(cas(l, 0, 1)); } ((void) 0)
#define pthread_mutex_unlock(l) { *l = 0; } ((void) 0)
#define pthread_attr_t int

int __VERIFIER_thread_create(pthread_t *, pthread_attr_t *, void * (*)(void *), void *);

int __VERIFIER_thread_join(pthread_t, void **);

#define pthread_exit(...) ((void) 0)

pthread_mutex_t __global_lock;
#define __VERIFIER_atomic_init() { __global_lock = 0 ;} ((void) 0)
#define __VERIFIER_atomic_begin() { pthread_mutex_lock(&__global_lock); } ((void) 0)
#define __VERIFIER_atomic_end() { pthread_mutex_unlock(&__global_lock); } ((void) 0)

#endif
