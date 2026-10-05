#ifndef __config_h_
#define __config_h_

#define NRXX_NAMESPACE_OPEN namespace NRXX {
#define NRXX_NAMESPACE_CLOSE }

#ifndef MAX_VECTOR_DIMENSION
#define MAX_VECTOR_DIMENSION 3
#endif // MAX_VECTOR_DIMENSION

#ifndef MAX_RKC_STEP
#define MAX_RKC_STEP 100
#endif // MAX_RKC_STEP

#ifndef MAX_ORDER
#define MAX_ORDER 50
#endif // MAX_ORDER

// 默认使用线性化
#if (!defined(LINEARIZATION) && !defined(NOLINEARIZATION))
#define LINEARIZATION
#endif

// 若要使用多线程，则需取消注释以下行以防投影时产生竞争条件
// #define MULTITHREAD_NRXX

// 若要谱方法做投影，则取消注释以下行
// #define SPECTRAL_PROJECTION

#endif // config.h
