/**
 * @file   config.h
 * @author HU Zhicheng <huzhicheng1986@gmail.com>
 * @date   Tue Sep 16 11:22:04 2014
 * 
 * @brief  一些共用的宏设置
 * 
 * SteadyState/ 目录下还是有几个 .cpp 文件会事先编译成 .o 再与
 * examples/ 的算例链接成最终的程序, 因而在 examples/ 下的各 cpp 中再指
 * 定具体的宏似乎并不合适. 所以建这个文件, 将用于编译时判别的宏定义等放
 * 到这里, 方便各例子的具体设置.
 */


#ifndef __SSP_CONFIG_H__
#define __SSP_CONFIG_H__


/**
 * 或者把这 V, W 循环放到模板中去? V 循环是在粗网格校正时递归调用一次多
 * 重网格算法. W 循环是递归调用两次. 所谓的 F 循环, 是把 W 循环粗网格校
 * 正时第二次调用自身的 W 循环替换为 V 循环, 也就是说粗网格校正时先调用
 * 一次自己 (F 循环), 再调用一次 V 循环, 据说这样兼顾两者的优点. F 循环,
 * 看起来有点像是每一步粗化到最粗网格, 然后使用 V cycle 的 full
 * multigrid.
 * 
 * 发现 F 循环在粗网格校正时不同的调用方式, 效果各有优势, 于是再细化一
 * 下. 其中 FVCYCLE 就是指上面提到的 F 循环, 即粗网格求解时第一次调用自
 * 身, 第二次 V 循环. 而 VFCYCLE 自然是反一反, 即第一次调用 V 循环, 第
 * 二次调用自身 (初步测试比 FVCYCLE 略快). VVCYCLE 是 WCYCLE 的简化版,
 * 粗网格求解两次都调用 V 循环, 其示意图更符合 W 形状, 不过通常讲得 W
 * 循环是 WCYCLE 那种形式的, 所以这里就用 VVCYCLE 表示这种简化版了. 最
 * 后一个可以看作是 FVCYCLE 和 VFCYCLE 的结合版本 (对称), 粗网格求解调
 * 用三次, 中间一次是调用自己, 前后则是 V 循环.
 * 
 */
enum CycleType { VCYCLE=1, WCYCLE, FVCYCLE, VFCYCLE, VVCYCLE, VFFVCYCLE, VFVVCYCLE, VVFVCYCLE };
/**
 * 多重矩方法的降阶策略也用一个枚举吧, 这样代码可能可以写得略有效率些.
 * 三个枚举值分别表示减半(下整), 减一, 减二 
 * 
 */
enum OrderReductionStrategy { ORSHALF=0, ORS1, ORS2 };
 
/**
 * g++ 不允许在未特化类模板时特化其成员模板函数, 在用重载函数方式实现类
 * 似特征时, 为了区别不同枚举值的重载函数, 引入如下的一个空结构体, 于是
 * 可以用 CycleTypeIdentity<VCYCLE> 这样的结构体作为函数参数类型. 这样
 * 重载的函数参数都有相同的特征, 可以再写一个模板函数作为统一接口调用这
 * 些重载函数.
 *
 * 简言之, 这个结构体可以使 CycleType 中的每个值都可以当做数据类型类使
 * 用.
 * 
 */
template <CycleType _cycle> struct CycleTypeIdentity { };
template <OrderReductionStrategy _ors> struct ORSIdentity { };

/**
 * 引入这个空的结构体与 CycleTypeIdentity<CycleType> 目的类似, 还是为了
 * 用重载函数的方式替代解决 "g++ 在未特化类模板时不允许特化单个成员函数
 * " 的问题.
 * 
 */
template <int _flag> struct ArtificialTypeIdentity { };

// 默认使用单层网格求解器
#if (!defined(MULTI_GRID_SOLVER) && !defined(SINGLE_GRID_SOLVER))
#define SINGLE_GRID_SOLVER
#endif

/**
 * 在计算半导体器件时, 以下两个宏关于材料的宏定义需要开启其中一个. 当然
 * 目前程序能算的还只是 Si 器件, 也就是 DEVICE_SILICON 宏.
 * 
 */
//#define DEVICE_GaAs
#define DEVICE_SILICON

/// 是否做重构
//#define RECONSTRUCT

#endif // __SSP_CONFIG_H__

/**
 * end of file
 * 
 */
