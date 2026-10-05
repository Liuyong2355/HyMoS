# 第三方组件说明

HyMoS 使用以下第三方组件，各组件保留原作者版权及各自的许可证。表中版本为本项目已验证的构建环境版本。

## NRxx 代码来源

感谢新加坡国立大学蔡振宁教授提供本项目使用的 NRxx 代码。相关源文件保留原有版权和作者署名。

## 随源码保留的组件

| 组件 | 用途及来源 | 许可证 |
| --- | --- | --- |
| Chameleon，Copyright © 2002–2004 René Nyffenegger | `src/common/chameleon.h`、`chameleon.cpp` 提供配置值的字符串和数值转换；在原实现上增加整数转换并调整接口。原始出处：[作者仓库](https://github.com/ReneNyffenegger/cpp-read-configuration-files)。两个文件均保留原版权声明并标记修改。 | [zlib 许可原文](licenses/CHAMELEON-ZLIB.txt) |

## 构建依赖

这些库通过系统开发包提供，未将其完整源码复制到本仓库。

| 组件 | 验证版本 | 用途 | 许可证及来源 |
| --- | --- | --- | --- |
| Boost | 1.83.0 | C++ 工具、序列化及线程相关头文件；构建配置链接 serialization、system、thread。 | [Boost Software License 1.0](licenses/BOOST-1.0.txt)；[上游](https://www.boost.org/LICENSE_1_0.txt) |
| Eigen | 3.4.0 | 矩阵及线性代数运算。 | [MPL 2.0](licenses/MPL-2.0.txt)；[文件级版权及许可清单](licenses/EIGEN-COPYRIGHT.txt) |
| FFTW | 3.3.10 | 矩方法中的谱投影相关代码及构建链接项。 | GPL 2.0 或更高版本；[版权声明](licenses/FFTW-COPYRIGHT.txt)、[GPL 2.0](licenses/GPL-2.0.txt)、[上游说明](https://www.fftw.org/faq/section1.html) |
| oneTBB | 2021.11.0 | 构建链接项。 | [Apache 2.0](licenses/ONETBB-APACHE-2.0.txt)；[上游对应版本](https://github.com/uxlfoundation/oneTBB/tree/v2021.11.0) |
| Open MPI | 4.1.6 | MPI 开发头文件及构建链接项。 | [上游版权及 BSD 类许可声明](licenses/OPENMPI-LICENSE.txt)；[上游对应版本](https://github.com/open-mpi/ompi/tree/v4.1.6) |

Eigen 的头文件参与编译，HyMoS 未修改这些文件。对应源码可从 [Eigen 3.4.0](https://gitlab.com/libeigen/eigen/-/tree/3.4.0) 获取；各文件适用的 MPL、可选双重许可及辅助文件 BSD 声明见上述版权清单。

## 运行环境

北太天元插件包包含插件入口和 HyMoS 核心库。北太天元 SDK 及 `libbex.so` 由用户安装的北太天元提供，其使用遵循厂商许可。

Linux 系统提供 C/C++ 数学及运行库（`libc`、`libm`、`libstdc++`、`libgcc_s`）和 OpenMP 运行库（`libgomp`）。其中 GCC 运行库适用 GPL 3.0 或更高版本及 GCC Runtime Library Exception 3.1，相关文本见 [GPL 3.0](licenses/GPL-3.0.txt) 和 [运行库例外](licenses/GCC-RUNTIME-EXCEPTION-3.1.txt)；glibc 适用 LGPL 2.1 或更高版本。上述系统运行库不随插件包分发。
