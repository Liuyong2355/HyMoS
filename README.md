# HyMoS-Baltamatica

**HyMoS（Hyperbolic Moment Solver）的北太天元插件**

HyMoS-Baltamatica 将 HyMoS 接入北太天元，采用快速迭代矩方法（FIM）求解一维稳态 Boltzmann–BGK 方程，支持 Couette、Fourier 和 Couette–Fourier 平行板流动。用户在北太天元图形界面中设置参数，调用插件完成求解、结果分析和数据导出。

本仓库面向 B3“自选开源库插件开发”。参赛说明见 [材料清单](submission/README.md)。

## 功能

- 配置壁面温度、切向速度、Knudsen 数、网格及矩阶数。
- 选择 FIM-1/2/3、单网格或非线性多重网格（NMG），设置 V/W/F 循环。
- 后台执行求解，查询状态、监视残差、安全停止任务。
- 返回密度、速度、温度、应力矩、热流和残差数组，供北太天元继续绘图及运算。
- 导出场量、残差历史、运行信息及分布展开数据，保存实际输入快照。

## 安装

已验证平台：**Ubuntu 24.04 x86_64（含 WSL2）和北太天元 2025 Linux 版**。

从 [Releases](https://github.com/Liuyong2355/HyMoS-Baltamatica/releases) 下载 `HyMoS-Baltamatica-1.0.2-linux-x86_64.tar.gz`，解压后将整个 `HyMoS-Baltamatica` 目录放入北太天元的 `plugins` 目录。默认位置为 `/opt/Baltamatica/plugins/HyMoS-Baltamatica`。

```text
plugins/HyMoS-Baltamatica/
├── main.so
├── lib/
├── ui/
├── assets/
└── examples/ParallelPlate1D/
```

包内包含 HyMoS 共享库、界面脚本、模板和算例；使用预编译包无需安装 Eigen、Boost、MPI、FFTW、oneTBB 开发包。北太天元提供宿主库，系统提供 C/C++ 及 OpenMP 运行库（Ubuntu 包 `libgomp1`）。

## 使用

在北太天元命令窗口执行：

```matlab
load_plugin('HyMoS-Baltamatica');
hymos_setup('1D');
```

1. 在“物理参数”页设置 `Kn`、两侧壁温和切向壁速；界面自动识别工况。
2. 在“网格与阶数”页设置 `Nx` 和 `ORDER`。均匀网格取 `Nx >= 8`；演示采用 `Nx = 128`、`ORDER = 7`。
3. 在“求解方法”页选择 FIM 和求解框架，计算线程数建议设为 1–4。使用 NMG 时继续设置循环类型和平滑参数。
4. 点击“检查参数”，然后“应用并关闭”。参数文件可从随包模板修改。
5. 运行、绘图并导出：

```matlab
hymos_run();
r = hymos_result(hymos_task);
plot(r.x, r.temperature);
files = hymos_export(hymos_task);
```

任务在后台求解；`hymos_run()` 同时显示前台监视。`Ctrl+C` 可结束前台监视，随后通过 `hymos_monitor(hymos_task)` 继续查看。`hymos_stop(hymos_task)` 请求在安全迭代边界停止计算。

详细说明：[技术文档](docs/HyMoS-Baltamatica-technical-guide.pdf) · [接口参考](docs/API.md) · [算例](examples/ParallelPlate1D/README_CN.md) · [演示文稿](submission/slides/HyMoS-Baltamatica-slides.pdf) · [演示视频](submission/HyMoS-Baltamatica-demo.mp4)。视频含作者配音和中文字幕，时长约 2 分 32 秒；另附 [SRT 字幕](submission/HyMoS-Baltamatica-demo.zh-CN.srt)。

## 源码构建

源码构建用于修改、复现和开发测试。需要北太天元 SDK、C++11 编译器、CMake，以及工程使用的开发依赖。发布包采用 GCC 9.5.0 构建。

```bash
sudo apt install build-essential cmake python3 libopenmpi-dev \
  libboost-serialization-dev libboost-system-dev libboost-thread-dev \
  libeigen3-dev libfftw3-dev libtbb-dev
bash scripts/build.sh
```

SDK 默认路径 `/opt/Baltamatica`。其他路径可通过脚本的 `--prefix` 选项指定，完整选项见 `bash scripts/build.sh --help`。

构建产物位于 `build-release/plugin/HyMoS-Baltamatica/`；发布和安装命令见 [构建说明](docs/BUILD.md)。

## 验证

```bash
ctest --test-dir build-release --output-on-failure
```

在北太天元中，以本仓库为当前目录运行宿主集成检查：

```matlab
load_plugin('HyMoS-Baltamatica');
addpath('tests/baltamatica');
hymos_smoke_tests('examples/ParallelPlate1D');
```

检查覆盖三类平行板算例的收敛、场量、导出，以及参数校验和任务控制。检查项目和测试入口见 [测试说明](tests/README.md)。

## 仓库结构

| 目录 | 内容 |
|---|---|
| `plugins/baltamatica/` | C++ SDK 接口、原生图形界面、默认模板 |
| `src/`、`NRxx/` | HyMoS 求解器、矩表示及任务管理 |
| `cases/`、`apps/` | 核心库构建使用的算例定义和 CLI |
| `examples/ParallelPlate1D/` | 三类插件算例 |
| `tests/` | 库回归基准、配置读取测试、宿主检查 |
| `docs/` | 技术文档、接口参考、构建说明 |
| `scripts/` | 构建、打包、安装工具 |
| `submission/` | 演示文稿、讲解视频、参赛说明 |

## 致谢

感谢新加坡国立大学蔡振宁教授提供 NRxx 代码。联系邮箱：[matcz@nus.edu.cn](mailto:matcz@nus.edu.cn)。

## 联系信息

南京航空航天大学：

| 姓名 | 邮箱 |
|---|---|
| 胡志成 | [huzhicheng@nuaa.edu.cn](mailto:huzhicheng@nuaa.edu.cn) |
| 刘盛琦 | [liushengqi@nuaa.edu.cn](mailto:liushengqi@nuaa.edu.cn) |
| 疏凌云 | [cieloudsly@nuaa.edu.cn](mailto:cieloudsly@nuaa.edu.cn) |
| 陈升光 | [sunnychen@nuaa.edu.cn](mailto:sunnychen@nuaa.edu.cn) |

## 许可

HyMoS-Baltamatica 采用 [MIT 许可证](LICENSE)，保留源码中的原作者版权及署名。第三方组件遵循各自的许可证，来源和许可见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。AI 辅助范围见 [AI_USAGE.md](AI_USAGE.md)。
