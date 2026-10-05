# 构建和发布

目标平台为 Ubuntu 24.04 x86_64，宿主为北太天元 2025 Linux 版。预编译包用户直接按根目录 README 安装。

## 开发依赖

安装北太天元及其 SDK；默认位置为 `/opt/Baltamatica`。项目使用 C++11、CMake、OpenMP，并查找 Eigen、Boost、Open MPI、FFTW、oneTBB 开发包。

```bash
sudo apt install build-essential cmake python3 libopenmpi-dev \
  libeigen3-dev libboost-serialization-dev libboost-system-dev \
  libboost-thread-dev libfftw3-dev libtbb-dev
bash scripts/build.sh --prefix /opt/Baltamatica --jobs 4
```

脚本优先使用环境变量 `CXX`，其次为 `g++-9`，再使用 `g++`；发布版本使用 GNU C++ 9.5.0。其他编译器可通过 `--cxx /path/to/g++` 指定。构建完成后自动运行 CTest，详情见[测试说明](../tests/README.md)。

## 产物

- `build-release/lib/libhymos.so.1.0.0`：HyMoS 计算共享库，ABI 名为 `libhymos.so.5`。
- `build-release/plugin/HyMoS/`：完整插件，入口遵循北太天元 SDK 的 `main.so` 命名。
- `build-release/bin/hymos`：命令行入口。

## 打包

```bash
bash scripts/package.sh --source
```

在 `dist/` 中生成二进制包、源码包及 SHA-256 校验文件。二进制包顶层目录为 `HyMoS/`，包含共享库、界面脚本、输入模板、示例、使用说明和许可证。SDK 库由北太天元提供，C/C++ 及 OpenMP 运行库由系统提供。实际依赖清单保存在包内 `RUNTIME_DEPENDENCIES.txt`。

## 安装

可直接将解压后的完整 `HyMoS/` 复制到北太天元的 `plugins/` 目录，也可使用随包脚本：

```bash
bash HyMoS/scripts/install.sh --prefix /opt/Baltamatica --dry-run
bash HyMoS/scripts/install.sh --prefix /opt/Baltamatica
```

从源码工作区安装时使用 `bash scripts/install.sh --prefix /opt/Baltamatica`。若通过 `--build-dir` 使用了其他构建目录，安装时以 `--bundle PATH/plugin/HyMoS` 指定对应插件包。`--plugin-dir PATH` 指定宿主插件父目录，脚本将在其中安装 `HyMoS/`。更新本工具已安装的版本时使用 `--replace`，原目录会保留为备份。安装到系统目录时需相应写入权限。

加载命令统一为 `load_plugin('HyMoS')`。更换已加载的插件文件后重启北太天元。

## 技术文档编译

文档源文件位于 `docs/technical/`，使用 XeLaTeX、ctex、TikZ 和 gbt7714。

```bash
cd docs/technical
latexmk -xelatex HyMoS.tex
```
