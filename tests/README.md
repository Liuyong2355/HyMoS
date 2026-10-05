# 测试说明

## 自动回归

`bash scripts/build.sh` 构建后自动运行两项 CTest，也可单独执行：

```bash
ctest --test-dir build-release --output-on-failure
```

- `config_file_contract`：INI 节和键、空白和注释、重复键、值中的等号、默认值、缺失项及格式错误提示。
- `couette_memory_and_cli`：内存结果和命令行输出同基准数据的一致性；FIM-1/2/3、NMG、输入快照、异常配置、迭代上限、异步结果及停止行为。历史数值基准保存在 `regression/nodealii/couette/runs/`。

测试入口为 `library/config_file.cpp`、`library/couette_memory.cpp` 和 `library/check_couette.py`。

## 北太天元集成检查

安装插件后，在北太天元中以仓库为当前目录运行：

```matlab
load_plugin('HyMoS-Baltamatica');
addpath('tests/baltamatica');
hymos_smoke_tests('examples/ParallelPlate1D');
```

脚本依次运行 Couette、Fourier、Couette–Fourier，检查收敛、结果尺寸、有限场量、密度和温度为正、壁速特征及导出文件；另验证最小网格 `Nx=8`，以及不合法网格和缺失输入文件的错误提示。模板保持原样，每个任务独立保存输入快照和结果。

成功时输出 `HYMOS_HOST_SMOKE_PASS`，失败时给出具体检查项。

## 界面检查

加载 `HyMoS-Baltamatica`，运行 `hymos_setup('1D')`，依次确认参数分类、模板加载、检查参数和应用；随后执行 `hymos_run()`，读取结果并绘图、导出。讲解视频展示这套流程。
