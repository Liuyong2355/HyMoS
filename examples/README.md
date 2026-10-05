# 算例

`ParallelPlate1D/` 包含插件使用的 Couette、Fourier、Couette–Fourier 配置。在北太天元参数窗口点击“加载参数”，选择对应文件即可，详见[平行板算例说明](ParallelPlate1D/README_CN.md)。

`Couette/`、`ShockStructure/`、`Cavity/` 保留核心库的命令行模板。运行 `bash scripts/build.sh` 后，命令行程序位于 `build-release/bin/hymos`。

在仓库根目录执行单次 Couette 计算：

```bash
bash examples/Couette/hymos_cli run couette examples/Couette/input.txt
```

算例启动脚本将结果写入仓库根目录的 `workspace_runs/<Case>_YYYYMMDD_HHMMSS/`。

交互运行时，先在一个终端启动服务：

```bash
cd examples/Couette
bash hymos_cli serve couette
```

随后在另一个终端进入同一目录，通过生成的控制脚本操作：

```bash
bash hymos_control.sh config input.txt
bash hymos_control.sh run
```
