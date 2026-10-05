# 一维平行板算例

在 HyMoS 参数窗口中点击“加载参数”，选择本目录的模板，然后修改参数、检查并应用。

| 文件 | 壁面条件 | 工况 |
|---|---|---|
| `couette.txt` | 同温，切向壁速不同 | Couette |
| `fourier.txt` | 壁温不同，切向壁速相同 | Fourier |
| `couette_fourier.txt` | 壁温、切向壁速均不同 | Couette–Fourier |

三份模板均采用 `Nx=128`、`ORDER=7`、`FIM=3`、单网格和 1 个线程。完整模板中保留 BGK 模型所需的 `Pr=1` 和均匀网格设置 `Type=0`。界面按物理参数、网格阶数、求解方法和多重网格分类编辑。

批量运行示例：

```matlab
check = hymos_validate('couette', 'examples/ParallelPlate1D/fourier.txt');
task = hymos_submit('couette', 'examples/ParallelPlate1D/fourier.txt');
s = hymos_wait(task, 60);
r = hymos_result(task);
```

这里的 `couette` 是一维平行板求解入口标识，实际工况由输入文件中的壁温、壁速确定。
