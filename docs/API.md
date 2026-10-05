# HyMoS-Baltamatica 接口参考

插件在北太天元中通过 `load_plugin('HyMoS-Baltamatica')` 加载。日常操作见[技术文档](HyMoS-Baltamatica-technical-guide.pdf)，本页用于脚本调用及开发查阅。配置和返回数组均使用无量纲量。

## 图形界面入口

| 调用 | 输入 | 输出及行为 |
|---|---|---|
| `info = hymos_info()` | 无 | 版本、支持的算例标识及能力说明结构体；可省略输出 |
| `hymos_setup('1D')` | 字符串 `1D`（也接受 `1d`、`couette`） | 打开参数窗口，无返回值；应用后将配置文本存为宿主工作区变量 `hymos_config_text` |
| `hymos_run()` | 无 | 提交已应用配置，将任务存为 `hymos_task`，并启动前台监视；无返回值 |

## 文件输入

| 调用 | 输入 | 输出 |
|---|---|---|
| `v = hymos_validate(case_id, path)` | 算例标识、完整 INI 文件的绝对或相对路径 | `valid`（逻辑值）、`case`、`message` |
| `task = hymos_submit(case_id, path)` | 算例标识、完整 INI 文件的绝对或相对路径 | 任务对象；提交后后台执行 |

`case_id` 使用 `couette`。Couette、Fourier、Couette–Fourier 共用该平行板求解入口，通过两侧壁温及切向壁速区分。完整输入以 `examples/ParallelPlate1D/` 的三个模板为准，参数说明见技术文档。

```matlab
v = hymos_validate('couette', 'examples/ParallelPlate1D/fourier.txt');
if ~v.valid
    error(char(v.message));
end
task = hymos_submit('couette', 'examples/ParallelPlate1D/fourier.txt');
s = hymos_wait(task, 60);
r = hymos_result(task);
files = hymos_export(task);
```

## 任务控制

| 调用 | 输入 | 输出及行为 |
|---|---|---|
| `s = hymos_status(task)` | 任务对象 | 当前状态结构体 |
| `s = hymos_wait(task, seconds)` | 任务对象、0–3600 秒实数 | 等待结束或超时，返回当时状态；0 为立即查询 |
| `s = hymos_stop(task)` | 任务对象 | 请求协作停止，返回当时状态；可省略输出 |
| `s = hymos_monitor(task, refresh_seconds)` | 任务对象、可选刷新间隔（0.05–60 秒，默认 0.5） | 前台打印进度，结束后返回状态；可省略输出 |

`hymos_monitor` 由随包 M 脚本提供，调用 `hymos_setup` 或 `hymos_run` 后自动加入搜索路径。直接编写脚本时，也可先将插件的 `ui/` 目录加入路径。

停止请求在迭代检查点生效，调用后可继续 `hymos_wait` 确认最终状态。监视期间按 Ctrl+C 可返回命令窗口，后台计算继续；再次调用 `hymos_monitor(task)` 可继续查看。

状态结构体字段：`state`、`iteration`、`residual`、`elapsed_seconds`、`message`、`workspace`、`resolved_config`、`tolerance`。

| `state` | 含义 |
|---|---|
| `pending`、`running` | 等待执行、正在计算 |
| `converged` | 残差达到容差 |
| `iteration_limit` | 达到迭代步数上限 |
| `stopped` | 已响应停止请求 |
| `invalid_config` | 配置校验失败，原因见 `message` |
| `busy` | 求解器已有活动任务 |
| `failed` | 求解过程中发生错误，原因见 `message` |

插件每次执行一个计算任务。任务对象复制后仍引用同一任务；释放最后一个引用或卸载插件时，会请求停止并等待线程退出。

## 结果读取

`r = hymos_result(task)` 返回已结束任务的结果；任务运行中会提示先等待。`invalid_config`、`busy`、`failed` 状态返回明确错误。

| 字段 | 尺寸 | 含义 |
|---|---|---|
| `x` | Nx × 1 | 网格单元中心坐标 |
| `density` | Nx × 1 | 密度 |
| `velocity` | Nx × 3 | 速度分量 u1、u2、u3；u2 为壁面切向 |
| `temperature` | Nx × 1 | 温度变量 θ |
| `stress_moments` | Nx × 3 | 二阶 Hermite 系数 f₂e₁、fₑ₁₊ₑ₂、f₂e₂ |
| `heat_flux` | Nx × 2 | 热流 q1、q2 |
| `residual` | K × 1 | 残差历史，K 为记录数 |
| `metadata` | 结构体 | 算例、状态、网格、矩阶数、线程数、Kn、容差、FIM、耗时及文件路径等 |

`metadata` 中的 `velocity_components`、`stress_definition`、`heat_flux_components`、`field_columns`、`nondimensional_definition` 给出数组分量及无量纲约定。

`files = hymos_export(task)` 导出结果并返回以下路径（可省略输出）：

| 字段 | 文件或目录 |
|---|---|
| `workspace` | 本次任务目录 |
| `fields` | `couette_fields.csv`，三类工况共用此文件名 |
| `residual` | `residual.csv` |
| `metadata` | `result_metadata.txt` |
| `distribution` | `DisSol1th.dat`，分布展开数据 |

任务目录还保存 `input.resolved.txt` 输入快照。默认根目录为 `/tmp/hymos_baltamatica_runs`；可在启动北太天元前通过环境变量 `HYMOS_RUNS_ROOT` 指定可写目录。

## 界面桥接接口

以下接口供随包界面调用，输入和输出规则也可用于自动化配置。

| 调用 | 输入及返回 |
|---|---|
| `v = hymos_validate_text(case_id, text)` | 校验完整 INI 文本，返回与 `hymos_validate` 相同的结构体 |
| `task = hymos_submit_text(case_id, text)` | 从完整 INI 文本提交任务，返回任务对象 |
| `value = hymos_ui_config(text, key)` | 从配置读取字符串值，`key` 格式为 `Section/Key`；缺失键返回空字符串 |
| `updated = hymos_ui_config(text, key, value)` | 修改或补入指定键，返回完整文本；重复键修改最后一项，值为单行字符串 |
| `path = hymos_ui_file(path)` | 返回规范化绝对路径；父目录须已存在，目标位于已加载插件目录之外 |
| `path = hymos_ui_file(path, text)` | 校验并原子保存完整配置，返回保存路径；路径规则同上 |

参数个数、类型、范围或任务状态不满足调用条件时，插件通过宿主错误机制给出原因；校验接口将配置问题放入 `valid` 和 `message`。
