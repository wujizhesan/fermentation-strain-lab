# Fermentation Strain Lab

公开数据驱动的菌株筛选与发酵建模研究原型，版本 0.1.1。

**主问题：能否利用菌株的早期观测，选出最终产物表现更好的菌株？**

本项目提供冻结数据、评价方案、计算程序、独立核验和一个统一入口。选株是主线；软测量和组学代谢建模使用其他数据，是两条独立研究线。当前只支持随包的固定公开案例，尚未提供任意自有数据训练接口。

## 四项计算与已有结果

| 模块 | 数据与评价 | 当前结论 |
| --- | --- | --- |
| `screening` | 酿酒酵母富马酸 ALE，9 株、3 个群体；18h OD，留一群体出 | 固定岭回归选中滴度 0.110767 g/L，随机期望 0.137973 g/L；岭回归低 19.72% |
| `softsensor` | Tamminen 2022 真菌发酵，20 个培养、2 轮；第 1 轮训练，第 2 轮测试 | CDW 当前信号 RMSE 3.240953 g/L，历史窗口 3.666398 g/L；历史未改善 |
| `metabolism` | 公开酵母 RNA/蛋白组，MARGE 原样与平衡修订各 12 情景 | 独立复核 LP 与通量区间；OptLang/GLPK 存在同名行覆盖，作者 CPLEX 路径未验证 |
| `external` | Yarrowia 3 株、大肠杆菌 4 株，24h→96h | 回顾性筛选规则回放；两个面板均未选中真实最优株，原 18h 岭回归尚无外部验证 |

外部案例中，Yarrowia 早期 OD 选中终点为面板最优的 98.07%，随机期望 99.27%；大肠杆菌早期产物规则为 98.29%，随机期望 75.09%。这些结果来自此前已查看的代表株面板，不是盲测，也未验证节省化验。

**证据边界：**没有同株数据贯通、前瞻选株或湿实验增产证据；严格论文复现仍缺五类输入；不能据此宣称模型优于论文或已实现菌株优化。通过软件测试表示实现符合当前评价方案。

## 下一阶段：前瞻验证准备

[ALE 选株模型前瞻验证方案](research/ale_prospective_validation_plan_20261015.md) 固定了外部测试、同株样本追踪与湿实验增产确认的门槛，并提供候选株、测量值及组学关联表头。该方案目前处于准备阶段：原公开终点滴度的具体采样时间尚未核实，因此严格外部验证尚未启动；没有新增预测、化验或增产结果。

## 无需安装依赖：查看参考结果

用 Python 3.13，在发行目录执行：

```text
python lab.py status
python lab.py report --output runs/reference
```

打开 `runs/reference/index.html`。此操作检查冻结参考报告、回执及展示文件，没有重训或重新求解。MARGE 随包参考只保留报告与展示所需的 34 个文件，不含大型 LP 快照；完整 LP 将在实际运行时生成。

## 安装与实际计算

Windows PowerShell：

```powershell
py -3.13 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-lab.txt
.\.venv\Scripts\python.exe lab.py prepare --download
.\.venv\Scripts\python.exe lab.py run --module screening --output runs/first
```

Linux/macOS：

```bash
python3.13 -m venv .venv
.venv/bin/python -m pip install -r requirements-lab.txt
.venv/bin/python lab.py prepare --download
.venv/bin/python lab.py run --module screening --output runs/first
```

将上面的 Python 命令替换成所建环境的实际路径，可执行：

```text
python lab.py run --module all --output runs/all
python lab.py verify --run-dir runs/all
python lab.py status --run-dir runs/all
```

也可分别指定 `softsensor`、`metabolism`、`external`。每次输出必须是尚不存在的新目录。成功运行后，`summary.json` 汇总本次指标，`lab_receipt.json` 绑定执行日志、模块回执与总览，`index.html` 可直接打开。MARGE 的完整 LP 与独立求解回执在 `core/metabolism/`。

`prepare --download` 从作者仓库、Zenodo、Ensembl、Europe PMC 和 PMC 公共存储获取原始输入，并核对冻结 SHA。MARGE 补材按 27 个成员核对，允许服务重新打包时间戳。公开服务若改变内容，命令会拒绝写入；此时应使用同版离线缓存，不能替换冻结输入后继续沿用原评价声明。选株快速运行只需要随包派生输入，可以跳过数据下载。

离线包已包含计算所需原始输入，可用 `python lab.py prepare` 检查。离线包中的第三方数据与模型保留原许可；GitHub 源码版通过下载入口获取它们。来源、许可记录、路径、SHA 和是否随包见 `RELEASE_SOURCES.json` 与 `THIRD_PARTY.md`。

## 目录与评价输入

```text
lab.py                  统一入口：status / report / prepare / run / verify
research/*protocol*     冻结评价方案与输入摘要
research/results/       参考结果和派生输入
data/raw/               公开原始输入缓存
tests/                  选株、软测量、外部案例与入口测试
gem_tests/              MARGE 来源、求解与展示测试
runs/                   用户新建的运行结果，不进入 Git
```

选株读取 `ale_chassis_audit/isolate_observations.csv` 和 `isolate_endpoints.csv`，使用菌株标识与群体标识对应早期观测和终点。这些 CSV 是从论文补材整理并已独立核验的冻结输入；本发行包不重新执行整个补材整理流程。新增菌株或修改评价时，需要另建协议与输入适配，不能直接改写这里的冻结文件。

软测量在原始化验点计分，不将插值点当作独立标签。MARGE 输入包括模型、基因映射、组学比值及来源审计；它预测允许的通量区间，不等同于实测代谢流量。外部案例独立读取作者数值工作簿，保留候选株、重复和排除依据。

## 测试与持续集成

完成全量运行后设置 `MARGE_RESULTS_ROOT` 为 `runs/all/core/metabolism` 的绝对路径，再执行：

```text
python -m unittest discover -s tests -v
python -m unittest discover -s gem_tests -v
```

Windows 可用 `$env:MARGE_RESULTS_ROOT=(Resolve-Path runs/all/core/metabolism).Path`，Linux 可用 `export MARGE_RESULTS_ROOT="$PWD/runs/all/core/metabolism"`。MARGE 的运行入口本身也会在新生成结果上执行 `gem_tests`。

GitHub Actions 配置 Windows/Ubuntu 的参考检查、快速选株及入口测试；手动触发并选择 `full_replay` 可下载输入并全量重跑。工作流是否在 GitHub 通过，以实际 Actions 记录为准。本地复算验收记录见 [VALIDATION.json](VALIDATION.json)。

依赖锁定在 `requirements-lab.txt` 及其包含的文件中。发行环境采用 Python 3.13；早期软测量归档记录的 Python/NumPy 版本不同，因此保留归档原件，并对本次结果重新核验，不声称原环境逐位一致。

## 来源与使用

主选株和 MARGE 的案例来自 [酵母 ALE 论文](https://doi.org/10.15252/msb.202110253)，软测量来自 [Tamminen 2022](https://doi.org/10.1186/s13568-022-01436-4) 及 [Zenodo 数据](https://zenodo.org/records/5846028)。外部案例见 `THIRD_PARTY.md`。请保留作者与来源归属。

本发行包没有为第三方论文、数据、模型或依赖重新授予许可证。本项目原创代码尚未指定对外开源许可证，公开仓库展示不等于授予再分发许可。
