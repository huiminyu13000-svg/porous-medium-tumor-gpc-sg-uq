# gPC：concave 完整实验

本包保留 concave 实验，对应论文 Test I(a)。原始代码为 `test1111.m` 的 concave 分支，原图数据为 `1127concave_max_vivo.mat`。

## 运行

在 MATLAB 中将当前目录设为本文件所在目录：

```matlab
run_experiment                  % 从初值计算 SG 和全部 SC 轨迹，三个时刻
run_experiment('cached_reference') % 计算 SG，复用原 SC 参考
run_experiment('reference_short')  % SG 和 SC 只计算第一个时刻
plot_saved                      % 直接绘制原保存数据
```

全计算会重新做初值投影、时间积分、均值/标准差计算和切片绘图；输出在
`results/`。原数据不被覆盖。直接画图只使用这一例的数据。

## 论文实验：Test I(a)

N. Jiang, L. Liu, H. Yu, *Efficient stochastic asymptotic-preserving scheme for tumor growth models with uncertain parameters*, JCAM 473 (2026), 116881，§4.1，式 (29)–(30)，图 1–4。

一维随机输入 z∼U[-1,1]；空间域 [-2.2,2.2]^2，dx=dy=0.1，dt=0.001，m=80。
论文采用 K=4 的 gPC-SG，并以 16 点 Gauss–Legendre stochastic collocation 为参考。
初始密度的 core 和 rim 分别为 0.15−0.1z 与 0.35+0.1z；
四个半径为 (0.4,0.5,0.6,0.7)(1+0.3z)，分段定义见式 (29)。
营养增长系数 G0=0.5(1−0.1z)。

图 1 比较不同时间的均值和标准差切片；图 2–3 展示 T=1 的统计场；
图 4 比较 T=0.5 的随机阶数误差。本文代码包只选取 concave 统计重构这一实验，
不包含整篇文章的全部参数扫描。

## 本包运行

`run_experiment` 使用包内初始样本和参考字段，`plot_saved` 使用已保存的原图数据。
实际运行输入记录在 `data/concave_compute.mat`；调用方式见上方“运行”。
`run_smoke` 在 T=0.01 执行 SG 与全部 16 个 SC 节点，再检查直接画图。
保存数据与此次短时计算分别输出，不用短时结果替代论文图。

## 文件与依赖

| 文件 | 用途 | 原名 |
|---|---|---|
| `run_experiment.m` | 初值投影、SG/SC 计算、统计、计算图 | 整理后的主入口 |
| `solver/solve_gpc_sg.m` | 耦合 gPC-SG 求解 | solver_2D_New_sg_couple |
| `solver/solve_sc.m` | 单节点 SC 求解 | solver_2D_collocation_couple |
| `solver/gpc_growth.m` | SG 增长项投影 | calculate_rhoGGG |
| `solver/sc_growth.m` | SC 增长项 | calculate_collocation_rhoGGG |
| `solver/lgwt.m` | Gauss–Legendre 节点和权重 | lgwt |
| `plot_saved.m` | 原图切片、场和已保存误差 | test1111 的画图部分 |
| `data/concave_compute.mat` | 初值、节点、基函数、参考统计 | 原 MAT 数值字段 |
| `data/gpc_concave_paper.mat` | 原始画图数据 | 1127concave_max_vivo.mat |

上述基础函数全部保留，调用关系为：主入口 → SG/SC solver → 增长项；
SG 增长项调用 lgwt。通量、矩阵乘法等局部函数仍在相应 solver 内。
本次只调整文件组织、命名、入口默认方式和部分中文说明，不新增数学模型。

## 论文结果与运行检查

论文报告均值、标准差与 SC 参考吻合，并展示随机阶数增加时误差下降。
本次只进行短时运行与画图检查；执行范围和结果见 `TESTS.md`。

GitHub 仓库根目录放置本 README、入口文件、`solver/` 和 `data/`。
不需要展平子目录；入口按自身所在位置定位数据和函数。`results/` 自动生成并已在 `.gitignore` 中排除。
