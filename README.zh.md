<p align="center"><img src="docs/logo/logo.png" alt="CP2Ksol" width="360"></p>

# CP2Ksol

[English](README.md) · [中文](README.zh.md)

开发者：[Leshen Zhang](https://github.com/leshenzhang)（UCLA Alexandrova 课题组）。

为 **CP2K 2026.2**（Quickstep, GPW）加入隐式溶剂、隐式电解质和恒电势（巨正则）电极。支持单点、几何优化和分子动力学，晶胞里可以有显式水，也可以没有。

代码是 CP2K 源码树里的纯 Fortran，不依赖外部库和内嵌 Python，不改变 MPI 的生命周期。

目录：[功能](#功能) · [安装](#安装) · [快速上手](#快速上手) · [模型怎么选](#模型怎么选) · [关键字](#关键字) · [推荐设置](#推荐设置) · [性能与精度](#性能与精度) · [输出](#输出) · [测试](#测试) · [已知限制](#已知限制) · [引用与许可](#引用与许可)

## 功能

- **隐式溶剂。** 两类空腔。`METHOD ISODENSITY` 把一个空腔放在电子密度的等值面上（https://doi.org/10.1063/1.4865107）。`METHOD SOLVENT_EXCLUDED` 用探针壳卷积构造分层的最近接近表面（溶剂、介电、离子、空化）（https://doi.org/10.1063/5.0176308）。
- **隐式电解质。** 线性化 Poisson–Boltzmann 屏蔽，G = 0 处有电势参考，所以体相电解质电势为零（https://doi.org/10.1063/1.5132354）；带离子尺寸的非线性 Poisson–Boltzmann 离子（https://doi.org/10.1103/physrevlett.79.435）；带偶极饱和的非线性介电（https://doi.org/10.1088/0965-0393/21/7/074005）。
- **显式水加隐式溶剂。** `&REGION` 把隐式溶剂限制在一个平板区域；`&REGION &FLEXIBLE` 让区域边界跟着显式水膜走；`&EXCLUSION` 在移动的显式原子周围关掉隐式介电和离子。三者都有解析力。
- **恒电势。** `&DFT &SCF &GCE` 让电子数取分数，并把功函数控制到目标电位。方案有 `CAPACITANCE`（电子数的割线更新）和 `DYNAMICS`（电子数作为带虚拟质量的动力学变量，可接 Nosé–Hoover 链）。电解质充当对电荷。电子数状态写进 restart 文件。
- **求解器。** `SOLVER NEWTON_CG`：对广义 Poisson 方程的真实残差做 Newton 迭代，内层预条件 CG，对自由能做线搜索。提速手段有 ψ 形式、两级预条件、电势历史、与 SCF 耦合的容差。
- **赝核空腔密度**，用于 GTH 赝势（`ADD_CORE_DENSITY`），空腔不会在原子内部打开。

## 安装

代码是 6 个新模块加 23 个被修改的 CP2K 文件，全在 `src/` 下（新文件列在 `src/CMakeLists.txt`）。从解压好的 CP2K 2026.2 源码树开始。

```bash
tar xjf cp2k-2026.2.tar.bz2
tools/install.sh cp2k-2026.2          # 把 src/ 拷到对应文件上，旧文件留作 *.orig
# 或者等价地:
cd cp2k-2026.2 && patch -p1 < ../patches/cp2k-2026.2.diff
```

`install.sh` 会拒绝 `src/cp2k_info.F` 不是 2026.2 的源码树。然后照常编译 CP2K（CMake，见 CP2K 的 `INSTALL.md`），用 `psmp` 版本（MPI 加 OpenMP）。

已经编译并运行过的工具链：

| 编译器与 MPI | 库 | 编译选项 |
|---|---|---|
| gcc 11.4, Intel MPI 2021.7.1 | MKL 2022.2.1, FFTW 3.3.10 | `-O3 -march=x86-64-v3 -mtune=generic` |
| gcc 13, OpenMPI 4.1 | OpenBLAS, build type GENERIC | `-O3 -march=x86-64-v3 -mtune=generic` |
| Cray MPICH | Cray 库 | `-O3 -march=znver3` |

在 Cray MPICH 上，变量未设置、或为 `CMA`、`XPMEM` 时，程序在 MPI 启动前把 `MPICH_SMP_SINGLE_COPY_MODE` 设为 `NONE`。单拷贝传输会在 ScaLAPACK 对角化里泄漏内存，丢掉约 4 GiB 之后各 rank 在 BLACS 里永远等待，长时间的恒电势运行会因此卡住。作业脚本里设的其他值会被保留。

并行布局属于作业设置，代码默认值不依赖它。其他布局的结果与 64 × 1 默认布局一致，|dE| ≤ 3.5e-10 Ha，max |dF| ≤ 8.8e-7 Ha/bohr。64 核节点可以从 32 个 rank × 2 线程开始，设 `OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close`，并在 `&DFT &MGRID &RS_GRID` 下写 `MAX_DISTRIBUTED_LEVEL 1`；16 核节点用默认的 16 × 1。

改了 `src/` 之后重新生成补丁：`tools/make_patch.sh <原版 cp2k-2026.2 根目录>`。

## 快速上手

例 1 和例 2 是 `tests/` 里某个算例的 `&SCCS`（和 `&SCF`）部分，晶胞、坐标、基组和其余输入照那个文件；例 3 是 203 原子膜的生产设置。三段输入都能通过 `cp2k.psmp --check` 解析。

### 1. 隐式溶剂里的分子

基础文件：`tests/regtest-sccs-lpb-1/H2O_isodensity_lpb.inp`（10 Å 盒子里的 H2O，PBE，`CUTOFF 400`）。把它的 `&SCCS` 块换成：

```
&SCCS on
  METHOD ISODENSITY
  RELATIVE_PERMITTIVITY 78.4
  GAMMA [mN/m] 18.47
&END SCCS
```

这是 ISO 模型。`ISODENSITY` 默认 `SOLVER NEWTON_CG`；`GAMMA` 是为分子溶剂化自由能拟合的空化表面张力。这个算例 18 步 SCF 收敛，能量约 −17.2305 Ha。

### 2. 恒电势下带电的金属表面

基础文件：`tests/regtest-sccs-lpb-3/Al111_gce_cap_lpb.inp`（Al(111) 2×2，3 层，Γ 点）。把它的 `&SCCS` 和 `&SCF` 块换成下面的内容。恒电势运行需要 `ADDED_MOS`、`&DIAGONALIZATION` 和 `&SMEAR METHOD FERMI_DIRAC`：

```
&SCCS on
  METHOD ISODENSITY
  RELATIVE_PERMITTIVITY 78.4
  DEBYE_LENGTH [angstrom] 3.0       ! 1 M 1:1 电解质, 充当电极的对电荷
&END SCCS
&SCF
  EPS_SCF 1.0E-6
  MAX_SCF 150
  ADDED_MOS 20
  &SMEAR on
    METHOD FERMI_DIRAC
    ELECTRONIC_TEMPERATURE [K] 300
  &END SMEAR
  &DIAGONALIZATION on
    ALGORITHM STANDARD
  &END DIAGONALIZATION
  &MIXING
    METHOD BROYDEN_MIXING
    ALPHA 0.2
    NBROYDEN 8
  &END MIXING
  &GCE on
    UPDATE_SCHEME CAPACITANCE
    UPDATE_EVERY SCF_CYCLE
    UPDATE_INTERVAL 5
    TARGET_POTENTIAL [eV] -0.6      ! U = -0.6 V vs SHE
    PHI_SHE [eV] 4.600
    FERMI_TOL [eV] 0.01
    MAX_OUTER_UPDATES 3
  &END GCE
&END SCF
```

电子数 N_e 不断调整，直到功函数达到目标。每个离子步打印一行 `GCE| trace`，记录达到的电位（见[输出](#输出)）；运行结束时 |U − U_T| 小于 0.01 V。

### 3. 恒电势下的显式水膜（SE 加显式水）

这是文中 203 原子 Cu(111)/水/K⁺ 膜所用的设置：4×4×4 的 Cu 平板（最下面 32 个原子固定），46 个水分子和一个 K⁺，晶胞 10.267 × 8.892 × 47.000 Å，Γ 点，PBE 加 DFT-D3(BJ)，`CUTOFF 400`，`REL_CUTOFF 40`，300 K 的 NVT（CSVR 恒温器），步长 1 fs，保持在 −0.600 V vs SHE，每 5 步更新一次电子数。种类有 `Cu`、`K`、`Ow`（水的 O）和 `Hw`（水的 H）。这里只给隐式溶剂、`&EXCLUSION`、`&REGION`/`&FLEXIBLE` 和恒电势几块；晶胞、坐标、`&MOTION` 和种类定义与普通 MD 输入相同。

```
&SCCS on
  RELATIVE_PERMITTIVITY 78.4
  GAMMA [mN/m] 0.0
  DELTA_RHO 2.0E-5
  DERIVATIVE_METHOD fft
  EPS_SCCS 1.0E-6
  EPS_SCF 0.5
  MAX_ITER 200
  SOLVER NEWTON_CG
  SOLVER_TOLERANCE 1.0E-8
  SOLVER_MAX_ITER 50
  DEBYE_LENGTH [angstrom] 3.0
  NONLINEAR_IONS F
  NONLINEAR_DIELECTRIC F
  SOLVER_GRID_LEVEL 1
  METHOD SOLVENT_EXCLUDED
  &SOLVENT_EXCLUDED
    RHO_CAVITY [angstrom^-3] 0.0158
    SWITCH_WIDTH 0.6
    PROBE_KERNEL_WIDTH [angstrom] 0.125
    SOLVENT_RADIUS [angstrom] 1.4
    ION_RADIUS [angstrom] 4.0
    DIELECTRIC_RADIUS [angstrom] 0.7
    CAVITATION_RADIUS [angstrom] 0.0
  &END SOLVENT_EXCLUDED
  &REGION
    AXIS 3
    WIDTH [angstrom] 0.5
    PROBE_ATOMS 65                  ! K+ 离子; 只用于报告
    PROBE_RADIUS [angstrom] 3.0
    &FLEXIBLE
      KIND Ow 1.0                   ! 只列溶剂的氧和离子
      KIND K 1.0
    &END FLEXIBLE
  &END REGION
  &EXCLUSION on
    KIND Ow 1.80 2.40               ! 种类, r_diel, r_ion, 单位 angstrom
    KIND Hw 1.10 1.70
    KIND K  2.00 2.60
    WIDTH [angstrom] 0.25
    PRINT_SUMMARY T
  &END EXCLUSION
  SCF_TOL_FACTOR 0.1
  SOLVER_STRIDE 1
  PSI_FORM T
&END SCCS
&SCF
  &GCE on
    UPDATE_SCHEME CAPACITANCE
    UPDATE_EVERY SCF_CYCLE
    UPDATE_INTERVAL 5               ! 每 5 个 MD 步更新一次 N_e
    TARGET_POTENTIAL [eV] -0.600    ! U = -0.6 V vs SHE
    PHI_SHE [eV] 4.600
    CAPACITANCE_INIT 0.4749
    FERMI_TOL [eV] 0.01
    FERMI_UPDATE_TOL [eV] 0.005
    MIN_STEPS_AFTER_UPDATE 12
    MAX_OUTER_UPDATES 6
    EXCESS_ELECTRONS 0.673537       ! 这一帧的起始 N_e
  &END GCE
  EPS_SCF 1.0E-6
  MAX_SCF 200
  SCF_GUESS ATOMIC
  ADDED_MOS 117
  &SMEAR on
    METHOD FERMI_DIRAC
    ELECTRONIC_TEMPERATURE [K] 300
  &END SMEAR
  &DIAGONALIZATION
    ALGORITHM STANDARD
  &END DIAGONALIZATION
  &MIXING on
    METHOD BROYDEN_MIXING
    ALPHA 0.2
    BETA 1.5
    NBROYDEN 8
  &END MIXING
&END SCF
```

溶剂是 `DIELECTRIC_RADIUS` 0.7 的 SE，加上 `&EXCLUSION` 和 `&REGION &FLEXIBLE`，即为显式水加排除球和原子跟随边界。`&FLEXIBLE` 列出水的氧和 K⁺，不列 Cu 原子。`EXCESS_ELECTRONS` 和 `CAPACITANCE_INIT` 是这一帧的起始值；换一个体系要用它自己的值，或者用默认值。完整的输入（含晶胞和坐标）在 `examples/film203_se/`。想用小体系试这几块，可以拿 `tests/regtest-sccs-lpb-5/K6H2O_dbg_se_exclusion_linear_all.inp`（K⁺ 加 6 个 H2O）作基础。

## 模型怎么选

| 模型 | 关键字 | 参数 |
|---|---|---|
| ISO | `METHOD ISODENSITY`，线性介电，线性 Poisson–Boltzmann | n_c = 3.00e-3 Å⁻³（`RHO_CAVITY`，默认值），σ = 0.6（`SWITCH_WIDTH`） |
| SE | `METHOD SOLVENT_EXCLUDED`，线性介电（`NONLINEAR_DIELECTRIC F`），`DIELECTRIC_RADIUS 0.7`，线性 Poisson–Boltzmann | n_c = 1.58e-2 Å⁻³，σ = 0.6，R_solv = 1.4 Å（`SOLVENT_RADIUS`），R_ion = 4.0 Å（`ION_RADIUS`），R_diel = 0.70 Å。晶胞含显式溶剂原子时，加 `&EXCLUSION` 和 `&REGION &FLEXIBLE`（排除球和原子跟随边界）；`&FLEXIBLE` 里只列溶剂的氧和离子。没有显式溶剂原子时就是普通的 SE |
| SE-NLD | `METHOD SOLVENT_EXCLUDED`，非线性介电（该方法的默认） | R_diel = 1.0 Å（默认值） |
| SE-NLPB | 在 SE-NLD 上加 `NONLINEAR_IONS T`、`ION_STERIC_RADIUS 4.0` | 带尺寸的非线性离子 |

所有模型的体相介电常数都是 ε_b = 78.4（`RELATIVE_PERMITTIVITY 78.4`）。

大多数情况一条规则就够：晶胞里没有显式溶剂用 ISO，有显式水用 SE（加 `&EXCLUSION` 和 `&REGION &FLEXIBLE`）。

- 没有显式溶剂的界面和分子：ISO。
- 有显式水的晶胞：SE 加 `&EXCLUSION` 和 `&REGION &FLEXIBLE`。`&FLEXIBLE` 里只列溶剂的氧和离子，不列电极原子。
- 孤立阳离子：ISO，`RHO_CAVITY [angstrom^-3] 0.0010`。SE 对孤立的小阳离子太强。
- SE-NLD 和 SE-NLPB 加入偶极饱和和有限尺寸离子。在 Ag(100) 上，线性模型给出平的 C(U)，SE-NLD 给出钟形，SE-NLPB 在零电荷电位处出现谷。

## 关键字

完整参考在 [`docs/keywords.zh.md`](docs/keywords.zh.md)（[English](docs/keywords.md)）：每个关键字的含义、默认值和单位，原名对照表，以及符号表。分节如下：

| 分节 | 控制什么 |
|---|---|
| `&DFT &SCCS` | 模型、求解器、电解质、非线性介电、赝核空腔、提速选项、诊断 |
| `&SCCS &ISODENSITY`、`&SCCS &SOLVENT_EXCLUDED` | 两种方法的空腔参数 |
| `&SCCS &REGION`、`&REGION &FLEXIBLE`、`&SCCS &EXCLUSION` | 隐式溶剂里的显式水 |
| `&DFT &SCF &GCE` | 恒电势：方案、目标、更新循环、restart 状态、电子数动力学 |

## 推荐设置

没写的都用默认值。298 K 的水：`RELATIVE_PERMITTIVITY 78.4`；1 M 1:1 电解质：`DEBYE_LENGTH [angstrom] 3.0`。有电解质时用 `&POISSON PERIODIC XYZ`、`POISSON_SOLVER PERIODIC`。

**电位标度。** U = W − Φ_SHE，W 是功函数，Φ_SHE = 4.600 eV，也是 `PHI_SHE` 的默认值。`TARGET_POTENTIAL [eV]` 里填 U，单位按 V vs SHE。

**表面张力。** 平板、吸附能、零电荷电位和电容：`GAMMA` 取 0（默认）。分子的溶剂化自由能：ISO 用 `GAMMA [mN/m] 18.47`，SE 用 `GAMMA [mN/m] 11.09`。

**恒电势 MD。** `UPDATE_SCHEME CAPACITANCE`、`UPDATE_EVERY SCF_CYCLE`、`UPDATE_INTERVAL 5`：每 5 个 MD 步更新一次电子数。代码默认值是 1（每步更新），推荐每 5 步，在速度和精度之间取得平衡。在 −0.600 V 的 203 原子 Cu/水/K⁺ 膜上（2 ps，步长 1 fs），每步更新时电位起伏 5.10 mV，每步 42.18 s；每 5 步更新时起伏 22.40 mV，每步 32.54 s，省时 22.9 %。`MAX_SCF` 计的是一个更新周期内的内层步数，60 就够。restart 文件带有全部状态，可续算的 `GEO_OPT` 或 MD 不需要额外关键字。

对从远离目标处出发的大晶胞单点，把下面几行加进 `&GCE` 块：

```
PIN_MIXER_G0 T
EARLY_UPDATE_DELTA 3.0E-3
EARLY_UPDATE_MIN_STEPS 6
EARLY_UPDATE_WF_STABLE [eV] 0.003
EARLY_UPDATE_MAX 4
EARLY_UPDATE_SECANT T
CAP_REFRESH_MIN [eV] 0.05
FINAL_TOL_FRAC 0.25
```

在 323 原子的 Cu/水/K⁺ 膜上（ISO，42 核），SCF 步数从 164 降到 92，墙钟时间从 862 s 降到 533 s。起点离目标 0.02 e 以内的晶胞只省 0–23 % 的步数，还可能更慢，所以先在一个点上比较。

**电子数动力学。** `UPDATE_SCHEME DYNAMICS`，`POTENTIOSTAT_PERIOD` 取远大于 10 个时间步，`POTENTIOSTAT_CHAIN_LENGTH 3` 做正则抽样，`POTENTIOSTAT_CAPACITANCE` 用静态 N_e 扫描测得。

**不推荐。** 电解质做对电荷时用 `UPDATE_SCHEME MU_MIXING`：电子数在每个 SCF 步都在变，某个瞬时值可能让离子没有可用体积，运行会因此停止。请用 `CAPACITANCE`。

## 性能与精度

| 量 | 结果 |
|---|---|
| 溶剂化自由能，FreeSolv，222 个分子，ISO | MAE 1.52 kcal/mol |
| 零电荷电位，九个 Cu、Ag、Au 晶面对实验，ISO，Φ_SHE = 4.600 eV | RMS 0.090 V |
| 连续介质求解，Pt(111) 3×3，3 层，ISO，1 M 电解质 | 每个 SCF 步 0.128 s |
| 恒电势 MD，203 原子膜，每步更新 / 每 5 步更新 | 每步 42.18 s，σ(U) 5.10 mV / 每步 32.54 s，σ(U) 22.40 mV |

## 输出

每个离子步一行记录电位：

```
GCE| trace  step <n>  N_e <x> e  W_f <y> eV  U <z> V vs SHE  U_target <t>  dU <d> V
GCE| omega  step <n>  E <e> Ha  n_exc <x> e  Omega <w> Ha  WF_target <t> eV
GCE| conserved  MD step <m>  Kin <k> Ha  E <e> Ha  dE <d> Ha  n_exc <x> e  Omega <w> Ha  Cons <c> Ha  Cons_Omega <co> Ha
GCE| dyn  step <n>  md_step <i>  N_e <x> e  dU <y> V  meas <T|F>  P_n <p>  Omega <o> Ha  Hextra <h> Ha  K_n <k>  E_th <t>
```

开了 `&GCE` 的运行每个离子步都打印 `trace` 和 `omega`，MD 运行打印 `conserved`，`DYNAMICS` 打印 `dyn`。`GCE| outer update`、`GCE| early update` 和 `GCE| carry` 标记 N_e 的更新；restart 之后会打印 `GCE| state read from the restart`。

`Omega = E + WF_target n_exc`，其中 `n_exc = N_e − N_input`，`WF_target = −mu_e`。N_e 变化时守恒量是 `Cons_Omega = Cons + WF_target n_exc`，不是 `Cons`。

## 测试

`tests/regtest-sccs-lpb-{1..7}` 沿用 CP2K 回归测试的目录布局，通过 `TEST_FILES.toml` 读取：H2O 上的能量层、有限差分力层（`RUN_TYPE DEBUG`）、Al(111) 上的恒电势层、保护性检查、排除掩模、柔性区域和 restart 状态。把它们拷到 CP2K 源码树的 `tests/QS/`，在 `tests/TEST_DIRS` 里加上 `QS/regtest-sccs-lpb-1` … `QS/regtest-sccs-lpb-7`，再照常运行 CP2K 回归测试器。

必须中止的算例、以及需要前一个算例 restart 文件的 restart 算例，没有列在 `TEST_FILES.toml` 里，按 `tests/README.md` 的说明手动运行；那份文件还列出各项容差和匹配器表达不了的判据（例如 Al(111) 算例的输出里 `grep -c "GCE| outer update"` 必须大于 0）。

`tests/fivepoint/` 是 Pt(111) 4×4×4 在 SE-NLPB 下的五点充电扫描：每个算例运行 `cp2k.psmp -i fp_x00.inp -o fp_x00.out`，然后 `python3 fit_fivepoint.py fp_*.out`。它复现 C = 1.1764 e/V = 17.70 µF/cm² 和 U_pzc = 0.8500 V vs SHE。


## 引用与许可

实现所依据的方法：

- 等密度空腔、赝核密度：https://doi.org/10.1063/1.4865107
- 隐式电解质，G = 0 电势参考的线性化 Poisson–Boltzmann：https://doi.org/10.1063/1.5132354
- 溶剂排除空腔、非线性介电、带尺寸的离子：https://doi.org/10.1063/5.0176308
- 非线性流体响应：https://doi.org/10.1088/0965-0393/21/7/074005
- 带尺寸修正的 Poisson–Boltzmann 方程：https://doi.org/10.1103/physrevlett.79.435
- 带尺寸修正的 Poisson–Boltzmann 方程的求解方案：https://doi.org/10.1021/acs.jctc.6b00435
- 自洽连续介质溶剂化，`METHOD ANDREUSSI`：https://doi.org/10.1063/1.3676407
- 广义 Poisson 与 Poisson–Boltzmann 求解器，ψ 形式：https://doi.org/10.1063/1.4939125
- CP2K 中的巨正则系综方法：https://doi.org/10.1021/acs.jctc.4c00671
- CP2K 中的溶剂感知 SCCS，`METHOD SAA_ANDREUSSI`：https://doi.org/10.1016/j.cpc.2025.109563
- Nosé–Hoover 链积分器：https://doi.org/10.1080/00268979600100761

CP2Ksol 由 Leshen Zhang（https://github.com/leshenzhang）开发和维护。使用本代码时，请引用本仓库（见 `CITATION.cff`），以及上面与你所用模型相关的论文。

GPL-3.0-or-later。CP2Ksol 是对 CP2K 的修改；CP2K 由其开发者按 GPL-2.0-or-later 发布，该许可允许按第 3 版或更新版本再分发。改自 CP2K 的文件保留上游文件头，本仓库中的修改和新增模块按 GPL-3.0-or-later 授权。见 `LICENSE`。
