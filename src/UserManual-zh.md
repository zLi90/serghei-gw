# SERGHEI 用户手册 — 地表-地下水动力、溶质输运与 WOFOST 作物生长耦合

> **范围**：本手册覆盖 SERGHEI 中的**地表水水动力**模块（[第 7 节](#7-地表水水动力swe)）、**地下水水动力**（Richards 方程地下水流）模块（[第 8 节](#8-地下水水动力地下水流)）、**溶质（反应性）输运**模块（地表水与地下水）与 **WOFOST 7.2 作物生长**模块，以及**触发式灌溉 / 触发式排水**地表边界类型（`extbc.input` 中的 `bctype` 13/14）。
>
> **参证测试工况**：`tests/t9-gaoyou-Ncycle-WOFOST/`（含氮循环的稻田）、`tests/t15-gaoyou-wofost-triggerIrrigation/`（基于规则的灌排管理，水深触发）、`tests/t16-gaoyou-wofost-triggerIrrigation/`（W1 浅湿露干制度，含根系层土壤水触发阶段）

---

## 目录

- [1. 地表水溶质输运](#1-地表水溶质输运)
  - [1.1 输运参数 rttransportsw.input](#11-输运参数-rttransportswinput)
  - [1.2 边界条件 rtswbc.input](#12-边界条件-rtswbcinput)
  - [1.3 源汇项 rtswss.input](#13-源汇项-rtswssinput)
  - [1.4 反应参数 reactionsw.input](#14-反应参数-reactionswinput)
- [2. 地下水溶质输运](#2-地下水溶质输运)
  - [2.1 输运参数 rttransportgw.input](#21-输运参数-rttransportgwinput)
  - [2.2 边界条件 rtgwbc.input](#22-边界条件-rtgwbcinput)
  - [2.3 反应参数 reactiongw.input](#23-反应参数-reactiongwinput)
  - [2.4 地下水源汇 rtgwss.input](#24-地下水源汇-rtgwssinput)
  - [2.5 初始条件文件](#25-初始条件文件)
- [3. WOFOST 作物生长模型](#3-wofost-作物生长模型)
  - [3.1 作物参数 cropparameter.input](#31-作物参数-cropparameterinput)
  - [3.2 气象数据 cropmeteo.input](#32-气象数据-cropmeteoinput)
  - [3.3 农业管理 agro.input](#33-农业管理-agroinput)
  - [3.4 WOFOST-SERGHEI 耦合接口](#34-wofost-serghei-耦合接口)
- [4. 触发式灌溉与排水边界](#4-触发式灌溉与排水边界)
  - [4.1 概述](#41-概述)
  - [4.2 边界定义 extbc.input](#42-边界定义-extbcinput)
  - [4.3 灌溉策略文件](#43-灌溉策略文件)
  - [4.4 排水策略文件](#44-排水策略文件)
  - [4.5 触发多边形](#45-触发多边形)
  - [4.6 运行行为与输出](#46-运行行为与输出)
- [5. 时间序列文件格式](#5-时间序列文件格式)
- [6. 输出文件与表头](#6-输出文件与表头)
- [7. 地表水水动力（SWE）](#7-地表水水动力swe)
  - [7.1 模拟控制参数 parameters.input](#71-模拟控制参数-parametersinput)
  - [7.2 地形 dem.input](#72-地形-deminput)
  - [7.3 初始状态与摩阻 sw.input](#73-初始状态与摩阻-swinput)
  - [7.4 外部边界 extbc.input](#74-外部边界-extbcinput)
  - [7.5 降雨 rainfall.input（可选）](#75-降雨-rainfallinput可选)
  - [7.6 入渗 infiltration.input（可选）](#76-入渗-infiltrationinput可选)
  - [7.7 风场 wind.input（可选）](#77-风场-windinput可选)
- [8. 地下水水动力（地下水流）](#8-地下水水动力地下水流)
  - [8.1 计算域与求解器 subsurface.input](#81-计算域与求解器-subsurfaceinput)
  - [8.2 土壤水力参数 vg.input](#82-土壤水力参数-vginput)
  - [8.3 土壤类型分布 soilID.input](#83-土壤类型分布-soilidinput)
  - [8.4 边界条件 gwbc.input](#84-边界条件-gwbcinput)
  - [8.5 源汇区 gwss.input](#85-源汇区-gwssinput)
  - [8.6 初始条件](#86-初始条件)
  - [8.7 地表与地下水耦合（小结）](#87-地表与地下水耦合小结)
- [附录 A 组分顺序约定](#附录-a-组分顺序约定)
- [附录 B 输入文件速查清单](#附录-b-输入文件速查清单)

---

## 1. 地表水溶质输运

### 1.1 输运参数 rttransportsw.input

**文件路径**：`<input_folder>/rttransportsw.input`
**源代码**：`src/RTInitSW.h` → `readRTFileSW()`

本文件定义溶质组分数、逐组分输运属性、初始条件以及地表水反应输运的全局求解器设置。

#### 文件格式

冒号分隔的键–值对。逐组分参数通过 `id` 选择器分组。以 `//` 开头的行为注释。

#### 全局参数

| 键 | 类型 | 必填 | 说明 |
|-----|------|----------|-------------|
| `n_mass` | int | **是** | 溶质组分（化学物种）总数。必须是文件中第一个与组分相关的键。 |
| `ReactionModule` | int | **是** | 反应模块开关：`0` = 关闭，`1` = 开启（需要 `reactionsw.input`）。 |
| `Advection_Scheme` | int | **是** | 对流输运的数值格式： |

| `Advection_Scheme` 取值 | 名称 | 说明 |
|:---:|------|-------------|
| 1 | 一阶迎风 | 一阶上游加权。数值耗散最大，最稳定。 |
| 2 | TVD Van Leer | 带Van Leer通量限制器的总变差减小格式。推荐用于大多数场景。 |
| 3 | TVD Superbee | 带 Superbee 通量限制器的总变差减小格式。锋面最锐利，但在间断附近可能产生振荡。 |

#### 逐组分参数

每个组分块以 `id : <组分序号>`（0 起）开始，后续所有键均作用于该组分，直至下一个 `id` 键。

| 键 | 类型 | 必填 | 说明 |
|-----|------|----------|-------------|
| `aq_mode` | int | **是** | 液相（水相）初始条件模式：`0` = 常数（`aq_val`），`1` = 空间文件（`aq_file`）。 |
| `aq_val` | real | 若 aq_mode=0 | 均匀初始浓度 [mg/L]。解析器不做其他浓度单位的换算。 |
| `aq_file` | string | 若 aq_mode=1 | 空间分布初始浓度的文件名（ESRI ASCII Grid 或简化格式；见[第 2.5 节](#25-初始条件文件)）。 |
| `diffusion_molecular` | real | **是** | 分子扩散系数 [m²/s]。常见范围：1.0×10⁻¹⁰ 至 1.0×10⁻⁸。 |
| `alpha_L` | real | **是** | 纵向弥散度 [m]。常见范围：0.01 至 100。 |
| `alpha_T` | real | **是** | 横向弥散度 [m]。常见范围：0.001 至 10，常取 `alpha_L` 的 1/10。 |
| `RainCon` | real | **是** | 该组分在降雨中的浓度 [mg/L]。若降雨不携带溶质则填 `0.0`。 |

#### 示例

```text
// --- Reactive Transport Parameters ---
n_mass : 5        // Total number of species (NH4, NO3, DO, DOC, DON)

// ========== Species 0: NH4 ==========
id : 0
aq_mode : 1
aq_file : rtinitialswNH4.input
RainCon : 0.5                           // [mg/L]
diffusion_molecular : 1.76e-9           // [m^2/s]
alpha_L : 0.5                           // [m]
alpha_T : 0.05                          // [m]

// ========== Species 1: NO3 ==========
id : 1
aq_mode : 1
aq_file : rtinitialswNO3.input
RainCon : 0.4                           // [mg/L]
diffusion_molecular : 1.76e-9           // [m^2/s]
alpha_L : 0.5                           // [m]
alpha_T : 0.05                          // [m]

// ========== Species 2: DO ==========
id : 2
aq_mode : 0
aq_val : 5.0                            // [mg/L]
RainCon : 4.5                           // [mg/L]
diffusion_molecular : 1.76e-9           // [m^2/s]
alpha_L : 0.5                           // [m]
alpha_T : 0.05                          // [m]

// ... species 3 (DOC), 4 (DON) follow the same pattern ...

ReactionModule : 1
Advection_Scheme : 2                    // TVD Van Leer
```

---

### 1.2 边界条件 rtswbc.input

**文件路径**：`<input_folder>/rtswbc.input`
**源代码**：`src/RTInitSW.h` → `readRTBCFileSW()`（解析与配对）、`src/BC.h` → `applyrtbc()`（门控解析与施加）、`src/RTFunctionSW.h`（每个输运半步同步灌溉/排水开关状态）

定义各组分在外边界上的溶质边界条件。边界几何（格胞索引）继承自水动力边界条件（`extbc.input`）。

> **配对规则**：`rtswbc.input` 的边界块按**顺序**与 `extbc.input` 的边界一一配对（第 k 个 `id:` 块 ↔ 第 k 个水动力边界），`id:` 后的名称本身不被解析。门控类型 5/6 依赖此配对读取水动力边界类型与运行时开关状态，**块顺序错位会导致门控失效**。

#### 边界条件类型

| 取值 | 常量名 | 说明 |
|:---:|------|-------------|
| 1 | `RTSW_BC_NONE` | 该边界不施加反应输运边界条件。 |
| 2 | `RTSW_BC_CON_CONST` | **定浓度**（Dirichlet）。边界格胞被固定在 `spec_X_bcval` 给定的浓度上。 |
| 3 | `RTSW_BC_ZEROGRAD` | **零梯度**（Neumann，自由出流）。边界浓度由内部外推（边界法向浓度梯度为零）。用于自由泄流出口。 |
| 4 | `RTSW_BC_CON_T` | **时变浓度**（Dirichlet 时间序列）。边界浓度按 `spec_X_bcfile` 指定的外部时间序列文件随时间变化。时间序列格式见[第 5 节](#5-时间序列文件格式)。 |
| 5 | `RTSW_BC_TRIGGERED_IRRIGATION` | **触发门控 Dirichlet**，与水动力 `bctype : 13` 配对。灌溉**开启**期间边界表现为类型 2/4（按 `spec_X_bcval` 定浓度，或若给了 `spec_X_bcfile` 时间序列则按其执行——`bcfile` 优先）；灌溉**关闭**期间回退为零梯度，边界浓度跟随内部水体，不强行进出任何溶质。质量收支中的对流入流项在开启期间已完整计入随水输入的溶质（输运核经幽灵单元把 Dirichlet 值平流进来），无需额外记账。 |
| 6 | `RTSW_BC_TRIGGERED_DRAINAGE` | **零梯度出流**，触发式排水边界（`bctype : 14`）的显式类型。与类型 3 完全相同：只有排水实际出流时溶质才离开；出口关闭期间不强行进出任何物质。 |

#### 文件格式

```text
bccount : <number_of_boundaries>

id : <boundary_name>            // 名称不被解析；第 k 个 id 块对应 extbc.input 的第 k 个边界
direction : <nx> <ny>           // 法向向量（会被归一化为单位长度）
rtbctype : <type>               // 边界级默认类型（可选；未写 spec_X_bctype 的组分继承此类型）

// 逐组分边界条件（X = 组分索引，0 基）
spec_X_bctype : <type>          // 边界类型（1–6，见上表）
spec_X_bcval  : <value>         // 定浓度 [mg/L]（用于类型 2 与 5）
spec_X_bcfile : <filename>      // 时间序列文件（用于类型 4，及类型 5 可选）
```

#### 交叉校验规则

- **类型 2**（定浓度）：使用 `spec_X_bcval`。**不要**指定 `spec_X_bcfile`。
- **类型 4**（时变）：使用 `spec_X_bcfile`。**不要**指定 `spec_X_bcval`。
- **类型 5**（触发式灌溉）：`spec_X_bcval`（开启期间定浓度）或 `spec_X_bcfile`（开启期间按时间序列；两者同给时时间序列优先，`bcval` 被静默忽略——请二选一）。对非触发式水动力边界，门控标志保持关闭，类型 5 退化为零梯度。
- **时间序列加载条件**：`spec_X_bcfile` **仅在 `spec_X_bctype` 为 4 或 5 时才被加载**。只写边界级 `rtbctype : 5` 而省略 `spec_X_bctype` 行时，时间序列被静默忽略，灌溉期退化为清水（Dirichlet 0）且不报错。配了 `bcfile` 的组分务必显式写 `spec_X_bctype : 5`。
- **触发式水动力边界上避免裸类型 4**：类型 4 不受门控，灌溉关闭期仍对边界施加 Dirichlet（无水流入也存在浓度源/汇，且该变化不计入边界通量记账）。肥水灌溉请使用类型 5。
- 解析器会拒绝不匹配的组合。

#### 示例 — 入流（时变）+ 出口（自由出流）

```text
bccount : 2

// --- Boundary 0: irrigation inlet ---
id : rtirrigation
direction : 0.0 1.0
spec_0_bctype : 4                        // NH4: time-varying
spec_0_bcfile : rtirrigationNH4.input
spec_1_bctype : 4                        // NO3: time-varying
spec_1_bcfile : rtirrigationNO3.input
spec_2_bctype : 4                        // DO: time-varying
spec_2_bcfile : rtirrigationDO.input
spec_3_bctype : 4                        // DOC: time-varying
spec_3_bcfile : rtirrigationDOC.input
spec_4_bctype : 4                        // DON: time-varying
spec_4_bcfile : rtirrigationDON.input

// --- Boundary 1: outlet ---
id : rtoutlet
direction : 0.0 1.0
spec_0_bctype : 3                        // NH4: free outflow
spec_1_bctype : 3                        // NO3: free outflow
spec_2_bctype : 3                        // DO: free outflow
spec_3_bctype : 3                        // DOC: free outflow
spec_4_bctype : 3                        // DON: free outflow
```

#### 示例 — 触发式灌溉/排水的门控溶质边界（类型 5/6，配水动力 `bctype : 13/14`）

与触发式水动力边界（[第 4 节](#4-触发式灌溉与排水边界)）配对时，溶质边界的推荐类型是 **5（灌溉）** 与 **6（排水）**。门控下时间序列的语义需特别说明：

> **`spec_X_bcfile` 时间序列在类型 5 下的含义**：它是**灌溉水自身浓度随日历时间（自模拟起算的秒数）变化的函数**，不是"何时灌水/施肥"的事件时刻表。何时进水、进多少水完全由水动力触发器在每个时间步动态决定；时间序列只回答"**若** t 时刻正在灌水，这些水的浓度是多少"。灌溉**开启**期间，边界按当前时刻从序列插值出的浓度施加 Dirichlet，实际入流通量 = Q(t)·c(t)——每次真实触发事件按其持续时段积分得到入流质量，是涌现量；灌溉**关闭**期间序列值不被读取，边界回退为零梯度。因此**无需预知触发时间**：把序列写成"施肥窗口"即可——窗口内触发灌水则进肥水，窗口外触发则进清水（序列值为 0）。

**情景 A — 肥水浓度随时间变化（时间序列）**

`rtswbc.input`（块顺序与 `extbc.input` 一致：第 1 块对应 `bctype : 13`，第 2 块对应 `bctype : 14`）：

```text
// --- 边界0：触发式灌溉（水动力 bctype : 13） ---
id : rtirrigation
direction : 0.0 1.0
rtbctype : 5                 // 门控：灌水=Dirichlet(查表)，停灌=零梯度；未列出的组分继承此类型

spec_0_bctype : 5            // NH4：必须显式写 5，否则下面的 bcfile 不会被加载
spec_0_bcfile : rtirrigationNH4.input
spec_1_bctype : 5            // NO3
spec_1_bcfile : rtirrigationNO3.input
spec_2_bctype : 5            // DO
spec_2_bcfile : rtirrigationDO.input
spec_3_bctype : 5            // DOC
spec_3_bcfile : rtirrigationDOC.input
spec_4_bctype : 5            // DON：清水（不给 bcval/bcfile → 灌水期 Dirichlet 0）

// --- 边界1：触发式排水（水动力 bctype : 14） ---
id : rtoutlet
direction : 0.0 1.0
rtbctype : 6                 // 恒零梯度出流；一个键即可，各组分自动继承（写 3 完全等价）
```

`rtirrigationNH4.input`（格式见[第 5 节](#5-时间序列文件格式)；数据文件中**不允许注释**）：

```text
ndata 14
0.0         0.0
489480.0    0.0
489600.0    0.768
497040.0    0.768
497200.0    0.0
835080.0    0.0
835200.0    0.639
842400.0    0.639
842600.0    0.0
1270440.0   0.0
1270740.0   1.348
1290540.0   1.348
1290600.0   0.0
10886400.0  0.0
```

上例定义了三个"施肥窗口"（0.768 / 0.639 / 1.348 mg/L，约 DAT 5.7 d、9.7 d、14.7 d；126 d = 10 886 400 s）。窗口内灌溉触发 → 进该浓度的肥水；窗口外触发 → 清水；全程停灌 → 零梯度、无任何溶质输入。

**情景 B — 灌溉水浓度恒定**

```text
// --- 边界0：触发式灌溉（水动力 bctype : 13） ---
id : rtirrigation
direction : 0.0 1.0
rtbctype : 5

spec_0_bctype : 5
spec_0_bcval : 0.768     // NH4  mg/L（按实测肥水水质填写）
spec_1_bctype : 5
spec_1_bcval : 0.10      // NO3
spec_2_bctype : 5
spec_2_bcval : 8.0       // DO（近饱和溶解氧量级）
spec_3_bctype : 5
spec_3_bcval : 2.0       // DOC
spec_4_bctype : 5        // DON：清水

// --- 边界1：触发式排水（水动力 bctype : 14） ---
id : rtoutlet
direction : 0.0 1.0
rtbctype : 6
```

整个生育期每次灌溉都按同一浓度进溶质；每次事件的入流质量仍由该次事件的流量与时长决定。

**注意事项**

1. 同一组分不要同时给 `bcfile` 与 `bcval`：时间序列优先，`bcval` 被静默忽略。
2. 配了 `bcfile` 的组分必须显式写 `spec_X_bctype : 5`（加载条件，见交叉校验规则）。
3. 时间序列首行建议放 `t = 0`；末行建议放到模拟总长——超出末点后**保持末值**，若末点留成非零浓度，会一直生效到模拟结束。
4. 验证方法：短工况运行后检查输入目录 `outRTSurfaceTimeSeries.out` —— 停灌期各组分 `_BCInflow[mg/s]` 应为 0，灌水期 `_BCInflow[mg/s] ≈ Q·c`，`_BalanceResidual[mg/s]` 全程闭合。

---

### 1.3 源汇项 rtswss.input

**文件路径**：`<input_folder>/rtswss.input`
**源代码**：`src/RTInitSW.h` → `readRTSSFileSW()`、`src/SourceSink.h`

定义向地表水域注入（或从中抽取）溶质质量的空间源汇区。每个源汇区由一个多边形确定其空间范围。

#### 全局源汇类型（`rtsstype`）

| `rtsstype` 取值 | 说明 |
|:---:|-------------|
| 3 | 面积质量通量注入（地表撒施，如施肥）[mg/m²/s]。代码自动把通量分配到格胞面积上。 |

#### 逐组分源汇类型（`spec_X_sstype`）

| 取值 | 说明 |
|:---:|-------------|
| 0 | **常数**。源汇速率取 `spec_X_ssval` 给定的固定值。 |
| 1 | **线性插值时间序列**。速率随时间按时间序列文件（`spec_X_ssfile`）变化，数据点间线性插值。 |
| 2 | **阶梯时间序列**。与 1 相同，但数据点间采用阶梯（分段常数）插值。用于施肥等瞬时事件。 |

#### 文件格式

```text
sscount : <number_of_sources>

id : <source_name>              // Source identifier
rtsstype : <global_type>        // Global source/sink type (3 = area flux)
polygon : <polygon_file>        // Polygon defining the spatial zone

// Per-species (X = species index)
spec_X_sstype : <type>          // 0=constant, 1=linear timeseries
spec_X_ssval  : <value>         // Area mass flux [mg/m^2/s] (if sstype=0)
spec_X_ssfile : <filename>      // Timeseries file (if sstype=1)
```

#### 示例 — 施肥

```text
sscount : 1

id : fertilizer
rtsstype : 3                            // Area mass flux [mg/m^2/s]
polygon : polygonfer.input              // Fertilizer application area

spec_0_sstype : 1                       // NH4: linear timeseries
spec_0_ssfile : rtfertNH4-DON-NH4.input

spec_1_sstype : 0                       // NO3: constant (no NO3 fertilizer)
spec_1_ssval  : 0

spec_2_sstype : 0                       // DO: no source
spec_2_ssval  : 0

spec_3_sstype : 0                       // DOC: no source
spec_3_ssval  : 0

spec_4_sstype : 1                       // DON: linear timeseries
spec_4_ssfile : rtfertDON-DON-NH4.input
```

> **注**：若找不到 `rtswss.input`，模型在不启用地表水溶质源汇的情况下继续运行（打印一条警告）。

---

### 1.4 反应参数 reactionsw.input

**文件路径**：`<input_folder>/reactionsw.input`
**源代码**：`src/RTInitSW.h` → `readRTReactionFile()`、`src/RTStateSW.h`

本文件**仅当** `rttransportsw.input` 中设置了 `ReactionModule : 1` 时才需要。它配置地表水的多组分 Monod 动力学反应模型，包括氮循环过程、复氧与环境校正因子。

#### 反应质量基准与固定组分顺序

激活的反应内核是固定零基索引的五组分体系：

| 索引 | 名称 | 浓度基准 |
|---:|---|---|
| 0 | NH4 | 化合物质量，`mg NH4/L`（摩尔质量约 18.04 g/mol） |
| 1 | NO3 | 化合物质量，`mg NO3/L`（摩尔质量约 62.00 g/mol） |
| 2 | DO | 氧质量，`mg O2/L` |
| 3 | DOC | 元素碳质量，`mg C/L` |
| 4 | DON | 元素氮质量，`mg N/L` |

这不是一个通用的 `mg N/L` 模型。输入读取器不做 `mg NH4-N/L` 到 `mg NH4/L`、或 `mg NO3-N/L` 到 `mg NO3/L` 的换算；初始条件、边界、源汇文件与反应参数都必须按上表基准填写。`SPECIES_ID` 与 `SPECIES_Name` 仅作信息记录，不会重排浓度数组。固定反应内核要求 `n_mass >= 5`；组分数少于五个时，反应步直接返回、不求解氮循环方程。

基准速率为 `R_min [mg N/L/s]`（DON 矿化）、`R_nit [mg NH4/L/s]`（NH4 硝化）与 `R_denit`/`R_aer [mg C/L/s]`（DOC 消耗）。源项应用以下硬编码质量比：每矿化 1 份 N 生成 NH4 用 `18.04/14.01 = 1.2877`；每氧化 1 份 NH4 生成 NO3 用 `62.00/18.04 = 3.4368`；每硝化 1 份 NH4 消耗 O2 用 `64.00/18.04 = 3.5477`；每矿化 1 份 N 消耗 O2 用 `32.00/14.01 = 2.2841`；每呼吸 1 份 C 消耗 O2 用 `32.00/12.01 = 2.6644`；每反硝化 1 份 C 生成 NO3 用 `(0.8*62.00)/12.01 ≈ 4.130`。其中 `0.8` 是硬编码的反硝化效率因子，不是输入参数。

**校验状态。** 解析器用 `-999` 哨兵值检查必需键是否存在，但一般不强制物理最小/最大值。下文范围是推荐率定范围，不是程序强制上限。

#### 全局配置

| 键 | 类型 | 必填 | 说明 |
|-----|------|----------|-------------|
| `ReactModel_Type` | int | **是** | 反应模型类型：`0` = 无反应，`1` = 多组分 Monod 动力学。 |
| `Nitrogen_Cycle_Simulation` | int | **是** | 氮循环开关：`0` = 关，`1` = 开。启用 NH4→NO3→N2 路径。 |
| `SPECIES_ID` | int 列表 | 否 | 信息性。组分索引列表（分号分隔），如 `0; 1; 2; 3; 4`。 |
| `SPECIES_Name` | string 列表 | 否 | 信息性。组分名称列表（分号分隔），如 `NH4; NO3; DO; DOC; DON`。 |

#### 一阶衰减速率

| 键 | 类型 | 单位 | 说明 |
|-----|------|-------|-------------|
| `Lambda_1` | real 列表 | [1/s] | 各组分在**液相**中的一阶衰减速率（分号分隔，每组分一个值）。常见范围 0 至 1×10⁻⁵。不衰减的组分填 `0.0`。 |
| `Lambda_2` | real 列表 | [1/s] | 各组分在**固相**中的一阶衰减速率。格式同 `Lambda_1`。地表水无固相，填 `0.0`。 |

#### 硝化（NH4⁺ → NO3⁻）

硝化细菌（Nitrosomonas）对铵的生物氧化。

| 键 | 类型 | 单位 | 常见范围 | 说明 |
|-----|------|-------|---------------|-------------|
| `Rate_Max_Nit` | real | [mg NH4/L/s] | 1×10⁻⁷ – 1×10⁻⁴ | NH4 化合物质量基准的最大硝化速率。 |
| `K_Monod_NH4` | real | [mg NH4/L] | 0.01 – 5.0 | 化合物质量基准的 NH4 半饱和浓度。 |
| `K_Monod_DO` | real | [mg O2/L] | 0.01 – 5.0 | DO 半饱和浓度。 |

#### 反硝化（NO3⁻ → N2）

缺氧条件下异养细菌把硝酸盐还原为氮气。这是淹水稻田中氮流失的主导途径。

| 键 | 类型 | 单位 | 常见范围 | 说明 |
|-----|------|-------|---------------|-------------|
| `Rate_Max_Denit` | real | [mg C/L/s] | 1×10⁻⁷ – 3×10⁻⁴ | DOC-C 基准的最大反硝化速率；硝酸盐消耗使用硬编码的 0.8 因子。 |
| `K_Monod_NO3` | real | [mg NO3/L] | 0.01 – 5.0 | 化合物质量基准的 NO3 半饱和浓度。 |
| `K_Monod_DOC` | real | [mg C/L] | 0.01 – 20.0 | DOC 半饱和浓度。 |
| `Ki_Inhib_DO` | real | [mg O2/L] | 0.01 – 5.0 | DO 抑制浓度；DO 越大反硝化越弱。 |

#### 有氧呼吸（DOC + O2 → CO2）

以氧为电子受体的溶解性有机碳异养消耗。

| 键 | 类型 | 单位 | 常见范围 | 说明 |
|-----|------|-------|---------------|-------------|
| `Rate_Max_Hetero` | real | [mg C/L/s] | 1×10⁻⁷ – 3×10⁻⁴ | DOC-C 基准的最大有氧呼吸速率。 |
| `K_Monod_DOC_Aerobic` | real | [mg C/L] | 0.01 – 20.0 | 元素碳基准的 DOC 半饱和浓度。 |
| `K_Monod_DO_Aerobic` | real | [mg O2/L] | 0.1 – 2.0 | DO 半饱和常数。 |

#### 矿化（DON → NH4⁺）

溶解性有机氮分解为铵。

| 键 | 类型 | 单位 | 常见范围 | 说明 |
|-----|------|-------|---------------|-------------|
| `Rate_Max_Min` | real | [mg N/L/s] | 1×10⁻⁸ – 1×10⁻⁴ | DON-N 基准的最大矿化速率。 |
| `K_Monod_DON` | real | [mg N/L] | 0.01 – 10.0 | 元素氮基准的 DON 半饱和浓度。 |
| `K_Monod_DO_Min` | real | [mg O2/L] | 0.01 – 5.0 | 矿化的 DO 半饱和浓度。 |

#### 地表水特有的源汇参数

这些参数控制大气复氧与沉积物–水溶质交换等地表水独有的过程。

| 键 | 类型 | 单位 | 常见范围 | 说明 |
|-----|------|-------|---------------|-------------|
| `Ka_sw` | real | [1/s] | 0 – 1×10⁻⁴ | 大气复氧速率常数。验证算例可填 `0.0` 关闭复氧。典型值 5×10⁻⁵。 |
| `DO_sat` | real | [mg O2/L] | 0 – 15 | 复氧使用的饱和溶解氧浓度。 |
| `DOC_eq_sw` | real | [mg C/L] | 2 – 20 | 地表水 DOC 平衡浓度（沉积物释放终点）。肥沃稻田约 10 mg/L。 |
| `DON_eq_sw` | real | [mg N/L] | 0.1 – 5.0 | 地表水 DON 平衡浓度（沉积物释放终点）。 |
| `K_rel_sw` | real | [1/s] | 0 – 1×10⁻⁵ | 沉积物–水溶质交换的释放速率常数。填 `0.0` 关闭。 |

#### 环境校正因子

地表水生物速率的温度校正。`Opt_pH` 目前仅为保持输入兼容性而保留。

| 键 | 类型 | 单位 | 常见范围 | 说明 |
|-----|------|-------|---------------|-------------|
| `Temp_Coeff_Theta` | real | [-] | 1.00 – 1.10 | 当前地表水乘子 `theta^(20-Opt_Temp)` 的底数。 |
| `Opt_Temp` | real | [°C] | 0 – 40 | 地表水乘子中的参考温度；当前 `T_env` 固定为 20 °C。 |
| `Opt_pH` | real | [-] | 0 – 14 | 被解析且必填，但当前反应方程未使用。 |

#### 示例

```text
# === Global ===
ReactModel_Type           : 1                  // Multispecies Monod kinetics
Nitrogen_Cycle_Simulation : 1                  // Enable nitrogen cycle
SPECIES_ID   : 0; 1; 2; 3; 4
SPECIES_Name : NH4; NO3; DO; DOC; DON

# === Decay ===
Lambda_1 : 0.0; 0.0; 0.0; 0.0; 0.0
Lambda_2 : 0.0; 0.0; 0.0; 0.0; 0.0

# === Nitrification ===
Rate_Max_Nit  : 1.0e-5        // [mg NH4/L/s]
K_Monod_NH4   : 0.5           // [mg NH4/L]
K_Monod_DO    : 1.5           // [mg O2/L]

# === Denitrification ===
Rate_Max_Denit : 5.0e-6       // [mg C/L/s]
K_Monod_NO3    : 0.5          // [mg NO3/L]
K_Monod_DOC    : 5.0          // [mg C/L]
Ki_Inhib_DO    : 1.0          // [mg O2/L]

# === Aerobic Respiration ===
Rate_Max_Hetero       : 5.0e-5  // [mg C/L/s]
K_Monod_DOC_Aerobic   : 5.0     // [mg C/L]
K_Monod_DO_Aerobic    : 1.0     // [mg O2/L]

# === Mineralization ===
Rate_Max_Min   : 5.0e-6        // [mg N/L/s]
K_Monod_DON    : 0.21015       // [mg N/L]
K_Monod_DO_Min : 0.48          // [mg O2/L]

# === Surface Water Source/Sink ===
Ka_sw      : 0.0               // [1/s] (disabled for validation)
DO_sat     : 8.1               // [mg/L]
DOC_eq_sw  : 10.0              // [mg C/L]
DON_eq_sw  : 0.5               // [mg N/L]
K_rel_sw   : 0.0               // [1/s]

# === Environmental Correction ===
Temp_Coeff_Theta : 1.047
Opt_Temp         : 25.0        // [°C]
Opt_pH           : 7.0
```

---

## 2. 地下水溶质输运

### 2.1 输运参数 rttransportgw.input

**文件路径**：`<input_folder>/rttransportgw.input`
**源代码**：`src/RTInitGW.h` → `readRTFile()`

格式与 `rttransportsw.input` 相同，另增加若干地下水专有键。主要差异：

输运与反应代码消费的所有液相浓度均使用 **mg/L**。解析器不做 kg/m3 或其他浓度单位的换算。分子扩散系数与弥散度按下列 SI 单位填写。

#### 新增全局参数

| 键 | 类型 | 必填 | 说明 |
|-----|------|----------|-------------|
| `rt_scheme` | int | **是** | 输运矩阵的数值求解器：`1` = Gauss-Seidel，`2` = Kokkos BiCGSTAB（推荐用于大型三维域）。 |
| `Up_Weighting_vplus` | real | **是** | 正速度方向的上游加权系数 [-]。范围 [0.5, 1.0]：`1.0` = 完全上游，`0.5` = 中心差分。 |
| `Up_Weighting_vminus` | real | **是** | 负速度方向的上游加权系数 [-]。范围 [0.0, 0.5]：`0.0` = 完全上游，`0.5` = 中心差分。 |

#### 新增逐组分参数

地下水输运在液相之外还包含**固相**（吸附浓度）。

| 键 | 类型 | 必填 | 说明 |
|-----|------|----------|-------------|
| `solid_mode` | int | 若 ReactionModule=1 | 固相初始条件模式：`0` = 常数（`solid_val`），`1` = 空间文件（`solid_file`）。反应关闭时不需要。 |
| `solid_val` | real | 若 solid_mode=0 | 均匀初始固相（吸附）浓度 [mg/kg]。 |
| `solid_file` | string | 若 solid_mode=1 | 空间分布固相初始浓度的文件名。 |
| `aq_phase` | int | 否 | 当前地下水读取器不解析此键。请使用 `aq_mode`；不要把 `aq_phase` 当作别名。 |

#### 示例

```text
n_mass : 5

// === Species 0: NH4 ===
id : 0
aq_mode : 1
aq_file : rtinitialgwNH4.input
solid_mode : 1
solid_file : rtinitialgwsolidNH4.input
diffusion_molecular : 1.76e-9
alpha_L : 0.5
alpha_T : 0.05

// === Species 2: DO ===
id : 2
aq_mode : 0
aq_val : 1.0                            // [mg/L]
solid_mode : 0
solid_val : 0.0                         // [mg/kg] (DO not adsorbed)
diffusion_molecular : 1.76e-9
alpha_L : 0.5
alpha_T : 0.05

// ... other species ...

ReactionModule : 1
rt_scheme : 2                            // BiCGSTAB solver
Up_Weighting_vplus : 1.0                 // Full upstream (positive)
Up_Weighting_vminus : 0.0               // Full upstream (negative)
```

---

### 2.2 边界条件 rtgwbc.input

**文件路径**：`<input_folder>/rtgwbc.input`
**源代码**：`src/RTBCGW.h`（47–55 行）、`src/RTInitGW.h`

定义三维地下水域的溶质边界条件。边界通过多边形文件与面方向指定。

#### 面方向代码

| 方向 | 面 | 说明 |
|:---:|------|-------------|
| 1 | +x | 东面 |
| 2 | −x | 西面 |
| 3 | +y | 北面 |
| 4 | −y | 南面 |
| 5 | +z | **顶面**（地表，向上） |
| 6 | −z | **底面**（含水层底板，向下） |

#### 边界条件类型

| 取值 | 常量名 | 说明 |
|:---:|------|-------------|
| 1 | `SUB_RT_BC_NOFLOW` | **无流量**。水与溶质均不穿越边界。 |
| 2 | `SUB_RT_BC_Dirichlet_CONST` | **定浓度**（Dirichlet）。取 `spec_X_bcval` 给定的固定浓度。 |
| 3 | `SUB_RT_BC_Neumann_CONST` | **定通量**（Neumann）。固定溶质通量 [mg/m²/s]。 |
| 4 | `SUB_RT_BC_Cauchy_CONST` | **混合（Cauchy）**。浓度/通量组合条件。 |
| 5 | `SUB_RT_BC_Dirichlet_T` | **时变 Dirichlet**。浓度来自时间序列文件。 |
| 6 | `SUB_RT_BC_Neumann_T` | **时变 Neumann**。通量来自时间序列文件。 |
| 7 | `SUB_RT_BC_Cauchy_T` | **时变 Cauchy**。混合条件来自时间序列文件。 |
| 8 | `SUB_RT_BC_SWE` | **地表水交换**。把地下水溶质输运与坡面流浓度耦合。用于地表，在表层与地下域之间交换溶质。 |
| 9 | `SUB_RT_BC_FD` | **自由排水**（零梯度出流）。溶质以与内部格胞相同的浓度离开域。用于含水层底板。 |

#### 示例 — 顶面（SWE 耦合）+ 底面（自由排水）

```text
bccount : 2

// --- Top boundary: coupled with surface water ---
id : rttop
direction : 6                            // ZMINUS: upper/surface face in the groundwater convention
polygon : polygontop.input
spec_0_bctype : 8                        // NH4: surface-water exchange
spec_1_bctype : 8                        // NO3: surface-water exchange
spec_2_bctype : 8                        // DO: surface-water exchange
spec_3_bctype : 8                        // DOC: surface-water exchange
spec_4_bctype : 8                        // DON: surface-water exchange

// --- Bottom boundary: free drainage ---
id : rtbottom
direction : 5                            // ZPLUS: lower/bottom face in the groundwater convention
polygon : polygontop.input
spec_0_bctype : 9                        // NH4: free drainage
spec_1_bctype : 9
spec_2_bctype : 9
spec_3_bctype : 9
spec_4_bctype : 9
```

---

### 2.3 反应参数 reactiongw.input

**文件路径**：`<input_folder>/reactiongw.input`
**源代码**：`src/RTInitGW.h` → `readRTReactionFile()`、`src/RTStateGW.h`

包含地表水反应参数的全部内容，另加**吸附**、**容重**与**地下水专有**的复氧/释放参数。仅当 `ReactionModule : 1` 时需要。

反应实现中液相浓度用 mg/L、固相浓度用 mg/kg、分配系数用 L/kg，容重在内部以 kg/L 计。`Bulk_Density` 是本文件中唯一的密度值：按 g/cm3 填写。数值上 1 g/cm3 = 1 kg/L，因此解析后无需换算。大于 10 的取值会被当作疑似 kg/m3 输入而拒绝。

#### 反应质量基准与求解器前提

地下水反应内核使用固定五组分顺序与混合质量基准：索引 0 `NH4` 为 `mg NH4/L`、索引 1 `NO3` 为 `mg NO3/L`、索引 2 `DO` 为 `mg O2/L`、索引 3 `DOC` 为 `mg C/L`、索引 4 `DON` 为 `mg N/L`。它不是通用的氮质量模型，输入读取器不做 `mg NH4-N/L` 或 `mg NO3-N/L` 到化合物质量的换算。`SPECIES_ID` 与 `SPECIES_Name` 不重排数组。固定索引反应求解器要求 `n_mass >= 5`，否则直接返回、不求解氮循环反应。

地下水基准速率为 `R_min [mg N/L/s]`、`R_nit [mg NH4/L/s]`、`R_denit`/`R_aer [mg C/L/s]`。源项使用硬编码质量比 `18.04/14.01`、`62.00/18.04`、`64.00/18.04`、`32.00/14.01`、`32.00/12.01` 与 `(0.8*62.00)/12.01`。地下水侧还会把各反应源项除以该组分的阻滞因子 `Rf`。复氧为 `Ka_DO * max(0, theta_s-theta) * max(0, DO_sat-DO)`，故孔隙被水充满时复氧减弱。

**校验状态。** 必需反应键用 `-999` 检查存在性；多数物理范围不被解析器强制。下文范围是推荐率定范围。`Bulk_Density` 是例外：必须按 `g/cm³` 提供，大于 `10` 的取值会被当作疑似 `kg/m³` 输入拒绝。对激活的五组分运行，每个逐组分列表都要给满五个值，且 `Sorption_Type` 只能取 `0` 到 `4`。

#### 新增全局参数

| 键 | 类型 | 单位 | 常见范围 | 说明 |
|-----|------|-------|---------------|-------------|
| `Bulk_Density` | real | [g/cm³] | 1.1 – 1.4 | 多孔介质容重（ρb）。输入必须为 g/cm³；模型内部以 kg/L 存储和使用。典型稻田耕层：1.30 g/cm³。 |

#### 吸附参数（逐组分，分号分隔）

| 键 | 类型 | 单位 | 说明 |
|-----|------|-------|-------------|
| `Sorption_Type` | int 列表 | [-] | 各组分的吸附模型：`0` = 无，`1` = 线性（`Kd`），`2` = Freundlich（`Kf`、`Nf`），`3` = Langmuir（`Alpha_D`、`Beta_D`），`4` = 动力学/非平衡（`Kd`、`Beta`）。 |
| `Kd` | real 列表 | [L/kg] | **线性**吸附的分配系数。C 用 mg/L、S 用 mg/kg 时，单位取 L/kg。稻田 NH4 典型 2.0–5.0。非吸附组分填 `0.0`。 |
| `Kf` | real 列表 | [(mg/kg)/(mg/L)^Nf] | Freundlich 常数。不用 Freundlich 时填 `0.0`。 |
| `Nf` | real 列表 | [-] | Freundlich 非线性指数。典型 0.5–1.2。不用时填 `1.0`。 |
| `Alpha_D` | real 列表 | [L/mg] | Langmuir 亲和常数，用于 `S = Alpha_D*Beta_D*C/(1+Alpha_D*C)`。不用 Langmuir 时填 `0.0`。 |
| `Beta_D` | real 列表 | [mg/kg] | Langmuir 最大吸附容量。不用时填 `0.0`。 |
| `Beta` | real 列表 | [1/s] | `Sorption_Type=4` 使用的动力学吸附速率常数。它不是 Langmuir 参数。 |

#### 吸附模型汇总

| `Sorption_Type` | 模型 | 所需参数 | 等温线方程 |
|:---:|-------|---------------------|-------------------|
| 0 | 无吸附 | — | — |
| 1 | 线性 | `Kd` | S = Kd × C |
| 2 | Freundlich | `Kf`、`Nf` | S = Kf × C^Nf |
| 3 | Langmuir | `Alpha_D`、`Beta_D` | S = (Alpha_D × Beta_D × C) / (1 + Alpha_D × C) |
| 4 | 动力学/非平衡 | `Kd`、`Beta` | 在输运装配中实现；`Beta` 是一阶质量传递速率 [1/s] |

其中 S = 吸附浓度 [mg/kg]，C = 溶解浓度 [mg/L]。

#### 根系溶质吸收 — `TSCF`

`TSCF` 是 `reactiongw.input` 中必填的逐组分列表，名称意为**蒸腾流浓度因子**（Transpiration Stream Concentration Factor）。它控制随根系蒸腾流带出地下水域的溶解态地下水浓度的比例，逐组分独立施加，并使用固定组分顺序 `NH4; NO3; DO; DOC; DON`。

| 键 | 类型 | 单位 | 推荐物理范围 | 解析状态 | 实现含义 |
|---|---|---|---:|---|---|
| `TSCF` | real 列表 | [-] | `0.0–1.0` | 仅检查存在性；数值范围不强制 | 乘在根系蒸腾溶质通量上的逐组分比例因子。 |

对一个格胞与组分，激活的源汇计算概念上为：

```text
solute_outflow_root = |transpGW| * cell_area * C_dissolved * 1000 * TSCF   [mg/s]
```

其中 `transpGW < 0` 表示向上根系蒸腾，`C_dissolved` 是当前溶解浓度 `[mg/L]`，因子 `1000` 把 `L` 换算为 `m³`，`TSCF` 无量纲。同一溶质质量同时从当地反应输运方程中扣除。因此根系吸收是地下水源汇过程，不是额外的 Monod 反应项。

推荐取值解释：

- `TSCF = 0.0`：蒸腾只移走水、不移走该组分的质量。除非有植物吸收或输运数据支持，这是 DO、DOC、DON 的合适中性值。
- `0.0 < TSCF < 1.0`：蒸腾流只携带当地溶解浓度的一部分。用于部分溶质排斥或组分依赖的吸收。
- `TSCF = 1.0`：蒸腾水以当地溶解浓度携带溶质；这是被动平流解释的常用上界。
- `TSCF > 1.0`：当前解析器数学上接受，会使蒸腾流浓度富集，但只应在有明确实验证据时使用并作为富集因子记录。
- `TSCF < 0.0`：解析器不拒绝，但激活的计算只在 `TSCF > 0` 时施加吸收，负值等效于停用；不要使用负值。

对五组分氮循环算例，一个保守的起始列表是：

```text
TSCF : 0.0; 0.0; 0.0; 0.0; 0.0   // NH4, NO3, DO, DOC, DON
```

若要显式表达植物吸氮，敏感性算例可给 NH4 与 NO3 赋非零值，例如 `TSCF : 0.3; 0.3; 0.0; 0.0; 0.0`。这是率定假设，不是通用水稻值。不要仅因为 WOFOST 计算根系吸水就推断 TSCF 非零：水吸收与溶质吸收是两回事，土壤蒸发实现为纯水量损失（`TSCF` 对蒸发无效）。

当前代码不强制物理 `0–1` 范围，用户须自行校验列表。每个已配置组分恰好一个值；对激活的五组分反应内核给五个值。相关实现在 `src/RTInitGW.h`、`src/RTStateGW.h`、`src/SourceSink.h` 与 `src/RTIntegrator.h`。

#### 地下水专有源汇参数

| 键 | 类型 | 单位 | 常见范围 | 说明 |
|-----|------|-------|---------------|-------------|
| `Ka_DO` | real | [1/s] | 0 – 1×10⁻⁵ | 地下水复氧速率常数。典型 1×10⁻⁵（≈ 0.5/3600）。 |
| `DOC_eq` | real | [mg C/L] | 2 – 20 | 沉积物释放的 DOC 平衡浓度。 |
| `DON_eq` | real | [mg N/L] | 0.1 – 5.0 | 沉积物释放的 DON 平衡浓度。 |
| `K_rel_DOC` | real | [1/s] | 0 – 1×10⁻⁵ | DOC 释放速率常数。填 `0.0` 关闭。 |
| `K_rel_DON` | real | [1/s] | 0 – 1×10⁻⁵ | DON 释放速率常数。填 `0.0` 关闭。 |

#### 地下水反应参数

以下参数按标量读取，`ReactModel_Type : 1` 时必填。推荐范围仅为率定指引；当前解析器只检查键存在性、不检查这些数值界限。

| 键 | 类型 | 单位 | 推荐范围 | 实现含义 |
|---|---|---|---:|---|
| `Rate_Max_Nit` | real | [mg NH4/L/s] | 1×10⁻⁷ 至 1×10⁻⁴ | 最大硝化基准速率；直接移除 NH4 化合物质量。 |
| `K_Monod_NH4` | real | [mg NH4/L] | 0.01 至 5.0 | NH4 半饱和浓度。 |
| `K_Monod_DO` | real | [mg O2/L] | 0.01 至 5.0 | 硝化的 DO 半饱和浓度。 |
| `Rate_Max_Denit` | real | [mg C/L/s] | 1×10⁻⁷ 至 3×10⁻⁴ | 最大反硝化基准速率；直接移除 DOC-C。 |
| `K_Monod_NO3` | real | [mg NO3/L] | 0.01 至 5.0 | 化合物质量基准的 NO3 半饱和浓度。 |
| `K_Monod_DOC` | real | [mg C/L] | 0.01 至 20.0 | 反硝化的 DOC 半饱和浓度。 |
| `Ki_Inhib_DO` | real | [mg O2/L] | 0.01 至 5.0 | `Ki/(Ki+DO)` 中的 DO 抑制浓度。 |
| `Rate_Max_Hetero` | real | [mg C/L/s] | 1×10⁻⁷ 至 3×10⁻⁴ | 最大有氧 DOC 呼吸速率。 |
| `K_Monod_DOC_Aerobic` | real | [mg C/L] | 0.01 至 20.0 | 有氧呼吸的 DOC 半饱和浓度。 |
| `K_Monod_DO_Aerobic` | real | [mg O2/L] | 0.01 至 5.0 | 有氧呼吸的 DO 半饱和浓度。 |
| `Rate_Max_Min` | real | [mg N/L/s] | 1×10⁻⁸ 至 1×10⁻⁴ | 最大矿化基准速率；直接移除 DON-N。 |
| `K_Monod_DON` | real | [mg N/L] | 0.01 至 10.0 | 元素氮基准的 DON 半饱和浓度。 |
| `K_Monod_DO_Min` | real | [mg O2/L] | 0.01 至 5.0 | 矿化的 DO 半饱和浓度。 |
| `Ka_DO` | real | [1/s] | 0 至 1×10⁻⁴ | 地下水复氧系数。源项额外乘以充气孔隙度。 |
| `DO_sat` | real | [mg/L] | 0 至 15 | 复氧使用的 DO 饱和浓度。 |
| `DOC_eq` | real | [mg C/L] | 0 至 50 | 向液相释放的 DOC 平衡浓度。 |
| `DON_eq` | real | [mg N/L] | 0 至 10 | 向液相释放的 DON 平衡浓度。 |
| `K_rel_DOC` | real | [1/s] | 0 至 1×10⁻⁴ | DOC 释放系数；为零关闭该释放源。 |
| `K_rel_DON` | real | [1/s] | 0 至 1×10⁻⁴ | DON 释放系数；为零关闭该释放源。 |
| `Temp_Coeff_Theta` | real | [-] | 1.00 至 1.10 | 被解析的环境温度系数，因地下水内核 `Temp_Corr` 固定为 1.0，当前不生效。 |
| `Opt_Temp` | real | [°C] | 0 至 40 | 被解析的参考温度，当前在地下水内核中不生效。 |
| `Opt_pH` | real | [-] | 0 至 14 | 被解析且必填，但当前反应方程未使用。 |

`Lambda_1` 与 `Lambda_2` 是逐组分一阶衰减系数 [1/s]，推荐范围 `0` 至 `1×10⁻⁵`；它们作用在输运装配中而非五组分反应源方程里。`Lambda_1` 作用于溶解浓度，`Lambda_2` 作用于吸附浓度。取零关闭该衰减途径。

#### 示例

```text
# === Global ===
ReactModel_Type           : 1
Bulk_Density              : 1.30                // input [g/cm^3], internal kg/L
Nitrogen_Cycle_Simulation : 1
SPECIES_ID   : 0; 1; 2; 3; 4
SPECIES_Name : NH4; NO3; DO; DOC; DON

# === Sorption ===
Sorption_Type : 1; 0; 0; 0; 1          // NH4: linear, DON: linear
Kd            : 3.0; 0.0; 0.0; 0.0; 0.5  // [L/kg]
Kf            : 0.0; 0.0; 0.0; 0.0; 0.0
Nf            : 1.0; 1.0; 1.0; 1.0; 1.0
Alpha_D       : 0.0; 0.0; 0.0; 0.0; 0.0
Beta_D        : 0.0; 0.0; 0.0; 0.0; 0.0
Beta          : 0.0; 0.0; 0.0; 0.0; 0.0

# === Decay ===
Lambda_1 : 0.0; 0.0; 0.0; 0.0; 0.0
Lambda_2 : 0.0; 0.0; 0.0; 0.0; 0.0

# === Nitrification ===
Rate_Max_Nit  : 5.0e-6                 // [mg NH4/L/s]
K_Monod_NH4   : 0.5                    // [mg NH4/L]
K_Monod_DO    : 0.5                    // [mg O2/L]

# === Denitrification ===
Rate_Max_Denit : 5.0e-5                // [mg C/L/s]
K_Monod_NO3    : 0.5                   // [mg NO3/L]
K_Monod_DOC    : 2.0                   // [mg C/L]
Ki_Inhib_DO    : 0.2                   // [mg O2/L]

# === Aerobic Respiration ===
Rate_Max_Hetero       : 1.0e-4         // [mg C/L/s]
K_Monod_DOC_Aerobic   : 5.0            // [mg C/L]
K_Monod_DO_Aerobic    : 0.5            // [mg O2/L]

# === Mineralization ===
Rate_Max_Min  : 3.0e-6                // [mg N/L/s]
K_Monod_DON   : 0.21015               // [mg N/L]
K_Monod_DO_Min: 0.48                  // [mg O2/L]

# === Reaeration / Release ===
Ka_DO      : 1.0e-5                    // [1/s]
DO_sat     : 8.1                       // [mg/L]
DOC_eq     : 20.0                      // [mg C/L]
DON_eq     : 2.0                       // [mg N/L]
K_rel_DOC  : 0.0                       // [1/s]
K_rel_DON  : 0.0                       // [1/s]

# === Environmental Correction ===
Temp_Coeff_Theta : 1.047
Opt_Temp         : 25.0
Opt_pH           : 7.0
```

---

### 2.4 地下水源汇 rtgwss.input

**文件路径**：`<input_folder>/rtgwss.input`
**源代码**：`src/RTInitGW.h`

定义地下水域中的溶质源汇区。格式与 `rtgwbc.input` 类似，采用基于多边形的空间选择。

#### 文件格式

```text
sscount : <number_of_sources>

id : <source_name>
rtsstype : <type>                // Global source/sink type
polygon : <polygon_file>         // Spatial extent polygon

// Per-species values or time series
spec_X_sstype : <type>
spec_X_ssval  : <value>
spec_X_ssfile : <timeseries_file>
```

#### 源汇类型（`sstype`）

| 取值 | 说明 |
|:---:|-------------|
| `rtsstype=1` | 在所选地下水格胞中的体积质量注入 [mg/L/s]。 |
| `rtsstype=3` | 面积质量通量注入 [mg/m²/s]，沿所选垂向柱厚度分配。 |

`rtsstype=0` 属于另一套 WOFOST 水源汇结构，不属于反应溶质 `RTGwSS` 计算。`rtsstype=2` 不被当前反应溶质注入分支处理，请勿在此使用。

#### 示例 — 用于 WOFOST 耦合的根系区

```text
sscount : 1
id : root1
rtsstype : 1                     // Reactive-solute volumetric injection
polygon : polygonroot.input      // Root zone polygon
spec_0_sstype : 0
spec_0_ssval : 0.0                // [mg/L/s]
```

根系吸水与土壤蒸发由地下水水量平衡状态与 WOFOST 耦合提供；不通过本反应溶质文件中的 `rtsstype=0` 配置。

---

### 2.5 初始条件文件

**源代码**：`src/RTInitSW.h` → `readRtICFileSW()`、`src/RTInitGW.h` → `readRtICFileGW()`

初始浓度场可用空间栅格给出。地表水接受灵活的 ESRI 风格表头解析器；地下水使用下述三记录标签式表头。

#### 文件命名约定

| 域 | 相 | 模式 | 示例 |
|--------|-------|---------|---------|
| 地表水 | 液相（水相） | `rtinitialsw<Name>.input` | `rtinitialswNH4.input` |
| 地下水 | 液相（水相） | `rtinitialgw<Name>.input` | `rtinitialgwNO3.input` |
| 地下水 | 固相（吸附） | `rtinitialgwsolid<Name>.input` | `rtinitialgwsolidNH4.input` |

#### 地下水格式：三条表头记录 + 三维场

```text
ncols <ncols>
nrows <nrows>
nodata_value <value>
<data values...>     // ncols × nrows × nz values, layer-major/row-major
```

地下水读取器把前三条表头记录当作"标签 + 数字"消费。请使用上述标签形式。`ncols` 与 `nrows` 必须等于全局水平域尺寸；数据段必须含 `nx_glob * ny_glob * nz_glob` 个值。没有 `nz` 表头字段。

#### 地表水格式：简化或 ESRI 风格二维栅格

```text
NCOLS <ncols>
NROWS <nrows>
XLLCORNER <x>         // Ignored by the model, logged only
YLLCORNER <y>         // Ignored by the model, logged only
CELLSIZE <size>       // Ignored by the model, logged only
NODATA_value <value>
<data values...>      // Row-major, nrows × ncols values
```

#### 要求

- 栅格尺寸必须与全局域**精确一致**：`NCOLS = nx_glob`、`NROWS = ny_glob`。
- 数据按行主序排列。
- Nodata 值被当作数值浓度读入；当前地下水读取器不做掩膜。
- 单位：液相 [mg/L]，固相 [mg/kg]。

---

## 3. WOFOST 作物生长模型

### 3.1 作物参数 cropparameter.input

**文件路径**：`<input_folder>/cropparameter.input`
**源代码**：`src/cropsrc/CropInit.h` → `readCropParameter()`、`src/cropsrc/CropState.h`

本文件包含全部 WOFOST 作物生长参数以及根系层土壤物理属性。文件采用冒号分隔的键–值对。表格参数用空格分隔的 x,y 数对给出。

#### 元数据

| 键 | 类型 | 说明 |
|-----|------|-------------|
| `CRPNAM` | string | 作物名称（信息性）。 |

`CROP_NO` 不被当前 `CropInit` 读取器解析。作物名称元数据只用 `CROPNAM`；所有数值作物参数使用下文列出的键。

#### 出苗参数

| 键 | 类型 | 单位 | 说明 |
|-----|------|-------|-------------|
| `TBASEM` | real | [°C] | 出苗下限温度。低于此值不出苗。 |
| `TEFFMX` | real | [°C] | 出苗有效上限温度。高于此值不再累积有效温度。 |
| `TSUMEM` | real | [°C·d] | 播种到出苗的温度和。 |

#### 物候参数

| 键 | 类型 | 单位 | 说明 |
|-----|------|-------|-------------|
| `IDSL` | int | [-] | 发育触发模式：`0` = 仅温度，`1` = 仅日长，`2` = 温度 + 日长。 |
| `DLO` | real | [hr] | 发育最适日长。`IDSL = 0` 时填 `-99`。 |
| `DLC` | real | [hr] | 临界日长（下限）。`IDSL = 0` 时填 `-99`。 |
| `TSUM1` | real | [°C·d] | 出苗到开花的温度和。 |
| `TSUM2` | real | [°C·d] | 开花到成熟的温度和。 |
| `DVSI` | real | [-] | 模拟起始时的初始发育阶段（DVS）。`0.0` = 出苗，`1.0` = 开花，`2.0` = 成熟。 |
| `DVSEND` | real | [-] | 生理成熟/收获时的发育阶段（通常 `2.0`）。 |

#### 初始条件

| 键 | 类型 | 单位 | 常见范围 | 说明 |
|-----|------|-------|---------------|-------------|
| `TDWI` | real | [kg/ha] | 10 – 200 | 出苗时初始总干物重。 |
| `LAIEM` | real | [ha/ha] | 0.0001 – 0.01 | 出苗时叶面积指数。 |
| `RGRLAI` | real | [ha/ha/d] | 0.001 – 0.05 | 指数生长期叶面积指数最大相对增量。 |

#### 绿面积参数

| 键 | 类型 | 单位 | 说明 |
|-----|------|-------|-------------|
| `SPAN` | real | [d] | 35 °C 下生长叶片的寿命。控制叶片衰老速率。 |
| `TBASE` | real | [°C] | 叶片老化的下限温度。 |
| `SPA` | real | [ha/kg] | 荚面积系数。无荚作物填 `0.0`。 |

#### 表格参数

所有表格用空格分隔的 x,y 数对给出。x 变量通常为发育阶段（DVS）或温度。模型在数据点间线性插值。

| 表 | x 变量 | y 变量 | 说明 |
|-------|-----------|-----------|-------------|
| `DTSMTB` | 平均温度 [°C] | 温度和日增量 [°C/d] | 定义物候的有效温度函数。 |
| `SLATB` | DVS [-] | 比叶面积 [ha/kg] | 单位叶干重的叶面积。控制 LAI 增长。 |
| `SSATB` | DVS [-] | 比茎面积 [ha/kg] | 单位茎干重的茎面积。常取 0。 |
| `KDIFTB` | DVS [-] | 漫射可见光消光系数 [-] | 控制冠层光截获。 |
| `EFFTB` | 日均温 [°C] | 单叶光能利用效率 [kg CO₂ ha⁻¹ hr⁻¹ / J m⁻² s⁻¹] | |
| `AMAXTB` | DVS [-] | 叶片最大 CO₂ 同化速率 [kg CO₂ ha⁻¹ hr⁻¹] | |
| `TMPFTB` | 平均温度 [°C] | AMAX 修正因子 [-] | 温度对光合的作用。范围 0–1。 |
| `TMNFTB` | 最低温度 [°C] | 总同化修正因子 [-] | 霜冻/损伤效应。范围 0–1。 |
| `RFSETB` | DVS [-] | 衰老修正因子 [-] | |

#### 转化效率参数

同化产物（CH₂O）向植物器官的干物质转化效率。

| 键 | 类型 | 单位 | 常见范围 | 说明 |
|-----|------|-------|---------------|-------------|
| `CVL` | real | [kg/kg] | 0.60 – 0.80 | 向**叶**的转化效率。 |
| `CVO` | real | [kg/kg] | 0.60 – 0.80 | 向**储存器官**（谷粒、块茎）的转化效率。 |
| `CVR` | real | [kg/kg] | 0.60 – 0.80 | 向**根**的转化效率。 |
| `CVS` | real | [kg/kg] | 0.60 – 0.80 | 向**茎**的转化效率。 |

#### 维持呼吸参数

| 键 | 类型 | 单位 | 说明 |
|-----|------|-------|-------------|
| `Q10` | real | [-] | 温度每升 10 °C 呼吸速率的相对增量。典型 2.0。 |
| `RML` | real | [kg CH₂O/kg/d] | **叶**的相对维持呼吸速率。 |
| `RMO` | real | [kg CH₂O/kg/d] | **储存器官**的相对维持呼吸速率。 |
| `RMR` | real | [kg CH₂O/kg/d] | **根**的相对维持呼吸速率。 |
| `RMS` | real | [kg CH₂O/kg/d] | **茎**的相对维持呼吸速率。 |

#### 分配表

定义总干物（按发育阶段）在各器官间的分配。

| 表 | 说明 | 约束 |
|-------|-------------|------------|
| `FRTB` | **总**干物向根的分配比例随 DVS 变化。 | 独立；余量进入地上部。 |
| `FLTB` | **地上**干物向叶的分配比例随 DVS 变化。 | 所有 DVS 处 FL + FS + FO = 1.0。 |
| `FSTB` | **地上**干物向茎的分配比例随 DVS 变化。 | |
| `FOTB` | **地上**干物向储存器官的分配比例随 DVS 变化。 | |

#### 死亡速率参数

| 键/表 | 类型 | 单位 | 说明 |
|-----------|------|-------|-------------|
| `PERDL` | real | [d⁻¹] | 水分胁迫导致的叶片最大相对死亡速率。 |
| `RDRRTB` | 表 | — | 根死亡速率随 DVS 变化 [d⁻¹]。 |
| `RDRSTB` | 表 | — | 茎死亡速率随 DVS 变化 [d⁻¹]。 |

#### 用水参数

| 键 | 类型 | 单位 | 说明 |
|-----|------|-------|-------------|
| `CFET` | real | [-] | 蒸腾速率修正因子。典型 `1.0`。 |
| `DEPNR` | real | [-] | 土壤水分消耗的作物组号（1–5）。值越大越耐旱。 |
| `IAIRDU` | int | [-] | 根系通气组织：`1` = 有（水稻等湿生作物），`0` = 无。 |
| `IOX` | int | [-] | 氧胁迫开关：`0` = 无氧胁迫，`1` = 氧胁迫激活。 |

#### 根系参数

| 键 | 类型 | 单位 | 常见范围 | 说明 |
|-----|------|-------|---------------|-------------|
| `RDI` | real | [cm] | 5 – 20 | 出苗时初始根深。 |
| `RRI` | real | [cm/d] | 0.5 – 2.0 | 根深最大日增量。 |
| `RDMCR` | real | [cm] | 30 – 150 | 最大根深（因作物而异）。水稻约 80 cm；小麦约 120 cm。 |

#### 养分参数（NPK）

##### 浓度界限

| 键 | 类型 | 单位 | 说明 |
|-----|------|-------|-------------|
| `NMINSO` | real | [kg N/kg] | 储存器官最低（残余）氮浓度。 |
| `NMINVE` | real | [kg N/kg] | 营养器官最低（残余）氮浓度。 |
| `NMAXSO` | real | [kg N/kg] | 储存器官最大氮浓度。 |
| `NMAXVE` | real | [kg N/kg] | DVS=0 时叶片最大氮浓度。 |
| `PMINSO` | real | [kg P/kg] | 储存器官最低（残余）磷浓度。 |
| `PMINVE` | real | [kg P/kg] | 营养器官最低（残余）磷浓度。 |
| `PMAXSO` | real | [kg P/kg] | 储存器官最大磷浓度。 |
| `PMAXVE` | real | [kg P/kg] | DVS=0 时叶片最大磷浓度。 |
| `KMINSO` | real | [kg K/kg] | 储存器官最低（残余）钾浓度。 |
| `KMINVE` | real | [kg K/kg] | 营养器官最低（残余）钾浓度。 |
| `KMAXSO` | real | [kg K/kg] | 储存器官最大钾浓度。 |
| `KMAXVE` | real | [kg K/kg] | DVS=0 时叶片最大钾浓度。 |

##### 吸收与转移

| 键 | 类型 | 单位 | 说明 |
|-----|------|-------|-------------|
| `RNUPTAKEMAX` | real | [kg N/ha/d] | 最大日吸氮速率。 |
| `RPUPTAKEMAX` | real | [kg P/ha/d] | 最大日吸磷速率。 |
| `RKUPTAKEMAX` | real | [kg K/ha/d] | 最大日吸钾速率。 |
| `DVS_NPK_STOP` | real | [-] | 高于此 DVS 后作物停止 NPK 吸收。 |
| `DVS_NPK_TRANSL` | real | [-] | 高于此 DVS 后 NPK 由营养器官向储存器官转移。 |

##### NPK 胁迫与分配

| 键 | 类型 | 单位 | 说明 |
|-----|------|-------|-------------|
| `NLAI_NPK` | real | [-] | NPK 胁迫下 LAI 降幅系数。 |
| `NSLA_NPK` | real | [-] | NPK 胁迫下 SLA 降幅系数。 |
| `NPART` | real | [-] | 氮胁迫对叶片分配影响的系数。 |
| `NLUE_NPK` | real | [-] | NPK 胁迫下 CO₂ 同化降幅系数。 |
| `NPK_TRANSLRT_FR` | real | [-] | 从根转移出的 NPK 比例。 |
| `RDRLV_NPK` | real | [d⁻¹] | NPK 胁迫导致的叶片最大相对死亡速率。 |

##### 器官浓度比例

| 键 | 说明 |
|-------------|-------------|
| `NMAXRT_FR` | 根最大氮浓度占叶最大氮浓度的比例。 |
| `NMAXST_FR` | 茎最大氮浓度占叶最大氮浓度的比例。 |
| `NCRIT_FR` | 临界氮浓度占最大氮浓度的比例。 |
| `NRESIDRT` | 根残余氮比例 [kg N/kg]。 |
| `PMAXRT_FR`、`PMAXST_FR`、`PCRIT_FR`、`PRESIDRT` | 磷的同上各量。 |
| `KMAXRT_FR`、`KMAXST_FR`、`KCRIT_FR`、`KRESIDRT` | 钾的同上各量。 |

##### 转移时间系数

| 键 | 类型 | 单位 | 说明 |
|-----|------|-------|-------------|
| `TCNT` | int | [d] | 氮转移时间系数。 |
| `TCPT` | int | [d] | 磷转移时间系数。 |
| `TCKT` | int | [d] | 钾转移时间系数。 |

##### NPK 浓度表

| 表 | x 变量 | y 变量 | 说明 |
|-------|-----------|-----------|-------------|
| `NMAXLV_TB` | DVS [-] | 叶片最大氮浓度 [kg N/kg] | |
| `PMAXLV_TB` | DVS [-] | 叶片最大磷浓度 [kg P/kg] | |
| `KMAXLV_TB` | DVS [-] | 叶片最大钾浓度 [kg K/kg] | |

##### 其他

| 键 | 类型 | 单位 | 说明 |
|-----|------|-------|-------------|
| `YZERO` | real | [kg/ha] | 零产量时营养器官的最大量。 |
| `NFIX` | real | [kg/kg] | 来自生物固氮（如豆科）的吸氮比例。 |

#### 根系层土壤物理属性

这些参数定义 WOFOST 用于水分胁迫计算的土壤持水属性。

| 键 | 类型 | 单位 | 常见范围 | 说明 |
|-----|------|-------|---------------|-------------|
| `SMW` | real | [cm³/cm³] | 0.05 – 0.15 | **凋萎点**土壤含水量（ψ ≈ −1500 cm 或 −150 kPa）。 |
| `SMFCF` | real | [cm³/cm³] | 0.20 – 0.40 | **田间持水量**土壤含水量（ψ ≈ −100 至 −340 cm）。南方粘土 0.35–0.40。 |
| `SM0` | real | [cm³/cm³] | 0.35 – 0.50 | **饱和**含水量（≈ 孔隙度 Φ）。 |
| `CRAIRC` | real | [cm³/cm³] | 0.04 – 0.10 | 通气**临界含气量**。低于此值时氧胁迫限制根系功能。 |

#### 示例

```text
// --- Meta ---
CRPNAM : Gaoyou Rice (European Rice ecotype)

// --- Emergence ---
TBASEM : 8.0          // [°C]
TEFFMX : 22.0         // [°C]
TSUMEM : 100.0        // [°C·d]

// --- Phenology ---
IDSL   : 0            // Temperature-driven only
TSUM1  : 1490.0       // [°C·d]
TSUM2  : 610.0        // [°C·d]
DVSI   : 0.33         // Initial DVS
DVSEND : 2.00

DTSMTB : 0.00 0.00 10.00 0.00 25.00 15.00

// --- Initial ---
TDWI   : 130.0        // [kg/ha]
LAIEM  : 0.0007       // [ha/ha]
RGRLAI : 0.0070       // [ha/ha/d]

// --- Green Area ---
SPAN   : 50.0         // [d]
TBASE  : 0.0          // [°C]
SPA    : 0.000

SLATB  : 0.00 0.0017 0.10 0.0017 0.58 0.0019 0.86 0.00200 1.10 0.0015 1.40 0.0015 2.00 0.0010
SSATB  : 0.00 0.0000 2.00 0.0000

// --- Assimilation ---
KDIFTB : 0.00 0.600 2.00 0.600
EFFTB  : 0.0 0.450 40.0 0.450
AMAXTB : 0.00 40.00 1.00 40.00 1.30 40.00 2.00 40.00
TMPFTB : 0.0 0.010 10.0 0.000 25.0 1.000 35.0 1.000 42.0 0.000
TMNFTB : 0.0 0.00 3.0 1.00

// --- Conversion ---
CVL : 0.754
CVO : 0.684
CVR : 0.754
CVS : 0.754

// --- Respiration ---
Q10  : 2.0
RML  : 0.0300
RMO  : 0.0150
RMR  : 0.0150
RMS  : 0.0200
RFSETB : 0.00 0.300 2.00 0.400

// --- Partitioning ---
FRTB : 0.000 0.200 0.500 0.200 0.800 0.150 1.000 0.150 1.100 0.000 2.000 0.000
FLTB : 0.000 0.400 0.580 0.520 0.800 0.380 0.850 0.000 0.900 0.000 1.000 0.000 2.000 0.000
FSTB : 0.000 0.600 0.580 0.480 0.800 0.620 0.850 0.400 0.900 0.000 1.000 0.000 2.000 0.000
FOTB : 0.000 0.000 0.580 0.000 0.800 0.000 0.850 0.600 0.900 1.000 1.000 1.000 2.000 1.000

// --- Death ---
PERDL  : 0.030
RDRRTB : 0.00 0.000 1.50 0.000 1.51 0.020 2.00 0.020
RDRSTB : 0.00 0.000 1.50 0.000 1.51 0.020 2.00 0.020

// --- Water ---
CFET   : 1.00
DEPNR  : 3.5
IAIRDU : 1
IOX    : 0

// --- Rooting ---
RDI    : 10.0         // [cm]
RRI    : 1.2          // [cm/d]
RDMCR  : 40.0         // [cm]

// --- Nutrients ---
NMINSO : 0.0043
NMINVE : 0.0021
NMAXSO : 0.0144
NMAXVE : 0.0600
// ... (all nutrient parameters) ...
RNUPTAKEMAX : 7.2    // [kg N/ha/d]
RPUPTAKEMAX : 0.8    // [kg P/ha/d]
RKUPTAKEMAX : 7.4    // [kg K/ha/d]

// --- Soil Physics ---
SMW    : 0.085        // [cm^3/cm^3]
SMFCF  : 0.27
SM0    : 0.46
CRAIRC : 0.19
```

---

### 3.2 气象数据 cropmeteo.input

**文件路径**：`<input_folder>/cropmeteo.input`
**源代码**：`src/cropsrc/MeteoState.h`、`src/cropsrc/MeteoInit.h`

为 WOFOST 作物模型提供逐日气象强迫。文件是空白分隔的数值表，至少九列；它不是 YAML，也不使用第 5 节的通用 `time value` 格式。

文件还必须在注释或纯文本行中包含位置元数据：`LAT`、`LON`、`ELEV`。可选的 `ANGSTA` 与 `ANGSTB` 覆盖默认 Angstrom 系数。

#### 格式

每时间步一行的时间序列。每行包含：

| 列 | 类型 | 单位 | 说明 |
|--------|------|-------|-------------|
| 1 | real | [d 或模型时间单位] | 时间/索引，存入 `MeteoState::time`；读取器不做换算。 |
| 2 | int/real | [-] | 日历年份，如 `2004`。 |
| 3 | int/real | [年内日序] | DOY，通常 1–366。 |
| 4 | real | [J/m²/d] | 逐日总太阳辐射（`IRRAD`）。 |
| 5 | real | [°C] | 逐日最低温度（`TMIN`）。 |
| 6 | real | [°C] | 逐日最高温度（`TMAX`）。 |
| 7 | real | [hPa] | 水汽压（`VAP`），不是 kPa。 |
| 8 | real | [m/s] | 2 m 高风速（`WIND`）。 |
| 9 | real | [mm/d] | 逐日降水（`RAIN`）。 |

读取器要求至少九个数值列。计算得到的开阔水面、裸土与参考作物蒸散分别以 `E0`、`ES0`、`ET0` 存储，内部计算结果由 mm/d 换算为 [cm/d]。

---

### 3.3 农业管理 agro.input

**文件路径**：`<input_folder>/agro.input`
**格式**：YAML

定义作物品种选择与播/收日程。

#### 键字段

| 键 | 类型 | 说明 |
|-----|------|-------------|
| `crop` | string | 作物类型名（如 `sugar-beet`、`rice`）。 |
| `variety` | string | 品种标识，须与作物 YAML 文件中的小节匹配。 |
| `crop_start_type` | string | 起始类型：`sowing`（播种）或 `emergence`（出苗）。 |
| `crop_start_date` | string | 作物起始日期（YYYY-MM-DD）。 |
| `crop_end_type` | string | 结束类型：`harvest`（收获）或 `maturity`（成熟）。 |
| `crop_end_date` | string | 作物结束日期（YYYY-MM-DD）。 |
| `max_duration` | int | 最长作物生长期 [d]。 |

#### 示例

```yaml
crop: sugar-beet
variety: sugar-beet-601
crop_start_type: emergence
crop_start_date: 2004-01-02
crop_end_type: harvest
crop_end_date: 2006-07-20
max_duration: 300
```

---

### 3.4 WOFOST-SERGHEI 耦合接口

**源代码**：`src/SourceSinkCrop.h`

WOFOST 作物生长模型通过根系层与 SERGHEI 地下水模型双向耦合。

#### 数据流

```
┌─────────────────────────────┐
│       SERGHEI（地下水）       │
│                             │
│  三维土壤含水量场             │──── 根系层土壤含水量 ────┐
│                             │    （[0, RD] 深度加权平均） │
└─────────────────────────────┘                          │
                                                         ▼
┌──────────────────────────────────────────────────────┐
│                  WOFOST 作物模型                       │
│                                                      │
│  输入：土壤含水量、气象数据                              │
│  输出：LAI、根深（RD）、生物量、DVS                      │
│                                                      │
│  内部过程：                                            │
│    - 物候发育（发育阶段）                                │
│    - CO₂ 同化与生长呼吸                                 │
│    - 生物量分配（叶、茎、根、储存器官）                    │
│    - 蒸腾与水分胁迫                                     │
│    - 养分（NPK）吸收与胁迫（当前未启用）                   │
└──────────────────────────────────────────────────────┘
        │                                    │
        │ LAI（叶面积指数）                    │ RD（根深）
        │                                    │
        ▼                                    ▼
┌──────────────────────────────────────────────────────┐
│       SERGHEI（地下水 + 地表水）                        │
│                                                      │
│  蒸散划分：                                            │
│    T_split = 1 - exp(-k_ext × LAI)    （蒸腾）         │
│    E_split = 1 - T_split              （土壤蒸发）      │
│                                                      │
│  根系吸水：                                            │
│    深度范围：[0, RD]                                   │
│    速率：由土壤含水量与潜在蒸散决定                       │
└──────────────────────────────────────────────────────┘
```

#### 耦合变量

| 变量 | 方向 | 单位 | 说明 | 源代码 |
|----------|-----------|-------|-------------|-------------|
| 根系层土壤含水量（RZSM） | 地下水 → WOFOST | [m³/m³] | [0, RD] 深度加权平均体积含水量。 | `SourceSinkCrop.h` 中的 `compute_root_zone_moisture()` |
| LAI | WOFOST → SERGHEI | [m²/m²] | 叶面积指数。范围 0–7。控制辐射截获。 | `SourceSinkCrop.h` 中的 `wofost_lai` |
| RD（根深） | WOFOST → SERGHEI | [cm] | 当前根深。范围 RDI 至 RDMCR。决定吸水深度。 | `SourceSinkCrop.h` 中的 `wofost_rd` |
| T_split | 计算 | [-] | 潜在蒸散中的蒸腾比例。 | `SourceSinkCrop.h` 中的 `compute_evap_partitioning()` |
| E_split | 计算 | [-] | 土壤蒸发比例。 | `E_split = 1 - T_split` |

#### 水分胁迫计算

WOFOST 用根系层土壤含水量计算水分胁迫因子：

```
RZSM = (Σ θ_i × Δz_i) / RD    （对 [0, RD] 内的层 i 求和）

胁迫因子取决于：
  - SMW（凋萎点）
  - SMFCF（田间持水量）
  - SM0（饱和含水量）
  - DEPNR（作物组号）
  - IAIRDU（通气组织，决定耐涝性）
```

- **干旱胁迫**：RZSM 接近 SMW 时，蒸腾被削减。
- **氧胁迫**：RZSM 接近 SM0（渍水）时根系功能受损（若 `IOX = 1`）。`IAIRDU = 1` 的水稻具有耐受性。

---

## 4. 触发式灌溉与排水边界

### 4.1 概述

**源代码**：`src/Parser.h` → `readExtBCFile()`（解析与校验）、`src/BC.h`（触发判定与水动力施加，`SWE_BC_TRIGGERED_IRRIGATION = 13`、`SWE_BC_TRIGGERED_DRAINAGE = 14`）、`src/TimeIntegrator.h`（每时间步调用 `ExtBC::apply()`）。

**参证测试工况**：`tests/t15-gaoyou-wofost-triggerIrrigation/`（水深触发灌溉）、`tests/t16-gaoyou-wofost-triggerIrrigation/`（W1 浅湿露干，含土壤水触发阶段）

`bctype : 13` 与 `bctype : 14` 实现面向稻田管理的、基于规则的自动控制地表水边界：

- **触发式灌溉**（`bctype : 13`，`SWE_BC_TRIGGERED_IRRIGATION`）：入流边界。管理策略表（按作物生育阶段）给出目标水深；当观测区平均水深低于该阶段下限时，按该阶段流量供水，直至达到上限。各阶段也可以（或同时）通过行级 `trig_type`/`theta_trig` 列以**根系层土壤水**触发（见[4.3](#43-灌溉策略文件)），例如浅湿露干灌溉制度。
- **触发式排水**（`bctype : 14`，`SWE_BC_TRIGGERED_DRAINAGE`）：出流边界。当观测区平均水深超过该阶段"排水开"阈值时排水——临界（自由）出流或固定泵流量——直至达到"排水关"阈值。

两类边界都在**每个水动力时间步**重新判定：从策略表选择阶段，把触发条件与用户定义观测区上的平均水深（以及土壤水）比较，据此施加或关闭边界流量。内置开/关滞回（见[4.3](#43-灌溉策略文件)与[4.4](#44-排水策略文件)）。

支持两种阶段选择模式：

- **DAT 模式**（`use_DVS : 0`）：按*移栽/播种后天数*选阶段，`DAT = floor((etime - crop_start_time)/86400) + 1`；当前代码中 `cropStartTime` 固定为 0，因此 `DAT` 实际等于自模拟起算的天数加一。
- **DVS 模式**（`use_DVS : 1`）：仅在 `CROP_GROWTH_MODEL` 构建中可用。按触发观测格胞上空间平均的 WOFOST 发育阶段（DVS）选阶段。平均 DVS 每个水动力步由 WOFOST 状态重新计算（启动时初始化为 `cropparameter.input` 的 `DVSI`，此后每个逐日作物步前进一次）。非作物构建中 `use_DVS : 1` 无法生效：首个时间步打印一次性警告（"DVS data are not connected …"）后灌溉/排水保持禁用（非致命错误）。注意阶段选择从 `DVSI` 起步：`DVSI : 0.33` 时 t=0 即匹配 0.10–0.35 的行，返青期行（0.00–0.10）永远不会被命中——请相应设置 `DVSI`。

### 4.2 边界定义 extbc.input

**文件路径**：`<input_folder>/extbc.input`
**源代码**：`src/Parser.h` → `readExtBCFile()`

#### 边界条件类型

| 取值 | 常量名 | 类别 | 说明 |
|:---:|------|:---:|-------------|
| 13 | `SWE_BC_TRIGGERED_IRRIGATION` | 入流 | 基于规则的灌溉入流口；水深（或土壤水）触发后施加流量。 |
| 14 | `SWE_BC_TRIGGERED_DRAINAGE` | 出流 | 基于规则的排水出口；触发后临界出流或定流量出流。 |

#### 逐边界键

| 键 | 类型 | 必填 | 说明 |
|-----|------|----------|-------------|
| `id` | string | **是** | 边界标识（必须位于该块其他键之前）。 |
| `bctype` | int | **是** | 触发式灌溉填 `13`，触发式排水填 `14`。 |
| `polygon` | string | **是** | 选择水进出域所在边界格胞的多边形文件（与其他水动力边界同格式；见[4.5](#45-触发多边形)）。 |
| `trigger_polygon` | string | 否* | 计算平均触发水深（DVS 模式下还有平均作物 DVS）的观测区。可为多边形文件名或关键字 `ALL_DOMAIN`。省略时使用该边界自身的 `polygon` 并打印警告。 |
| `direction` | real ×2 | **是** | 边界单位法向（内部归一化），如 `0.0 1.0`。灌溉沿 `+direction` 方向施加；排水沿外法向。 |
| `bcvals` | real ×1 | **是** | 恰好一个值：灌溉的**全局触发模式**（见[4.3](#43-灌溉策略文件)）；策略文件中的行级 `trig_type` 列可逐行覆盖。排水边界解析但不使用；为一致性填 `1`。 |
| `irrigation_strategy` | string | 若 bctype=13 | 灌溉策略表文件名（见[4.3](#43-灌溉策略文件)）。`bctype : 13` 必填，缺失即致命错误。 |
| `drainage_strategy` | string | 若 bctype=14 | 排水策略表文件名（见[4.4](#44-排水策略文件)）。`bctype : 14` 必填，缺失即致命错误。 |

#### 解析与校验规则

读取器执行以下校验（除特别注明外均为致命错误）：

- 触发式灌溉/排水的 `bcvals` 必须**恰好 1 个**值（触发模式）。
- `trigger_polygon` 可为 `ALL_DOMAIN`（不分大小写）或一个多边形文件：**1 个点**（最近格胞）、**≥ 3 个点**（多边形区域）或 0 个点（全部格胞）；其他情况拒绝。
- 策略文件必须定义 `use_DVS`（0 或 1）与正的 `ndata`，其后恰好 `ndata` 行数据；策略文件内的未知键会被拒绝。
- 策略行必须满足 DAT 与 DVS（及水深）范围的 `max ≥ min`，且各行在 DAT 与 DVS 上**有序且不重叠**。
- `irrigation_enabled : 1` 的灌溉行必须 `Qirr ≠ 0`；`drainage_enabled : 1` 且 `drain_mode : 2` 的排水行必须 `Qdrain ≠ 0`。
- 排水行要求 `drain_on_h ≥ drain_off_h` 且 `drain_mode ∈ {1, 2}`。

#### 示例 — `tests/t15-gaoyou-wofost-triggerIrrigation/extbc.input`

```text
bccount : 2

id : trigger_drainage
bctype : 14
polygon : polygonoutlet.input
trigger_polygon : ALL_DOMAIN
direction : 0.0 1.0
bcvals : 1
drainage_strategy : test2.input

id : trigger_irrigation
bctype : 13
polygon : polygoninlet.input
trigger_polygon : ALL_DOMAIN
direction : 0.0 1.0
bcvals : 1
irrigation_strategy : test.input

//bcvals(0)：全局触发模式，`1` 为只用地表水深，`2` 为只用根系层土壤水（策略行需提供 theta_trig），`3` 为任一满足触发，`4` 为同时满足触发；策略行的 trig_type 列可逐阶段覆盖
```

> **注**：触发模式 2–4 即[4.3](#43-灌溉策略文件)所述的**根系层土壤水触发**。要求 `CROP_GROWTH_MODEL` 与 `SERGHEI_RE_MODEL` 同时启用的构建，且有效模式为 2–4 的每个策略行都必须提供 `theta_trig`（否则启动报错）。未开启这些宏的构建中，模式 2–4 与任何 `theta_*` 列都会在启动时被拒绝。

### 4.3 灌溉策略文件

**引用自**：`extbc.input` 的 `irrigation_strategy` 键（文件必须与 `extbc.input` 同目录）
**源代码**：`src/Parser.h` → `readExtBCFile()`（读取）；`src/BC.h` → `selectIrrigationStage()`、`checkIrrigationTrigger()`、`applyTriggeredIrrigationStage()`

#### 文件格式

```text
use_DVS : <0|1>
ndata : <行数>
// dat_min dat_max dvs_min dvs_max h_min h_max gw_depth_trigger Qirr irrigation_enabled drainage_enabled [trig_type] [theta_trig] [theta_stop]
<dat_min> <dat_max> <dvs_min> <dvs_max> <h_min> <h_max> <gw_depth_trigger> <Qirr> <irrigation_enabled> <drainage_enabled> [<trig_type>] [<theta_trig>] [<theta_stop>]
... （共 ndata 行）
```

`//` 注释行与空行被跳过。说明列顺序的注释表头行是惯例、并非必需。

#### 键与列

| 键/列 | 类型 | 单位 | 说明 |
|-----|------|-------|-------------|
| `use_DVS` | int | - | 阶段选择模式：`0` = 按 DAT，`1` = 按平均作物 DVS（仅作物构建）。 |
| `ndata` | int | - | 其后策略行数。 |
| `dat_min`、`dat_max` | real | [d] | 该阶段 DAT 窗口（闭区间）。`use_DVS : 0` 时使用。 |
| `dvs_min`、`dvs_max` | real | [-] | 该阶段平均 DVS 窗口（闭区间）。`use_DVS : 1` 时使用。 |
| `h_min` | real | [m] | 水深触发：关闭状态下平均观测水深 < `h_min` 时灌溉**开启**。 |
| `h_max` | real | [m] | 水深停止：开启状态下平均观测水深 ≥ `h_max` 时灌溉**停止**。若 `h_max < h_min`，内部抬升为 `h_min`。 |
| `gw_depth_trigger` | real | [m] | 地下水埋深触发阈值（保留）。格式必填但**尚未接通**；填大哑值如 `999.0`。 |
| `Qirr` | real | [m³/s] | 触发开启期间施加的灌溉流量。`irrigation_enabled : 1` 时必须非零。 |
| `irrigation_enabled` | int (0/1) | - | 该阶段是否允许灌溉（`0` 禁止本阶段灌溉）。 |
| `drainage_enabled` | int (0/1) | - | 信息性列，记录配对排水策略在该阶段是否允许排水。灌溉边界本身**不使用**此列；排水边界由其自己的策略文件控制。 |
| `trig_type` | int (0–4) | 否（默认 0） | 追加在 `drainage_enabled` 之后的**行级触发类型**（旧 10 列文件继续有效）：`0` = 沿用全局 `bcvals(0)` 模式；`1` = 只用地表水深；`2` = 只用根系层土壤水；`3` = 任一满足即开；`4` = 须同时满足。 |
| `theta_trig` | real | 若模式 ≥ 2 | 土壤水触发阈值，以**相对饱和度分数**给出，∈(0,1]（如 `0.70` = 70% θ_sat）。触发格胞平均相对饱和度 S = Σ(wc·dz)/Σ(θ_s·dz) 低于该值时灌溉开启。有效模式为 2–4 时必填（否则启动报错）。 |
| `theta_stop` | real | 否（默认 −1） | 可选提前停阈值（相对饱和度分数，∈[theta_trig,1]）：S ≥ `theta_stop` 后灌溉停止（例如降雨已补足土壤墒情）。`-1` 关闭提前停；任何情况下 `h_max` 仍为硬上限。 |

#### 触发模式（`bcvals(0)` 与行级 `trig_type`）

某阶段的有效模式 = 行内给出的 `trig_type`，否则用全局 `bcvals(0)`：

| 取值 | 开启条件（关闭状态下） |
|:---:|---------|
| 1 | 平均水深 < `h_min`（地表水深触发）。 |
| 2 | 根系层相对饱和度 S < `theta_trig`（土壤水触发；该行必须提供 `theta_trig`）。 |
| 3 | 任一条件满足（水深**或**土壤水）。 |
| 4 | 两条件同时满足（水深**且**土壤水）。 |
| 其他 | 运行时警告拒绝；回退为水深触发。 |

开启状态下，**四种模式共用同一条持续规则**：`平均水深 < h_max` **且** S < `theta_stop` 同时成立时灌溉继续（土壤水子句仅在 `theta_stop > 0` 时生效）。`h_max` 恒为硬上限——即使土壤仍然干旱，任何阶段都不会越过它继续灌水。

S 在**每个水动力步**用触发格胞上的当前地下水含水量场 `wc` 重新计算：逐柱在 WOFOST 根系层 [0, RD] 内做深度加权，并以逐层饱和含水量归一（θ_s 取自 van Genuchten 土壤参数表 `vg.input`）。异步耦合下 S 相对地表水步最多滞后一个地下水子循环。模式 2–4 要求 `CROP_GROWTH_MODEL` 与 `SERGHEI_RE_MODEL` 同时启用的构建；其他构建在启动时拒绝此类行。

#### 施加的水动力学

触发开启期间，边界先沿边界断面把水面摊平（`flattenWaterSurface`，同 `SWE_BC_Q_T`），再把总流量 `Qirr` 分配到边界格胞（`distributeDischarge`）：

- 湿边界格胞按水深加权获得单宽流量，`hu = Q·(h/Σh)/ds · nx`，`hv = ... · ny`。
- 若整条边界干涸，先把各格胞抬升至临界水深 `hc = (q²/g)^(1/3)`（不低于 `hmin`），再施加单宽流量 `q = Q/(N·ds)`，避免启动时动量归零。
- 触发关闭（或阶段禁用/缺位）时，边界 `hu`、`hv` 置零。

#### 示例 — 灌溉策略（`test.input`）

```text
use_DVS : 0
ndata : 7
// dat_min dat_max dvs_min dvs_max h_min h_max gw_depth_trigger Qirr irrigation_enabled drainage_enabled
1 7 0.00 0.10 0.030 0.070 999.0 0.003 1 0
8 25 0.10 0.35 0.020 0.060 999.0 0.003 1 0
26 40 0.35 0.55 0.000 0.000 999.0 0.000 0 1
41 65 0.55 0.85 0.030 0.070 999.0 0.003 1 0
66 85 0.85 1.05 0.020 0.060 999.0 0.003 1 0
86 110 1.05 1.65 0.000 0.030 999.0 0.002 1 0
111 125 1.65 2.00 0.000 0.000 999.0 0.000 0 1
```

该示例在营养生长期（阶段 1–5）保持 2–7 cm 水层，抽穗前把目标水深降到 0–3 cm（阶段 6，Qirr 降为 0.002 m³/s），并在末期排水阶段（第 111–125 天，`irrigation_enabled : 0`）禁灌。

#### 示例 — W1 浅湿露干全生育期策略（`tests/t16-gaoyou-wofost-triggerIrrigation/test.input`）

```text
use_DVS : 0
ndata : 7
// dat_min dat_max dvs_min dvs_max h_min h_max gw_depth_trigger Qirr irrigation_enabled drainage_enabled trig_type theta_trig theta_stop
1 7 0.00 0.10 0.005 0.020 999.0 0.003 1 0 1 -1 -1
8 25 0.10 0.35 0.000 0.010 999.0 0.003 1 0 2 0.70 -1
26 40 0.35 0.55 0.000 0.000 999.0 0.000 0 1
41 65 0.55 0.85 0.010 0.020 999.0 0.003 1 0 1 -1 -1
66 85 0.85 1.05 0.005 0.015 999.0 0.003 1 0 1 -1 -1
86 110 1.05 1.65 0.005 0.010 999.0 0.002 1 0 1 -1 -1
111 125 1.65 2.00 0.000 0.000 999.0 0.000 0 1
```

逐阶段说明：返青期（第 1–7 天）用水深触发维持 5–20 mm 水层；**分蘖前期（第 8–25 天）切换为土壤水触发**（`trig_type : 2`、`theta_trig : 0.70`）——根系层相对饱和度 S 降到 70% θ_sat 以下时灌溉开启，灌至 `h_max` = 10 mm（模式 2 下 `h_min` 不参与开启判定，但 `h_max` 仍是硬停）；第 26–40 天为中期排水晒田阶段，禁灌；拔节孕穗、抽穗扬花与乳熟期回到浅水深触发（10→20、5→15、5→10 mm）；第 111–125 天为黄熟落干，禁灌。策略的排水一半（阶段 3 与 7 的排干）由配对的排水策略文件承担。三个可选列在任意行都可整体省略（上面禁灌行保持了旧 10 列形式）；土壤水触发行的 `theta_stop` 可把 `-1` 换成如 `0.95`，在降雨补足墒情后提前停灌。注意策略解析器**不剥离**行内 `//` 注释，因此行的最后一列之后不能再有任何内容。

---

### 4.4 排水策略文件

**引用自**：`extbc.input` 的 `drainage_strategy` 键（文件必须与 `extbc.input` 同目录）
**源代码**：`src/Parser.h` → `readExtBCFile()`（读取）；`src/BC.h` → `selectDrainageStage()`、`checkDrainageTrigger()`、`applyTriggeredDrainageStage()`

#### 文件格式

```text
use_DVS : <0|1>
ndata : <行数>
// dat_min dat_max dvs_min dvs_max drain_on_h drain_off_h drain_mode Qdrain drainage_enabled
<dat_min> <dat_max> <dvs_min> <dvs_max> <drain_on_h> <drain_off_h> <drain_mode> <Qdrain> <drainage_enabled>
... （共 ndata 行）
```

#### 键与列

| 键/列 | 类型 | 单位 | 说明 |
|-----|------|-------|-------------|
| `use_DVS` | int | - | 阶段选择模式：`0` = 按 DAT，`1` = 按平均作物 DVS（仅作物构建）。 |
| `ndata` | int | - | 其后策略行数。 |
| `dat_min`、`dat_max` | real | [d] | 该阶段 DAT 窗口（闭区间）。 |
| `dvs_min`、`dvs_max` | real | [-] | 该阶段平均 DVS 窗口（闭区间）。 |
| `drain_on_h` | real | [m] | 关闭状态下平均观测水深 > `drain_on_h` 时排水**开启**。必须 ≥ `drain_off_h`。 |
| `drain_off_h` | real | [m] | 开启状态下平均观测水深 ≤ `drain_off_h` 时排水**停止**。 |
| `drain_mode` | int | - | 出流模式：`1` = 临界/自由出流（每个 `h ≥ hmin` 的边界格胞以 Froude=1 泄流，`|u| = √(gh)`）；`2` = 固定泵流量 `Qdrain`，按水深比例分配到湿边界格胞。 |
| `Qdrain` | real | [m³/s] | `drain_mode : 2` 的泵流量。`drain_mode : 1` 不使用（填任意值，如 `0.005`）。`drainage_enabled : 1` **且** `drain_mode : 2` 时必须非零。 |
| `drainage_enabled` | int (0/1) | - | 该阶段是否允许排水。 |

触发关闭（或阶段禁用/缺位）时，边界 `hu`、`hv` 置零。

#### 示例 — 排水策略（`test2.input`）

```text
use_DVS : 0
ndata : 7
// dat_min dat_max dvs_min dvs_max drain_on_h drain_off_h drain_mode Qdrain drainage_enabled
1 7 0.00 0.10 0.090 0.040 1 0.005 1
8 25 0.10 0.35 0.090 0.040 1 0.005 1
26 40 0.35 0.55 0.005 0.000 2 0.005 1
41 65 0.55 0.85 0.080 0.060 1 0.005 1
66 85 0.85 1.05 0.080 0.040 1 0.005 1
86 110 1.05 1.65 0.040 0.010 2 0.005 1
111 125 1.65 2.00 0.005 0.000 2 0.005 1
```

该示例在早期阶段田面水深超过 8–9 cm 时即排水（临界出流），在中期晒田与末期用固定泵流量把田面排至近干（`drain_mode : 2`，`Qdrain = 0.005` m³/s 排至约 0 水深），并在抽穗前放宽阈值。

### 4.5 触发多边形

**引用自**：`extbc.input` 的 `trigger_polygon` 键
**源代码**：`src/Parser.h`（文件读取）；`src/BC.h` → `find_trigger_cells()`

计算平均水深（与平均 DVS）的观测区由 `trigger_polygon` 选定：

| 选项 | 格式 | 含义 |
|--------|--------|---------|
| `ALL_DOMAIN` | 关键字（不分大小写） | 域内全部非 nodata 格胞均为观测格胞。 |
| 多边形文件，≥ 3 个点 | `NPOINTS <n>` 随后 n 行 `x y` 顶点 | 中心位于多边形内的格胞。 |
| 多边形文件，1 个点 | `NPOINTS 1` 随后一行 `x y` | 中心距该点最近的那个格胞。 |

多边形文件采用与其他边界/源汇多边形相同的 `NPOINTS` 表头格式，例如：

```text
NPOINTS 4
27.875 -0.125
27.875 0.625
28.625 0.625
28.625 -0.125
```

若一个观测格胞都找不到，打印警告且触发永久关闭。若省略 `trigger_polygon`，该边界自身的 `polygon` 被复用为观测区（一次性警告）。

### 4.6 运行行为与输出

- **判定频率。** 触发式边界在每个水动力时间步重新判定（`TimeIntegrator.h` 中的 `ExtBC::apply`），先于 halo 交换；流量随后像其他边界一样被积分。
- **流量记账。** `bctype : 13` 归类为*入流*边界、`bctype : 14` 归类为*出流*边界，其通量出现在 `domainTimeSeries.out` 的 `BoundaryInflow`/`BoundaryOutflow` 以及逐边界 `BoundaryFlow_<i>` / `BoundaryAccum_<i>` 列（下标 `i` = `extbc.input` 中的边界顺序，0 起；入流记为正，出流记为负）。
- **与溶质输运的相互作用。** 若启用地表反应输运而 `rtswbc.input` 未对触发式边界给出逐组分边界条件，则应用门控回退（各一次性警告）：
  - 触发式灌溉 → `RTSW_BC_TRIGGERED_IRRIGATION`：灌溉**开启**期间为浓度 `0.0` 的 Dirichlet（清水），**关闭**期间为零梯度（边界浓度跟随内部；无水进入时不强行施加溶质）。开启期间随水输入的溶质已由质量收支的对流入流项计入。
  - 触发式排水 → `RTSW_BC_TRIGGERED_DRAINAGE`（零梯度出流，浓度取内部值；只有排水实际出流时溶质才离开）。
  要覆盖回退，可在 `rtswbc.input` 中给出 `spec_X_bctype` 条目（见[第 1.2 节](#12-边界条件-rtswbcinput)）。含肥灌溉水建议使用门控类型 `spec_X_bctype : 5` 配 `spec_X_bcval` 或 `spec_X_bcfile` 时间序列，使浓度边界跟随灌溉启停循环。普通类型 2/3/4 仍被接受（无条件施加，与触发状态无关）。
- **与 WOFOST 的耦合。** `CROP_GROWTH_MODEL` 构建中使用作物感知重载：支持 `use_DVS : 1`（按触发格胞平均作物 DVS 选阶段），并为土壤水触发（模式 2–4）提供根深（RD）。DVS 与 RD 取自 WOFOST 状态；不需要额外输入文件。

---

## 5. 时间序列文件格式

若干输入文件对时变数据使用统一的时间序列格式，包括边界条件时间序列（`spec_X_bcfile`）、源汇时间序列（`spec_X_ssfile`）与强迫数据。

#### 格式

```text
ndata <N>
<time_1> <value_1>
<time_2> <value_2>
...
<time_N> <value_N>
```

| 字段 | 类型 | 单位 | 说明 |
|-------|------|-------|-------------|
| `ndata` | int | — | 数据点数。必须位于首行、任何数据之前。 |
| `time` | real | [s] | 自模拟起算的秒数。必须单调递增。 |
| `value` | real | 视场景 | 数据值。单位视场景而定（见下）。 |

#### 场景单位

| 文件类型 | 数值单位 | 示例 |
|-----------|-------------|---------|
| 边界条件时间序列（`spec_X_bcfile`） | [mg/L] | `rtirrigationNH4.input` — 灌溉水 NH4 浓度 |
| 源汇时间序列（`spec_X_ssfile`） | [mg/m²/s] | `rtfertNH4-DON-NH4.input` — NH4 施肥面积质量通量 |
| 蒸散时间序列 | [m/s] | `et.input` — 潜在蒸散速率 |
| 降雨时间序列 | [m/s] | `rainfall.input` — 降雨强度 |

#### 插值

- **边界条件时间序列**（`spec_X_bcfile`，类型 4/5）：恒为线性插值。
- **sstype = 1**（线性）：相邻数据点间线性插值。
- **sstype = 2**（阶梯）：分段常数（阶梯函数）。数值保持到下一个时间点。
- **阶梯效果**也可在线性插值下实现：相邻数据点取相同值即可。
- **端点行为**：t 超出最后一个数据点后**保持末值**（若末点为非零浓度/通量，会一直生效到模拟结束，末行建议放到模拟总长）；t 早于首点时按第一段直线插值（可能向后外推，首行建议放 t = 0）。

#### 示例 — 施肥（`rtfertNH4-DON-NH4.input`）

```text
ndata 126
0.00000e+00 0.00000e+00
561600.00    1.31250e+00
838800.00    1.12500e+00
1080000.00   0.00000e+00
3240000.00   1.68750e+00
3499200.00   0.00000e+00
...
95040000.0   0.00000e+00
```

（第 2–6 行依次对应：DAT 6.5 d NH4 脉冲开始、9.7 d 施用、12.5 d 脉冲结束、37.5 d 第二次施用、40.5 d 结束；末行为模拟结束。）

> **注 1**：`ndata` 关键字通过先读一个 token 再读计数的方式解析。解析器忽略标签、读取随后的数字。
> **注 2**：数据文件必须是**纯数字两列**——解析器用原始数据流读取 token，**不允许注释行或行内注释**（`rtswbc.input` 等键值文件的注释规则不适用于时间序列数据文件）。

---

## 6. 输出文件与表头

输出文件写入配置的输出目录。文本时间序列为空白分隔列；首行为表头，其后各行为科学计数法数值。`spec0`、`spec1` 等是 0 基组分索引，不是 `SPECIES_Name` 中的字符串。

### 6.1 反应输运文本输出

#### `RTSubsurfaceTimeSeries.out`

初始化时与每个反应输运输出时刻各写一行。确切表头取决于运行配置：

| 表头模式 | 单位 | 含义 |
|---|---|---|
| `Time[s]` | s | 模拟时间。 |
| `specX_LiquidMass[mg]` | mg | 地下水格胞中组分 X 的溶解态总质量。 |
| `specX_SolidMass[mg]` | mg | 吸附态总质量；组分 X 有吸附时出现。 |
| `specX_TotalMass[mg]` | mg | 溶解态 + 固态质量。 |
| `specX_StorageChange[mg/s]` | mg/s | 该输运步的储量变化率。 |
| `specX_MassExch[mg/s]` | mg/s | 输运积分器使用的质量交换诊断量。 |
| `specX_ExchangeIntoGW[mg/s]` | mg/s | 地表水向地下水的交换贡献。 |
| `specX_BCInflow[mg/s]`、`specX_BCOutflow[mg/s]` | mg/s | 反应输运边界总入流与总出流。 |
| `specX_RTBC..._Inflow[mg/s]`、`..._Outflow[mg/s]` | mg/s | 每个已配置地下水 RT 边界的入/出流分项。表头编码了边界索引、位置、方向与类型。 |
| `specX_SSInflow[mg/s]`、`specX_SSOutflow[mg/s]` | mg/s | 进入/离开地下水溶质系统的源汇质量。 |
| `specX_Reaction[mg/s]` | mg/s | 净反应贡献；反应模块激活时出现。 |
| `specX_BalanceResidual[mg/s]` | mg/s | 绝对质量平衡残差。 |
| `specX_BalanceResidualRelThroughput[-]` | - | 残差除以吞吐量。 |

当 `SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS=1` 时，会在最终残差列之前插入附加列 `BalanceRHS`、`BalanceResidualRel`、`TransportHalf1Residual`、`TransportHalf3Residual` 与 `NegativeClippedMassRate`。其单位编码在程序输出的表头中。

#### `RTSurfaceTimeSeries.out`

每个组分包含 `SurfaceMass[mg]`、`StorageChange[mg/s]`、`BCInflow[mg/s]`、`BCOutflow[mg/s]`、`SSInflow[mg/s]`、`SSOutflow[mg/s]` 及最后的平衡残差列。启用地表–地下交换时出现 `MassExch[mg/s]`，启用反应时出现 `Reaction[mg/s]`，启用降雨强迫时出现 `RainInflow[mg/s]`。可选的质量守恒指标列遵循上述同一编译期开关。

### 6.2 水动力与耦合文本输出

| 文件 | 表头字段与含义 |
|---|---|
| `domainTimeSeries.out` | `Time` [s]、`SurfaceVolume` [m3]、可选的边界入/出流与累计体积 [m3/s 或 m3]、可选的逐边界 `BoundaryFlow_i` 与 `BoundaryAccum_i`、可选 `RainFlux`/`RainAccum`、`SurfaceEvapFlux`/`SurfaceEvapAccum`、入渗 `InfFlux`/`InfAccum`，以及可选 `SurfaceSolidVolume`。实际列取决于边界、降雨、入渗与编译期泥沙设置。 |
| `SubsurfaceTimeSeries.out` | `Time` [s]、`SubSurfaceVolume [m3]`、`ExchangeRate [m3/s]`、`BoundaryInflow`、`BoundaryOutflow`、`Source/SinkInflow [m3/s]`、`Source/SinkOutflow [m3/s]`、`Source/SinkOutflow_RootTransp[m3/s]` 与 `Source/SinkOutflow_SoilEvap[m3/s]`。 |
| `RootZoneWaterContentTimeSeries.out` | 时间列随后为每个已配置根系层耦合地下水源汇对象一个 `RootZoneWaterContent_Mean [-]` 值。 |

### 6.3 NetCDF 输出与变量表头

| 文件 | 维度 | 主要变量与单位 |
|---|---|---|
| `output_subsurface.nc` | `t,z,y,x` | `t` [s]、坐标 `x,y,z` [m]、`hd` [m]、`wc` [-]、`qx,qy,qz` [m3/s 或模型通量约定]；启用 ET 时还有 `evapGW` 与 `transpGW`。 |
| `output_transport.nc` | `t,n,z,y,x` | `t` [s]、`n` 0 基组分索引、`x,y,z` [m]、`c` [mg/L]、`c_solid` [mg/kg]。调试构建可能增加 `c_difxz`、`c_difzx`、`c_difxx`、`c_difzz`。 |
| `output_transportSW.nc` | `t,n,y,x` | `t` [s]、`n` 0 基组分索引、`x,y` [m]、`csw` [mg/L]。 |
| `crop.nc` | `t,y,x` | `t` [s]、`x,y` [m]、`DVS` [-]、`LAI` [m2/m2]、`RD` [cm]、`TAGP/TWSO/TWLV/TWST/TWRT` [kg/ha]、`RFTRA` [-] 与 `RZSM` [cm3/cm3]。 |

输运 NetCDF 写出器在每个输出索引存储当前状态。当前源码中变量属性有限，上表单位即模型约定，后处理时应照此使用。

---

## 7. 地表水水动力（SWE）

地表水模块在基于 DEM 的均匀网格上求解二维浅水方程。水平离散（`nx` × `ny` 个单元、单元尺寸、左下角坐标）**唯一地**由 `dem.input` 定义；其余所有栅格输入（糙率、初始条件、溶质初始场）必须带有完全一致的头信息。触发式灌溉/排水边界类型（`bctype` 13/14）的专用输入见[第 4 节](#4-触发式灌溉与排水边界)；溶质输入见[第 1 节](#1-地表水溶质输运)。

最小地表水动力配置只需 `parameters.input`、`dem.input` 与 `sw.input`。`extbc.input`、`rainfall.input`、`infiltration.input` 与 `wind.input` 均为可选，可直接省略。

### 7.1 模拟控制参数 parameters.input

**文件路径**：`<input_folder>/parameters.input`
**源码**：`src/Parser.h` → `readParamsFile()`

冒号分隔的 `key : value` 文件，允许 `//` 开头的整行注释与行尾 `//` 注释。**未识别的关键字会使运行中止**，所有必填关键字必须出现。

#### 参数

| 键 | 类型 | 必填 | 说明 |
|-----|------|------|------|
| `simLength` | 实数 | **是** | 模拟总时长 [s]。 |
| `endTime` | 实数 | 否 | 显式停止时刻 [s]；给出时运行在 `endTime` 停止（末步相应缩短）。 |
| `cfl` | 实数 | **是** | 自适应水动力时间步的目标 Courant 数（如 0.5）。 |
| `parNx` | 整数 | **是** | x 方向 MPI 进程数（须与总进程数一致，`parNx × parNy`）。 |
| `parNy` | 整数 | **是** | y 方向 MPI 进程数。 |
| `outFreq` | 实数 | **是** | 场（NetCDF/VTK）输出的间隔 [s]。 |
| `obsFreq` | 实数 | **是** | 时间序列输出的间隔 [s]。 |
| `nScreen` | 整数 | 否 | 控制台进度打印的时间步间隔（默认 100）。 |
| `outFormat` | 字符串 | **是** | 输出格式：`NETCDF`、`VTK` 或 `BIN`。 |
| `BCtype` | 字符串 | **是** | 四个域边（未被 `extbc.input` 多边形覆盖处）的默认边界：`REFLECTIVE`、`TRANSMISSIVE` 或 `PERIODIC`。 |
| `writeRain`、`writeLandUse`、`writeSoilMap`、`writeRoughness`、`writeInfParameters` | 整数 | 否 | 置 1 时将相应场写入输出文件。 |

#### 示例 — `tests/t9-gaoyou-Ncycle-WOFOST/parameters.input`

```text
simLength  : 86400  // How long to simulate (s)
cfl        : 0.5      // CFL value to use
outFreq    : 3600      // Frequency for output in seconds
obsFreq    : 3600     // Frequency for time series output
nScreen    : 60
parNx      : 2       // Number of parallel tasks in the x-direction
parNy      : 2        // Number of parallel tasks in the y-direction
outFormat  : NETCDF   // Output format: NETCDF or VTK
BCtype     : REFLECTIVE // outer boundaries: TRANSMISSIVE, REFLECTIVE or PERIODIC
```

### 7.2 地形 dem.input

**文件路径**：`<input_folder>/dem.input`
**源码**：`src/Parser.h` → `readHeaderDEMFile()` / `readDEMFile()`；`src/rasterTools.h` → `readRasterField()`

ESRI-ASCII 风格栅格。六条头记录（每行一个 `标签 数值` 对——**只读取数值部分，标签不校验**），随后为高程数据：

| 记录 | 含义 |
|---|---|
| 1 | 列数 `nx` |
| 2 | 行数 `ny` |
| 3 | 左下角 x 坐标 `xll` [m] |
| 4 | 左下角 y 坐标 `yll` [m] |
| 5 | 单元尺寸 `dx` [m]（正方形、均匀单元） |
| 6 | NODATA 值 |

排序与校验规则：

- **第一行数据对应最北端（y 最大）的行**；行自北向南排列，行内数值自 x 小到大（标准 ESRI ASCII 顺序）。
- 等于 NODATA 值的单元被钝化（不参与计算）。
- 所有经 `readRasterField()` 读入的其它栅格（`roughness.input`、`hini.input`、`uini.input`、`vini.input`、溶质初始场）重复同样六条头记录；头信息（nx、ny、xll、yll、dx）必须与 `dem.input` **完全一致**，否则运行中止。

#### 示例 — `tests/t9-gaoyou-Ncycle-WOFOST/dem.input`

```text
n_cols 61
n_rows 221
xllcorner 0.0
yllcorber 0.0
cellsize 0.5
nodata_value -9999
7.310 7.310 7.309 ...   （每行 61 个值，共 221 行）
```

（标签中的笔误——如上面的 `yllcorber`——可以容忍，因为解析是按位置的。）

### 7.3 初始状态与摩阻 sw.input

**文件路径**：`<input_folder>/sw.input`
**源码**：`src/Parser.h` → `readSWFile()`

| 键 | 类型 | 必填 | 说明 |
|-----|------|------|------|
| `initialMode` | 字符串 | **是** | 地表初始条件：`dry`（干域）、`h`（恒定水深）、`h+z`（恒定自由水面高程；h = max(wse − z, 0)）或 `file`（读取 `hini.input`、`uini.input`、`vini.input` 栅格）。小写。 |
| `initialValue` | 实数 | 模式为 h 或 h+z 时 | 模式 `h`（水深 [m]）或 `h+z`（自由水面高程 [m]）的常数值。 |
| `friction` | 字符串 | **是** | 摩阻模型：`none`、`manning`、`darcyweisbach` 或 `chezy`。必须与编译期 `SERGHEI_FRICTION_MODEL` 一致（本构建为 Manning）。 |
| `roughness` | 字符串 | friction ≠ none 时 | 糙率系数数值（全域统一），或 `file` 读取 `roughness.input`（与 DEM 同头的栅格），或 `landuse` 由土地利用图推求。 |
| `dryDepth` | 实数 | **是** | 干单元处理的水深阈值 [m]，低于该值速度归零，如 0.001。 |

#### 示例 — `tests/t9-gaoyou-Ncycle-WOFOST/sw.input`

```text
initialMode :  h+z      // dry, h, h+z, file
initialValue : 7.4     // depth or water surface elevation constant value

friction : manning       // none, manning
roughness : 0.15         // roughness coefficient value, or set to 'file' to read from file

dryDepth : 0.001
```

### 7.4 外部边界 extbc.input

**文件路径**：`<input_folder>/extbc.input`
**源码**：`src/Parser.h` → `readExtBCFile()`；边界施加见 `src/BC.h`

以多边形方式定义外部（域边）地表边界。文件缺失时，全域各边采用 `parameters.input` 中的全局 `BCtype`。文件以 `bccount`（边界块数）开头，随后每个边界一块；每块以 `id` 开始：

| 键 | 类型 | 说明 |
|-----|------|------|
| `bccount` | 整数 | 边界块数（必须为首键）。 |
| `id` | 字符串 | 边界标识（开启一个新块）。 |
| `bctype` | 整数 | 边界类型码（见下表）。 |
| `direction` | 实数 ×2 | 外法向分量 `(normalx, normaly)`，如 `0.0 1.0`；向量在内部归一化，并决定给定通量的方向。 |
| `polygon` | 字符串 | 选取边界单元的二维多边形文件（`NPOINTS n` + `x y` 行）。 |
| `bcvals` | 实数（若干） | 常数边界值；含义随 `bctype` 而定（见下表）。 |
| `hydrograph` | 字符串 | 类型 10/11/12 的时间序列文件（格式见[第 5 节](#5-时间序列文件格式)；纯数字、不允许注释）。 |
| `trigger_polygon`、`irrigation_strategy`、`drainage_strategy` | 字符串 | 仅用于 `bctype` 13/14 — 见[第 4 节](#4-触发式灌溉与排水边界)。 |

未识别的关键字会使运行中止；块数必须等于 `bccount`；每个多边形至少须选中一个单元。

#### 边界条件类型

| 码 | 名称 | 取值来源 | 说明 |
|------|------|--------------|------|
| 1–3 | `PERIODIC` / `REFLECTIVE` / `TRANSMISSIVE` | — | 全域边界类型；仅可用于 `parameters.input` 的 `BCtype`，不能作为 `extbc` 块类型。 |
| 5 | `CRITICAL` | — | 临界流出（沿外法向 Froude = 1）。 |
| 6 | `H_CONST` | `bcvals : h` | 恒定水深 [m]。 |
| 7 | `Q_CONST` | `bcvals : <占位> qx qy` | 恒定单宽流量 [m²/s]：**第 2、3 个**值给出 x/y 分量；若给定流量为超临界，边界水深抬升至临界水深。 |
| 8 | `WSE_CONST` | `bcvals : wse` | 恒定自由水面高程 [m]。 |
| 9 | `FREE_OUTFLOW` | — | 零梯度自由出流。 |
| 10 | `HZ_T_INLET` | `hydrograph` = WSE(t) [m] | 入口水位过程；边界水深 = max(wse(t) − z, 0)。 |
| 11 | `HZ_T_OUTLET` | `hydrograph` = WSE(t) [m] | 出口水位过程；允许反向流动（如潮汐入流）。 |
| 12 | `Q_T` | `hydrograph` = Q(t) [m³/s] | 流量过程，分配到边界各单元；正值沿 `direction` 方向作用。 |
| 13 | `TRIGGERED_IRRIGATION` | 策略文件 | 基于规则的灌溉 — 见[第 4 节](#4-触发式灌溉与排水边界)。 |
| 14 | `TRIGGERED_DRAINAGE` | 策略文件 | 基于规则的排水 — 见[第 4 节](#4-触发式灌溉与排水边界)。 |

#### 示例 — `tests/t9-gaoyou-Ncycle-WOFOST/extbc.input`（灌溉入口与出口，均为流量过程）

```text
bccount : 2
id : irrigation
bctype : 12
direction : 0.0 1.0
polygon : polygoninlet.input
hydrograph: inlet_data_rate-2-input.input
id : outlet
bctype : 12
direction : 0.0 -1.0
polygon : polygonoutlet.input
hydrograph: outlet_data_rate-2-input.input
```

### 7.5 降雨 rainfall.input（可选）

**源码**：`src/Parser.h` → `readRainfallFile()`

文件缺失则无降雨。头信息：

```text
timeUnits s          // s 或 h
rainUnits mm/s       // mm/s 或 mm/h
TimeSeriesNp <np>    // 时间点数
nx <nx>              // 降雨场的 x 分区数
ny <ny>              // 降雨场的 y 分区数
```

随后为 `np` 行：`时间`（换算为 [s]）+ `nx × ny` 个降雨值（自动换算为 [m/s]）。`nx = ny = 1` 时降雨全域均匀；分布式降雨场按 x 最快排序（分区 `k = j*nx + i` 对应 DEM 网格），每个分区占表中的一列。

### 7.6 入渗 infiltration.input（可选）

**源码**：`src/Parser.h` → `readInfiltrationFile()`

文件缺失时地表不透水（`INF_NONE`）。**地表-地下水耦合运行（本构建）请保持 `infmodel : none`** ——此时入渗由[第 8.4 节](#84-边界条件-gwbcinput)的交换边界驱动，若再开地表入渗模型会重复扣除水量。

| 键 | 类型 | 说明 |
|-----|------|------|
| `infmodel` | 字符串 | `none`、`constant`、`horton` 或 `greenampt`。 |
| `spatialInput` | 字符串 | 参数的空间分布：`soilmap`、`landuse`、`netcdf` 或 `raster`。 |
| `infiltrationclasses` | 整数 | 参数类别数（标签 0 恒为不透水）。 |
| `rate` | 实数 | constant 模型：各类入渗能力 [mm/s]（自动换算为 m/s）。 |
| `k`、`fc`、`f0` | 实数 | Horton 模型：各类的衰减率与初始/最终能力 [mm/s]。 |
| `ks`、`psi`、`dtheta` | 实数 | Green-Ampt 模型：饱和导水率 [m/s]、吸力水头 [m]、含水率亏缺 [-]。 |

### 7.7 风场 wind.input（可选）

**源码**：`src/Parser.h` → `readWindFile()`

三条头记录——数据点数、风应力系数 `CwT`、起算风应力的最小水深 `hwmin`——随后为 `t vx vy` 行（风速分量 [m/s]）。文件缺失则无风驱动。

---

## 8. 地下水水动力（地下水流）

地下水流模块在地形贴合网格上求解混合形式的三维 Richards 方程：水平方向复用 DEM 网格（`nx` × `ny`），垂向在地形之下悬挂 `ndepth` 层。模块以 `-DSERGHEI_SWE_GW=ON` 编译进来（宏 `SERGHEI_RE_MODEL`；本构建已启用），且**始终激活**：`subsurface.input` 为必需文件。其输出（`SubsurfaceTimeSeries.out`、`output_subsurface.nc`）已在[第 6 节](#6-输出文件与表头)列出。

| 文件 | 作用 |
|---|---|
| `subsurface.input` | 垂向离散、土壤数、初始模式、求解器与耦合开关（**必需**） |
| `vg.input` | 各土壤类型的 van Genuchten–Mualem 参数 |
| `soilID.input` | 三维土壤类型索引场 |
| `gwbc.input` | 地下水边界条件（可选；未覆盖面默认不透水） |
| `gwss.input` | 地下水源汇区：蒸散发、补给、排水、暗管排水（可选） |
| `head.input` / `theta.input` / `wt.input` | `initialMode` = 2/3/4 时的初始条件场 |

### 8.1 计算域与求解器 subsurface.input

**文件路径**：`<input_folder>/subsurface.input`
**源码**：`src/GwInit.h` → `readGwFile()`

冒号分隔的 `key : value` 文件，带 `//` 注释。所有必填键缺失会使运行中止；未知键被**静默忽略**（例如旧测试用例中出现的 `BCtype` 行并无作用）。

| 键 | 类型 | 必填 | 说明 |
|-----|------|------|------|
| `ndepth` | 整数 | **是** | 垂向单元数（`nz`）。 |
| `height` | 实数 | **是** | 计算域厚度 [m]；`dz_multiplier : 1.0` 时 `dz = height / ndepth`。 |
| `dz_multiplier` | 实数 | **是** | 垂向拉伸比；`1.0` = 均匀 `dz`；`> 1` = 几何增长的 `dz`（实际厚度将超过 `height`）。 |
| `dz_base` | 实数 | multiplier > 1 时 | 拉伸网格最细（顶部）单元的厚度 [m]。 |
| `nSoilID` | 整数 | **是** | 土壤类型数 = `vg.input` 中每个键的条目数。 |
| `initialMode` | 整数 | **是** | 初始条件：`1` 饱和、`2` `head.input`、`3` `theta.input`、`4` `wt.input` — 见[第 8.6 节](#86-初始条件)。 |
| `dt_init` | 实数 | **是** | 地下水初始（也是最小）时间步 [s]。 |
| `dt_max` | 实数 | **是** | 地下水最大时间步 [s]。 |
| `gw_scheme` | 整数 | **是** | 非线性格式：`1` = 修正 Picard（PCA）、`2` = Picard 迭代。 |
| `cg_iter` | 整数 | 否 | CG 最大迭代数（默认 10⁶）。 |
| `cg_tol` | 实数 | 否 | CG 收敛容差（默认 1e-8）。 |
| `aev` | 实数 | **是** | 修正 van Genuchten 模型的进气值 [m]，用于平滑饱和附近的导水率/持水曲线；`0.0` 关闭平滑。 |
| `async` | 整数 | **是** | 地表-地下水耦合模式：`1` = 异步（地下水模块以自己的自适应 `dt` 在 `dt_init` 与 `dt_max` 之间推进，追上地表即交换；地表步长以地下水步长为上限）；`0` = 同步（两模块每步取相同的较小 `dt`）。 |

#### 示例 — `tests/t9-gaoyou-Ncycle-WOFOST/subsurface.input`

```text
BCtype  : REFLECTIVE  // legacy, ignored
height : 0.4          // Thickness of the subsurface domain
dz_multiplier : 1.0   // Non-uniform dz
dz_base : 0.05
ndepth  : 8           // Number of cells in z-direction
nSoilID : 4           // Number of soil types
initialMode : 1       // saturated(1), file-h(2), file-theta(3), file-wt(4)
dt_init : 0.01        // initial (also minimum) dt for subsurface
dt_max  : 100.0       // maximum allowable dt
gw_scheme : 1         // numerical schemes, 1:PCA, 2:Picard
async : 1             // Asynchronous coupling
aev : 0.0
```

### 8.2 土壤水力参数 vg.input

**文件路径**：`<input_folder>/vg.input`
**源码**：`src/GwInit.h` → `readVGParameters()`

每个 van Genuchten–Mualem 参数一个键；`nSoilID` 类土壤的取值以**分号分隔列表**给出（无空格）：

| 键 | 单位 | 说明 |
|-----|-------|------|
| `Ks` | [m/s] | 饱和导水率。 |
| `Phi` | [-] | 孔隙度（可与 `ThetaS` 相等）。 |
| `ThetaS` | [-] | 饱和含水率。 |
| `ThetaR` | [-] | 残余含水率。 |
| `n` | [-] | van Genuchten `n`（> 1）；`m = 1 − 1/n` 由内部计算。 |
| `alpha` | [1/m] | van Genuchten `alpha`。 |

各单元的土壤索引（来自 `soilID.input`）选取对应的列。**请始终提供 `vg.input`**——文件缺失时的回退路径依赖不再解析的旧标量键。

#### 示例 — `tests/t9-gaoyou-Ncycle-WOFOST/vg.input`

```text
alpha   : 0.56;0.68;1.69;1.93
n       : 1.65;1.61;1.33;1.31
Ks      : 4.86e-6;2.94e-6;5.61e-7;4.88e-7
Phi     : 0.51;0.47;0.43;0.39
ThetaR  : 0.048;0.046;0.036;0.035
ThetaS  : 0.51;0.47;0.43;0.39
```

### 8.3 土壤类型分布 soilID.input

**文件路径**：`<input_folder>/soilID.input`
**源码**：`src/GwInit.h` → `readSoilID()`

首记录 `n_soil <N>`（必须等于 `nSoilID`），随后为 `nx × ny × nz` 个整数（0 基土壤索引，对应 `vg.input` 的列）。排序为 **x 最快、次 y、次 z**：每层内按 DEM 自北向南的行主序；层自上而下。`nSoilID : 1` 时数据块（乃至文件本身）可省略——所有单元使用土壤类型 0。

### 8.4 边界条件 gwbc.input

**文件路径**：`<input_folder>/gwbc.input`
**源码**：`src/GwInit.h` → `readGwBCFile()`；边界施加见 `src/GwBC.h`

可选文件；未被边界块覆盖的地下水面为不透水（零通量）。结构与 `extbc.input` 相同：`bccount` 开头，随后以 `id` 开始的块。**`direction` 必须写在 `bcvals`/`bcfile` 之前**（它决定边界数据数组的长度）。

| 键 | 类型 | 说明 |
|-----|------|------|
| `bccount` | 整数 | 边界块数（首键）。 |
| `id` | 字符串 | 边界标识（开启一个块）。 |
| `direction` | 整数 | 边界面：`1` +x（东）、`2` −x（西）、`3` +y（北）、`4` −y（南）、`5` 底面（计算域基底）、`6` 顶面（地表）。 |
| `bctype` | 整数 | 类型码（见下表）。 |
| `polygon` | 字符串 | 二维多边形文件；多边形投影到所选面。 |
| `bcvals` | 实数 | 广播到所有边界面的常数；含义随类型而定。 |
| `timeseries` | 字符串 | 类型 5/6/7 的时间序列文件。 |
| `bcfile` | 字符串 | 类型 2/3/4 的空间分布取值文件。 |

#### 边界条件类型

| 码 | 名称 | 取值 | 说明 |
|------|------|-------|------|
| 1 | `NOFLOW` | — | 不透水面。 |
| 2 | `H_CONST` | `bcvals : h` [m] | 恒定压力水头（Dirichlet）。 |
| 3 | `Q_CONST` | `bcvals : q` [m/s] | 恒定 Darcy 通量（Neumann）。 |
| 4 | `WT_CONST` | `bcvals : wt` [m] | 恒定潜水面高程；水头 = wt − z。 |
| 5 | `H_T` | `timeseries` h(t) [m] | 压力水头时间序列（线性插值）。 |
| 6 | `Q_T` | `timeseries` q(t) [m/s] | 通量时间序列。 |
| 7 | `WT_T` | `timeseries` wt(t) [m] | 潜水面时间序列；支持逐单元列。 |
| 8 | `SWE` | — | **地表-地下水交换**；`direction` 必须为 `6`。幽灵单元水头取地表水头；区分三种交换状态（不交换 / 入渗 / 出渗，出渗以饱和导水率 `Ks` 封顶）。耦合运行的标准顶边界。 |
| 9 | `FD` | — | 自由排水（单位重力梯度，q = −k）；`direction` 必须为 `5`。标准底边界。 |

时间序列文件（类型 5/6/7）以 `ndata <np> <ncell>` 开头；每行为 `t value`（`ncell = 1`）或 `t v1 … v_ncell`（逐单元取值，用于 `WT_T`）。`bcfile` 文件以 `ndata <nx> <ny> <nz>` 开头，随后为取值，维数须与所选面匹配（`direction` 1/2：`ny × nz`；3/4：`nx × nz`；5/6：`nx × ny`）。

#### 示例 — `tests/t9-gaoyou-Ncycle-WOFOST/gwbc.input`（顶部 SWE 交换 + 底部自由排水；侧向恒定水头已注释）

```text
bccount : 2
id : surface
direction : 6
bctype : 8
bcvals : 0.0
polygon : polygontop.input
id : bottom
direction : 5
bctype : 9
polygon : polygontop.input

//id : left
//direction : 2
//bctype : 2
//bcvals : 0
//polygon : polygonleft.input
//bcfile : headleft.input
```

### 8.5 源汇区 gwss.input

**文件路径**：`<input_folder>/gwss.input`
**源码**：`src/GwInit.h` → `readGwSSFile()`；施加见 `src/SourceSink.h`

可选文件（缺失则无地下水源汇；如需停用可改名，`tests/t15` 即如此处理）。结构：`sscount` 开头，随后以 `id` 开始的块：

| 键 | 类型 | 说明 |
|-----|------|------|
| `sscount` | 整数 | 源汇块数（首键）。 |
| `id` | 字符串 | 标识（开启一个块）。 |
| `sstype` | 整数 | 类型码（见下表）。 |
| `polygon` | 字符串 | **三维**多边形文件：`NPOINTS n` + `x y z` 行；`z` 在地表以下取负值（如 `-0.3` = 地表下 0.3 m）。4 顶点 = 足印 + 底部深度（该深度以上整柱生效）；8 顶点 = 两个深度之间的闭合长方体（用于暗管排水）。 |
| `ssvals` | 实数 | 常数值（类型 2：固定水头 [m]）。 |
| `timeseries` | 字符串 | 类型 0：ET 文件（格式见下）；类型 1：`ndata` + `t q` 文件。 |
| `TileDrainageFile` | 字符串 | 类型 3：暗管参数文件（见下）。 |
| `shapecorrection` | 实数 | 类型 2：导水率修正系数 `Cpipe` [-]。 |
| `direction` | 整数 | 类型 2：排水特殊化的轴向（1–4；默认 z）。 |

#### 源汇类型

| 码 | 类型 | 说明 |
|------|------|------|
| 0 | 蒸散发 | 根系层水分摄取 — ET 文件格式见下。 |
| 1 | 通量 | 体积源汇速率 [1/s]（单位时间含水率变化；正 = 补给），取自 `timeseries` 文件。 |
| 2 | 水头排水 | 向固定水头 `ssvals` 的内部排水，导水率修正系数 `shapecorrection`。 |
| 3 | 暗管排水 | 暗管汇项 `Γ·h`（仅在压力水头为正处激活），参数来自 `TileDrainageFile`。 |

#### ET 时间序列文件（`sstype : 0` 的 `timeseries`）

```text
ndata <N>                     // 数据行数
f <0.5>                       // E/T 划分的 Beer 定律消光系数
stress <h1> <h2> <h3> <h4>    // Feddes 应力折减阈值，压力水头 [m]
rootx <xs> <ys> <zs>          // 根系分布参数（已解析；见注）
rootp <px> <py> <pz>          // 根系分布参数（已解析；见注）
<t> <ET0> <LAI> <RD>          // N 行：时间 [s]、潜在蒸散发 [m/s]、LAI [-]、根系深度 [m]
```

`E_pot = ET0·exp(−f·LAI)`，`T_pot = ET0 − E_pot`；蒸腾按 Feddes 应力函数 α(h)（`h1…h4` 之间）与归一化根系密度形状在根系层内分配。**`CROP_GROWTH_MODEL` 构建（本仓库）中 LAI 与 RD 每步取自 WOFOST，本文件只提供 ET0**（混合耦合）；此时 LAI/RD 两列为占位。当前根系密度形状使用固定参数（`zs = 0.2·RD`、`pz = 1`）——`rootx`/`rootp` 记录为兼容性而读取，不影响当前分布。潜在土壤蒸发同时发布给地表模块（有积水处蒸发积水）。

#### 暗管参数文件（`TileDrainageFile`，如 `drains.input`）

| 键 | 单位 | 说明 |
|-----|-------|------|
| `diameter` | [m] | 暗管有效直径 `d_e`。 |
| `scale_factor` | [-] | 网格尺寸效应的折减系数（常取 1、2 或 4）。 |

无量纲阻力修正 `Cd` 由 `d_e`、`scale_factor` 与网格尺寸按 HYDRUS 公式在内部计算。

#### 示例 — `tests/t9-gaoyou-Ncycle-WOFOST/gwss.input`

```text
sscount : 1
id : root1
sstype : 0
ssvals : 0.0
polygon : polygonroot.input
timeseries : et.input

//id : TileDrains
//sstype : 3
//polygon : polygondrains.input
//TileDrainageFile : drains.input
```

### 8.6 初始条件

**源码**：`src/GwInit.h` → `setGwState()` / `readGwICFile()`

由 `subsurface.input` 的 `initialMode` 选择：

| 模式 | 来源 | 说明 |
|------|--------|------|
| `1` | — | 饱和：`wc = ThetaS`；压力水头延续地表水柱（`h = 地表水深 + 地表以下深度`）。 |
| `2` | `head.input` | 三维压力水头场 [m]。 |
| `3` | `theta.input` | 三维含水率场 [-]；水头经 van Genuchten 反算。 |
| `4` | `wt.input` | 二维潜水面高程 [m]；逐单元 `h = wt − z`；NODATA 单元置为饱和。 |

文件格式（`head.input` / `theta.input` / `wt.input`）：三条头记录（`nx`、`ny`、NODATA 值——标签不校验），随后为取值，排序 **x 最快、次 y、次 z**（二维：`nx × ny`，与 DEM 相同的自北向南行序；三维：每层一个 `nx × ny` 块，自上而下）。维数必须与 DEM 一致，且模式 2/3/4 下文件**必须存在**。

### 8.7 地表与地下水耦合（小结）

- 采用 `SWE` 顶边界（[第 8.4 节](#84-边界条件-gwbcinput)）时，每个地下水步依据地表与地下水头计算交换通量；出渗以饱和导水率封顶，交换水量直接施加于地表水深。交换速率以 `ExchangeRate` 列于 `SubsurfaceTimeSeries.out`（见[第 6.2 节](#62-水动力与耦合文本输出)）。
- `async : 1` 解耦两模块时间步（推荐；所有参证工况均如此）；`async : 0` 强制同步。
- 本构建中降雨恒为地表驱动；蒸发经[第 8.5 节](#85-源汇区-gwssinput)的 ET 区在地表积水与顶层土壤之间分配。

---

## 附录 A 组分顺序约定

氮循环模型使用固定组分顺序。组分索引（0 基）必须在所有输入文件中一致。

| 索引 | 名称 | 化学式 | 在氮循环中的作用 |
|:---:|------|-----------------|----------------------|
| 0 | NH4（铵） | NH₄⁺ | 由矿化产生；被硝化消耗 |
| 1 | NO3（硝酸盐） | NO₃⁻ | 由硝化产生；被反硝化消耗 |
| 2 | DO（溶解氧） | O₂ (aq) | 被硝化与有氧呼吸消耗；由复氧补充 |
| 3 | DOC（溶解性有机碳） | C (org) | 反硝化的电子供体；被有氧呼吸消耗；由沉积物释放 |
| 4 | DON（溶解性有机氮） | N (org) | 由 DOC/DON 释放产生；被矿化为 NH4 消耗 |

### 氮循环路径

```
                    Mineralization（矿化）
    DON ──────────────────────────→ NH4
                                      │
                        Nitrification │ （硝化，需 O2）
                                      ▼
    O2 ←── Reaeration（复氧）      NO3
                                      │
                      Denitrification │ （反硝化，消耗 DOC，受 O2 抑制）
                                      ▼
                                     N2 (lost，流失)

    DOC + O2 ──→ CO2    （Aerobic Respiration，有氧呼吸）
    DOC ←── Sediment Release（沉积物释放）
    O2  ←── Atmospheric Reaeration（大气复氧）
```

---

## 附录 B 输入文件速查清单

下表列出耦合模型（水动力、溶质输运、作物生长）的全部输入文件。运行模拟前请确认输入目录中所需文件齐全。

### 地表水水动力

| 文件 | 必需 | 用途 |
|------|----------|---------|
| `parameters.input` | **是** | 模拟控制：时长、CFL、MPI 布局、输出频率/格式（[第 7.1 节](#71-模拟控制参数-parametersinput)） |
| `dem.input` | **是** | 地形栅格；定义水平网格（[第 7.2 节](#72-地形-deminput)） |
| `sw.input` | **是** | 初始状态、摩阻、糙率（[第 7.3 节](#73-初始状态与摩阻-swinput)） |
| `extbc.input` | 否* | 地表外边界，含 `bctype : 13`/`14` 触发式灌溉/排水（[第 4 节](#4-触发式灌溉与排水边界)与[第 7.4 节](#74-外部边界-extbcinput)；*使用全局 `BCtype` 以外的边界时必需） |
| `roughness.input`、`hini.input`、`uini.input`、`vini.input` | 由 sw.input 引用时 | 空间糙率/初始状态栅格（[第 7.3 节](#73-初始状态与摩阻-swinput)） |
| `rainfall.input` | 否 | 降雨强迫（[第 7.5 节](#75-降雨-rainfallinput可选)） |
| `infiltration.input` | 否 | 地表入渗模型——地表-地下水耦合运行请保持 `none`（[第 7.6 节](#76-入渗-infiltrationinput可选)） |
| `wind.input` | 否 | 风驱动（[第 7.7 节](#77-风场-windinput可选)） |

### 地表水溶质输运

| 文件 | 必需 | 用途 |
|------|----------|---------|
| `rttransportsw.input` | **是** | 组分数、输运属性、初始条件、求解器设置 |
| `rtswbc.input` | 否* | 溶质边界条件（*存在外部边界时必需） |
| `rtswss.input` | 否 | 溶质源汇区（施肥等） |
| `reactionsw.input` | 若 ReactionModule=1 | 反应动力学、氮循环参数 |
| `rtinitialsw*.input` | 若 aq_mode=1 | 各组分空间初始浓度栅格 |

### 地下水水动力

| 文件 | 必需 | 用途 |
|------|----------|---------|
| `subsurface.input` | **是** | 离散、求解器、耦合开关、初始模式（[第 8.1 节](#81-计算域与求解器-subsurfaceinput)；本构建中该模块恒激活） |
| `vg.input` | **是** | 各土壤类型 van Genuchten–Mualem 参数（[第 8.2 节](#82-土壤水力参数-vginput)） |
| `soilID.input` | 若 nSoilID > 1 | 三维土壤类型索引场（[第 8.3 节](#83-土壤类型分布-soilidinput)） |
| `gwbc.input` | 否 | 地下水边界：顶部 SWE 交换、底部 FD、水头/通量（[第 8.4 节](#84-边界条件-gwbcinput)） |
| `gwss.input` | 否 | 蒸散发/补给/排水/暗管源汇区（[第 8.5 节](#85-源汇区-gwssinput)） |
| `head.input` / `theta.input` / `wt.input` | 若 initialMode = 2/3/4 | 初始条件场（[第 8.6 节](#86-初始条件)） |
| `et.input` | 若 gwss sstype = 0 | 含 Feddes/根系参数的 ET0 序列（[第 8.5 节](#85-源汇区-gwssinput)） |
| `drains.input` | 若 gwss sstype = 3 | 暗管参数（[第 8.5 节](#85-源汇区-gwssinput)） |

### 地下水溶质输运

| 文件 | 必需 | 用途 |
|------|----------|---------|
| `rttransportgw.input` | **是** | 组分数、输运属性、初始条件、求解器设置 |
| `rtgwbc.input` | 否* | 溶质边界条件 |
| `rtgwss.input` | 否 | 溶质源汇区（根系吸收等） |
| `reactiongw.input` | 若 ReactionModule=1 | 反应动力学、吸附、氮循环 |
| `rtinitialgw*.input` | 若 aq_mode=1 | 空间液相初始浓度 |
| `rtinitialgwsolid*.input` | 若 solid_mode=1 | 空间固相（吸附）初始浓度 |

### WOFOST 作物生长

| 文件 | 必需 | 用途 |
|------|----------|---------|
| `cropparameter.input` | **是** | 作物参数、土壤属性、养分参数 |
| `cropmeteo.input` | **是** | 逐日气象强迫数据 |
| `agro.input` | **是** | 作物品种、播/收日期 |

### 支持文件

| 文件 | 必需 | 用途 |
|------|----------|---------|
| `extbc.input` | 使用地表边界时 | 水动力外边界条件（[第 7.4 节](#74-外部边界-extbcinput)）；`bctype : 13`/`14` 配置触发式灌溉/排水（见[第 4 节](#4-触发式灌溉与排水边界)） |
| 灌溉/排水策略文件（如 `test.input`、`test2.input`） | 若 bctype=13/14 | 触发式灌溉/排水的阶段表（见[第 4 节](#4-触发式灌溉与排水边界)） |
| `polygon*.input` | 按引用 | 边界/源汇空间区与触发观测区的多边形顶点坐标 |
| 时间序列文件（如 `rtfert*.input`） | 按引用 | 时变边界/源数据 |

### 主要输出文件

| 文件 | 产生条件 | 用途 |
|------|---------------|---------|
| `domainTimeSeries.out` | 启用地表水动力 | 域尺度地表水体积、边界流量、降雨、蒸发与入渗诊断。 |
| `SubsurfaceTimeSeries.out` | 启用地下水水动力 | 域尺度地下水体积、交换、边界流量与源汇诊断。 |
| `RTSurfaceTimeSeries.out` | 启用地表溶质输运 | 逐组分地表溶质质量与质量平衡诊断。 |
| `RTSubsurfaceTimeSeries.out` | 启用地下水溶质输运 | 逐组分液/固/总质量、边界/源汇/反应项与质量平衡诊断。 |
| `RootZoneWaterContentTimeSeries.out` | 启用根系层耦合 | 已配置耦合对象的根系层含水量。 |
| `output_transportSW.nc` | 启用地表溶质输运 | 时变二维溶解浓度场 `csw[t,n,y,x]`。 |
| `output_transport.nc` | 启用地下水溶质输运 | 时变三维溶解与吸附浓度场 `c[t,n,z,y,x]` 与 `c_solid[t,n,z,y,x]`。 |
| `crop.nc` | 启用 WOFOST | 时变作物发育、生物量、根深、胁迫与根系层含水量场。 |

完整文本表头与 NetCDF 变量清单见[第 6 节](#6-输出文件与表头)。可选调试与质量守恒列由编译期开关控制，因此并非每次运行都出现。

---

*用户手册完 — 地表-地下水动力、溶质输运与 WOFOST 作物生长*
