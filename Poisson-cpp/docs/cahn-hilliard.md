# Cahn–Hilliard 求解器:方程、时空离散方案与参数说明

本文档详细介绍 `cahn_hilliard` 算例:它求解什么方程、用什么边界条件、整套的**空间(间断
Galerkin)**与**时间(稳定化半隐 IMEX)**离散方案,以及 JSON 配置文件里**所有可调参数的含义与
取值建议**。代码位于 `src/cahn_hilliard/`,复用了仓库已有的 2D DG 脚手架(`src/common`、`src/ipdg`)。

> 速览见 [`../README.md`](../README.md) 第 6 节;本文是完整版。

---

## 目录

1. [方程与物理背景](#1-方程与物理背景)
2. [边界条件](#2-边界条件)
3. [空间离散:间断 Galerkin(SIPG)](#3-空间离散间断-galerkinsipg)
4. [时间离散:一阶 / 二阶 IMEX](#4-时间离散一阶--二阶-imex)
5. [完整算法](#5-完整算法)
6. [初值与可视化](#6-初值与可视化)
7. [可调参数(JSON 配置)](#7-可调参数json-配置)
8. [运行时诊断与预期行为](#8-运行时诊断与预期行为)
9. [收敛阶测试(MMS,`cahn_hilliard_convergence`)](#9-收敛阶测试mmscahn_hilliard_convergence)
10. [构建与运行](#10-构建与运行)
11. [代码结构](#11-代码结构)
12. [参考文献](#12-参考文献)

---

## 1. 方程与物理背景

Cahn–Hilliard(CH)方程描述二元混合物的**相分离**:一个守恒的序参量(相场)$c(\mathbf x,t)$
(例如两组分的浓度差,取值约在 $[-1,1]$,$\pm1$ 表示两个纯相)在自由能驱动下自发分离成两相,
并随时间**粗化**(coarsening,小畴融合成大畴)。

### 1.1 自由能泛函

Ginzburg–Landau 型自由能:

$$
\mathcal E[c]=\int_\Omega\Big(\underbrace{\tfrac{\varepsilon^2}{2}\,|\nabla c|^2}_{\text{界面/梯度能}}
+\underbrace{F(c)}_{\text{体相势}}\Big)\,\mathrm d\mathbf x,
\qquad F(c)=\tfrac14\,(c^2-1)^2 .
$$

- $F$ 是**双井势**,两个极小在 $c=\pm1$(两个纯相),驱动 $c$ 趋向 $\pm1$;
- 梯度项惩罚 $c$ 的剧烈变化,使两相之间形成**有限厚度的扩散界面**,厚度 $\sim\varepsilon$;
- 两项竞争决定平衡界面廓线 $c(x)=\tanh\!\big(x/(\sqrt2\,\varepsilon)\big)$。

势的导数(后面要反复用到):

$$
F'(c)=c^3-c,\qquad F''(c)=3c^2-1 .
$$

### 1.2 化学势与演化方程

**化学势**是自由能对 $c$ 的变分导数:

$$
\mu=\frac{\delta\mathcal E}{\delta c}=F'(c)-\varepsilon^2\Delta c .
$$

CH 方程是自由能在 $H^{-1}$ 度量下的**守恒型梯度流**(质量守恒的 Allen–Cahn 对应物):

$$
\boxed{\ \frac{\partial c}{\partial t}=\nabla\!\cdot\!\big(M\,\nabla\mu\big)=M\,\Delta\mu,\qquad
\mu=F'(c)-\varepsilon^2\Delta c\ }
$$

其中 $M>0$ 是**迁移率**(本代码取常数)。把 $\mu$ 代入即得四阶非线性抛物方程
$\partial_t c=M\,\Delta\big(c^3-c-\varepsilon^2\Delta c\big)$。

### 1.3 两条结构性质:质量守恒与能量耗散

在无通量边界(见 §2)下,CH 有两条本质性质,后面的离散格式会**逐条保持**:

- **质量守恒**:$\dfrac{\mathrm d}{\mathrm dt}\displaystyle\int_\Omega c
  =\int_\Omega M\Delta\mu=\oint_{\partial\Omega}M\,\partial_{\mathbf n}\mu=0.$
- **能量耗散**:$\dfrac{\mathrm d}{\mathrm dt}\mathcal E[c]
  =\displaystyle\int_\Omega\mu\,\partial_t c
  =\int_\Omega\mu\,M\Delta\mu=-M\!\int_\Omega|\nabla\mu|^2\le0.$

即总质量恒定、总自由能单调不增——这正是相分离 + 粗化的热力学根源。

### 1.4 混合(拆分)形式

直接离散四阶算子需要 $C^1$ 元或 $C^0$ 内罚(本仓库的双调和求解器走的就是后者)。CH 更常用
**混合 / 拆分形式**:把 $\mu$ 作为独立未知量,把一个四阶方程拆成两个二阶方程

$$
\partial_t c=M\,\Delta\mu,\qquad \mu=F'(c)-\varepsilon^2\Delta c .
$$

这样每个方程只含 $\Delta$,可以**直接复用已有的二阶 SIPG 拉普拉斯算子**(见 §3),无需 $C^1$ 元。
代价是未知量翻倍($c$ 与 $\mu$),得到一个 $2\times2$ 分块系统。

---

## 2. 边界条件

### 2.1 无通量(齐次 Neumann)

本算例在单位正方形 $\Omega=[0,1]^2$ 上,采用物理上最常用的**无通量(no-flux)**边界条件:

$$
\partial_{\mathbf n}\mu=0\quad(\text{无质量通量穿过边界}),\qquad
\partial_{\mathbf n}c=0\quad(\text{界面垂直于边界}),
$$

$\mathbf n$ 为外法向。$\partial_{\mathbf n}\mu=0$ 保证 §1.3 的质量守恒;两个齐次 Neumann 条件一起
保证能量耗散。直观上:盒子是封闭的,物质既不进也不出。

### 2.2 在弱形式中是"自然"的

对二阶算子 $-\Delta$ 做分部积分会产生边界项 $\oint_{\partial\Omega}(\partial_{\mathbf n}\cdot)\,(\cdot)$。
齐次 Neumann 条件让该边界项**直接为零**,因此在弱形式里**不需要任何边界罚项或 Nitsche 项**——
这类边界条件称为"自然边界条件"。具体到 SIPG:内罚只施加在**内部边**上,边界边不加任何项(见 §3.2)。
这与本仓库 Poisson IPDG(强 Dirichlet,内部边内罚)在装配上完全一致,可直接复用。

> 若改成 Dirichlet 边界(给定 $c$ 或 $\mu$),则需在边界边加 Nitsche/罚项,**且质量守恒会被破坏**
> ——本算例刻意保持无通量,正是为了精确守恒(见 §4.4 的代数前提 $A\mathbf 1=0$)。

---

## 3. 空间离散:间断 Galerkin(SIPG)

### 3.1 DG 空间

把 $\Omega$ 剖分成三角形网格 $\mathcal T_h$(`Mesh::getMesh` 生成 $[0,1]^2$ 上的结构化三角网,
每边 $N$ 个单元,$h=1/N$)。**间断**有限元空间为逐单元的 $P_k$ 多项式、单元间不要求连续:

$$
V_h=\{\,v\in L^2(\Omega):v|_K\in P_k(K)\ \ \forall K\in\mathcal T_h\,\}.
$$

自由度逐单元独立编号(`FEM::getDOF`),基函数取单元上的节点 Lagrange 基。$c_h,\mu_h\in V_h$。

### 3.2 $-\Delta$ 的对称内罚双线性型(SIPG)

对间断函数,$-\Delta$ 的对称内罚(Symmetric Interior Penalty)双线性型为

$$
a(u,v)=\sum_{K}\int_K\nabla u\cdot\nabla v
-\sum_{e\in\Gamma_{\mathrm{int}}}\int_e\Big(\{\!\{\partial_{\mathbf n}u\}\!\}[\![v]\!]
+\{\!\{\partial_{\mathbf n}v\}\!\}[\![u]\!]\Big)
+\sum_{e\in\Gamma_{\mathrm{int}}}\frac{\sigma}{h_e}\int_e[\![u]\!]\,[\![v]\!],
$$

其中 $\{\!\{\cdot\}\!\}$、$[\![\cdot]\!]$ 是边上的平均与跳量,$\sigma$ 是内罚系数,$h_e$ 是边长。

- 三项依次是:单元内能量、保证相容性/对称性的**consistency + adjoint-consistency**项、保证稳定性的
  **惩罚**项。$\sigma$ 足够大时 $a(\cdot,\cdot)$ 强制(coercive)。
- **关键**:求和只跑**内部边** $\Gamma_{\mathrm{int}}$,边界边不出现——这正对应 §2 的自然无通量条件。

记其矩阵为 $A$。代码里 $A=$ `assembleK_Poi2D`(单元内 $\int\nabla\phi_i\!\cdot\!\nabla\phi_j$)
$+$ `assembleIP_Poi2D(...,\sigma,\beta{=}1)`(内部边的内罚)。`beta=1` 即对称 SIPG,使 $A=A^\top$。
$A$ 对称、半正定,且 $a(u,v)\approx(-\Delta u,v)$。

### 3.3 关键性质:$A\mathbf 1=0$

把全局常数函数代入:单元内 $\nabla(\text{const})=0$,内部边上跳量 $[\![\text{const}]\!]=0$、
$\partial_{\mathbf n}(\text{const})=0$,故 $a(\mathbf 1,v)=0\ \forall v$,即

$$
A\,\mathbf 1=0 .
$$

常数在 $A$ 的零空间里。这是后面**离散质量精确守恒**的唯一代数前提(§4.4)。程序在装配后会断言
$\lVert A\mathbf 1\rVert_\infty\approx0$,实测 $\sim 7\times10^{-16}$。

### 3.4 质量矩阵与非线性项

- **DG 质量阵** $M_h$:$(M_h)_{ij}=\sum_K\int_K\phi_i\phi_j$,逐单元块对角、对称正定
  (`assembleMass_DG2D`)。
- **非线性载荷** $b^n$:$b^n_i=\sum_K\int_K F'(c^n_h)\,\phi_i$,$F'(c)=c^3-c$,显式装配
  (`assembleNonlinearCH`)。被积函数次数为 $4\,\text{ord}$($c_h^3$ 是 $3\,\text{ord}$,乘检验函数 $\text{ord}$),
  用 $4\,\text{ord}+2$ 阶求积。

### 3.5 半离散系统

在 $V_h$ 上检验,用 $a(\cdot,\cdot)\approx(-\Delta\cdot,\cdot)$(注意 $(\Delta\mu,v)=-a(\mu,v)$):

$$
\big(\partial_t c_h,v\big)+M\,a(\mu_h,v)=0,\qquad
\big(\mu_h,w\big)=\big(F'(c_h),w\big)+\varepsilon^2\,a(c_h,w),\qquad\forall v,w\in V_h.
$$

矩阵形式($\mathbf c,\boldsymbol\mu$ 为系数向量):

$$
M_h\,\dot{\mathbf c}+M\,A\,\boldsymbol\mu=0,\qquad
M_h\,\boldsymbol\mu-\varepsilon^2 A\,\mathbf c=\mathbf b(\mathbf c).
$$

### 3.6 复用了哪些已有部件

| 部件 | 来源 | 说明 |
|---|---|---|
| 网格 / 边-单元关系 | `Mesh`(`src/common`) | `getMesh`、`getEdge2Side` |
| $P_k$ 基、自由度、求积 | `FEM`、`Quadrature`(`src/common`) | DG 节点基 |
| 刚度阵 $K$ | `assembleK_Poi2D`(`src/ipdg`) | $\sum_K\int\nabla\phi\!\cdot\!\nabla\phi$ |
| 内罚阵 $P$ | `assembleIP_Poi2D`(`src/ipdg`) | 内部边 SIPG,`beta=1` |
| 质量阵 $M_h$ / 非线性 $b$ / 能量 / 栅格化 | `CahnHilliard.{h,cpp}`(**新增**) | CH 专用 |

---

## 4. 时间离散:一阶 / 二阶 IMEX

四阶项使显式格式有极苛刻的步长限制 $\tau\lesssim h^4/(M\varepsilon^2)$;全隐式又要每步解非线性方程。
本代码用**稳定化半隐(linearly stabilized semi-implicit / 凸分裂)IMEX** 格式:线性的四阶部分
**隐式**、非线性 $F'(c)$ **显式**,加稳定化项保持线性且换取稳定性。配置项 `time_order` 选择阶数:

- **`time_order=1`**(默认,§4.1–4.5):一阶稳定化后向欧拉型,**无条件能量稳定**、最稳健;
- **`time_order=2`**(§4.6):二阶 **SBDF2**(BDF2 + 非线性项二阶外推),时间精度更高。

两者的系统矩阵都与时间无关、只分解一次,且都精确守恒质量(依赖 $A\mathbf 1=0$)。

### 4.1 一阶格式

给定 $c^n$,求 $(c^{n+1},\mu^{n+1})$:

$$
\frac{c^{n+1}-c^n}{\tau}=M\,\Delta\mu^{n+1},\qquad
\mu^{n+1}=-\varepsilon^2\Delta c^{n+1}+S\,(c^{n+1}-c^n)+F'(c^n).
$$

- $-\varepsilon^2\Delta c^{n+1}$ 隐式(数值稳定);$F'(c^n)$ 显式(避免每步非线性求解,系统保持**线性**);
- $S(c^{n+1}-c^n)$ 是稳定化项,$S\to S$ 越大越稳但一阶误差略增。这相当于把 $F$ 做凸分裂
  $F=F_{\text{凸}}-F_{\text{凹}}$ 的线性化版本。

### 4.2 分块线性系统(只分解一次)

在 $V_h$ 上离散($(\Delta\mu,v)=-a(\mu,v)$,$(-\varepsilon^2\Delta c,w)=\varepsilon^2 a(c,w)$),得未知
$\mathbf x=[\mathbf c^{n+1};\boldsymbol\mu^{n+1}]$ 的 $2\times2$ 分块系统:

$$
\underbrace{\begin{bmatrix}\frac1\tau M_h & M\,A\\[4pt] -(\varepsilon^2 A+S M_h) & M_h\end{bmatrix}}_{=:J}
\begin{bmatrix}\mathbf c^{n+1}\\[2pt]\boldsymbol\mu^{n+1}\end{bmatrix}
=\begin{bmatrix}\frac1\tau M_h\,\mathbf c^n\\[4pt]\mathbf b^n-S M_h\,\mathbf c^n\end{bmatrix},
\qquad b^n_i=\int F'(c^n_h)\,\phi_i .
$$

矩阵 $J$ 由 $M_h,A,\tau,M,\varepsilon^2,S$ 组成,**全部与时间无关**。因此程序只对 $J$ 做**一次**
`SparseLU` 分解(系统非对称),之后每步只需:装配 $\mathbf b^n$、构造右端、一次回代。逐帧出片很快。

### 4.3 能量稳定性($S\ge L/2$)

记 $\delta=c^{n+1}-c^n$。用 $\mu^{n+1}$ 检验第一式、用 $\delta$ 检验第二式,配合对称 $A$ 的极化恒等式
与带余项的 Taylor 展开,可得**离散自由能**
$E_h(\mathbf c)=\tfrac{\varepsilon^2}{2}\mathbf c^\top A\mathbf c+\sum_K\int_K F(c_h)$ 的逐步差分:

$$
E_h(\mathbf c^{n+1})-E_h(\mathbf c^n)\le
-\tau M\lVert\mu^{n+1}\rVert_A^2
-\tfrac{\varepsilon^2}{2}\,\delta^\top A\,\delta
-\big(S-\tfrac{L}{2}\big)\lVert\delta\rVert_{M_h}^2 .
$$

只要

$$
\boxed{\,S\ \ge\ \tfrac{L}{2},\qquad L=\max|F''(c)|=\max|3c^2-1|\,}
$$

右端三项全 $\le0$,从而 $E_h$ **单调不增**(对**原始**离散能量成立,不是修正能量)。

- 这里 $L$ 是 $F''$ 在 $c$ 实际取值范围上的上界。四阶势的 $F''$ 无界,故严格无条件性需对 $\|c\|_\infty$
  有先验界 $M_\infty$:$S\ge(3M_\infty^2-1)/2$,即 $S$ 能覆盖到 $|c|\le\sqrt{(2S+1)/3}$。
- **$S=2$** 覆盖 $|c|\le\sqrt{5/3}\approx1.29$,足以应付 spinodal 初期的轻微过冲(实测峰值 $|c|\approx1.04$)。
  若过冲更大,程序会打印告警(检测 $3\max|c|^2-1>2S$),此时把 $S$ 调大(如 $S=3$ 覆盖 $|c|\le1.53$)。
- $S$ 与 $\tau,h,\varepsilon,M$ 都**无关**。

### 4.4 离散质量守恒(精确)

用全 1 向量左乘分块系统的第一行:

$$
\tfrac1\tau\mathbf 1^\top M_h\mathbf c^{n+1}+M\,(\mathbf 1^\top A)\boldsymbol\mu^{n+1}
=\tfrac1\tau\mathbf 1^\top M_h\mathbf c^{n}.
$$

由 $A=A^\top$ 与 $A\mathbf 1=0$(§3.3)得 $\mathbf 1^\top A=0$,中间项对**任意** $\boldsymbol\mu^{n+1}$ 消失,于是

$$
\mathbf 1^\top M_h\,\mathbf c^{n+1}=\mathbf 1^\top M_h\,\mathbf c^{n}\quad(\text{即}\ \textstyle\int_\Omega c_h\ \text{逐步精确守恒}).
$$

这只用到第一式,与非线性项 $\mathbf b^n$、$S$、$\varepsilon^2$、以及非线性是否显式**都无关**。实测漂移
$\sim10^{-14}$(纯属 LU 回代的浮点舍入)。**前提**仅是 $A\mathbf 1=0$,即必须保持无通量(不给 $c$ 加边界项)。

### 4.5 适定性

$M_h$ 对称正定,$A$ 对称半正定且 $A\mathbf 1=0$,$\tau,\varepsilon^2,S,M>0$。对 Schur 补
$S_c=\tfrac1\tau M_h+M\varepsilon^2 A M_h^{-1}A+M S\,A$(三项均半正定,且 $\tfrac1\tau M_h$ 严格锚定常数模)
做能量估计可得 $J\mathbf x=0\Rightarrow\mathbf x=0$,即 $J$ 非奇异。$A$ 的纯 Neumann 零空间**不会**传到 $J$,
因为时间导数项 $\tfrac1\tau M_h$ 控制了常数模。实测 $\mathrm{cond}(J)$ 在 $\tau$ 很小时约 $10^5$,随 $\tau$ 线性增长,
`SparseLU` 完全够用。

> 以上三条性质(能量稳定、质量守恒、适定)均由独立的多智能体核验工作流逐项推导确认,并与运行时诊断吻合。

### 4.6 二阶格式:SBDF2(`time_order=2`)

二阶半隐 BDF(SBDF2)= **BDF2** 时间导数 + 非线性项的**二阶外推** + 二阶相容的稳定化。
记外推 $c^\ast=2c^n-c^{n-1}$,给定 $c^n,c^{n-1}$ 求 $(c^{n+1},\mu^{n+1})$:

$$
\frac{3c^{n+1}-4c^n+c^{n-1}}{2\tau}=M\,\Delta\mu^{n+1},\qquad
\mu^{n+1}=-\varepsilon^2\Delta c^{n+1}+F'(c^\ast)+S\,(c^{n+1}-2c^n+c^{n-1}).
$$

分块系统(未知 $[\mathbf c^{n+1};\boldsymbol\mu^{n+1}]$):

$$
J_2=\begin{bmatrix}\tfrac{3}{2\tau}M_h & M\,A\\[3pt] -(\varepsilon^2 A+S M_h) & M_h\end{bmatrix},\qquad
\text{rhs}=\begin{bmatrix}\tfrac{1}{2\tau}M_h(4\mathbf c^n-\mathbf c^{n-1})\\[3pt]
\mathbf b(\mathbf c^\ast)-S M_h\,\mathbf c^\ast\end{bmatrix},\quad b(\mathbf c^\ast)_i=\int F'(c^\ast_h)\phi_i .
$$

- 与一阶 $J$ 相比**只有 $(1,1)$ 块系数变了**($\tfrac1\tau\!\to\!\tfrac{3}{2\tau}$);其余块、所有符号都一样。$J_2$ 仍与
  时间无关,只分解一次。
- **二阶精度**:BDF2 的时间导数是二阶;$c^\ast$ 是 $c^{n+1}$ 的二阶外推,故 $F'(c^\ast)=F'(c^{n+1})+O(\tau^2)$;
  稳定化项 $c^{n+1}-2c^n+c^{n-1}=O(\tau^2)$ 是相容扰动,不降阶。
- **启动**:第一步($c^{-1}$ 不存在)用一阶格式自举;单步一阶的解误差为 $O(\tau^2)$,且 BDF2 零稳定
  (特征根 $\{1,\tfrac13\}$),不破坏全局二阶。
- **质量守恒**同样精确:对第一行左乘 $\mathbf 1^\top$ 得 $3m^{n+1}=4m^n-m^{n-1}$($m^k=\mathbf 1^\top M_h\mathbf c^k$),
  自举给 $m^1=m^0$,归纳得 $m^n=m^0$。实测漂移 $\sim10^{-14}$。
- **能量稳定性提醒**:一阶的 $S\ge L/2$ 无条件**原始**能量递减结论**不**原样适用于 SBDF2(SBDF2 通常只在
  步长限制下、或对某个修正能量稳定)。因此 `time_order=2` 路径**不**断言原始 $E_h$ 单调;演示时用较小 $\tau$。

> SBDF2 的分块系统、自举、阶数与质量守恒均经独立多智能体工作流逐项核验,并由 §9 的收敛测试经验证实
> (时间二阶、空间 $k{+}1$ 阶)。

### 4.7 自适应时间步(电影用,`adaptive=true`)

电影里早期(自发相分离)变化极快、后期(粗化)很慢。固定 $\tau$ 要么早期欠采样("一闪而过")、要么后期白白多算。**自适应步长**按"每步相对变化量"自动调 $\tau$:

$$r_n=\frac{\|c^{n+1}-c^n\|_{M_h}}{\|c^{n+1}\|_{M_h}},\qquad
\tau_{n+1}=\mathrm{clamp}\Big(\tau_n\cdot\mathrm{clamp}\big(0.9\,\tfrac{r_{\mathrm{tol}}}{r_n},\,0.5,\,2\big),\ \tau_{\min},\ \tau_{\max}\Big).$$

变化大 → $\tau$ 自动变小;变化小 → $\tau$ 自动变大,使每步变化量稳定在目标 $r_{\mathrm{tol}}$(`rel_change`)附近。

- **正好解决"前几帧太快"**:每 `save_every` 步存一帧;因为每步变化量 ≈ $r_{\mathrm{tol}}$,"每 K 步一帧" ≈ **每帧等量变化**——早期 $\tau$ 小、步密 → 出帧密(慢放);后期 $\tau$ 大、步疏 → 出帧疏(快进略过粗化)。
- **效率**:$\tau$ 不连续变,而是**snap 到 $\tau_{\min}\cdot 2^{j}$ 的离散档**;每档的分块矩阵**只分解一次并缓存**,整段运行只发生 $O(\log_2(\tau_{\max}/\tau_{\min}))\approx7\text{–}8$ 次分解,保持直接解的高效。
- **格式**:自适应用**一阶后向欧拉**(无条件能量稳定;任意 $\tau$ 都精确守恒质量,故变步长不破坏这两条不变量,实测漂移 $\sim10^{-14}$)。要二阶请用固定步(变步 BDF2 未实现)。
- 停止于物理时间 `t_end`(而非固定步数),`max_steps` 为安全上限。

---

## 5. 完整算法

```
输入:JSON 配置(见 §7)
1.  生成网格 Mesh::getMesh(1/N);建立 DG 空间 FEM(ord);自由度 elem2dof;内部边 getEdge2Side
2.  装配 M_h = assembleMass_DG2D
            A  = assembleK_Poi2D + assembleIP_Poi2D(sigma, beta=1)
    断言 ||A·1||_inf ≈ 0
3.  组装常数分块矩阵 J = [[ (1/τ)M_h ,  M·A      ],
                         [ -(ε²A+S·M_h), M_h    ]]
    一次性 SparseLU 分解 J
4.  初值 c ← cbar + amp·U(-1,1)(逐自由度,seed 可复现);记录 mass0 = 1ᵀM_h c
5.  for n = 1 .. n_steps:
        b ← assembleNonlinearCH(c)             // F'(cⁿ) = cⁿ³ - cⁿ
        rhs ← [ (1/τ) M_h c ;  b - S M_h c ]
        x ← J 回代 rhs
        c ← x.head                              // cⁿ⁺¹
        每 save_every 步:写一帧 PPM、打印 E、mass-drift、max|c|
6.  ffmpeg 把 ch_frames/*.ppm 合成 MP4
```

每步成本 = 一次非线性装配(便宜)+ 一次稀疏 LU 回代;$J$ 的分解只做一次。

---

## 6. 初值与可视化

### 6.1 随机 spinodal 初值

逐自由度取 $c^0=\bar c+\text{amp}\cdot U(-1,1)$(`initSpinodal`,用固定 `seed` 可复现)。
$\bar c=0$ 给对称 50/50 的双连通相分离图样;`amp` 越大初始越随机/对比越强。

### 6.2 栅格化与配色

`writeFramePPM` 把 DG 场 $c_h$ 扫描到 $N_{\text{pix}}\times N_{\text{pix}}$ 像素:逐三角形求其像素包围盒,
对每个像素算重心坐标、用顶点值线性插值得到 $c$(对 P1 精确,P2 亦足够清晰),再过**发散
coolwarm 配色**(蓝 $=-1$,白 $=0$,红 $=+1$),写成二进制 **P6 PPM**。图像行从上(大 $y$)到下。

### 6.3 逐帧归一化(`normalize`)

初始随机场幅度只有 $\sim$`amp`(如 0.1),若用固定配色范围 $[-1,1]$,一开始几乎全白、看不出结构。
开启 `normalize=true` 后,每帧按 $m=\max|c|$ 对称缩放到 $[-1,1]$ 再配色(即用 $c_{\min}=-m,c_{\max}=m$):
$c=0$ 仍是中性白,但微弱的初始涨落会被拉满对比。**代价**是各帧颜色标度不同(早期 $\pm0.1$ 即满色,
后期 $\pm1$ 满色)。若要严格统一标度做定量比较,设 `normalize=false` 并用固定的 `cmin/cmax`。

### 6.4 合成视频

```bash
ffmpeg -y -framerate 25 -i ch_frames/frame_%05d.ppm \
       -c:v libx264 -pix_fmt yuv420p -crf 18 cahn_hilliard.mp4
```

`yuv420p` 要求宽高为偶数,故 `Npix` 必须是偶数(奇数会被程序自动 +1 并提示)。

---

## 7. 可调参数(JSON 配置)

参数全部写在 JSON 文件里,**改完无需重新编译**。程序按以下顺序定位配置:

1. 命令行给出的路径:`./build/cahn_hilliard my_config.json`;
2. 否则当前目录下的 `ch_config.json`(若存在);
3. 否则使用内置默认值。

JSON 支持 `//` 行注释与 `/* */` 块注释;任何缺省的键都回退到内置默认值。仓库根目录(`Poisson-cpp/`)
附带的 `ch_config.json` 是一个 **512×512、初值更随机、逐帧归一化**的较大算例。

### 7.1 参数表

| 键 | 含义 | 类型 | `ch_config.json` | 内置默认 | 约束 / 建议 |
|---|---|---|---|---|---|
| `ord` | 多项式次数 $P_k$ | int | 1 | 1 | $\ge1$;渲染用顶点值,P1 最快 |
| `N` | 每边网格数,$h=1/N$ | int | 128 | 64 | $\ge1$;越大越精细也越慢($\sim N^2$ 自由度) |
| `eps` | 界面厚度 $\varepsilon$ | double | 0.015 | 0.02 | $>0$;越小畴越细越多 |
| `sigma` | SIPG 内罚 | number 或 `"auto"` | `"auto"` | `"auto"` | `"auto"`$=3\,\text{ord}(\text{ord}+1)$;太小会失稳 |
| `beta` | SIPG 变体 | double | 1.0 | 1.0 | **保持 1**(对称,能量稳定) |
| `mob` | 迁移率 $M$ | double | 1.0 | 1.0 | $>0$;只缩放时间尺度 |
| `S` | 稳定化常数 | double | 2.0 | 2.0 | $\ge L/2$;2 安全到 $|c|\le1.29$ |
| `adaptive` | 自适应时间步(§4.7) | bool | **true** | false | 早期 $\tau$ 小、后期大;**让前几帧慢放** |
| `t_end` | 自适应:终止物理时间 | double | 0.30 | — | `adaptive=true` 时的停止条件(0→取 `n_steps`·`tau`) |
| `tau_min` / `tau_max` | 自适应:$\tau$ 上下界 | double | 2e-6 / 2e-4 | — | $\tau$ snap 到 $\tau_{\min}\!\cdot\!2^{j}$ 离散档,各档分解一次缓存 |
| `rel_change` | 自适应:每步目标相对变化 | double | 2.5e-3 | — | 控制节奏;越小越慢越精细 |
| `max_steps` | 自适应:步数安全上限 | int | 200000 | — | 防止配置失误跑不停 |
| `time_order` | 时间精度阶(固定步) | int | 1 | 1 | **1**=稳定化 BE;**2**=SBDF2(§4.6)。**自适应模式恒用 1** |
| `tau` | 固定步长 $\tau$ | double | 2.5e-5 | 2.5e-5 | 仅 `adaptive=false` 用;无 CFL,过大会过度光滑 |
| `n_steps` | 固定模式总步数 | int | 12000 | 16000 | 仅 `adaptive=false` 用;$T=$ `n_steps`$\cdot\tau$ |
| `save_every` | 出帧间隔(步) | int | 40 | 80 | $\ge1$;自适应下"每 K 步一帧"≈每帧等量变化 |
| `cbar` | 初值均值 $\bar c$ | double | 0.0 | 0.0 | 0 给对称相分离;非 0 偏向某相 |
| `amp` | 初值噪声幅度 | double | 0.10 | 0.05 | 远小于 1;越大初始越随机 |
| `seed` | 随机种子 | int | 20240601 | 20240601 | 改变得到不同随机实现 |
| `Npix` | 每帧像素 | int | 512 | 256 | **偶数**(奇数自动 +1) |
| `normalize` | 逐帧归一化到 $[-1,1]$ | bool | true | true | 见 §6.3 |
| `cmin`/`cmax` | 固定配色范围 | double | -1 / 1 | -1 / 1 | 仅 `normalize=false` 用 |
| `frames_dir` | PPM 输出目录 | string | `"ch_frames"` | `"ch_frames"` | 启动时清空其中旧 `*.ppm` |

### 7.2 取值建议与相互约束

- **网格分辨率 vs. $\varepsilon$**:90–10 扩散界面宽度 $W\approx4.16\,\varepsilon$;为分辨界面需约 4–6 个单元
  跨过它,即 $h\lesssim0.83\,\varepsilon$(P1)。例:$\varepsilon=0.02\Rightarrow N\gtrsim48$;
  $\varepsilon=0.015\Rightarrow N\gtrsim90$(取 128);$\varepsilon=0.01\Rightarrow N\gtrsim96\text{–}128$。
  **$\varepsilon$ 越小,要求的 $N$ 越大、粗化越慢**,总开销迅速上升。
- **时间步 $\tau$**:格式无条件能量稳定(没有 CFL 限制),$\tau$ 只受精度/过度光滑约束。经验
  $\tau\sim0.1\,h^2$;$\tau\gtrsim10^{-3}$ 会抹平初期快速相分离。$M$ 增大等价于加快时间,可相应减小 $\tau$。
- **总时间 $T$**:粗化遵循 LSW 标度 $L(t)\sim(M\,\sigma_{\mathrm{st}}\,t)^{1/3}$,
  $\sigma_{\mathrm{st}}=\tfrac{2\sqrt2}{3}\varepsilon$。要看到明显粗化,$\varepsilon=0.02$ 需 $T\approx0.3\text{–}0.5$。
- **稳定化 $S$**:默认 2;若日志出现 `[WARN energy-stability]`(即 $3\max|c|^2-1>2S$),调大 $S$。
- **性能**:自由度 $\approx N^2\cdot(\text{ord}{+}1)(\text{ord}{+}2)$;分块系统是其 2 倍。`SparseLU` 一次分解 +
  每步回代;$N=64$(P1)约 6 ms/步,$N=128$ 约 33 ms/步。

### 7.3 几个现成配方

- **轻量验证版(内置默认)**:`N=64, eps=0.02, tau=2.5e-5, n_steps=16000, Npix=256, amp=0.05`,约 100 秒。
- **较大演示版(附带 `ch_config.json`)**:`N=128, eps=0.015, Npix=512, amp=0.10, normalize=true,
  n_steps=12000`。
- **更锐界面(更贵)**:`eps=0.01, N=128, tau=1e-5, T≈0.6`。
- **更高阶**:`ord=2, N=48, sigma=18`(`"auto"` 会自动给 18)。

---

## 8. 运行时诊断与预期行为

每 `save_every` 步打印一行,例如:

```
step    240  t=0.00600  frame=   3  E=0.249279  mass_drift=-5.69e-18  max|c|=0.2077
```

- **`E`(离散自由能)应单调不增**——这是格式能量稳定性的运行时验证;
- **`mass_drift = 1ᵀM_h cⁿ − 1ᵀM_h c⁰` 应在机器精度量级**($\sim10^{-14}$);
- **`max|c|`** 从 $\approx$`amp` 增长到略大于 1(相分离趋向 $\pm1$);若 $3\max|c|^2-1>2S$ 会追加
  `[WARN energy-stability]`,提示调大 `S`。

典型演化:① 初始白噪声 → ② 自发相分离出双连通的红/蓝畴(界面变锐)→ ③ 缓慢粗化(小畴融合、
界面总长减少),对应 `E` 持续下降。

---

## 9. 收敛阶测试(MMS,`cahn_hilliard_convergence`)

`cahn_hilliard_convergence` 是一个**不产生视频**的可执行文件,用**制造解(method of manufactured
solutions, MMS)**定量验证时空收敛阶。取精确解 $c_e=\cos(\pi x)\cos(\pi y)\,e^{-t}$(它自动满足无通量
$\partial_{\mathbf n}c_e=\partial_{\mathbf n}\mu_e=0$,故内部边 SIPG 相容、收敛阶不被边界污染),只在 $c$ 方程
加一个源项 $g=\partial_t c_e-M\Delta\mu_e$(见 `ExactSolutionCH`;$g$ 可分离为 $g=a(t)\,g_1+a(t)^3 g_3$,
两块空间载荷只预计算一次)。初值取 $c_e(\cdot,0)$ 的 $L^2$ 投影,误差用 $14$ 阶求积对 $c_e(\cdot,T)$ 算
$L^2$、$H^1$ 范数。

- **空间扫描**:固定极小 $\tau$ 并用 SBDF2(令时间误差远低于空间误差),加密 $h$。预期 $L^2\sim h^{k+1}$、$H^1\sim h^k$。
- **时间扫描**:固定细的高阶空间($P_3,N{=}32$),加密 $\tau$,分别跑 `time_order` 1 与 2。阶数用**逐次差分
  Richardson** $\;\log_2\!\big(\lVert c_\tau-c_{\tau/2}\rVert_{M_h}/\lVert c_{\tau/2}-c_{\tau/4}\rVert_{M_h}\big)$
  读出(抵消公共的空间误差),直接误差表作交叉验证。

运行(约 3 分钟,Release 构建):

```bash
./build/cahn_hilliard_convergence
```

本机实测(节选):

```
空间 (SBDF2, tau=1e-5):           时间 (P3, N=32, Richardson):
  P1:  L2 rate -> 1.98  H1 -> 1.00    time_order=1 (BE)   : -> 1  (direct 0.94, Rich. 0.84, 渐近趋 1)
  P2:  L2 rate -> 3.07  H1 -> 2.00    time_order=2 (SBDF2): -> 2  (Richardson 1.95)
  P3:  L2 rate -> 4.01  H1 -> 3.00
                (即 L2 阶 = k+1, H1 阶 = k)
```

与理论完全一致:$P_k$ 元给 $L^2$ 阶 $k{+}1$、$H^1$ 阶 $k$;时间一阶 / 二阶分别给 $1$、$2$ 阶。
(时间扫描里**直接**误差在最细 $\tau$ 会触及固定空间误差的"地板"而变缓——这正是用 Richardson 逐次差分
作主指标的原因,它抵消了公共空间误差。)

---

## 10. 构建与运行

```bash
# 在 Poisson-cpp/ 目录下
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/homebrew
cmake --build build -j

# 演示:运行(自动读取 ch_config.json)+ 合成视频
./build/cahn_hilliard
ffmpeg -y -framerate 25 -i ch_frames/frame_%05d.ppm \
       -c:v libx264 -pix_fmt yuv420p -crf 18 cahn_hilliard.mp4

# 收敛阶测试(无视频)
./build/cahn_hilliard_convergence
```

依赖:Eigen(求解)、ffmpeg(仅演示合成视频用)。Eigen 3.4 与 5.0 均可,详见 [`../README.md`](../README.md) §2–3。

---

## 11. 代码结构

```
src/cahn_hilliard/
├── CahnHilliard.h / .cpp     质量阵、非线性载荷 F'(c)、能量/质量、spinodal 初值、PPM 栅格化、
│                             标量载荷、以及 CHIntegrator(一阶 BE / 二阶 SBDF2 时间推进)
├── Json.h                    无依赖的小型 JSON 读取器(支持 // 与 /* */ 注释)
├── ExactSolutionCH.h         收敛测试的制造解 c_e、梯度、源项 g(可分离)
├── ch_main.cpp               → cahn_hilliard:读配置 → 装配 → 时间推进 → 逐帧出图 + 诊断
└── ch_convergence_main.cpp   → cahn_hilliard_convergence:MMS 时空收敛阶测试(无视频)
```

复用:`src/common`(Mesh/FEM/Quadrature)、`src/ipdg`(`assembleK_Poi2D`、`assembleIP_Poi2D`)。

| 函数 / 类 | 文件 | 作用 |
|---|---|---|
| `assembleMass_DG2D` | `CahnHilliard.cpp` | DG 质量阵 $M_h$ |
| `assembleNonlinearCH` | `CahnHilliard.cpp` | 显式非线性载荷 $b_i=\int F'(c_h)\phi_i$ |
| `assembleScalarLoad` | `CahnHilliard.cpp` | 通用载荷 $\int f(x,y)\phi_i$(MMS 源项) |
| `computeEnergyCH` / `computeMassCH` | `CahnHilliard.cpp` | 离散能量 $E_h$ / 总质量 $\mathbf 1^\top M_h\mathbf c$ |
| `initSpinodal` | `CahnHilliard.cpp` | 随机 spinodal 初值 |
| `writeFramePPM` | `CahnHilliard.cpp` | 栅格化 + 配色 + 写 PPM |
| `CHIntegrator` | `CahnHilliard.{h,cpp}` | 固定步分块系统装配/分解 + 一阶/二阶单步推进(含 SBDF2 自举) |
| `CHAdaptiveStepper` | `CahnHilliard.{h,cpp}` | 自适应一阶 BE:按相对变化调 τ,离散档缓存分解(§4.7) |
| `cfgjson::parse` / `readFile` | `Json.h` | 读取 JSON 配置 |
| `ExactSolutionCH` | `ExactSolutionCH.h` | 制造解、源项 |
| `main` | `ch_main.cpp` / `ch_convergence_main.cpp` | 演示 / 收敛测试驱动 |

---

## 12. 参考文献

1. J. W. Cahn, J. E. Hilliard, *Free Energy of a Nonuniform System. I. Interfacial Free Energy*,
   J. Chem. Phys. 28 (1958) 258–267.
2. D. J. Eyre, *Unconditionally Gradient Stable Time Marching the Cahn–Hilliard Equation*,
   MRS Proc. 529 (1998) 39–46.(凸分裂)
3. J. Shen, X. Yang, *Numerical approximations of Allen–Cahn and Cahn–Hilliard equations*,
   DCDS-A 28 (2010) 1669–1691.(稳定化半隐格式与能量稳定性)
4. U. M. Ascher, S. J. Ruuth, B. T. R. Wetton, *Implicit-Explicit Methods for Time-Dependent
   PDEs*, SIAM J. Numer. Anal. 32 (1995) 797–823.(SBDF / 半隐 BDF IMEX)
5. D. N. Arnold, F. Brezzi, B. Cockburn, L. D. Marini, *Unified Analysis of DG Methods for Elliptic
   Problems*, SIAM J. Numer. Anal. 39 (2002) 1749–1779.(SIPG)
6. G. Wells, E. Kuhl, K. Garikipati, *A discontinuous Galerkin method for the Cahn–Hilliard equation*,
   J. Comput. Phys. 218 (2006) 860–877.(CH 的 DG 离散)
