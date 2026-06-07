# 不可压 Navier–Stokes 求解器:圆柱绕流、网格自动生成与时空离散方案

本文档详细介绍 `navier_stokes_cylinder` 算例:它求解什么方程、用什么计算区域与边界条件、如何**按
单一离散尺度 $h$ 自动生成高质量三角网格**,以及整套**空间(间断 Galerkin,任意阶 $\mathbb{dP}_k$)**
与**时间(高阶分裂 IMEX,二阶)**离散方案,最后给出 JSON 配置里所有可调参数和收敛性验证。
代码位于 `src/navier_stokes/`,复用了仓库已有的 2D DG 脚手架(`src/common`、`src/ipdg`)。

> 速览见 [`../README.md`](../README.md) 第 7 节;本文是完整版。

---

## 目录

1. [方程与物理背景](#1-方程与物理背景)
2. [计算区域与边界条件](#2-计算区域与边界条件)
3. [网格自动生成(DistMesh 式)](#3-网格自动生成distmesh-式)
4. [空间离散:间断 Galerkin](#4-空间离散间断-galerkin)
5. [时间离散:高阶分裂 IMEX(二阶)](#5-时间离散高阶分裂-imex二阶)
6. [完整算法](#6-完整算法)
7. [初值、可视化与受力系数](#7-初值可视化与受力系数)
8. [可调参数(JSON 配置)](#8-可调参数json-配置)
9. [收敛阶测试(Taylor–Green)](#9-收敛阶测试taylorgreen)
10. [构建与运行](#10-构建与运行)
11. [代码结构](#11-代码结构)
12. [参考文献](#12-参考文献)

---

## 1. 方程与物理背景

求解二维**不可压 Navier–Stokes 方程**(无量纲形式):

$$
\frac{\partial \mathbf u}{\partial t} + (\mathbf u\!\cdot\!\nabla)\,\mathbf u
= -\nabla p + \nu\,\Delta \mathbf u,
\qquad \nabla\!\cdot\!\mathbf u = 0,
\qquad \mathbf u=(u,v),
$$

其中 $\mathbf u$ 是速度、$p$ 是压力(已除以密度 $\rho=1$)、$\nu$ 是运动黏度。以来流速度 $U_\infty$
和圆柱直径 $D=2r$ 为特征量,Reynolds 数为

$$
\mathrm{Re}=\frac{U_\infty D}{\nu}.
$$

**算例:圆柱绕流(von Kármán 涡街)。** 均匀来流绕过一个圆柱,当 $\mathrm{Re}\gtrsim 47$ 时,
圆柱后方的剪切层失稳,交替地从上、下两侧脱落旋向相反的涡,形成经典的**卡门涡街**。
在 $\mathrm{Re}=100$ 时是稳定的周期性脱涡,可由以下无量纲量刻画:

- **Strouhal 数** $\mathrm{St}=fD/U_\infty$($f$ 为脱涡频率),文献值约 $0.16\text{–}0.17$;
- **阻力系数** $C_D=\dfrac{F_x}{\tfrac12\rho U_\infty^2 D}$,均值约 $1.3\text{–}1.4$(无界、二维);
- **升力系数** $C_L=\dfrac{F_y}{\tfrac12\rho U_\infty^2 D}$,随脱涡周期振荡,幅值约 $\pm0.3$。

本求解器在较粗网格、滑移侧壁的有限计算域上给出 $C_D\approx 1.5$、$\mathrm{St}\approx 0.18$,与文献
量级一致(侧壁阻塞与网格分辨率会使数值略偏高)。

---

## 2. 计算区域与边界条件

计算域是一个**矩形挖去一个圆盘**(圆柱截面):

$$
\Omega=\big([x_a,x_b]\times[y_a,y_b]\big)\setminus \overline{B(\mathbf c,r)}.
$$

四类边界、四种边界条件(用户要求的"经典入流-出流 + 两侧可滑移"):

| 边界 | 物理 | 速度 $u$ | 速度 $v$ | 压力 $p$ |
|---|---|---|---|---|
| 左 $x=x_a$ | **入流** | Dirichlet $u=U_\infty$ | Dirichlet $v=0$ | 高阶 Neumann(见 §5) |
| 右 $x=x_b$ | **出流**(do-nothing) | Neumann $\partial_n u=0$ | Neumann $\partial_n v=0$ | 高阶 Neumann |
| 上/下 $y=y_a,y_b$ | **可滑移**(对称面) | Neumann $\partial_n u=0$ | Dirichlet $v=0$ | Neumann $\partial_n p=0$ |
| 圆柱面 | **无滑移**(壁) | Dirichlet $u=0$ | Dirichlet $v=0$ | 高阶 Neumann |

要点:

- **可滑移侧壁**即对称(自由滑移)条件:法向速度为零($v=0$)、切向无应力($\partial_n u=0$);
  因此上下壁相当于对称面,圆柱处于近似无界来流中。
- **无滑移圆柱**是唯一产生尾流/涡街的来源(用户的入流-出流 + 滑移壁设定下,若圆柱也滑移则没有涡街)。
- **出流"do-nothing"**:速度自然边界 + 压力锚定 $p=0$,是分裂格式中最稳健的出流处理。
- 每条**速度分量**的边界类型可不同(滑移壁上 $u$ 是 Neumann、$v$ 是 Dirichlet),所以代码对 $u$、$v$
  各自组装一套带 Nitsche 边界的 Helmholtz 算子 $A_u$、$A_v$,压力另有一套 $A_p$。

边界条件在代码里用三个**逐边数组** `bcU`、`bcV`、`bcP`(每条网格边一个整型标记:内部 / Dirichlet /
Neumann)表达,加上随时间/空间给值的回调 `velDir`、`velAcc`(边界速度的时间导数,用于压力 BC)、
`presDir`。这样同一个积分器既能跑圆柱绕流,也能跑收敛测试的全 Dirichlet 盒子(§9)。

---

## 3. 网格自动生成(DistMesh 式)

用户要求"**根据空间离散尺度 $h$ 自动生成高质量三角网格**"。`src/navier_stokes/MeshGen.{h,cpp}`
从零实现了一个 **DistMesh**(Persson–Strang 2004)式的非结构网格生成器,只需给一个目标尺度 $h$。

### 3.1 几何的符号距离函数

矩形挖圆盘用符号距离函数(SDF)表示,流体域内为负:

$$
d(\mathbf x)=\max\!\big(d_{\text{rect}}(\mathbf x),\,-d_{\text{circ}}(\mathbf x)\big),\quad
d_{\text{circ}}=\lVert\mathbf x-\mathbf c\rVert-r,
$$

即"在矩形内 **且** 在圆外"。$\max(A,-B)$ 是布尔差(rectangle $\setminus$ circle)的标准 SDF 组合。

### 3.2 尺度场(网格加密)

尺度按到圆柱面的距离**渐变**:圆柱面附近 $\approx h$,远场放大到 $h\cdot\texttt{far\_ratio}$:

$$
f_h(\mathbf x)=\min\!\Big(h\cdot\texttt{far\_ratio},\; h\,\big(1+\tfrac{\texttt{grade}}{h}\,
\max(0,\lVert\mathbf x-\mathbf c\rVert-r)\big)\Big).
$$

于是涡街所在的近尾流被细密分辨,远场用大单元省算力。

### 3.3 算法(力平衡 + Delaunay)

1. **固定点**:矩形四角 + 圆周上 $\lceil 2\pi r/h\rceil$ 个等分点(让圆周被节点精确分辨)。
2. **初始布点**:六角点阵填充包围盒,按 SDF 保留域内点,再按 $1/f_h^2$ **拒绝采样**实现加密分布。
3. **迭代力平衡**:把网格的每条边看作只排斥的弹簧,期望长度 $L_0 \propto f_h$;
   合力推动节点,出域的点沿 SDF 梯度投影回边界;固定点不动。每隔若干步用 **Delaunay** 重剖分,
   按重心是否在域内筛掉圆内/域外三角形。收敛(内部点位移足够小)或到达上限即停。
4. **Delaunay 三角化**:`delaunayTriangulate` 是从零实现的 **Bowyer–Watson** 增量算法(超三角形 +
   外接圆判定 + 空腔重连),输出 CCW 定向三角形。

### 3.4 质量

生成器输出三角形数、节点数、**最小/平均内角**和面积范围。典型结果:**最小内角 $\sim 30^\circ$、
平均 $60^\circ$**——这正是高质量等边趋向的三角网格。圆周被节点精确分辨,矩形边界平直。

> 该生成器是 C++ 版独有的:仓库其它算例用结构化方形/六边形/圆盘网格,本算例需要"矩形挖洞 + 加密"
> 的非结构网格,故单独实现。`classifyEdges` 把每条边界边按位置判为入流/出流/壁/圆柱,供边界条件使用。

---

## 4. 空间离散:间断 Galerkin

$u$、$v$、$p$ 都在同一个**逐元 $\mathbb{dP}_k$ 间断空间** $V_h$ 上(`FEM::getDOF` 给出间断编号),
$k$ 任意。所有算子在 $V_h$ 上组装一次、整个时间推进里不变,因此**只分解一次**。

### 4.1 质量矩阵与黏性/压力的 SIPG 算子

- **DG 质量阵** $M$:逐元块对角、SPD(`assembleScalarMassDG`)。
- **SIPG $-\Delta$**:复用脚手架的体刚度 `assembleK_Poi2D` + 内部罚 `assembleIP_Poi2D`(只在内部边),
  得到无边界项的对称内罚算子;**Dirichlet 边界**用 Nitsche 法弱加(`assembleNitscheDirichlet`):

$$
a^{\text{bd}}_e(u,v)=\int_e\Big(-(\nabla u\!\cdot\!\mathbf n)\,v-\beta(\nabla v\!\cdot\!\mathbf n)\,u
+\tfrac{\sigma}{h_e}\,uv\Big)\,\mathrm ds,\qquad \beta=1\ (\text{对称 SIPG}).
$$

  按各自的 Dirichlet 边集合得到三套算子:

$$
A_u=K+\mathrm{IP}+B_D(\Gamma^u_{\!D}),\quad
A_v=K+\mathrm{IP}+B_D(\Gamma^v_{\!D}),\quad
A_p=K+\mathrm{IP}+B_D(\{\text{出流}\}).
$$

  其中 $\Gamma^u_D=\{\text{入流},\text{圆柱}\}$,$\Gamma^v_D=\{\text{入流},\text{圆柱},\text{壁}\}$。
  罚参数 $\sigma=\texttt{sigma\_fac}\cdot(k+1)^2$ 保证强制性(coercivity),三套算子皆 SPD。

### 4.2 弱一阶导算子 $G_x,G_y$:同时充当散度与梯度

用中心通量的 DG 弱导数算子(`assembleWeakGrad`),$(G_m p)_i=\int_\Omega(\partial_m p)\,\phi_i$:

$$
(G_m)_{ij}=\underbrace{-\int_K \phi_j\,\partial_m\phi_i}_{\text{体}}
+\underbrace{\int_{e_{\text{int}}}\{\phi_j\}\,n_m\,[\![\phi_i]\!]}_{\text{内部边,中心通量}}
+\underbrace{\int_{e_{\text{bd}}}\phi_j\,n_m\,\phi_i}_{\text{边界,单侧迹}} .
$$

这对矩阵**一物两用**:

- **散度**:$\;\mathrm{div}_h(\mathbf w)=G_x\,w_x+G_y\,w_y$(压力方程右端要用);
- **梯度载荷**:$\;(G_x p,\,G_y p)=\int(\nabla p)\phi$(投影修正要用)。

二者共用同一对矩阵 $\Rightarrow$ 离散散度与离散梯度互为伴随,投影自洽。可验证 $G_m\mathbf 1=0$(常数无导数)。

### 4.3 对流项:显式 Lax–Friedrichs 通量

对流项 $\mathbf N(\mathbf u)=\nabla\!\cdot\!(\mathbf u\otimes\mathbf u)$(散度/守恒形)显式处理,
弱形式载荷 $c_m=\int N_m\,\phi$:体项 $-\int_K(u_m\mathbf u)\!\cdot\!\nabla\phi$,面上用**局部
Lax–Friedrichs** 数值通量

$$
\widehat{F}_m\!\cdot\!\mathbf n=\tfrac12\big(F_m^-+F_m^+\big)\!\cdot\!\mathbf n
+\tfrac12\,\lambda\,(u_m^- - u_m^+),\qquad
\lambda=\max\big(|\mathbf u^-\!\cdot\!\mathbf n|,\ |\mathbf u^+\!\cdot\!\mathbf n|\big),
$$

$F_m=u_m\mathbf u$。边界处:入流/圆柱用给定 $\mathbf u_b$ 作外部态,出流用内迹(透射),滑移壁因
$\mathbf u\!\cdot\!\mathbf n=0$ 故对流通量取零。显式对流带来 CFL 限制 $\Delta t\lesssim
h/\big((2k+1)U\big)$。

---

## 5. 时间离散:高阶分裂 IMEX(二阶)

采用经典的**刚性稳定高阶分裂格式**(Karniadakis–Israeli–Orszag 1991)的二阶 BDF2/EX2 骨架,但压力
不再使用由中间速度散度得到的投影 PPE;默认使用**直接 Pressure-Poisson(PPE) 压力**加
**weak grad-div 稳定化**。记
$\gamma_0=\tfrac32$,历史 $\widehat{\mathbf u}^{\,*}=2\mathbf u^n-\tfrac12\mathbf u^{n-1}$,
外插 $\mathbf N^*=2\mathbf N^n-\mathbf N^{n-1}$(第一步用一阶 BDF1/EX1 自启动)。

每个时间步三子步:

**(1) 显式中间速度**(对流):

$$
\widehat{\mathbf u}=\big(2\mathbf u^n-\tfrac12\mathbf u^{n-1}\big)-\Delta t\,\mathbf N^*.
$$

**(2) 直接压力 Poisson(PPE)**:对动量方程取散度并用 $\nabla\!\cdot\!\mathbf u=0$,得

$$
\Delta p^{n+1}=-\nabla\!\cdot\!\mathbf N^*,
$$

离散为 $A_p\,p=G_x(M^{-1}c_x^*)+G_y(M^{-1}c_y^*)-\delta(G_xu^*+G_yv^*)+(\text{边界})$,
其中 `ppe_div_damping` $=\delta$ 默认取 10。这个 PPE 散度阻尼来自 PPE 重构里的
$\delta\nabla\!\cdot u$ 项,不改变连续不可压解,但能在离散层面耗散散度,且不耦合速度两分量。旧的投影压力
`pressure_mode="projection"` 仍可复现,但会把 $\nabla\!\cdot\!\widehat u/\Delta t$ 的离散散度误差放大进压力,
在 KIO Neumann/混合边界下容易出现约一阶到 $3/2$ 阶的时间降阶。
**高阶压力 Neumann 边界条件**(KIO)是保证时间二阶的关键——在速度 Dirichlet 边上,把动量方程投影到
法向:

$$
\frac{\partial p}{\partial \mathbf n}=\mathbf n\!\cdot\!\big(-\mathbf N^*+\nu\,\Delta\mathbf u^*-\mathbf a_b\big),
\qquad \mathbf a_b=\partial_t\mathbf u_b,
$$

其中黏性项默认用旋转形式 $-\nu\nabla\times\omega$ 代替 $\nu\Delta u$;对连续散度自由速度二者等价,
但离散速度不完全散度自由,旋转形式避免把 $\nabla(\nabla\!\cdot u)$ 泄漏到压力边界。滑移壁仍取
$\partial_n p=0$(对称)。圆柱、入流、出流均可用 KIO/Gresho-Sani 法向动量 Neumann;纯 Neumann 压力系统
用一个很小的质量矩阵 gauge 固定常数。

**(3) 压力梯度 + 黏性 Helmholtz**:默认保持快速的逐分量 Helmholtz 回代:

$$
\left[\frac{\gamma_0}{\Delta t}M+\nu A_{u/v}\right]\mathbf u^{n+1}
=\frac{1}{\Delta t}M\,\widehat{\mathbf u}^{\,*}-\mathbf c^*-G_{x/y}\,p^{n+1}
+\nu\,(\text{Nitsche Dirichlet 载荷}),
$$

`grad_div>0` 时可额外加入体积分 $(\nabla\!\cdot u,\nabla\!\cdot v)$ 稳定化,这会把速度改成 u-v
耦合 SPD 系统;默认 `grad_div=0`,用 PPE 散度阻尼代替,避免求解速度退化。$\mathbf c^*=2\mathbf c^n-\mathbf c^{n-1}$
是对流载荷外插。等价地满足
$(\gamma_0\mathbf u^{n+1}-\widehat{\mathbf u}^{\,*})/\Delta t+\mathbf N^*=-\nabla p^{n+1}+\nu\Delta\mathbf u^{n+1}$,
正是 BDF2/EX2 动量方程。

三套系统矩阵 $H_u=\tfrac{\gamma_0}{\Delta t}M+\nu A_u$、$H_v$、$A_p$ 都与时间无关、对称正定,用
`SimplicialLDLT` **各分解一次**;之后每步只做几次回代 + 一次显式对流装配(主要开销)。

---

## 6. 完整算法

给定 $\mathbf u^n,\mathbf u^{n-1}$ 和上一拍对流载荷 $\mathbf c^{n-1}$:

```
1.  c^n = 显式 LF 对流载荷(u^n)                 // assembleConvection
2.  c* = 2 c^n - c^{n-1}    (首步用 c^n)
    û* = 2 u^n - ½ u^{n-1}  (首步用 u^n)
3.  ûhat = û* - Δt · M⁻¹ c*                       // 中间速度(场)
4.  b_p = pNeumannHO(u^n,u^{n-1}) + pDirLift + Gx M⁻¹ c*_x + Gy M⁻¹ c*_y
          - δ(Gx u*_x + Gy u*_y)                 // δ=ppe_div_damping
5.  解 A_p p = b_p                                 // 压力(回代)
6.  rhs = (1/Δt) M û* - c* - G p + ν·velDirLift
    默认解 H_u/H_v;若 grad_div>0 则解 coupled [u,v] 系统
7.  更新历史:u^{n-1}←u^n, u^n←u^{n+1}, c^{n-1}←c^n
```

矩阵 $M,A_u,A_v,A_p,G_x,G_y$ 在构造时组装并分解;$M^{-1}$ 用 $M$ 的 LDLT 实现(块对角,廉价)。

---

## 7. 初值、可视化与受力系数

- **初值**:均匀来流 $\mathbf u=(U_\infty,0)$,叠加一个近尾流的小幅**横向扰动** $v$(幅值 `perturb`)
  以打破对称、**触发脱涡**;脱涡自持后扰动细节无关紧要。两者用 $L^2$ 投影到 DG 空间。
- **可视化**:逐帧把**涡量** $\omega=\partial_x v-\partial_y u$(由 $M\omega=G_x v-G_y u$ 投影得到)
  或速度大小 $|\mathbf u|$ 栅格化为 PPM。栅格化在每个像素处**按真正的 $\mathbb{dP}_k$ 基函数**求值
  (而非顶点线性插值),所以高阶场渲染光滑。圆内/域外像素留白(`inDomain` 掩膜)。用 ffmpeg 合成 MP4。
- **受力系数**:在圆柱边上积分应力 $\boldsymbol\sigma\!\cdot\!\mathbf n$,
  $\boldsymbol\sigma=-p\mathbf I+\nu(\nabla\mathbf u+\nabla\mathbf u^\top)$,$\mathbf n$ 为**圆柱外法向**
  (= 流体单元外法向取负),得 $\mathbf F=(F_x,F_y)$,再无量纲化为 $C_D,C_L$。
- **Strouhal 数**:记录 $C_L(t)$,取后半段升力上行过零点的平均周期 $T_s$,$\mathrm{St}=D/(U_\infty T_s)$。
  力的时间序列写入 `ns_forces.csv`,终端打印均值 $C_D$、$\mathrm{St}$ 与 $C_L$ 振荡区间。

---

## 8. 可调参数(JSON 配置)

参数全写在 JSON 里(改完无需重编译):优先读命令行路径,否则读当前目录 `ns_config.json`,再否则用内置默认。

| 键 | 含义 | 默认 |
|---|---|---|
| `ord` | $\mathbb{dP}_k$ 多项式阶 $k$($\ge1$,推荐 2 或 3) | `2` |
| `h` | **圆柱面上的目标单元尺度**(网格据此自动加密) | `0.07` |
| `far_ratio` | 远场单元尺度 / `h` | `7` |
| `grade` | 网格渐变率(越小尾流越细) | `0.13` |
| `radius`,`cx`,`cy` | 圆柱半径与圆心 | `0.5, 0, 0` |
| `xa,xb,ya,yb` | 矩形域(上下为滑移壁) | `-5,20,-7,7` |
| `Re` | Reynolds 数 $U_\infty D/\nu$ | `100` |
| `Uinf` | 来流速度 | `1.0` |
| `perturb` | 触发脱涡的初始横向扰动幅值 | `0.35` |
| `cfl` | 由 CFL 定步长 $\Delta t=\texttt{cfl}\cdot h/((2k+1)U_\infty)$ | `0.5` |
| `dt` | 时间步(>0 则直接用,覆盖 CFL) | `0`(=自动) |
| `t_end` | 终止时间(Re=100 时 $\sim$20 个脱涡周期) | `120` |
| `n_frames` | 影片帧数(导出 `save_every`) | `240` |
| `sigma_fac` | SIPG/Nitsche 罚参 $=\texttt{sigma\_fac}\,(k+1)^2$ | `8` |
| `field` | 渲染量:`"vorticity"` 或 `"speed"` | `vorticity` |
| `vort_clip` | 涡量色标范围 $\pm$ | `3.0` |
| `Wpix` | 帧宽(像素;高由渲染窗纵横比定) | `1280` |
| `render_xa..ryb` | 渲染窗口 | `-3,18,-5,5` |
| `frames_dir` | 帧输出目录 | `ns_frames` |

---

## 9. 收敛阶测试(Taylor–Green)

第二个可执行 `navier_stokes_convergence`(`ns_convergence_main.cpp`,无视频)用 **Taylor–Green
衰减涡**——不可压 NS 的一个**解析解**——在单位正方形上验证时空收敛阶:

$$
u=-\cos x\sin y\,e^{-2\nu t},\quad v=\sin x\cos y\,e^{-2\nu t},\quad
p=-\tfrac14(\cos 2x+\cos 2y)\,e^{-4\nu t}.
$$

四壁均给速度 Dirichlet(由解析解),压力在右壁锚 Dirichlet、其余壁用 §5 的高阶 Neumann。预期
速度 $L^2$ 误差**空间 $\sim h^{k+1}$、时间 $\sim\Delta t^2$**。实测(本仓库验证):

```
[Spatial]  dt=2.5e-4 固定
  P1: rate 2.77, 2.58, 2.15   → k+1 = 2
  P2: rate 3.19, 3.18, 3.05   → k+1 = 3
  P3: rate 4.43, 4.16         → k+1 = 4
[Temporal] P3,N=16 自收敛(Richardson,消去空间误差)
  ||u(8e-3)-u(4e-3)|| → 4e-3 比较给出 rate 1.94  ≈ 2
```

(时间用**自收敛**差分:同一网格上空间误差精确抵消,得纯时间阶;更小 $\Delta t$ 的速率受线性解
舍入地板限制而略降。)空间阶 $k+1$、时间阶 2 均被证实。

---

## 10. 构建与运行

构建见 [`../README.md`](../README.md) 第 3 节(需 Eigen3 与 CMake;ffmpeg 仅用于合成视频)。

```bash
# 圆柱绕流影片:自动读取 ns_config.json,帧写入 ./ns_frames/,再合成 MP4
./build/navier_stokes_cylinder                 # 或 ./build/navier_stokes_cylinder my.json
ffmpeg -y -framerate 25 -i ns_frames/frame_%05d.ppm \
       -c:v libx264 -pix_fmt yuv420p -crf 18 cylinder_vortex.mp4

# Taylor–Green 时空收敛阶测试(无视频):验证空间 k+1 阶、时间 2 阶
./build/navier_stokes_convergence
```

视频与帧图为生成产物,已在 `.gitignore` 中忽略。想要更细的尾流:减小 `h`/`grade`、增大 `far_ratio`
的同时延长域,但算量随之增加;显式对流的 CFL 会自动收紧 $\Delta t$。

---

## 11. 代码结构

```
src/navier_stokes/
├── MeshGen.{h,cpp}          DistMesh 式网格生成(SDF + 力平衡 + Bowyer–Watson Delaunay)、
│                            网格质量统计、边界边分类(入流/出流/壁/圆柱)
├── NavierStokes.{h,cpp}     DG 质量阵、Nitsche 边界、弱导算子 Gx/Gy、LF 对流、
│                            高阶压力 Neumann、BDF2/EX2 分裂积分器 NSIntegrator、
│                            涡量/受力、高阶 PPM 栅格化
├── ns_main.cpp              → navier_stokes_cylinder(圆柱绕流影片 + 力/Strouhal)
└── ns_convergence_main.cpp  → navier_stokes_convergence(Taylor–Green 收敛阶)
```

CMake:静态库 `navier_stokes`(复用 `dg_assembly`、`poisson_common`)+ 两个可执行。
压力/黏性/质量系统都对称正定,用 Eigen 内置 `SimplicialLDLT`(无需 SuiteSparse/CHOLMOD)。

---

## 12. 参考文献

1. G. E. Karniadakis, M. Israeli, S. A. Orszag, *High-order splitting methods for the incompressible
   Navier–Stokes equations*, J. Comput. Phys. **97** (1991) 414–443.(高阶分裂格式与压力边界条件)
2. D. L. Brown, R. Cortez, M. L. Minion, *Accurate projection methods for the incompressible
   Navier–Stokes equations*, J. Comput. Phys. **168** (2001) 464–499.(投影压力边界层与二阶压力修正)
3. R. R. Rosales, B. Seibold, D. Shirokoff, D. Zhou, *High-order finite element methods for a pressure
   Poisson equation reformulation of the Navier–Stokes equations*, J. Comput. Phys. **431** (2021) 110099.
   (PPE 重构、压力 Neumann 与散度阻尼)
4. M. Piatkowski, S. Müthing, P. Bastian, *A stable and high-order accurate discontinuous Galerkin based
   splitting method for the incompressible Navier–Stokes equations*, J. Comput. Phys. **356** (2018) 220–239.
   (DG 分裂格式、rotational pressure correction 与离散散度)
5. P.-O. Persson, G. Strang, *A simple mesh generator in MATLAB*, SIAM Review **46** (2004) 329–345.
   (DistMesh)
6. A. Bowyer, *Computing Dirichlet tessellations*, Comput. J. **24** (1981) 162–166;
   D. F. Watson, 同卷 167–172.(增量 Delaunay)
7. B. Cockburn, G. Kanschat, D. Schötzau, *A locally conservative LDG method for the incompressible
   Navier–Stokes equations*, Math. Comp. **74** (2005) 1067–1095.(DG 不可压流)
8. M. Schäfer, S. Turek, *Benchmark computations of laminar flow around a cylinder*, in
   *Flow Simulation with High-Performance Computers II*, Vieweg (1996) 547–566.(圆柱绕流基准)
