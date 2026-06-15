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
9. [圆柱+柔性尾巴 FSI](#9-圆柱柔性尾巴-fsinavier_stokes_filament)
10. [圆柱+蝌蚪漂移](#10-圆柱蝌蚪漂移)
11. [圆形碗+汤勺搅拌起涡](#11-圆形碗--汤勺搅拌起涡navier_stokes_spoon)
12. [收敛阶测试(Taylor–Green)](#12-收敛阶测试taylorgreen)
13. [构建与运行](#13-构建与运行)
14. [代码结构](#14-代码结构)
15. [参考文献](#15-参考文献)

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
| `frames_dir` | 颜色场视频的帧输出目录 | `ns_frames` |
| `outputs` | 视频组合(数组) | `["vorticity","flow"]` |
| `n_particles` | 粒子数(每个粒子是一条丝带) | `1500` |
| `trail_len` | 每条丝带保留的历史采样数 | `28` |
| `trail_stride` | 多少步记录一次轨迹(>1 让丝带更长) | `2` |
| `bg_dim` | 背景明度 0=纯黑 1=全亮 viridis | `0.12` |
| `flow_frames_dir` | 粒子帧输出目录 | `ns_flow_frames` |
| `particle_seed` | 粒子撒点 RNG 种子 | `12345` |

### 选择要输出的视频

`outputs` 是一个字符串数组,每项独立产生一份帧目录 + 一行 ffmpeg 命令,合法值:

- `"vorticity"` — 经典的 cool-warm 涡量场(原版视频)
- `"speed"` — viridis 配色的速度幅值 $\lvert\mathbf u\rvert$
- `"flow"` — **示踪粒子视频**:大量白色小点被 $(u,v)$ 实时平流,直观展示流体怎么走

常见组合:

```json
"outputs": ["vorticity"]            // 只要老视频
"outputs": ["flow"]                 // 只要直观流动演示
"outputs": ["vorticity", "flow"]    // 两个都要(默认)
"outputs": ["speed", "flow"]        // 速度幅值 + 粒子
```

> 兼容性:若 `outputs` 字段缺省,会回退到旧的 `field` / `render_flow` 语义。

> **示踪粒子视频**:粒子在每个时间步用 RK2 推进,离开渲染窗或被吸入圆柱后会在入口栅栏处
> 重新撒入,形成稳定的"烟线",绕过圆柱后卷入交替的 Karman 涡街。运行末尾会同时打印
> 所有通道对应的 ffmpeg 合成命令。

---

## 9. 圆柱+柔性尾巴 FSI(`navier_stokes_filament`)

把上面的圆柱算例加上一条**柔性尾巴**(flexible filament / splitter plate),
跑出**von Kármán 涡街锁频驱动的尾巴拍打**.

```bash
./build/navier_stokes_filament examples/ns_filament_config.json
```

### 9.1 数学模型

- **流体**:与圆柱算例完全一致(DG dP_k + BDF2/EX2 + 直接 PPE).圆柱网格不动.
- **结构**:1-D 离散 Cosserat 弹性杆(Bergou 等, SIGGRAPH 2008),
  $N+1$ 个节点 $\mathbf X_i$ 通过拉伸势 $V_s = \tfrac12 EA \sum_i (\ell_i-\ell_{0,i})^2/\ell_{0,i}$
  与弯曲势 $V_b = \tfrac12 EI \sum_i \theta_i^2/\bar\ell_i$ 耦合.根部 $\mathbf X_0,\mathbf X_1$
  钳支在圆柱后驻点,锁住位置和切向.时间推进用 **Newmark-β** ($\beta=1/4,\gamma=1/2$)
  + 牛顿迭代.

- **耦合**:
  - **`coupling: "oneway"`(默认)**:线性 Stokes 阻力
    $\mathbf F_k = -c_{\mathrm{drag}}(\mathbf V_k - \mathbf u_h(\mathbf X_k))\Delta s_k$
    通过**隐式**写到 Newmark 切线矩阵里(无条件稳定);流体不受反作用.对于"演示
    流场怎么把尾巴吹动"的目的,这是最稳健的选择.
  - **`coupling: "twoway"`(实验)**:Uhlmann 2005 的直接强制 IB,
    $\mathbf F_k = \alpha (\mathbf V_k - \mathbf u_h(\mathbf X_k))$
    既加给流体也以反向加给结构,满足牛顿第三定律.物理自洽,但分离式
    partitioned 实现对**轻杆** ($m^*<1$) 有 added-mass 不稳定(Causin-Gerbeau-Nobile 2005),
    需要把 `fil_mass` 调到 ≥ 1 并降 `ib_alpha`.

- **传输核**:不用 Peskin 离散 δ,直接用 DG 基函数本身做核.网格→杆是
  `meshToRod`(每个杆节点处的 $\phi_j(\mathbf X_k)$ 内积),杆→网格是其离散伴随
  `rodToMesh`,**离散一致性误差 ≈ 机器精度**(单元测试给出 2.3e-16).

### 9.2 数值验证

| 项 | 值 | 备注 |
|---|---|---|
| 杆孤立振动周期(`cosserat_test`) | 1.66% | vs 解析 Euler-Bernoulli 第一阶 |
| IB 转移核伴随恒等(`cosserat_test`) | 2.3e-16 | 机器精度 |
| 默认 flutter ($K_B=0.005$) 尖端幅度 | $A/D \approx 0.45$ | 大幅拍打 |
| Strouhal St_CL | 0.159 | 锁频在圆柱脱涡 |
| mean CD | 1.346 | 比裸圆柱 1.41 略低 (splitter-plate 减阻) |

### 9.3 主要 JSON 参数

继承 `ns_config.json` 的全部流体参数,新增:

| 键 | 含义 | 默认 |
|---|---|---|
| `fil_N` | 杆段数 | `60` |
| `fil_L` | 长度 / D | `2.0` |
| `fil_mass` | $m^* = \rho_s h / (\rho_f D)$ | `0.5` |
| `K_B` | $EI / (\rho_f U^2 D^3)$ | `0.005`(柔,大幅 flutter) |
| `K_S` | $EA / (\rho_f U^2 D)$ | `1e3`(准不可伸) |
| `fil_damp` | 结构 Rayleigh-mass 阻尼 | `0.02` |
| `fil_kick` | 初始横向半正弦扰动幅度 | `0.05` |
| `coupling` | `"oneway"` | `"twoway"` | `"oneway"` |
| `c_drag` | 单向阻力系数(只 oneway) | `8.0` |
| `ib_alpha` | 直接强制系数(只 twoway) | `0.5` |

### 9.4 已知局限

- 极刚杆 ($K_B \gg 0.1$) 时 Newton 求解器(用解析弯曲梯度+数值弯曲 Hessian)
  对初始大冲击不收敛,会发散.对于"流体推动柔性尾巴"的**演示意图**这没影响—
  那种刚度区间本来就该几乎不动.要做硬杆极限需要写解析 6×6 弯曲 Hessian 或换 BFGS.
- `coupling: "twoway"` 在轻杆下不稳定,文献已知;实现里留了 hook,但默认关掉.

---

## 10. 圆柱+蝌蚪漂移

有两个独立算例:

```bash
# 刚性尾巴:保留 36c0dba 的三自由度刚体漂移模型
./build/navier_stokes_tadpole examples/ns_tadpole_config.json

# 弹性尾巴:Cosserat 杆 + FE immersed-boundary 反作用力
./build/navier_stokes_tadpole_elastic examples/ns_tadpole_elastic_config.json
```

### 10.1 刚性尾巴(`navier_stokes_tadpole`)

这是 `36c0dba` 版本的保留算例:圆形头部 + 直线刚性尾巴组成一个三自由度刚体
$(x_h,y_h,\theta)$. 受力来自头部/尾巴采样点的线性 Stokes 阻力、可选 wake-refuge 力
$-k\nabla|\mathbf u|^2$ 和上游自推进力. 蝌蚪只单向采样流体,不反作用到流场.

### 10.2 弹性尾巴(`navier_stokes_tadpole_elastic`)

弹性尾巴是一个有限质量 Cosserat/Kirchhoff 杆,前两个节点钳在移动头部后缘. 与固定圆柱
filament 不同,这个钳支点每步同步**位置和速度**:

$$
\mathbf X_0=\mathbf x_h-r_h\mathbf e_h,\quad
\mathbf X_1=\mathbf X_0+\ell_0(-\mathbf e_h),\quad
\dot{\mathbf X}_i=\mathbf v_h+\omega_h\,\mathbf e_z\times(\mathbf X_i-\mathbf x_h).
$$

流固耦合采用文献中的 immersed-boundary 思路,但使用本代码的 FE/DG 传输算子:

- `meshToRod`:在杆节点 $\mathbf X_k$ 处插值 DG 速度 $\mathbf u_h(\mathbf X_k)$.
- `rodToMesh`:把节点力乘弧长权重 $\Delta s_k$ 后用伴随算子散布回 DG RHS.
- 流体反作用力为 relaxed direct forcing
  $$\mathbf F^{fluid}_k = {\alpha\over\Delta t}
      \left(\mathbf V_k-\mathbf u_h(\mathbf X_k)\right).$$
- 杆本身用隐式 drag
  $$\mathbf F^{rod}_k = -c_t
      \left(\mathbf V_k-\mathbf u_h(\mathbf X_k)\right)\Delta s_k$$
  加到 Newmark 方程中,所以尾巴感受真实局部流速,不再混合到头部速度.
- Newmark 步之后从钳支端向自由端做一次弧长投影,把每段投回 rest length. 这把尾巴
  保持为近不可伸长的弹性边界,避免轻尾巴在分区 IB 力下被轴向拉长,而弯曲自由度仍保留.

为抑制轻尾巴分区 added-mass 抖动,对 IB 力做幅值截断(`tad_ibForceCap`)和 1-2-1 沿尾平滑
(`tad_ibSmooth`). 当尾巴节点碰到圆柱 Dirichlet 接触区时,下一步的 IB 力会跳过该节点,
避免"几何投影 + IB 反力"自激.
弹性算例还对头部角速度设置 `tad_maxOmega`,并用较强 `tad_dampAng`,避免小圆头转动惯量把
采样力矩噪声放大成尾根快速闪动. `tad_maxOmegaAccel` 进一步限制相邻时间步的角速度变化,
避免姿态速度在数值上逐步翻符号.
`tadpole_elastic_diagnostics.csv` 会额外写出 `tailRmsCurv`、`tailMaxCurv` 和
`tailRoughness`;其中 roughness 是相邻离散曲率差的 RMS,用于排查节点级高频锯齿.

### 10.3 细网格稳定性(远场加密后的再调参)

最初的配置(`h=0.14, far_ratio=9`)在蝌蚪出生点附近的远场单元边长 ~0.9–1.4,比整只
蝌蚪还大:散布的 IB 反力在涡量图上画出单元状的不规则纹影,且全程在 t≈22 发散(NaN).
把远场整体加密到 `h=0.12, far_ratio=2`(全域单元 ≤0.24,~11k 三角形,最小角 37°,
全程 250 s 约 40 分钟)消除了纹影,但依次暴露出五层各自独立的失稳机制,逐一定参解决:

1. **IB 增益随网格失稳**:`F=(α/Δt)(V−u)` 的反馈增益正比于 α/(局部单元面积).
   粗网格上调好的 `α=0.15` 在细网格上 t≈6 即节点级发散;α≤0.08 配合 `ibSmooth`
   稳定(示例用 α=0.06 + 3 次平滑).另外头部对流体是透明的(只有单向阻力采样),
   尾根却钳在头上随头运动,根部 (V−u)≈O(1) 的集中力偶极会在头尾衔接处画出单元
   尺度的细碎纹影:`tad_ibRootRamp`(默认 6,示例 10)把反力沿根部若干节点线性
   渐变,且 1-2-1 平滑包含第一个受力节点,消除该端点突变.
2. **流体散度模式**:细网格上裸圆柱流(完全无蝌蚪)也会在涡街饱和时(t≈30,与 dt
   无关)爆出全域单元尺度的散度棋盘格;`ppe_div_damping` 需从 30 提至 60(其自身
   稳定窗约 60–150,300 会立刻发散).这与历史上 12→30 的修正是同一种病.
3. **头部翻滚 → 尾巴受压屈曲**:头部力矩只来自圆周采样阻力差,没有任何航向恢复
   机制;细网格解析出的锐利剪切层会像轴承一样持续扭转头部,直到蝌蚪尾巴朝前
   "倒飞",近不可伸长的软尾巴受轴向压缩载荷,远超欧拉临界值而屈曲(对 KB 在
   0.01–0.1 全区间均成立).`tad_dampAng=150 + tad_maxOmega=0.03` 等效于主动航向
   保持,根治.
4. **接触锤击**:头部以 ~1.1 的速度撞圆柱时,一步内的速度反射会猛拽尾巴钳支,
   在贴边单元注入 IB 力尖峰击穿 PPE.`ns_tadpole_elastic_main.cpp` 中加入了
   **软接近带**(厚 3·rHead,允许的法向接近速度随间隙线性收窄到零),把一次锤击化
   为几十步轻推;刹车减速度必须低于尾巴的动力屈曲极限 ~K_B π²/(ρ_t L³),这也是
   带宽取 3·rHead 与 `tad_maxSpeed=1.2` 的依据.`restitution=0`(不回弹).
5. **回流泡内的分区耦合失稳**:停靠后尾巴整体浸在回流泡(流速慢、被剪切层围住、
   流体响应强)中,双向 IB 在该区域无条件地指数发散(α、KB、阻尼如何调都只挪时
   刻不除根).新参数 `tad_ibWallMargin`(默认 0.9)在距圆柱面或任一外边界该距离内
   不散布 IB 反力(IB 直接施力与边界上的 Dirichlet 数据本质上互相打架),局部退化
   为单向耦合,杆照常感受流场并隐式求解.

注意:恒定 −x 推力 + 纯阻力头部的模型下,圆柱正后方对侧滑是**不稳定平衡**(推力
的切向分量把任何偏离放大),且钳支不向头部回传尾巴的拉力(无第三定律配对),所以
长时间后蝌蚪会绕过圆柱肩部、悬停在圆柱上方缓慢漂移——这是该力学模型的真实动力学,
不是数值伪影.若要真正"驻留尾迹",需要给头部加航向/位置闭环或把钳支反力回传头部.

详细文献记录见 [`elastic-tail-literature.md`](elastic-tail-literature.md).

---

## 11. 圆形碗 + 汤勺搅拌起涡(`navier_stokes_spoon`)

把一只**圆形汤碗**(整圈无滑移壁)装满静止的"汤",一把**汤勺**沿碗心绕一段圆弧
**搅动一下**后抬出;勺面是宽边迎流的桨叶,**在勺的两侧各卷起一个涡**,两个反号涡
配成一个**自推进涡偶极子(vortex dipole)**离开勺子,在碗内自行平移、沿弯壁滑掠并
逐步耗散。这正是俗话"搅一下汤、两边起旋涡"的流体力学版本。

```bash
./build/navier_stokes_spoon                 # 或 ./build/navier_stokes_spoon my.json
ffmpeg -y -framerate 25 -i ns_spoon_frames/frame_%05d.ppm \
       -c:v libx264 -pix_fmt yuv420p -crf 18 spoon_vortex.mp4
```

### 11.1 物理图像与文献定位

迎流平板/桨叶上,两个面各自长出**异号**的边界层涡量;流体绕过桨叶两端时把边界层涡量
卷成两个反向旋转的螺旋核——一端一个,符号相反(绕左端与绕右端的回转方向相反)。这两个
**起动涡(starting vortices)** 合起来就是一个自推进涡偶极子(Lamb–Chaplygin 偶极子的黏性
版本)。勺子**停下**时,异号的**停止涡**从同样的两端脱落并并入偶极子。脱体后两核互相
诱导,使偶极子整体平移:对强度 $\pm\Gamma$、间距 $d$ 的点涡对,平移速度
$U_{\text{dip}}=\Gamma/(2\pi d)$。在**封闭圆碗**里,偶极子撞上无滑移弯壁时,壁面成为**涡量
源**:生出异号边界层涡、分离卷成**次级涡**,与入射主涡重新配对而**反弹/绕行**——经典的
偶极子–壁碰撞。圆形(相对方形)壁对碗心无净力矩,故总角动量守恒良好,无自发 spin-up。

> 这与**茶叶悖论**(Einstein:搅拌停止后底部摩擦驱动的二次 Ekman 环流把茶叶聚到杯心)
> 是**两回事**——本算例刻画的是搅拌**当下**勺子尾流里脱落的涡,而非停转后的二次流。

最贴近的经典对照是**法向运动的脉冲起动平板**(Xu & Nitsche 2015;Taneda & Honji 1971):
迎流平板从静止脉冲起动后甩出一对反号涡。无量纲参数遵循该文献(详见 §11.4)。

### 11.2 计算区域与压力规范

区域是单个圆盘 $\Omega=\{r<R\}$,`MeshGen` 里新增 `BowlGeom`(圆盘 SDF)+ `generateBowlMesh`
(同一套 DistMesh 力平衡 + Bowyer–Watson,只是把"矩形挖圆"换成圆盘,逃逸点沿半径投回边沿)。
网格**按勺子轨迹自适应加密**:尺度场在勺子扫过的**环带** $r\in[\,r_{\rm pivot}-a-m_{\rm in},\,
r_{\rm pivot}+a+m_{\rm out}\,]$ 内取细尺度 $h_{\rm fine}$(外侧裕量 $m_{\rm out}$ 覆盖脱体偶极子外移
走廊),离开环带按 `grade` 渐变粗到 $h_{\rm fine}\!\cdot\!$`far_ratio`(静止的碗心与近壁带)——
与圆柱算例"按到圆柱距离渐变加密"是同一机制,只是把"到圆柱的距离"换成"到勺子环带的距离"。
(NS 分裂积分器的算子常系数、只分解一次,故这里用的是**轨迹自适应的静态加密**,而非每步重分网格
的动态 AMR。)
**整圈边界都是静止无滑移壁**:速度 $\mathbf u=\mathbf 0$(Nitsche Dirichlet),压力用 §5 的高阶
Neumann。没有入流/出流,于是压力 PPE 是**纯 Neumann——奇异**(压力只定到相差一个常数)。
解法是把**碗底一条边**(远离搅拌区)钉为 Dirichlet $p=0$ 锁定规范常数;这条小边对压力梯度
(进而对速度)无影响,只去掉零空间。其余与圆柱算例同:Direct-PPE + `ppe_div_damping`。

### 11.3 汤勺模型与浸入边界施力

勺面(浸入水平面的截面)是一片刚体桨叶,用一片**拉格朗日标记点云**填充(间距 $\Delta s\approx h$),
形状可选**椭圆**或**月牙(crescent / scoop)**:月牙是半径 $R_c=a/\sin\alpha$、角半宽 $\alpha$、带厚 $2b$ 的
圆弧带,其**凹口(开口)朝 $+\mathbf e_t$ 即运动方向**(曲率中心在前方,弧带在其后侧)——像一把朝运动
方向"舀"汤的勺子;两个尖角在 $\pm a$ 径向。因桨叶整体随 $\varphi$ 刚性旋转,**开口始终对着(旋转中的)
运动方向**。整片随**给定的搅拌律**绕碗心枢轴 $P$ 刚性旋转:

$$
\mathbf X_k(t)=P+R(\varphi(t))\,\mathbf r^0_k,\qquad
\mathbf V_k(t)=\omega(t)\,\mathbf e_z\times\bigl(\mathbf X_k-P\bigr),
$$

其中角速度取**二次(抛物)廓线**——从零起、到零止:

$$
\omega(\tau)=\omega_{\max}\,4\tau(1-\tau),\quad \tau=\frac{t-t_{\rm enter}}{t_{\rm stroke}}\in[0,1],\qquad
\omega_{\max}=\frac{3}{2}\frac{\Delta\varphi}{t_{\rm stroke}},
$$

其积分 $\varphi(\tau)=\omega_{\max}t_{\rm stroke}\!\left(2\tau^2-\tfrac43\tau^3\right)$ 恰好扫过总转角 $\Delta\varphi$。
$\omega$ 在搅动两端均为零(无速度跳变,不激出数值脉冲),峰值在中点。搅动结束后抬勺、撤约束、
不再绘制("把勺子提出汤面"),之后让流场**自由弛豫一段较长时间**再观察。默认**转两圈**:
$\Delta\varphi=720^\circ$($=4\pi$);$t_{\rm stroke}$ 按比例取 $21$,使峰值角速度仍是已验证稳定的
$\omega_{\max}\approx0.9$(桨叶速度 $\sim0.57$)——**要搅更多圈/更久,就让 $t_{\rm stroke}\propto\Delta\varphi$、
保持峰值速度(进而流动能量/CFL)不变;切忌只加大扫角而不加长时间(那会把桨叶加速、冲破对流 CFL)**。

**浸入边界:显隐混合(半隐)约束,而非显式直接力。** 早期版本用显式 direct forcing
$\mathbf F_k=\tfrac{\alpha}{\Delta t}(\mathbf V_k-\mathbf u_h)\mathrm dA$ 把力加到右端项;但该反馈增益 $\alpha$
经一致质量阵放回再隐式求解,存在与**阶数/网格/Re 都相关**的稳定上限——细网格、长快行程下它在
搅动峰附近**超指数发散**(与 §10.3 弹性尾巴同病)。现改用**把无滑移当作约束、隐式求解约束力**的
做法(Taira & Colonius 2007 的 IB 投影法;Kallemov 2016 的刚体约束;Goza & Colonius 2017 的强耦合)。
不再有"增益",因此对约束强度**无条件稳定**。

把无滑移 $\mathbf u_h(\mathbf X_k)=\mathbf V_k$ 当作和压力同地位的**拉格朗日乘子**约束,附加到隐式黏性
(Helmholtz)求解上,得到鞍点系统(逐分量)

$$
\begin{bmatrix} H & I^{\!\top}\\ I & 0\end{bmatrix}
\begin{bmatrix}\mathbf u^{n+1}\\ \mathbf f\end{bmatrix}=
\begin{bmatrix}\mathbf b\\ \mathbf V\end{bmatrix},\qquad
H=\tfrac{\gamma_0}{\Delta t}M+\nu A,
$$

其中 $I$ 是把 DG 场插值到标记的算子($=$`meshToRod`),$I^{\!\top}$ 是其基转置散布($=$`rodToMesh`,
此处不带面权重以保证对称),$\mathbf f$ 是标记上的约束力(刚体反力)。块 LU 消元给出

$$
\underbrace{(I\,H^{-1}I^{\!\top}+\varepsilon\mathbf I)}_{=\,G,\ \text{SPD}}\,\mathbf h=I\,\mathbf u^\ast-\mathbf V,\qquad
\mathbf u^{n+1}=\mathbf u^\ast-H^{-1}(I^{\!\top}\mathbf h),
$$

$\mathbf u^\ast=H^{-1}\mathbf b$ 是不带约束的黏性解。**Schur 矩阵 $G$ 只有"标记数 × 标记数"那么小**,且关键在于:
对它做 CG 时每次矩阵-向量乘 $G\mathbf p=I\,H^{-1}(I^{\!\top}\mathbf p)$ 只需**一次回代**——**复用已分解好的
`luHu_`/`luHv_`**(随勺子移动每步重建的只是这点积,本身不重新分解 $H$)。本仓库实测 CG 每步 **~8–10 次
迭代**收敛、约束残差 $\max|\mathbf u_h(\mathbf X_k)-\mathbf V_k|\sim10^{-4}$。代码:`NSIntegrator::setIBConstraint`
(每步在新标记位置武装约束)+ `step()` 内的 Schur 校正(`applyIBConstraint_`);驱动里 `spoon.active(t)`
在搅动期开约束、抬勺后撤除(`clearIBConstraint`)。$\varepsilon$ 是一个极小的 Tikhonov 正则,只在标记过密、
$I$ 行近线性相关时改善 $G$ 的条件数(这就是"显隐混合"里"隐"之外的那点"显")。

**标记间距取 $\Delta s\approx h$(每单元约一个标记)**:太密($\lesssim h/2$)使 $I$ 行近重复、$G$ 病态,
约束只能满足到几个百分点并在桨叶上画出条纹;太疏($\gtrsim 2h$)则封不住桨叶而漏流(Goza/Kallemov
建议 $1.5\!-\!2h$)。

**还剩的稳定约束(与 IB 无关)**:(1) 对流仍是**显式** EX2,$\Delta t$ 必须满足对流 CFL——精确约束让
桨叶边界层很锐、局部流速可超过桨叶速度,故取 `cfl≈0.3`(比绕流算例更保守);(2) 约束在压力投影**之后**
施加,会注入一点散度,由 `ppe_div_damping` 兜底。把多圈搅动放慢($t_{\rm stroke}=21$,峰值
$\omega\approx0.9$、桨叶速度 $\sim0.57$)把流动能量压到这两条都能从容满足的范围(enstrophy 峰 $\sim20$,
**全程无条件稳定**)。两圈一搅 + 长弛豫约 $2.9\times10^4$ 步、零 NaN、$\max|\!\int\!\omega|\sim10^{-4}$;
注意**抬勺后约束撤除,纯 NS 步比带约束的 IB 步快好几倍**,所以后段长观察很便宜。

### 11.4 推荐参数(对齐起动平板 / 受限偶极子文献)

以**桨叶长 $L=2a$** 和**搅动峰值速度 $U=\omega_{\max}r_{\text{pivot}}$** 定义 $\mathrm{Re}=UL/\nu$:

| 量 | 默认 | 区间 | 依据 |
|---|---|---|---|
| $\mathrm{Re}=UL/\nu$ | **500** | 250–1000 | Xu–Nitsche 基准:干净对称的层流涡对、剪切层可解析;$\gtrsim$2500 出现次级失稳需更细网格 |
| 桨叶长/碗半径 $L/R$ | **0.26** | 0.15–0.35 | 偶极子需 $\gtrsim$3–4 个桨叶长的净水程才能脱体自行平移 |
| 搅动弧 $\Delta\varphi$ | **720°**(两圈) | 45°–几圈 | 多圈把汤完整旋起;`tStroke` 须随之 $\propto\Delta\varphi$ 加长以保持峰值速度 |
| 桨叶形状 | **月牙**(开口朝运动方向) | 椭圆 / 月牙 | `spoon_shape`;月牙像朝运动方向"舀"汤的勺子 |
| 行程比(行程/$L$) | **~24**(两圈) | 2–4 单偶极子 / >4 持续脱涡 | $\gg$ 成形数 $F\approx4$(Gharib):沿多圈**持续脱涡** + 把整碗汤**旋起**(真实搅汤);抬勺后一组对转涡缓慢弛豫 |
| 速度廓线 | **二次(抛物)** $4\tau(1-\tau)$ | — | 从零起、到零止,峰值在中点;无速度跳变,免奇异压力尖峰 |
| 浸入边界 | **半隐约束**(Schur,复用 $H$ 分解) | — | 无增益、对约束强度无条件稳定(见 §11.3);取代了原显式直接力 |
| 标记间距 $\Delta s$ | $\approx h$(每单元约一个) | $h$–$2h$ | 太密 → $G$ 病态、桨叶起条纹;太疏 → 漏流(Goza/Kallemov 取 $1.5\!-\!2h$) |
| 对流 CFL `cfl` | **0.30** | 0.25–0.40 | 精确约束使边界层锐、局部流速高,比绕流算例更保守 |
| Schur 正则 $\varepsilon$ / 容差 | $\varepsilon{=}10^{-6}$,`ib_tol`${=}3\!\times\!10^{-3}$ | — | $\varepsilon$ 仅作病态保护;CG ~8–10 次迭代收敛,残差 $\sim10^{-4}$ |

### 11.5 验证与诊断

`spoon_diagnostics.csv` 每步记录 $t,\varphi,\omega,g,\mathbf F_{\text{stir}},\;\text{KE},\;\text{enstrophy}=\!\int\!\omega^2,\;\text{circulation}=\!\int\!\omega$。
定性上应看到:(1) 起步时勺两侧两个**异号**起动涡卷起(关于搅动轴近镜像对称);(2) 长行程下
沿扫掠弧**持续脱涡** + 碗内建立**大尺度环流**(行程比 $>$ 成形数),收尾再甩一对停止涡;
(3) 脱体涡**自平移**、沿碗弯壁**滑掠/反弹**并生出壁面次级涡。定量上:

- **总环量 $\int\omega\,\mathrm dA\equiv0$**(无滑移封闭域 $\oint\mathbf u\cdot\mathrm d\boldsymbol\ell=0$,
  Kelvin/对称)——实测全程 $|\!\int\!\omega|\lesssim 3\times10^{-5}$,为机器/求积噪声,是干净的守恒检验;
- **enstrophy** 在搅动末($t\approx t_{\text{stroke}}$)达峰后随黏性扩散单调下降,**KE** 缓慢衰减——
  与脉冲注入后自由衰减的图像一致。

默认配置(`spoon_config.json`:$\mathbb{dP}_2$、自适应网格 ~10.8k 三角形、65k DOF/分量、Re$=500$、
$\Delta\varphi=720^\circ$(两圈,月牙桨叶)、$t_{\rm stroke}=21$、$t_{\rm end}=42$、440 帧 × 涡量/示踪两路)整跑约半小时量级,全程**零 NaN**,CG 每步
~8–10 次迭代、约束残差 $\sim10^{-4}$、enstrophy 峰 $\sim20$。注:为把每步的 Schur 回代(复用 $H$ 分解)
压到可接受耗时,默认用 $\mathbb{dP}_2$ + 勺周自适应细化(而非全域 $\mathbb{dP}_3$);两者在该尺度的视觉
平滑度相当。

---

## 12. 收敛阶测试(Taylor–Green)

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

## 13. 构建与运行

构建见 [`../README.md`](../README.md) 第 3 节(需 Eigen3 与 CMake;ffmpeg 仅用于合成视频)。

```bash
# 圆柱绕流影片:自动读取 ns_config.json,帧写入 ./ns_frames/,再合成 MP4
./build/navier_stokes_cylinder                 # 或 ./build/navier_stokes_cylinder my.json
ffmpeg -y -framerate 25 -i ns_frames/frame_%05d.ppm \
       -c:v libx264 -pix_fmt yuv420p -crf 18 cylinder_vortex.mp4

# Taylor–Green 时空收敛阶测试(无视频):验证空间 k+1 阶、时间 2 阶
./build/navier_stokes_convergence

# 蝌蚪刚性/弹性尾巴漂移
./build/navier_stokes_tadpole examples/ns_tadpole_config.json
./build/navier_stokes_tadpole_elastic examples/ns_tadpole_elastic_config.json
```

视频与帧图为生成产物,已在 `.gitignore` 中忽略。想要更细的尾流:减小 `h`/`grade`、增大 `far_ratio`
的同时延长域,但算量随之增加;显式对流的 CFL 会自动收紧 $\Delta t$。

---

## 14. 代码结构

```
src/navier_stokes/
├── MeshGen.{h,cpp}          DistMesh 式网格生成(SDF + 力平衡 + Bowyer–Watson Delaunay)、
│                            网格质量统计、边界边分类(矩形挖圆 + 圆形碗 BowlGeom)
├── NavierStokes.{h,cpp}     DG 质量阵、Nitsche 边界、弱导算子 Gx/Gy、LF 对流、
│                            高阶压力 Neumann、BDF2/EX2 分裂积分器 NSIntegrator、
│                            **半隐浸入边界约束**(setIBConstraint:Schur/CG 复用 H 分解)、
│                            涡量/受力、高阶 PPM 栅格化
├── CosseratFilament.{h,cpp} 离散 Cosserat/Kirchhoff 杆 + Newmark 隐式步进
├── IBCoupler.{h,cpp}        DG↔杆的 FE immersed-boundary 插值/力散布伴随算子
├── Tadpole.{h,cpp}          36c0dba 刚性尾巴三自由度漂移模型
├── ElasticTadpole.{h,cpp}   移动头部钳支的弹性尾巴模型
├── Spoon.{h,cpp}            给定运动的刚性勺面桨叶(标记云 + 二次搅拌律;提供标记位置/速度给半隐约束)
├── ns_main.cpp              → navier_stokes_cylinder(圆柱绕流影片 + 力/Strouhal)
├── ns_filament_main.cpp     → navier_stokes_filament(圆柱+柔性 filament)
├── ns_tadpole_main.cpp      → navier_stokes_tadpole(刚性尾巴)
├── ns_tadpole_elastic_main.cpp → navier_stokes_tadpole_elastic(弹性 IB 尾巴)
├── ns_spoon_main.cpp        → navier_stokes_spoon(圆形碗+汤勺搅拌起涡)
└── ns_convergence_main.cpp  → navier_stokes_convergence(Taylor–Green 收敛阶)
```

CMake:静态库 `navier_stokes`(复用 `dg_assembly`、`poisson_common`)+ 各算例可执行
(`navier_stokes_cylinder` / `_filament` / `_tadpole` / `_tadpole_elastic` / `_spoon` / `_convergence`)。
压力/黏性/质量系统都对称正定,用 Eigen 内置 `SimplicialLDLT`(无需 SuiteSparse/CHOLMOD)。

---

## 15. 参考文献

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
9. C. S. Peskin, *The immersed boundary method*, Acta Numerica **11** (2002) 479–517.
   (沉浸边界方法综述)
10. L. Zhu, C. S. Peskin, *Simulation of a flapping flexible filament in a flowing soap film by the
    immersed boundary method*, J. Comput. Phys. **179** (2002) 452–468.(弹性 filament IB)
11. F.-B. Tian, H. Luo, L. Zhu, J. C. Liao, X.-Y. Lu, *An efficient immersed boundary-lattice
    Boltzmann method for the hydrodynamic interaction of elastic filaments*, J. Comput. Phys.
    **230** (2011) 7266–7283.(圆柱尾流中的弹性 filament / Karman gait)
12. J. C. Liao, D. N. Beal, G. V. Lauder, M. S. Triantafyllou, *The Karman gait: novel body
    kinematics of rainbow trout swimming in a vortex street*, J. Exp. Biol. **206** (2003) 1059–1073.
13. L. Heltai, F. Costanzo, *Variational implementation of immersed finite element methods*,
    Comput. Methods Appl. Mech. Engrg. **229–232** (2012) 110–127.(FE-IB 变分耦合与伴随传输)

**§11 汤勺搅拌起涡 / 受限涡偶极子(本算例新增):**

14. L. Xu, M. Nitsche, *Numerical study of viscous starting flow past a flat plate*, J. Fluid Mech.
    **776** (2015) 223–249.(法向脉冲起动平板甩出反号涡对——本算例桨叶的直接对照;Re∈[250,2000])
15. S. Taneda, H. Honji, *Unsteady flow past a flat plate normal to the direction of motion*,
    J. Phys. Soc. Japan **30** (1971) 262–272.(迎流平板起动涡对的经典实验)
16. H. J. H. Clercx, C.-H. Bruneau, *The normal and oblique collision of a dipole with a no-slip
    boundary*, Computers & Fluids **35** (2006) 245–279.(偶极子–无滑移壁碰撞的标准基准)
17. Yu. D. Afanasyev, *Formation of vortex dipoles*, Phys. Fluids **18** (2006) 037103.
    (脉冲/连续射流头部成偶极子,$U_{\text{dip}}/U_{\text{jet}}\approx0.5$)
18. H. J. H. Clercx, G. J. F. van Heijst, *Two-dimensional Navier–Stokes turbulence in bounded
    domains*, Applied Mechanics Reviews **62** (2009) 020802.(受限二维涡动力学综述:圆vs方、壁面涡量生成)
19. K. Schneider, M. Farge, *Decaying two-dimensional turbulence in a circular container*,
    Phys. Rev. Lett. **95** (2005) 244502.(圆形"碗":无滑移壁作为涡量源、终态轴对称单极)
20. M. Gharib, E. Rambod, K. Shariff, *A universal time scale for vortex ring formation*,
    J. Fluid Mech. **360** (1998) 121–140.(成形数 $F\approx4$:单股相干涡脱体的行程比上界)
21. G. J. F. van Heijst, J. B. Flór, *Dipole formation and collisions in a stratified fluid*,
    Nature **340** (1989) 212–215.(自推进涡偶极子范式)
22. M. Uhlmann, *An immersed boundary method with direct forcing for the simulation of particulate
    flows*, J. Comput. Phys. **209** (2005) 448–476.(运动刚体 direct forcing,标记间距 $\Delta s\approx h$)
23. P. Angot, C.-H. Bruneau, P. Fabrie, *A penalization method to take into account obstacles in
    incompressible viscous flows*, Numer. Math. **81** (1999) 497–520.(体积罚:运动实心桨叶的等价视角)

**§11.3 半隐(显隐混合)浸入边界约束(本算例采用):**

24. K. Taira, T. Colonius, *The immersed boundary method: a projection approach*, J. Comput. Phys.
    **225** (2007) 2118–2137.(把 IB 力当压力式拉格朗日乘子,分式步/Schur 求解——本法的范本)
25. B. Kallemov, A. P. S. Bhalla, B. E. Griffith, A. Donev, *An immersed boundary method for rigid
    bodies*, Comm. Appl. Math. Comput. Sci. **11** (2016) 79–141.(刚体约束的鞍点形式 + Schur/mobility
    预条件;隐式黏性 + 显式对流;标记间距与条件数)
26. A. Goza, T. Colonius, *A strongly-coupled immersed-boundary formulation for thin deforming
    surfaces*, J. Comput. Phys. **336** (2017) 401–411.(块 LU 把迭代限制在标记维子系统、复用流体分解;
    薄体;BiCGSTAB 解 Schul 2–8 次收敛)
27. E. P. Newren, A. L. Fogelson, R. D. Guy, R. M. Kirby, *A comparison of implicit solvers for the
    immersed boundary equations*, Comput. Methods Appl. Mech. Engrg. **197** (2008) 2290–2304.
    (隐式 IB 线性系统的若干求解器/预条件比较)

> 上述 §11.3 文献的 PDF 已下载到 [`literature/semi-implicit-ib/`](literature/semi-implicit-ib/)
> (含 Schneider 2015 体积罚综述、Engels 2015 运动体体积罚、Kallemov/Goza/Newren 等),供本算例的
> 显隐混合 IB 方案设计参考。
