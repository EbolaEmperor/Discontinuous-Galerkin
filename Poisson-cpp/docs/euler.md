# DG + 二阶 IMEX 求解二维可压缩 Euler 方程

本文档描述 `Poisson-cpp` 中新增的二维可压缩 Euler 方程求解器：间断 Galerkin（DG, `dP_k`）
空间离散 + 二阶 IMEX（隐显）Runge–Kutta 时间离散，并给出时空收敛阶验证与双马赫反射（Double
Mach Reflection, DMR）算例。

源码位于 `src/euler/`：
- `EulerDG.{h,cpp}` —— 静态库 `euler`：守恒变量 DG 装配、Rusanov 数值通量、ARS(2,2,2) IMEX
  积分器、Persson–Peraire 人工粘性、Zhang–Shu 保正限制器、PPM 渲染。
- `euler_convergence_main.cpp` —— 可执行 `euler_convergence`：等熵涡时空收敛阶验证（无视频）。
- `euler_dmr_main.cpp` —— 可执行 `euler_dmr`：双马赫反射影片 + 局部放大 + 数值纹影。

---

## 1. 控制方程

二维可压缩 Euler 方程（守恒形式）

$$ \partial_t U + \nabla\cdot F(U) = 0,\qquad U=(\rho,\ \rho u,\ \rho v,\ E)^\top, $$

$$ F_x=\begin{pmatrix}\rho u\\ \rho u^2+p\\ \rho uv\\ (E+p)u\end{pmatrix},\quad
   F_y=\begin{pmatrix}\rho v\\ \rho uv\\ \rho v^2+p\\ (E+p)v\end{pmatrix},\quad
   p=(\gamma-1)\Big(E-\tfrac12\rho(u^2+v^2)\Big),\ \gamma=1.4. $$

声速 $c=\sqrt{\gamma p/\rho}$。

## 2. DG 空间离散

每个守恒分量在破碎多项式空间 $V_h=\{v:\ v|_K\in P_k(K)\}$（节点 Lagrange 基，逐单元
独立自由度）中展开。弱形式：对每个单元 $K$、检验函数 $\phi$，

$$ \int_K \partial_t U\,\phi
   = \int_K F(U):\nabla\phi\,\mathrm dx \;-\; \oint_{\partial K} \hat H(U^-,U^+,n)\,\phi\,\mathrm ds. $$

界面数值通量默认采用 **HLLC**（Harten–Lax–van Leer–Contact）近似黎曼求解器，亦可切换为
**Rusanov / 局部 Lax–Friedrichs**：

$$ \hat H_{\text{Rus}}(U^-,U^+,n)=\tfrac12\big(F_n(U^-)+F_n(U^+)\big)-\tfrac12\,\lambda\,(U^+-U^-),\quad
   \lambda=\max\big(|u^-\!\cdot n|+c^-,\ |u^+\!\cdot n|+c^+\big). $$

> **为什么默认 HLLC**：Rusanov 用 $|u\cdot n|+c$ 对**所有**波(含接触波)统一耗散,会严重抹平
> **接触间断/滑移线**——而双马赫反射尾部涡街正是滑移线的 Kelvin–Helmholtz 失稳。HLLC 显式还原
> 接触波(波速 $S_*$),对滑移线几乎不加耗散,涡街才能清晰卷起。单步开销仅比 Rusanov 多约 10%,
> 且在光滑区与精确通量一致(实测 HLLC 让 $\mathbb{dP}_2$ 的 $L^2$ 阶从 ~2.6 提升到干净的 ~3.0)。

边界条件全部通过"鬼状态"（exterior state）经同一数值通量弱施加（见 §5）。

质量矩阵逐单元块对角，且由仿射映射有 $M_K=\mathrm{area}_K\cdot \hat M$（$\hat M$ 为参考单元
质量阵），因此显式更新只需一次参考逆质量阵 $\hat M^{-1}$ 按 $1/\mathrm{area}_K$ 缩放，开销极低。

## 3. 二阶 IMEX 时间离散：ARS(2,2,2)

把半离散系统写成

$$ M\,\dot U = \underbrace{F(U)}_{\text{显式：无粘通量}} \;-\; \underbrace{A(\varepsilon)\,U}_{\text{隐式：人工粘性扩散}} , $$

其中 $A(\varepsilon)$ 是 §4 的人工粘性 SIPG 算子。采用 L-稳定、二阶的
**ARS(2,2,2)**（Ascher–Ruuth–Spiteri 1997）加性 Runge–Kutta：无粘（非刚性）通量显式处理，
人工粘性（刚性抛物项）隐式处理，从而隐式扩散不再限制时间步。

记 $\gamma=1-\tfrac{\sqrt2}{2}$，$\delta=-\tfrac{1}{\sqrt2}$。每步只需 **两次通量求值 + 两次
隐式解**（两次隐式解共用同一矩阵 $K=M+\gamma\,\mathrm dt\,A$）：

```
f1 = F(U^n)
解 K U2 = M U^n + γ dt f1                          (stage 2, 节点 c=γ)
f2 = F(U2)
解 K U3 = M U^n + dt(δ f1 + (1−δ) f2) + (1−γ) dt (−A U2)   (stage 3, 节点 c=1)
U^{n+1} = U3                                        (两个表都 stiffly-accurate)
```

当人工粘性 $\varepsilon\equiv 0$（光滑流，收敛性测试）时 $A=0$，$K=M$（块对角，直接用
逆质量阵），整个格式退化为稳定的二阶显式 RK，稳定函数 $R(z)=1+z+z^2/2$，配合 Rusanov 通量
的迎风耗散对 DG 对流是稳定的。

隐式解：$K=M+\gamma\,\mathrm dt\,A$ 对称正定，用对角预条件的共轭梯度（CG）求解；$\varepsilon$
随流动演化变化，CG 无需重新分解（$M$ 占主导，迭代数很少）。$A$ 每 `av_refresh` 步重装一次。

## 4. 激波捕捉

**Persson–Peraire 亚单元人工粘性**（隐式处理 = IMEX 的"隐"）：用模态衰减传感器

$$ S_K=\frac{\|u-\Pi_{k-1}u\|^2_{L^2(K)}}{\|u\|^2_{L^2(K)}},\quad s_K=\log_{10}S_K, $$

经正弦斜坡映射到逐单元常数粘性 $\varepsilon_K$（阈值 $s_0=-4\log_{10}k$，幅值
$\varepsilon_0=c_{AV}(h_K/k)(|u|+c)_K$）。**指示量默认用压力**（`av_indicator=1`）：压力在激波处
跳变、但在接触间断/滑移线处连续，因此人工粘性只作用在激波上、**不抹平滑移线的 Kelvin–Helmholtz
卷起**（这正是双马赫反射尾部涡街所需）。粘性算子 $-\nabla\cdot(\varepsilon\nabla U)$ 用
$\{\varepsilon\}$-加权对称内罚（SIPG）离散，对每个守恒分量同样作用。

> ⚠️ **强制系数（coercivity）**：$\{\varepsilon\}$-加权 SIPG 必须有足够大的内罚才能保证 $A$ 半正定
> （否则 $A$ 出现负特征值 → 隐式项变成"反扩散" → 加大粘性反而更快爆掉）。本实现需要
> $\sigma=\sigma_{ip}(k+1)^2$ 中 $\sigma_{ip}\approx40$（而非常规的 $2\!-\!4$），DMR 默认即取 40。
> 由于激波单元集合随流动变化、$A$ 的稀疏结构随之改变，每次刷新都用 `SimplicialLDLT::compute`
> 重新分析+分解（复用旧符号结构会悄悄给出错误分解）。

**Zhang–Shu 保正限制器**：每个 RK 级后，把多项式向其单元均值压缩
$U\leftarrow \bar U+\theta(U-\bar U)$，$\theta\in[0,1]$ 取到使所有体/边求积点上
$\rho\ge\varepsilon_\rho,\ p\ge\varepsilon_p$（压力沿线段为凹函数，闭式解二次方程）。两个要点：
(1) 求积点集必须取**两个边定向**（dir=0 与 dir=1 是边上不同的物理点，通量两者都会用到）；
(2) 密度与压力限制器**复合**：$\theta=\theta_\rho\cdot\theta_p$（$\theta_p$ 是在密度限制后线段
$\bar U\to \bar U+\theta_\rho(U_h-\bar U)$ 上量得的，故应相乘而非取 $\min$）。该限制器精确保持单元
均值（守恒），并在初值投影（不连续初值的 Gibbs 过冲）上先调用一次以保证首步通量良定义。

## 5. 边界条件（鬼状态）

| 边界 | 处理 |
|------|------|
| Dirichlet 入流 | exterior = 给定状态 |
| 超声速出流 | exterior = 内部状态（零梯度） |
| 反射/滑移壁 | exterior = 内部状态但法向动量取反 $m\mapsto m-2(m\!\cdot n)n$ |
| 远场/精确解 | exterior = 精确解（收敛性测试用） |

## 6. 时空收敛阶验证（等熵涡）

可执行 `euler_convergence`。等熵（Shu）涡是 Euler 的精确光滑解：定常涡叠加在均匀来流
$(\rho_\infty,u_\infty,v_\infty,p_\infty)=(1,1,1,1)$，强度 $\beta=5$，随来流平移不变形。域
$[0,10]^2$，精确解作为所有边界的 Rusanov 鬼状态（等价于周期边界）。

**空间阶**（固定 $\mathrm dt=5\times10^{-4}$ 使时间误差可忽略，$T=1$，加密网格）：

| $k$ | $L^2(\rho)$ 收敛阶（HLLC 通量） |
|-----|--------------------------|
| 1 | 1.96, 2.01, 2.02（→ $k+1=2$） |
| 2 | 2.97, 3.02, 3.04（→ $k+1=3$，HLLC 低耗散消除了 LF 在偶次的前渐近偏低） |
| 3 | 3.87, 3.89, 4.00（→ $k+1=4$） |

**时间阶**（固定 `dP3` 网格，逐次减半 $\mathrm dt$，Richardson 逐差）：**2.00, 2.00, 2.00**
—— ARS(2,2,2) 严格二阶。

## 7. 双马赫反射（DMR）

可执行 `euler_dmr`（参数见 `dmr_config.json`）。Woodward–Colella 标准设置：域 $[0,4]\times[0,1]$，
$\gamma=1.4$，马赫 10 激波与 $x$ 轴成 $60^\circ$，在 $x=1/6$ 触壁，$T=0.2$。

- 前方（未扰）状态 $(\rho,u,v,p)=(1.4,0,0,1)$。
- 后方（激波后）状态 $(\rho,u,v,p)=(8,\ 8.25\cos30^\circ,\ -8.25\sin30^\circ,\ 116.5)
  =(8,\ 7.14471,\ -4.125,\ 116.5)$。（由 $M_s=10$ 正激波 Rankine–Hugoniot 严格导出。）
- 初始激波线 $x=1/6+y/\sqrt3$，线后为后方状态。
- 边界：左 = 后方入流（Dirichlet）；右 = 零梯度出流；下 $x<1/6$ 为后方入流、$x\ge1/6$ 为
  反射壁；上 = 随动激波 $x_s(t)=1/6+(1+20t)/\sqrt3$，左后右前。

输出**三段影片**(PPM→ffmpeg)：全域密度、滑移线区域**局部放大·密度**、同窗口**局部放大·数值纹影**
(灰度,最能看清涡街)，外加最终静帧。**数值纹影逐像素**用真实 $\mathbb{dP}_k$ 多项式梯度
$|\nabla\rho|$ 着色 $s=\exp(-\beta|\nabla\rho|/\max|\nabla\rho|)$，并按渲染窗口自适应对比度(放大窗口
里凸显涡街而非主激波)；密度场同样逐像素用 $\mathbb{dP}_k$ 基函数采样(非顶点线性插值)。

**默认配置**(`dmr_config.json`)：$\mathbb{dP}_2$、$n_y=120$（$\approx11.5$ 万三角形、有效分辨率
$\approx360$）、HLLC、$c_{AV}=2$、$\sigma_{ip}=20$、`av_refresh=5`、$T=0.2$。

**性能(纯效率优化,输出逐位不变)**:
- **按单元并行**(`std::thread`):残差(单元中心、无写冲突的 race-free 装配)、限制器、块对角质量
  apply,均按单元并行。
- **分区隐式解**:人工粘性是**局部**的(只在激波单元 ε>0),故 $K=M+\gamma\,\mathrm dt\,A$ 在"活跃"
  (激波∪邻居,本例仅 ~0.7% 自由度)与其余之间**严格块对角**——其余 99.3% 退化为并行块逆质量,只对
  极小的活跃子块做 Cholesky 分解/回代。隐式解 43→6.9 ms、分解刷新 15→2.6 ms(实测与整体分解解
  相对误差 ~1e-16,即**精确**而非近似)。
- **纹影渲染**:整数幂查表替代逐像素 `std::pow`;矩阵-向量乘合并($A\!\cdot\!U_2$ 一次稀疏乘、
  $M\!\cdot\!U$ 用并行块对角)。

合计:单步 147→74 ms、含渲染端到端约 **0.076 s/步**,`ny=120` 全程约 **14–15 分钟**(优化前 ~35 分钟,
未优化串行 ~100 分钟;18 核)。

> 运行后用打印的 ffmpeg 命令合成 `dmr.mp4` / `dmr_zoom.mp4` / `dmr_schlieren_zoom.mp4`。

## 8. 运行

```bash
cmake --build build -j --target euler_convergence euler_dmr euler_dmr_amr
./build/euler_convergence                  # 时空收敛阶（约 2 分钟）
./build/euler_dmr     [dmr_config.json]     # 均匀网格 DMR 影片 + 放大 + 纹影
./build/euler_dmr_amr [dmr_amr_config.json] # 自适应网格 DMR（见 §9），约 1/8 单元
ffmpeg -y -framerate 25 -i dmr_frames/frame_%05d.ppm     -c:v libx264 -pix_fmt yuv420p -crf 16 dmr.mp4
ffmpeg -y -framerate 25 -i dmr_amr_frames/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 dmr_amr.mp4
ffmpeg -y -framerate 25 -i dmr_amr_frames_mesh/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 dmr_amr_mesh.mp4
```

构建依赖见 [poisson-cpp-build-env]，与现有求解器一致（Eigen3 + CMake；CG 求解无需 CHOLMOD）。

---

## 9. 双马赫反射的自适应网格加密（h-AMR）

`euler_dmr_amr`（源 `AdaptiveMesh.{h,cpp}` + `euler_dmr_amr_main.cpp`）在**同一套** DG 数值（HLLC
通量、Persson–Peraire 人工粘性、Zhang–Shu 保正限制器、ARS(2,2,2) IMEX、分区隐式解、渲染）之上，
把网格变为**自适应**：只在激波 / 滑移线附近加密、在光滑区粗化，从而用**远少于均匀细网格**的单元
（约 $1/8$）达到与均匀 $n_y{=}150$ 同量级的效果。控制方程、通量、限制器、边界条件与 §1–§5 完全相同。

### 9.1 为什么用**协调**网格（最新顶点二分，NVB）

DG 的面通量装配 `Mesh::getEdge2Side` 依赖**协调网格**（每条内部边恰被两个三角形共享）；一旦出现
**悬挂结点**，该内部边会只看到一侧而被误判为边界，悄无声息地施加错误的鬼状态通量。为彻底复用已调好的
数值，本实现选择**最新顶点二分（Newest-Vertex Bisection, NVB）**——它在任意加密 / 粗化下都保持网格
协调，故求解器（残差、人工粘性、限制器、IMEX）**一字不改**，每次重剖后仅重建一次（廉价）即可。

**森林结构**：每个三角形是二分森林的一个结点，存为 $(v_0{=}\text{peak}, v_1, v_2)$，**加密边**取
$(v_1,v_2)$（即 peak 对面的边）。二分在该边中点 $m$ 处把三角形劈成两个共享 peak→$m$ 段的子三角形：
$$\text{child}_0=(m,v_0,v_1),\qquad \text{child}_1=(m,v_2,v_0),$$
两子的新加密边各取 peak 旁的一条"腿"——这是标准 NVB 子规则。

**基网格标号**：`makeRectMesh` 把每个方格按对角线劈成两个直角三角形；把**两个三角形的加密边都设为
该共享对角线**（peak = 对角线外的那个顶点）。于是每条网格边要么是其两侧三角形**共同的**加密边
（对角线），要么**都不是**（横竖网格线）——这正是 Stevenson 的"相容标号"条件，保证协调闭合（closure）
递归**有限终止**，且任一方格的首次二分就是局部的红色 $1\!\to\!4$ 细分。本基网格属**单一相似类**（全为
直角等腰三角形，长宽比恒为 $\sqrt2$），故加密任意层数都不退化。

**协调闭合（加密）**：要二分 $T$ 的加密边，必须**同时**二分其对面邻居 $F$。若 $F$ 的加密边不是这条共享
边，则**先递归加密 $F$**——NVB 保证 $F$ 二分后，紧邻该共享边的子三角形会把它变成**自己的**加密边，于是
有界递归后该边"相容可分"，两侧在**同一中点**一起二分——**永不产生悬挂结点**。

**粗化**：中点 $m$ 可移除，当且仅当**所有**含 $m$ 的当前叶子都以 $m$ 为 peak（即 $m$ 之下没有更细的悬挂）
**且**它们都被标记粗化；此时把 $m$ 周围**所有父对一起**反二分（父恢复为叶子）。这是二分的精确逆，保持协调。

### 9.2 解的传递：加密**精确**、粗化**守恒**

dP$_k$ 系数挂在叶子上。加密时把父多项式**精确限制**到每个子三角形（$k$ 次多项式在子三角形上仍是 $k$ 次，
故零误差、严格守恒）；粗化时把子多项式 $L^2$ 投影回父空间（含常数 → 积分守恒）。二者都化为**常数预计算
矩阵** $T_\text{restrict}[c]$、$T_\text{coarsen}[c]$，因为子→父重心坐标映射 $B_c$ 与面积比（二分恒为 $1/2$）
是固定的：
$$U_\text{child}=T_\text{restrict}[c]\,U_\text{parent},\qquad
  U_\text{parent}=\tfrac12\!\sum_{c}T_\text{coarsen}[c]\,U_{\text{child}_c}.$$
**自检**（`AMR_SELFTEST=1`）：一次 1 次多项式场经"全加密→全粗化"后，$\int\rho$ 漂移 $\sim10^{-15}$
（守恒）、$\int\rho^2$ 漂移 $\sim10^{-14}$（精确）。在真实 Mach-10 DMR 中（含激波处保正限制器），每次重剖
传递的守恒量漂移 $\lesssim10^{-13}$——即传递**严格守恒**。

### 9.3 加密指示子（关键：要同时抓住滑移线）

逐单元指示子
$$\eta_K=\frac{h_K\,\lVert\nabla\rho_h\rVert_{\infty,K}}{\bar\rho_K}\ \ \big(\sim\text{单元内相对密度变化}\big),
  \quad\text{并与相邻单元的密度跳量取 max}.$$
$\eta_K>\theta_\text{ref}$ 加密、$\eta_K<\theta_\text{crs}$ 粗化（$\theta_\text{crs}<\theta_\text{ref}$，
**迟滞**防抖），并对加密集做 $N$ 层**缓冲膨胀**，使运动激波在两次重剖间跑不出细网格区。
**为何用密度梯度而非压力**：DMR 的滑移线是**接触间断**——压力连续但**密度有跳变**（涡街正是密度图里看
得见的特征）；§4 的人工粘性传感器有意用**压力**以避开滑移线（不抹平涡街），而**加密**恰恰要抓住它，故
指示子用密度。这样激波（$\rho,p$ 都跳）与滑移线（仅 $\rho$ 跳）都被加密。在激波处 $\eta_K\approx\Delta\rho/\rho$
与 $h$ 无关，故会一路加密到 `max_gen` 上限；光滑区 $\eta_K\sim h\to0$ 自动粗化。

### 9.4 CFL 控制（加密时缩步长）

显式部分要求 $\mathrm dt\le\text{cfl}\cdot h_\min/((2k{+}1)\lambda_\max)$。每次重剖后用当前**最细**单元
重算全局步长，$h_\min=\min_K 2\,\text{area}_K/\text{longest edge}_K$（取**最长边上的高**，对二分产生的长宽比
稳健，而非最长边本身）。加密使 $h_\min$ 变小 → $\mathrm dt$ 自动变小；粗化则放大。重剖只在**整步之间**进行。
要点：因为缓冲区使单元在激波**到达前**就已加密，重剖当下重算的 $\mathrm dt$ 始终对应当前最细单元，故全程
有效 CFL 稳定在 $\approx0.45$（实测峰值）。

### 9.5 结果

$\mathbb{dP}_2$、base $n_y{=}50$、`max_gen=3`（激波处最细 $h\approx1/141\approx$ 均匀 $n_y{=}150$）、
$T{=}0.30$、域 $[0,4.5]\times[0,1]$：自适应网格约 **2.6–4 万**三角形随激波推进**动态平衡**（加密 ≈ 粗化），
而同等最细分辨率的**均匀**网格需约 **20 万**三角形——**单元数减少约 $5$–$8$ 倍**（每步代价同比下降），
最终密度场 / 纹影与 §7 的均匀 $n_y{=}150$ 结果**同量级**（激波清晰、三波点、滑移线 Kelvin–Helmholtz 涡街）。
输出在 §7 三段静帧之外，多一段 **网格叠加影片**（`dmr_amr_frames_mesh/`，把自适应三角网叠在密度场上）以
直观展示加密随激波 / 滑移线的自适应过程。

**配置**（`dmr_amr_config.json`，新增项）：`base_ny`、`max_gen`、`remesh_every`（重剖间隔步数）、
`buffer_layers`、`th_ref` / `th_crs`（加密/粗化阈值，迟滞）、`init_passes`（初值激波的初始加密遍数）。

**实现要点 / 踩过的坑**：
- `std::map<pair,array<int,2>>` 的新条目被**值初始化为 `{0,0}`**——会被误当作"已被叶子 0 占用"。边邻接
  表必须用 `{-1,-1}` 显式初始化（否则首次二分即误报"一边 >2 叶子"）。
- 粗化反二分父三角形后，必须把旧的两个子三角形标记为 `alive=false`，否则它们仍 `leaf()==true` 成为
  **游离叶子**，与复活的父三角形重叠 → 重建网格时同一条边出现 3 个叶子。
- 加密的协调闭合递归里，邻居 $F$ 的闭合可能**绕回**经由 $T$ 的"腿"边把 $T$ 自身也二分了；外层调用恢复
  时必须**重新检查 $T$ 是否还是叶子**，否则会对已二分的 $T$ 再 `split` 一次（重复二分）。
- 初值在**协调加密后的网格上重新投影解析 IC**（而非从粗网格传递），避免粗基网格把初始间断抹平。
- `EulerDG` 持有 `mesh/elem2dof/edge/...` 的引用：重剖时必须**先析构旧求解器**再改这些驱动器持有的对象、
  再构造新求解器（顺序），避免悬垂引用；每次重剖边界标号按几何重建（`x=1/6` 下边界入流/壁面分界由 BC
  回调按积分点 $x$ 现场判断，故只需把下边界整体标 `BOTTOM`）。

### 9.6 装配优化（`AMR_PROFILE=1` 实测驱动）

驱动器内置轻量计时（`AMR_PROFILE=1` 打印 step / remesh / render 三段及 remesh 子项）。实测每步计算约占
**75%**、重剖装配约 **23%**、渲染约 **2%**。基于实测对**装配（重剖）**做了 5 项优化（均保持数值不变——
`euler_convergence` 阶仍 2.0/3.0/4.0、AMR 自检仍精确守恒、涡街不变），把重剖从 17.5s 压到 6.1s（≈2.9×，
t=0.04 基准）。全量 DMR 端到端确认：**763.7s → 599.9s = 1.27×**（其中 ~1.18× 是纯装配提速，~1.08× 来自
形心指示子少 ~11% 单元、质量已验证不变）：

1. **指示子是重剖最大单项（5.7s！）**：原来每单元算 4 个点的多项式梯度（`std::pow`+分配）。改为**只在形心
   算一次**（面上的激波由相邻单元密度跳量项兜底）并**预存参考 d/dλ**（`grad = Dlam_t·(dphiC·ρ_block)`）→ 0.6s。
2. **`EulerDG::Mblk_`（块对角质量稀疏阵）每次构造都建，却只在 `EULERCHK` 调试路径用**（生产里 M 按单元用
   `Mref` 施加）→ 用 `getenv("EULERCHK")` 门控，EulerDG 构造 3.9s→1.6s。（参考量表与 NT 无关、很便宜，缓存
   它们收益甚微——真正的开销是 O(NT) 的 Mblk 装配。）
3. **`Mesh::getEdge2Side` 改排序法**（(min,max,elem,side) 记录排序后分组）替代 `std::map` → 输出逐位一致、
   约 2.5× 快（所有 DG 求解器共用）。
4. **`unordered_map`** 替换森林顶点重映射（buildMesh 1.2s→0.2s）与 EulerDG 的 `ee` 边索引表（打包 int64 键）。
5. **`FEM` 构造加 `withHessian` 开关**，DG 跳过只给 Argyris 用的 Voigt-Hessian `R`（FEM 0.6s→0.3s）；
   守恒诊断 `conservedTotals` 改为 `check_conservation` 选项（默认关，仅验证时开）。

**没动的**：隐式 AV 解已是最优（分区直接解 + 活跃小块缓存 Cholesky，CG 慢 ~6×）——换求解器无收益；每步计算
（占 88%）已高度优化。更大的潜在收益是**增量重剖**（每次仅 ~1% 单元变化却重建整个求解器），但侵入大、需重验，
暂未实现。

### 9.7 局部时间步进（LTS，实验特性，`use_lts`，默认关）

`EulerDG::stepLTS` 实现 **n-level 递归通量寄存器子循环**（Berger–Colella reflux + Krivodonova /
Constantinescu–Sandu）。按 dt 把单元分类：class $c$ 的步长 $2^c\,\mathrm dt_{\min}$（最细 class 0 走
$\mathrm dt_{\min}$、最粗走 $2^{L-1}\mathrm dt_{\min}$），类间做 **2:1 平衡**；递归
`ltsIntegrate(level,\mathrm dt)` 先推进本级一步、再把更细层级子循环两次（$\mathrm dt/2$），然后在每个
类界面做 **reflux** $U\mathrel{+}=M^{-1}(\Phi_f-\Phi_c)$。因 HLLC 守恒（$\hat H(a,b,n)=-\hat H(b,a,-n)$），
两侧时间积分通量一致 → **质量/能量机器精度守恒**（graded 网格 2 级与 3 级均实测漂移 $\sim10^{-14}$，
`LTS_CONS=1`）；收敛阶与 DMR 结果不变；Mach-10 全程稳定。

**紧凑逐类数组（关键提速）**：LTS 子步在**打包数组**（尺寸 $n_c\cdot\text{locDof}$，行 $c_i\cdot\text{locDof}+i
\leftrightarrow$ 单元 `cells[ci]`、自由度 $i$）上运算，每步开销 $O(\text{class})$ 而非 $O(n_{\rm Dof})$。组件：
`g2c_`（全局单元→紧凑下标，每子步设/清、平时全 $-1$）；`inviscidResidualMaskedCompact`（同级邻居从紧凑
`Uc` 经 g2c_、跨级从 `ltsSnap_`）；紧凑质量/质量逆/限制器；`limitCellCore`（Zhang-Shu 压缩抽成对
locDof×4 块运算，全局与紧凑限制器共用、数值不变——收敛阶仍 2/3/4）；AV 活跃块在 $n_a$ 空间求解：`Aaa_`=
$A$ 限制到活跃×活跃（refreshImplicit 中建）做紧凑 $A\,U_2$，活跃 Cholesky 解小的 $n_a\times4$ 右端。DG 自由度
块连续（$\text{elem2dof}(t,i)=t\cdot\text{locDof}+i$），全局↔紧凑经 $/\text{locDof}$ 映射。

**提速（干净背靠背，max_gen=4 DMR；紧凑重构后）**：global 16.9s｜2-level 16.0s（1.06×）｜
**3-level 11.8s（1.43×，确定性，246 宏步 vs global 979）**。`LTSPROF` 每宏步 `integrate` 由 47ms→38.5ms。
DMR 密度场与 global-dt 逐像素均差 0.012/255（无像素差 >10）——激波/三波点/滑移线一致。

**5 个硬坑（n-level / 紧凑特有）**：① 跨级 ghost 用**宏步起点快照**（不能 live `U_` 前向读未来态，激波处不稳）；
② **AV 限制在 class 0**（refreshImplicit 前把非 class-0 的 $\varepsilon$ 清零）；③ 临时数组**零初始化**（紧凑后
天然满足，每行都写）；④ **AV 每宏步刷新**（宏步含 $2^{L-1}$ 个细步，按 av_refresh 冻结会让激波跑出粘性带）；
⑤ **AV 活跃集并非 class 0 的子集**——SIPG AV 把 class-0 激波单元耦合到其 class-1 邻居（面罚用 class-0 侧的
$\varepsilon>0$），故 `activeDofs_` 含 class-1 自由度；全 nDof 旧路径天然容忍（其右端为 0、结果丢弃，但仍在
$K_{aa}$ 内耦合并贡献 $A\,U_2$）；紧凑路径其紧凑行记为 $-1$（越界即 NaN，光滑 `LTS_CONS` 因 $n_a{=}0$ 测不出），
须把它们活跃右端取 0、保留其解供 $A\,U_2$、把 $-(1{-}g)\,dt\,(A U_2)$ 喂进第 2 级活跃右端的越级行、并丢弃回写。
每宏步重新分类（追激波）、class 0 膨胀 2 圈。
