# DG + 二阶 IMEX 求解二维可压缩 Euler 方程

本文档描述 `Poisson-cpp` 中的二维可压缩 Euler 方程求解器：间断 Galerkin（DG, `dP_k`）
空间离散 + 二阶 IMEX（隐显）Runge–Kutta 时间离散，并给出时空收敛阶验证与双马赫反射（Double
Mach Reflection, DMR）算例。

源码位于 `src/euler/`：
- `EulerDG.{h,cpp}` —— 静态库 `euler`：守恒变量 DG 装配、HLLC / Rusanov 数值通量、ARS(2,2,2) IMEX
  积分器、Persson–Peraire 人工粘性、Zhang–Shu 保正限制器、PPM 渲染。
- `euler_convergence_main.cpp` —— `euler_convergence`：等熵涡时空收敛阶验证（无视频）。
- `euler_dmr_main.cpp` —— `euler_dmr`：均匀网格双马赫反射影片 + 局部放大 + 数值纹影。
- `euler_dmr_amr_main.cpp` + `AdaptiveMesh.{h,cpp}` —— `euler_dmr_amr`：自适应网格双马赫反射（§9）。

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
逆质量阵），整个格式退化为稳定的二阶显式 RK。

## 4. 激波捕捉

**Persson–Peraire 亚单元人工粘性**（隐式处理 = IMEX 的"隐"）：用模态衰减传感器

$$ S_K=\frac{\|u-\Pi_{k-1}u\|^2_{L^2(K)}}{\|u\|^2_{L^2(K)}},\quad s_K=\log_{10}S_K, $$

经正弦斜坡映射到逐单元常数粘性 $\varepsilon_K$（阈值 $s_0=-4\log_{10}k$，幅值
$\varepsilon_0=c_{AV}(h_K/k)(|u|+c)_K$）。**指示量默认用压力**（`av_indicator=1`）：压力在激波处
跳变、但在接触间断/滑移线处连续，因此人工粘性只作用在激波上、**不抹平滑移线的 Kelvin–Helmholtz
卷起**（这正是双马赫反射尾部涡街所需）。粘性算子 $-\nabla\cdot(\varepsilon\nabla U)$ 用
$\{\varepsilon\}$-加权对称内罚（SIPG）离散；内罚需足够大以保证 $A$ 半正定（`sigma_ip` 默认 20）。

**Zhang–Shu 保正限制器**：每个 RK 级后，把多项式向其单元均值压缩
$U\leftarrow \bar U+\theta(U-\bar U)$，$\theta\in[0,1]$ 取到使所有体 / 边求积点上
$\rho\ge\varepsilon_\rho,\ p\ge\varepsilon_p$（压力沿线段为凹函数，闭式解二次方程）。该限制器
**精确保持单元均值**（守恒），光滑区不激活（无精度损失），并在初值投影上先调用一次以保证首步
通量良定义。

## 5. 边界条件（鬼状态）

| 边界 | 处理 |
|------|------|
| Dirichlet 入流 | exterior = 给定状态 |
| 超声速出流 | exterior = 内部状态（零梯度） |
| 反射 / 滑移壁 | exterior = 内部状态但法向动量取反 $m\mapsto m-2(m\!\cdot n)n$ |
| 远场 / 精确解 | exterior = 精确解（收敛性测试用） |

## 6. 时空收敛阶验证（等熵涡）

可执行 `euler_convergence`。等熵（Shu）涡是 Euler 的精确光滑解：定常涡叠加在均匀来流
$(\rho_\infty,u_\infty,v_\infty,p_\infty)=(1,1,1,1)$，强度 $\beta=5$，随来流平移不变形。域
$[0,10]^2$，精确解作为所有边界的鬼状态（等价于周期边界）。

**空间阶**（固定 $\mathrm dt=5\times10^{-4}$ 使时间误差可忽略，$T=1$，加密网格）：

| $k$ | $L^2(\rho)$ 收敛阶（HLLC 通量） |
|-----|--------------------------|
| 1 | 1.96, 2.01, 2.02（→ $k+1=2$） |
| 2 | 2.97, 3.02, 3.04（→ $k+1=3$） |
| 3 | 3.87, 3.89, 4.00（→ $k+1=4$） |

**时间阶**（固定 `dP3` 网格，逐次减半 $\mathrm dt$，Richardson 逐差）：**2.00, 2.00, 2.00**
—— ARS(2,2,2) 严格二阶。

## 7. 双马赫反射（均匀网格）

可执行 `euler_dmr`（参数见 `dmr_config.json`）。Woodward–Colella 标准设置：域 $[0,4]\times[0,1]$，
$\gamma=1.4$，马赫 10 激波与 $x$ 轴成 $60^\circ$，在 $x=1/6$ 触壁，$T=0.2$。

- 前方（未扰）状态 $(\rho,u,v,p)=(1.4,0,0,1)$。
- 后方（激波后）状态 $(\rho,u,v,p)=(8,\ 7.14471,\ -4.125,\ 116.5)$（$M_s=10$ 正激波 Rankine–Hugoniot）。
- 初始激波线 $x=1/6+y/\sqrt3$，线后为后方状态。
- 边界：左 = 后方入流（Dirichlet）；右 = 零梯度出流；下 $x<1/6$ 为后方入流、$x\ge1/6$ 为
  反射壁；上 = 随动激波 $x_s(t)=1/6+(1+20t)/\sqrt3$，左后右前。

输出**三段影片**（PPM→ffmpeg）：全域密度、滑移线区域**局部放大·密度**、同窗口**局部放大·数值纹影**
（灰度，最能看清涡街），外加最终静帧。密度与纹影均逐像素用真实 $\mathbb{dP}_k$ 多项式采样
（非顶点线性插值），纹影按渲染窗口自适应对比度以凸显涡街。

**默认配置**（`dmr_config.json`）：$\mathbb{dP}_2$、$n_y=120$（$\approx11.5$ 万三角形）、HLLC、
$c_{AV}=2$、$\sigma_{ip}=20$、`av_refresh=5`、$T=0.2$。

**性能**：残差、限制器、块对角质量均**按单元多线程并行**。人工粘性是**局部**的（只在激波单元
$\varepsilon>0$），故 $K=M+\gamma\,\mathrm dt\,A$ 在"激波活跃区"（本例仅 ~0.7% 自由度）之外**严格
块对角**——只对极小的活跃子块做 Cholesky **直接解**（精确，非迭代近似），其余退化为并行块逆质量。
`ny=120` 全程约 **14–15 分钟**（18 核）。

## 8. 运行

```bash
cmake --build build -j --target euler_convergence euler_dmr euler_dmr_amr
./build/euler_convergence                  # 时空收敛阶（约 2 分钟）
./build/euler_dmr     [dmr_config.json]     # 均匀网格 DMR 影片 + 放大 + 纹影
./build/euler_dmr_amr [dmr_amr_config.json] # 自适应网格 DMR（见 §9）
ffmpeg -y -framerate 25 -i dmr_frames/frame_%05d.ppm     -c:v libx264 -pix_fmt yuv420p -crf 16 dmr.mp4
ffmpeg -y -framerate 25 -i dmr_amr_frames/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 dmr_amr.mp4
ffmpeg -y -framerate 25 -i dmr_amr_frames_mesh/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 dmr_amr_mesh.mp4
```

构建依赖见 [poisson-cpp-build-env]，与现有求解器一致（Eigen3 + CMake）。

---

## 9. 双马赫反射的自适应网格加密（h-AMR）

`euler_dmr_amr` 在**同一套** DG 数值（§1–§5：HLLC、人工粘性、保正限制器、ARS IMEX、分区隐式解、
渲染）之上把网格变为**自适应**：只在激波 / 滑移线附近加密、光滑区粗化，用远少于均匀细网格的单元
达到同等效果。

**自适应机制**。采用**最新顶点二分（Newest-Vertex Bisection, NVB）**——一种在任意加密 / 粗化下
都保持网格**协调**（无悬挂结点）的二分加密，故求解器一字不改、每次重剖后仅重建一次。
解的传递：加密时父多项式**精确限制**到子单元（$k$ 次多项式在子三角形上仍是 $k$ 次 → 零误差、守恒），
粗化时 $L^2$ 投影回父空间（守恒）；实测每次重剖守恒量漂移 $\lesssim10^{-13}$。
**加密指示子用密度梯度**（$\eta_K=h_K\lVert\nabla\rho_h\rVert_\infty/\bar\rho_K$，与相邻单元密度跳量
取 max，迟滞防抖 + $N$ 层缓冲膨胀）：DMR 的滑移线是**接触间断**——压力连续、但密度跳变，故密度指示子
同时抓住激波与滑移线涡街，而 §4 的人工粘性传感器有意用**压力**以避开滑移线、不抹平涡街。
**CFL 控制**：每次重剖后按当前最细单元重算全局步长（$\mathrm dt\le\text{cfl}\cdot h_\min/((2k{+}1)\lambda_\max)$），
加密自动缩步、粗化放大；缓冲区使激波到达前已加密。生产 `cfl=0.70`（有效 CFL $\approx\text{cfl}/1.15=0.61$）；
扫描到 $T=0.30$：$\text{cfl}\le0.7$ 稳定，$\ge0.8$ 有**晚期爆破**风险（湍流可在一个重剖窗口内加速 $>1.15\times$ →
局部 CFL 尖峰），`0.7` 是验证过稳定 + 质量一致的最高值，要更稳可用 `0.6`。

**配置**（`dmr_amr_config.json`）：`base_ny`、`max_gen`（最细单元 = 基网格经 `max_gen` 次二分）、
`remesh_every`（重剖间隔步数）、`buffer_layers`、`th_ref` / `th_crs`（加密 / 粗化阈值，迟滞）、
`init_passes`（初值激波的初始加密遍数）。

**结果**。生产配置 `base_ny=25, max_gen=5`（最细 $h_\min\approx1/200$）、$T=0.30$、域 $[0,4.5]\times[0,1]$：
自适应网格约 **0.9–1.8 万**三角形随激波推进**动态平衡**（加密 ≈ 粗化），而同等最细分辨率的**均匀**网格需
约 20 万——最终密度场 / 纹影与均匀细网格**同量级**（激波清晰、三波点、滑移线 Kelvin–Helmholtz 涡街）。
输出在 §7 静帧之外多一段**网格叠加影片**（`dmr_amr_frames_mesh/`，把自适应三角网叠在密度场上）。

> **配置调优**：在**同样最细分辨率**下，"**粗基网格 + 深加密**"（`base_ny=25, max_gen=5`）比"细基网格 +
> 浅加密"（`base_ny=50, max_gen=3`）单元数更少——光滑远场几乎不花代价。实测全程 **603.8s → 407.3s
> （1.48×）**，最终帧效果一致（密度场逐像素均差 0.89/255）。这是 DMR 上最有效的提速手段。

**局部时间步进（LTS，`use_lts`，实验特性、默认关）**。`stepLTS` 实现守恒的 n-level 通量寄存器
子循环，让粗单元走更大步长。**机器精度守恒、收敛阶与 DMR 结果均不变**，但**不建议用于 DMR**：
发展后的双马赫反射中激波 + 涡街使细单元占比达 50%+，时间步类跨度坍塌，加上每宏步的额外开销，
实测反比 global-dt **慢约 1.8×**。LTS 仅在"细网格区域**小而固定**"的问题（如局部几何细节、稳定边界层）
上才划算，故对 DMR 默认关闭。

---

## 10. 更多可压缩 Euler 算例

在双马赫反射之外，新增 2 个经典算例，共用同一套 **h-AMR 场景驱动器**
`src/euler/euler_amr_scene.h`（`AMRScene` + `runAMRScene`）——把第 9 节的整套自适应机制
（NVB 自适应、守恒/精确传递、密度梯度+跨单元跳变指示子、CFL 控制、HLLC、Persson–Peraire
人工粘性、Zhang–Shu 保正、帧/网格叠加影片）封装为复用引擎，每个算例只提供四个钩子：初值
原始量场 `primIC(x,y)`、鬼状态边界 `bc`、边界标记 `tagEdge(mx,my)`、可选渲染掩膜 `inDomain`。

| 可执行 / 源文件 | 设置 | 说明 |
|---|---|---|
| `euler_riemann` | 单位方格四象限常状态，交点 `cross`（默认 0.8），**弱远场边界**，`config` 选 3/6/12 | 二维黎曼问题；config 3 中心蘑菇喷流、config 6 四臂剪切涡片、config 12 螺旋接触 |
| `euler_shock_bubble` | 域 $[0,2.5]\times[0,1]$，环境 $(\rho,u,p){=}(1,0,1)$，马赫 $1.22$ 激波，圆形密度跳变，**特征型右出流** | 单 $\gamma$ 代理。默认**重气泡** $\rho_b/\rho{=}3$（聚焦+射流，鲁棒）；**轻气泡** $0.138$（涡环）声速高、晚期成近真空细丝，较刚 |

**人工粘性传感器（关键、实测）**。双马赫反射是**激波主导**，用**压力传感器**（`av_indicator=1`）保滑移线锐利。
但这两个算例都是**接触/剪切主导**（二维黎曼的滑移线、激波-气泡的密度界面）：Zhang–Shu 保正只保正、
不抑制**保持正值**的高阶振荡，而本求解器**无斜率限制器**（`tvbLimit` 为空操作），故欠分辨的剪切层会
**无界增长**（$\rho\to10^5$–$10^6$ 后 NaN）。对策：

- **二维黎曼**：**密度传感器** + **低 Persson 阈值** `av_s0=-3`（及早起粘）+ `av_c=2`、`cfl=0.2`。三种配置均稳定且物理。
  （单换 Rusanov 通量**无效**——及早起粘的阈值才是关键。）
- **激波-气泡**：密度传感器，但要**温和**（默认 Persson 阈值、`av_c≈1.2`）。过激的低 `av_s0` 会在接触附近过度起粘，
  使 $\{\varepsilon\}$-加权 SIPG 在 `sigma_ip=20` 下丧失强制性（**反扩散→更快发散**）。默认用**重气泡**（声速低、稳健）。
  另有 `cfl_rho_floor`：仅给 dt 的波速估计设密度地板，防止虚假近真空节点把 $|u|+c$ 冲爆、压崩 dt（不碰物理）。

**二维黎曼边界条件**。自相似结构随 $x/t$ 增长，居中交点（`cross=0.5`）只能跑到标准短时（config 3 约 t=0.3）
就会冲出 $[0,1]^2$。默认改用 **Schulz-Rinne 离心交点 `cross=0.8`**：主结构往左下大区域发展，固定域也能干净
跑到 **t=0.8**（HOCUS-BVD / PyClaw 做法）。边界用**弱远场鬼态**：恒取该侧象限常状态、经 HLLC 逐特征处理
（超声速出流忽略鬼态、亚声速混入远场、超声速入流全施加），既修掉超声速入流壁（左/下）的伪侵入，也钉住
亚声速剪切出流壁（右）的入射声学特征，消除长时间的伪卷起。

**激波-气泡右出流**。入射激波在 $t\approx(x_b-x_{s0})/(M_s c)$ 抵右壁后，后激波出流是**亚声速**（$u\sim0.4<c\sim1.26$）；
裸零梯度（`return Uin`）在亚声速出流处反射 → 向左生长的低密度**稀疏波黑楔**最终吞域 NaN。固定背压也不行：远场
压力激波前是环境（1.0）、激波后是后激波（1.57），定值会在另一态注入伪波。正解是**特征型非反射出流**：出射不变量
（熵 $s$、切向速度、$R^+{=}u_n{+}2c/(\gamma{-}1)$）取自内部，入射 $R^-{=}u_n{-}2c/(\gamma{-}1)$ 取自环境远场，再重构鬼态。
它自适应环境/后激波两态、钉住 $R^-$ 不让压力漂移，从而**固定域 $[0,2.5]$ 干净跑到 t=2.4、无黑楔**。

运行命令见 README §8.1；参数见 `examples/{riemann,shock_bubble}_config.json`。
