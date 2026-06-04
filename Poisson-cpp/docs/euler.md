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
cmake --build build -j --target euler_convergence euler_dmr
./build/euler_convergence                 # 时空收敛阶（约 2 分钟）
./build/euler_dmr [dmr_config.json]        # DMR 影片 + 放大 + 纹影
ffmpeg -y -framerate 25 -i dmr_frames/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 dmr.mp4
```

构建依赖见 [poisson-cpp-build-env]，与现有求解器一致（Eigen3 + CMake；CG 求解无需 CHOLMOD）。
