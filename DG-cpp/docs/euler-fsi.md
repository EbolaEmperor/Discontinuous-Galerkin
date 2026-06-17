# Euler-ALE 流固耦合基础设施

本文档描述 `src/fsi/` 中新增的可压缩 Euler 流固耦合基础设施。当前版本包含一维和结构化二维
保守 ALE 核心：移动网格 Euler、移动壁面通量、弹簧活塞固体、收敛阶验证，以及 shock tube /
shock channel + spring piston 可视化算例。它的目标是把 FSI 的守恒接口、网格运动接口、验证方式、
能量账本和出图流水线固定下来；现有非结构二维 `EulerDG` 也已经接入 direct moving-wall flux hook，
后续再把同一套 ALE 几何推进到非结构二维 DG、刚体 / 弹性壁面算例。

源码位于 `src/fsi/`：

- `EulerALE1D.{h,cpp}`：静态库 `euler_fsi`，包含保守 ALE 有限体积核心、边界条件、移动壁面通量、
  弹簧活塞模型、误差计算和 PPM 渲染。
- `EulerALE2D.{h,cpp}`：结构化二维保守 ALE 有限体积核心，含 x/y ALE 通量、固定/移动滑移壁、
  二维误差计算、密度 PPM 渲染和带活塞结构层的 PPM 渲染。
- `FSIDiagnostics.{h,cpp}`：气体壁面做功、阻尼耗散、气体能量漂移与耦合总能量漂移账本，供 1D/2D
  活塞算例复用。
- `../euler/EulerDG.{h,cpp}`：非结构二维 DG 侧的 ALE 法向通量、exact moving-wall flux、
  direct boundary-flux residual / time-step hook，供后续非结构 ALE 网格运动复用。
- `euler_fsi_convergence_main.cpp`：`euler_fsi_convergence`，验证 GCL / 均匀流保持、移动网格空间二阶、
  固定网格时间二阶。
- `euler_fsi_flux_test_main.cpp`：`euler_fsi_flux_test`，直接验证 1D/2D ALE 等状态通量、移动壁面通量、
  非结构 DG direct boundary-flux 路径与 `FSIEnergyBudget` 公式。
- `euler_fsi_convergence2d_main.cpp`：`euler_fsi_convergence2d`，验证二维 GCL、移动网格空间二阶、
  固定网格时间二阶。
- `euler_fsi_piston_main.cpp`：`euler_fsi_piston`，shock tube + spring-mounted piston 演示，输出帧图、
  诊断 CSV、最终静帧，并打印 ffmpeg 合成命令。
- `euler_fsi_piston2d_main.cpp`：`euler_fsi_piston2d`，二维 shock channel + spring piston 演示，
  输出二维密度 + 活塞结构帧、诊断 CSV 和最终静帧。
- `euler_fsi_piston2d_showcase_main.cpp`：`euler_fsi_piston2d_showcase`，强激波 + 轻质长行程弹簧活塞
  展示算例，使用高对比渲染、活塞轨迹和可见弹簧，专门用于视觉演示。
- `euler_fsi_vortex_showcase_main.cpp`：`euler_fsi_vortex_showcase`，强激波穿过正弦扰动密度界面后与
  移动弹簧活塞反射波耦合，触发 Richtmyer-Meshkov / Kelvin-Helmholtz 风格涡卷；渲染中红/青分别表示
  正/负涡量，适合展示“流固耦合 + 涡结构”的视觉效果。

---

## 1. 控制方程与 ALE 通量

一维 Euler 方程写成守恒形式

$$ U_t + F(U)_x = 0,\qquad
   U=(\rho,\rho u,E)^T,\quad
   F=(\rho u,\rho u^2+p,(E+p)u)^T, $$

其中

$$ p=(\gamma-1)\left(E-\frac12\rho u^2\right),\qquad \gamma=1.4. $$

在移动网格界面速度为 $w$ 的 ALE 框架下，界面通量为

$$ F_{\mathrm{ALE}}(U;w)=F(U)-wU. $$

内部界面使用 ALE Rusanov 通量

$$ \hat F =
\frac12\left(F(U_L)-wU_L+F(U_R)-wU_R\right)
-\frac12\lambda (U_R-U_L), $$

其中

$$ \lambda=\max(|u_L-w|+c_L,\ |u_R-w|+c_R). $$

有限体积未知量存的是**移动单元上的守恒积分**，不是固定体积平均值：

$$ Q_i(t)=\int_{x_{i-1/2}(t)}^{x_{i+1/2}(t)}U(x,t)\,dx. $$

所以时间更新直接用

$$ \frac{dQ_i}{dt}=\hat F_{i-1/2}-\hat F_{i+1/2}. $$

这一点是后续做 FSI 的关键：网格压缩 / 拉伸不会凭空产生质量、动量和能量。

二维结构化版本使用同样的思想。守恒变量为

$$ U=(\rho,\rho u,\rho v,E)^T, $$

移动竖直面使用 $F_x(U)-w_xU$，移动水平面使用 $F_y(U)-w_yU$；未知量仍是移动网格单元上的
守恒积分。当前二维实现先支持结构化四边形通道 / 周期盒，用来验证二维 GCL、空间/时间阶和移动壁面
能量账本，作为后续迁入非结构 `EulerDG` 的保守 ALE 参考实现。

## 2. 移动壁面与固体侧

壁面速度为 $V$ 时，滑移不可穿透条件给出壁面通量

$$ \hat F_{\mathrm{wall}}=(0,\ p,\ pV)^T. $$

二维任意法向版本在 `EulerDG` 和 `EulerALE2D` 中统一为

$$ \hat F_{\mathrm{wall}}=(0,\ pn_x,\ pn_y,\ p\mathbf V_w\cdot\mathbf n)^T, $$

其中 $\mathbf n$ 是流体域外法向，$\mathbf V_w$ 是壁面速度。这个通量是 ALE 壁面通量，不依赖
静止网格 ghost state，因此可以直接作为 FSI 的压力载荷和气体做功来源。

当前固体侧实现为单自由度弹簧活塞

$$ m\ddot x = A(p_{\mathrm{gas}}-p_{\infty})
-k(x-x_0)-c\dot x. $$

`SpringPiston::advanceExplicit` 先由上一时刻壁面压力推进活塞，再把旧位置到新位置之间的仿射网格
运动传给 Euler-ALE 求解器。这个耦合是显式分区格式，适合作为第一版 FSI infra 的可读、可验证基线；
强激波 / 轻质量结构下后续可以替换为子迭代或 Robin-Neumann 耦合。

## 3. 网格与边界接口

`EulerALE1D::step` 接收两个函数对象：

```cpp
using GridFn = std::function<Grid1D(double)>;
using BoundaryFn = std::function<BoundaryCondition1D(double)>;
```

- `GridFn(t)` 给出当前时刻所有界面位置 `xFace` 和界面速度 `wFace`。
- `BoundaryFn(t)` 给出左右边界类型，可选周期、透射、固定壁、移动壁、Dirichlet。

已有网格生成器：

- `Grid1D::affinePiston(n, length, velocity)`：左端固定、右端随活塞运动的仿射网格。
- `Grid1D::oscillatory(n, length, t, amplitude, omega)`：用于 GCL / 移动网格收敛验证的周期振荡网格。
- `Grid2D::affinePiston(nx, ny, lengthX, lengthY, velocityX)`：二维通道，右端随活塞运动。
- `Grid2D::oscillatoryX(nx, ny, lengthX, lengthY, t, amplitude, omega)`：二维周期盒，x 方向振荡移动网格。

空间重构可选：

- `SlopeMode::Unlimited`：光滑收敛测试用，得到二阶空间精度。
- `SlopeMode::Minmod`：含激波算例默认，抑制 Gibbs 振荡。
- `SlopeMode::FirstOrder`：调试 / 极端鲁棒模式。

时间推进使用二阶 SSP-RK 型 predictor-corrector。第一阶段在 $t^n$ 网格上算残差，第二阶段在
$t^{n+1}$ 网格上算残差，并更新到新网格。

非结构二维 `EulerDG` 保留原有 ghost-state API：

```cpp
using ExteriorStateFn =
    std::function<Vector4d(double, double, double, const Vector4d&, double, double, int)>;
```

同时新增 direct numerical boundary flux API：

```cpp
struct BoundaryFluxFn;
bool EulerDG::stepWithBoundaryFlux(double dt, double tEnd,
                                   const BoundaryFluxFn& boundaryFlux);
void EulerDG::inviscidResidualWithBoundaryFlux(const MatrixXd& U, double t,
                                               const BoundaryFluxFn& boundaryFlux,
                                               MatrixXd& R) const;
```

移动壁面可直接使用 `movingWallBoundaryFlux(...)`，避免把 ALE 壁面通量硬塞成静止网格 ghost state。
低层已有 `aleNormalFlux` / `aleRusanov` / `movingWallFlux`，下一步只需把 stage mesh、扫掠法向速度和
随时间变化的单元面积接入 `EulerDG` 的体积分、面残差和质量矩阵。

## 4. 验证程序

构建并运行：

```bash
cmake --build build -j --target euler_fsi_convergence
./build/euler_fsi_convergence
cmake --build build -j --target euler_fsi_convergence2d
./build/euler_fsi_convergence2d
ctest --test-dir build -L fsi --output-on-failure
```

验证包含三部分：

1. **GCL / 均匀流保持**：在周期振荡移动网格上推进均匀流，使用端点位置给出的离散扫掠速度，
   primitive drift 应保持在舍入误差量级。
2. **移动网格空间收敛阶**：周期熵波 / 接触波精确解，移动网格上测 `rho` 的 $L^2$ 误差，预期二阶。
3. **固定网格时间收敛阶**：固定细网格上逐次减半时间步，用 successive differences 验证二阶。

`euler_fsi_flux_test` 验证低层 ALE 通量、移动壁面通量、非结构 DG direct boundary-flux 残差/步进路径
和 FSI 能量账本公式。
`euler_fsi_convergence` 和 `euler_fsi_convergence2d` 带有自动门槛：GCL、空间误差、空间阶和时间阶
任一失败都会返回非零码。`euler_fsi_piston` 与 `euler_fsi_piston2d` 也会检查正性、质量漂移、
气体-活塞壁面做功、耦合能量账本、活塞位移、帧数与最终静帧。五者已注册为 CTest 的 `fsi` 标签测试。

一次实测输出：

```text
[1] Geometric conservation law / uniform-flow preservation
  n=  40  max primitive drift=7.772e-16
  n=  80  max primitive drift=1.110e-15
  n= 160  max primitive drift=1.554e-15
  n= 320  max primitive drift=1.998e-15

[2] Spatial convergence on a moving mesh
  n=  50  h=0.02000  L2(rho)=1.368e-04  rate=0.00
  n= 100  h=0.01000  L2(rho)=3.006e-05  rate=2.19
  n= 200  h=0.00500  L2(rho)=7.502e-06  rate=2.00
  n= 400  h=0.00250  L2(rho)=1.909e-06  rate=1.97

[3] Temporal convergence on a fixed mesh (successive differences)
  dt=1.000e-03 -> 5.000e-04  ||rho(dt)-rho(dt/2)||_L2=7.194e-08  rate=0.00
  dt=5.000e-04 -> 2.500e-04  ||rho(dt)-rho(dt/2)||_L2=1.798e-08  rate=2.00
  dt=2.500e-04 -> 1.250e-04  ||rho(dt)-rho(dt/2)||_L2=4.496e-09  rate=2.00
  dt=1.250e-04 -> 6.250e-05  ||rho(dt)-rho(dt/2)||_L2=1.124e-09  rate=2.00

Verification PASS
```

二维验证的一次实测输出：

```text
[1] 2D GCL / uniform-flow preservation
  nx=  24 ny=  18  max primitive drift=1.110e-15
  nx=  48 ny=  36  max primitive drift=1.443e-15
  nx=  96 ny=  72  max primitive drift=2.109e-15

[2] 2D spatial convergence on a moving mesh
  nx=  24 ny=  18  L2(rho)=5.945e-04  rate=0.00
  nx=  48 ny=  36  L2(rho)=8.362e-05  rate=2.83
  nx=  96 ny=  72  L2(rho)=1.377e-05  rate=2.60
  nx= 192 ny= 144  L2(rho)=2.851e-06  rate=2.27

[3] 2D temporal convergence on a fixed mesh
  dt=1.000e-03 -> 5.000e-04  rate=0.00
  dt=5.000e-04 -> 2.500e-04  rate=2.00
  dt=2.500e-04 -> 1.250e-04  rate=2.00

Verification PASS
```

## 5. Shock Tube + Spring Piston 算例

运行：

```bash
cmake --build build -j --target euler_fsi_piston
./build/euler_fsi_piston
ffmpeg -y -framerate 30 -i out/euler_fsi_piston_frames/frame_%05d.ppm \
       -c:v libx264 -pix_fmt yuv420p -crf 18 out/euler_fsi_piston.mp4

cmake --build build -j --target euler_fsi_piston2d
./build/euler_fsi_piston2d
ffmpeg -y -framerate 30 -i out/euler_fsi_piston2d_frames/frame_%05d.ppm \
       -c:v libx264 -pix_fmt yuv420p -crf 18 out/euler_fsi_piston2d.mp4

cmake --build build -j --target euler_fsi_piston2d_showcase
./build/euler_fsi_piston2d_showcase
ffmpeg -y -framerate 30 -i out/euler_fsi_piston2d_showcase_frames/frame_%05d.ppm \
       -c:v libx264 -pix_fmt yuv420p -crf 18 out/euler_fsi_piston2d_showcase.mp4

cmake --build build -j --target euler_fsi_vortex_showcase
./build/euler_fsi_vortex_showcase
ffmpeg -y -framerate 30 -i out/euler_fsi_vortex_showcase_frames/frame_%05d.ppm \
       -c:v libx264 -pix_fmt yuv420p -crf 18 out/euler_fsi_vortex_showcase.mp4
```

默认设置：

- 区域左端固定壁，右端为弹簧活塞。
- 初值为高压腔 shock-tube：左侧 `p=8`，右侧 `p=1`。
- 活塞参数：`mass=7.5`、`stiffness=180`、`damping=6`、`externalPressure=1`。
- `n=420`，`t_end=0.62`，输出 180 个 PPM 帧。
- 二维通道默认 `nx=240, ny=48`，通道高度 `Ly=0.25`，`t_end=0.50`，输出 150 个 PPM 帧。
- 视觉展示版默认 `nx=360, ny=72`，`Ly=0.32`，左侧高压 `p=28`，轻质活塞
  `mass=0.75`、`stiffness=10`、`damping=0.85`，`t_end=0.95`，输出 240 个 PPM 帧。
- 涡量展示版默认 `nx=420, ny=126`，`Ly=0.36`，左侧高压 `p=28`，中部为正弦扰动重/轻气体界面，
  右侧为弹簧活塞；`t_end=0.76`，输出 210 个 PPM 帧。

输出文件：

- `out/euler_fsi_piston_frames/frame_*.ppm`：逐帧密度图。
- `out/euler_fsi_piston_final.ppm`：最终静帧。
- `out/euler_fsi_piston_diagnostics.csv`：时间、活塞位置/速度、壁面压力、最小密度/压力、守恒积分、
  活塞动能、弹簧势能、外压势能、气体壁面做功、阻尼耗散、能量账本漂移。
- `out/euler_fsi_piston.mp4`：ffmpeg 合成后的视频。
- `out/euler_fsi_piston2d_frames/frame_*.ppm`、`out/euler_fsi_piston2d_final.ppm`、
  `out/euler_fsi_piston2d_diagnostics.csv`、`out/euler_fsi_piston2d.mp4`：二维通道对应产物；
  画面使用固定视窗，右侧会显示移动活塞板、外侧结构区和弹簧参考。
- `out/euler_fsi_piston2d_showcase_frames/frame_*.ppm`、`out/euler_fsi_piston2d_showcase_final.ppm`、
  `out/euler_fsi_piston2d_showcase_diagnostics.csv`、`out/euler_fsi_piston2d_showcase.mp4`：
  展示版产物；画面使用高对比密度/压力混合 colormap，活塞行程约 `0.25`，可见轨迹和弹簧压缩。
- `out/euler_fsi_vortex_showcase_frames/frame_*.ppm`、`out/euler_fsi_vortex_showcase_final.ppm`、
  `out/euler_fsi_vortex_showcase_diagnostics.csv`、`out/euler_fsi_vortex_showcase.mp4`：
  涡量展示版产物；画面叠加密度/压力与符号涡量，诊断 CSV 记录 `max_abs_vorticity`、活塞响应和能量账本。

诊断 CSV 中应检查：

- `rho_min > 0` 且 `p_min > 0`，确认没有非物理状态。
- `piston_x` 和 `piston_v` 随冲击波响应，确认结构确实参与耦合。
- `mass` 在封闭移动域内保持稳定，确认 ALE 积分形式没有漏质量。
- `gas_work_drift` 和 `coupled_energy_drift` 保持小量，确认移动壁面通量、活塞做功、外压与阻尼耗散在同一账本上。

默认实测末态为 `piston_x≈1.03698`，`rho_min≈0.3151`，`p_min≈1.5890`，质量漂移约 `5.55e-16`；
全程最大 `gas_work_drift≈1.97e-04`，最大 `coupled_energy_drift≈1.44e-04`。
二维通道末态为 `piston_x≈1.00548`，`rho_min≈0.3339`，`p_min≈1.7235`，质量漂移约 `1.69e-15`；
全程最大 `gas_work_drift≈1.95e-05`，最大 `coupled_energy_drift≈1.75e-05`。

## 6. 下一步扩展

这一版已经把一维和结构化二维 FSI infra 的最小闭环做完整，并且给非结构二维 `EulerDG` 接上了
可验证的 ALE / moving-wall 通量入口。剩余最大的工程块是把非结构 DG 的几何也做成随 stage 变化。
建议后续按下面顺序扩展：

1. 在现有二维 `EulerDG` 中加入 stage mesh、ALE 体积分与界面速度项，先做非结构 DG 的 moving-wall MMS / GCL。
2. 把结构化 `Grid2D` 的边界速度和压力载荷接口推广成通用边界网格 / 边界段接口。
3. 加刚体 6-DOF 或 2D 截面刚体模型，做 shock-loaded rigid cylinder / piston channel。
4. 加结构有限元或 Cosserat beam，做 shock-loaded flexible panel / panel flutter。
5. 引入强耦合子迭代、Robin-Neumann 或 Aitken 松弛，应对轻结构 added-mass 问题。

### 6.1 现有二维 `EulerDG` 的接入点审计

当前二维可压缩 Euler 求解器位于 `src/euler/EulerDG.{h,cpp}`，边界条件通过 `ExteriorStateFn`
鬼状态回调进入。为 moving-wall / FSI，已经在 `EulerDG.h` 中抽出：

- `slipWallExterior(Uin,nx,ny)`：固定滑移壁，供现有 DMR、RMI、shock-bubble、corner diffraction 共用。
- `movingSlipWallExterior(Uin,nx,ny,wallU,wallV)`：移动滑移壁，镜像相对壁面速度并按原压力重算能量。
- `aleNormalFlux(U,nx,ny,w_n)` 与 `aleRusanov(UL,UR,nx,ny,w_n)`：任意法向 ALE 通量与 Rusanov 谱半径。
- `movingWallFlux(Uin,nx,ny,wallU,wallV)`：exact ALE moving-wall flux。
- `BoundaryFluxFn`、`movingWallBoundaryFlux(...)`、`EulerDG::inviscidResidualWithBoundaryFlux`、
  `EulerDG::stepWithBoundaryFlux`：绕过 ghost state 的边界数值通量入口。

`euler_fsi_flux_test` 已验证：等状态 ALE 通量退化为 `F_n-w_nU`，移动壁面通量满足
$(0,pn_x,pn_y,p\mathbf V_w\cdot\mathbf n)$，direct boundary-flux residual / step 对非结构 DG
常量流保持到舍入误差。

真正的二维 ALE 还需要改下面四处：

- **残差装配**：`EulerDG::inviscidResidual` 当前用静态 `mesh_.node` 与静态边长/法向。ALE 版需要在
  每个 RK stage 传入 stage mesh 与 swept face normal velocity，并保证内边一份扫掠通量两侧反号使用。
- **质量矩阵与几何**：当前 `M_e=area_e*Mref`，`applyMassInverse` 与 IMEX stage 都假设面积不变。
  ALE 版需要把 `area_e(t)` 作为 stage 几何的一部分，或把未知量改成移动单元守恒积分，像 1D 版一样
  先满足 GCL。
- **ALE-HLLC**：低层已经有 ALE-Rusanov；HLLC-ALE 可以在非结构 Rusanov-GCL 验证通过后再接。
- **结构耦合**：边界压力载荷应从同一条 moving-wall 通量积分得到，传给 rigid body / panel / beam；
  结构返回边界位移和速度，再由 mesh-motion 模块生成下一步体网格。

二维第一组验证建议沿用一维顺序：固定自由流的 exact-GCL 到舍入误差、二维移动网格等熵涡空间阶、
固定网格时间阶，最后再上 shock-loaded rigid cylinder 或 flexible panel。
