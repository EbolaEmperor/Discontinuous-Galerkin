# DG-cpp

二维(间断)Galerkin 有限元求解器的 C++ 实现,

---

## 1. 基础算例

均在区域 $\Omega \subset \mathbb{R}^2$ 上求解,$f$ 由制造解反推。

| 可执行文件 | 方程 | 边界条件 | 方法 | 默认区域 |
|---|---|---|---|---|
| `poisson_pk` | $-\Delta u = f$ | $u = g \ \text{on}\ \partial\Omega$(强/本质 Dirichlet) | 协调 Lagrange $P_k$ | 正六边形 |
| `poisson_ipdg` | $-\Delta u = f$ | $u = g \ \text{on}\ \partial\Omega$(强 Dirichlet;内罚仅作用于内部边) | 内罚间断 Galerkin(IPDG) | 正六边形 |
| `biharmonic_ipcg` | $\Delta^2 u = f \quad (k \ge 2)$ | 两种边界条件,由 `bc_type` 切换:**clamped**(固支)$u=g_1,\ \partial_{\mathbf{n}}u=g_2$($u$ 强 Dirichlet,$\partial_{\mathbf{n}}u$ 用 Nitsche 弱施加);**simply supported**(简支)$u=g,\ \partial_{\mathbf{nn}}u=h$($u$ 强 Dirichlet,$\partial_{\mathbf{nn}}u$ 作为自然边界条件进右端) | C⁰ 内罚 Galerkin(C⁰-IPCG) | 单位正方形 |
| `biharmonic_argyris` | $\Delta^2 u = f$ | **clamped**(固支)$u=g_1,\ \partial_{\mathbf{n}}u=g_2$,全部强施加(边界顶点除纯二阶法向导 $\partial_{\mathbf{nn}}u$ 外全部钉死,角点全钉) | **Argyris C¹ 协调元**(分片 $P_5$,21 自由度/三角形,真 $H^2$ 协调,纯能量 $\int D^2u:D^2v$,无内罚) | 单位正方形 |
| `poisson_mixed_hdg` | $\boldsymbol{\sigma} = \nabla u,\quad -\nabla\!\cdot\!\boldsymbol{\sigma} = f$ | $u = g \ \text{on}\ \partial\Omega$(弱 Dirichlet) | HDG 混合元 $(\mathrm{vec}DP_k,\, DP_k)$,**整体求解鞍点系统**(SparseLU) | 单位正方形 |
| `poisson_hdg_hybrid` | $\boldsymbol{\sigma} = \nabla u,\quad -\nabla\!\cdot\!\boldsymbol{\sigma} = f$ | $u = g \ \text{on}\ \partial\Omega$(弱 Dirichlet,边界 trace 取 $g$ 的 $L^2$ 投影) | **HDG 杂交化**:引入 trace $\lambda=\hat u\in P_k(e)$,逐单元静态凝聚 $(\boldsymbol\sigma,u)$,**仅对 $\lambda$ 求解对称正定 Schur 补系统**(CHOLMOD) | 单位正方形 |

- **边界数据全部取自制造解**:$g = u_{\text{exact}}\big|_{\partial\Omega}$,双调和 clamped 的法向导数 $g_2 = \nabla u_{\text{exact}}\cdot\mathbf{n}$,
  simply supported 的法向二阶导 $h = \mathbf{n}^{\!\top}(\nabla^2 u_{\text{exact}})\,\mathbf{n}$。
  因此四个程序本质上都是**非齐次 Dirichlet** 类问题的收敛性验证;目前 C++ 版未实现 Neumann 边界
  (MATLAB 版的 1D 混合求解器才支持逐端选 D/N)。
- **simply supported 的弱形式**:仍用 $\sum_K\int_K \nabla^2u:\nabla^2v$ 能量,分部积分在边界自然产生
  $+\oint_{\partial\Omega}(\partial_{\mathbf{nn}}u)(\partial_{\mathbf n}v)$。简支下 $\partial_{\mathbf n}v$ 自由、$\partial_{\mathbf{nn}}u=h$ 已知,
  该项直接进载荷向量;边界**不加**任何罚/Nitsche LHS(内部边的 C⁰ 内罚照旧)。k 阶元 H¹ 收敛阶为 k。
- **IPDG / C⁰-IPCG** 支持三种内罚变体:**SIPG**($\beta=1$)、**NIPG**($\beta=-1$)、**IIPG**($\beta=0$)。
- **mixed HDG** 在解出 $(\boldsymbol{\sigma}_h, u_h)$ 后,逐单元求解局部 Poisson 问题把 $u_h$ 投影到
  $DP_{k+1}$ 空间得到 $u_h^{\ast}$,用于展示超收敛(误差阶比 $u_h$ 高一阶)。
- **hybrid HDG**(`poisson_hdg_hybrid`)是真正的**杂交化**实现,与 `poisson_mixed_hdg`(直接组装并求解
  整体鞍点系统)的区别在于:它把通量变量取为 $\boldsymbol q=-\nabla u$(满足 $\nabla\!\cdot\!\boldsymbol q=f$),
  得到**对称**的单元局部解算子
  $\mathcal A=\begin{bmatrix}A_{qq}&-D\\ D^{\top}&H_{uu}\end{bmatrix}$;再用数值通量
  $\widehat{\boldsymbol q}\!\cdot\!\boldsymbol n=\boldsymbol q\!\cdot\!\boldsymbol n+\tau(u-\lambda)$ 与
  通量守恒/传输条件,逐单元**静态凝聚**掉 $(\boldsymbol q,u)$,得到只含骨架 trace $\lambda$ 的全局
  **对称正定** Schur 补系统
  $\mathbb A\,\lambda=\boldsymbol b,\ \ \mathbb A=\sum_K\big(M_{\lambda\lambda}-S\,\mathcal A^{-1}R\big)$,
  用 **CHOLMOD**(`CholmodSupernodalLLT`)求解;随后逐单元回代恢复 $(\boldsymbol q_h,u_h)$,取
  $\boldsymbol\sigma_h=\nabla u=-\boldsymbol q_h$。边界面上 $\lambda$ 直接取 $g$ 在该面 $P_k$ 上的 $L^2$ 投影
  并消去。全局未知量从 $(\mathrm{vec}DP_k,DP_k)$ 的体积自由度降到只剩骨架自由度,数量大幅下降。
- **收敛阶**(制造解,$\tau=O(1)$):$\boldsymbol\sigma_h$ 与 $u_h$ 的 $L^2$ 阶均为 $k+1$;局部后处理
  $u_h^\ast\in DP_{k+1}$ 的 $L^2$ 阶为 $k+2$(超收敛)。已在 $k=1,2,3,4$ 上数值验证。可执行文件接受可选参数
  `poisson_hdg_hybrid [k] [tau] [Nref] [solver]`(默认 `2 1.0 5 3`)。
- 制造解为 $u(x,y) = \sin\big(\pi(x-c)\big)\,\sin\big(\pi(y-c)\big)$,中心 $c$ 由 `ExactSolution sol(0.3)`
  给定(对应 MATLAB 里的 `sinsin(0.3)`)。

---

## 2. 依赖

| 依赖 | 必需? | 说明 |
|---|---|---|
| CMake ≥ 3.10 | ✅ | 构建系统 |
| C++17 编译器 | ✅ | Apple Clang / Clang / GCC 均可 |
| [Eigen3](https://eigen.tuxfamily.org) ≥ 3.3 | ✅ | 头文件库,提供稠密/稀疏矩阵与线性求解器 |
| OpenMP | ⛳ 可选 | 用于 DG 装配的并行;找不到时自动退化为单线程 |
| SuiteSparse / CHOLMOD | ⛳ 可选 | 提供更快的稀疏 Cholesky 直接解(`CholmodSupernodalLLT`);找不到时回退到 Eigen 自带的 `SimplicialLLT` |

### macOS(Homebrew)

```bash
brew install cmake eigen
brew install suite-sparse   # 可选:启用 CHOLMOD
brew install libomp         # 可选:启用 OpenMP 并行装配
```

CMake 会在 `/opt/homebrew/opt/suite-sparse` 查找 CHOLMOD(见 `CMakeLists.txt`),
Homebrew 默认路径即可被自动发现。

### Linux(Debian/Ubuntu)

```bash
sudo apt install cmake g++ libeigen3-dev
sudo apt install libsuitesparse-dev   # 可选:CHOLMOD
# OpenMP 通常随 g++ 提供,无需额外安装
```

### Windows(w64devkit / MinGW-w64,本仓库已配置)

本机使用 [w64devkit](https://github.com/skeeto/w64devkit)(GCC 14,msvcrt)作为编译器。Eigen 与
SuiteSparse/CHOLMOD 通过**预编译依赖**放在 `external/`(已 git-ignore,不入库):

- `external/eigen3/` — Eigen 3.4.0 头文件(纯头文件库)。
- `external/ssprefix/mingw64/` — 取自 MSYS2 的 **MinGW 预编译** SuiteSparse + OpenBLAS 及其运行期
  DLL(与 w64devkit 同为 msvcrt ABI,CHOLMOD 是 C 接口,链接无碍)。`CMakeLists.txt` 会优先在此
  前缀中查找 `cholmod.h` 与各 `*.dll.a` 导入库;构建后用 `POST_BUILD` 把所需运行期 DLL
  (`libcholmod / libopenblas / libgfortran-5 / libgomp-1 / …`)拷到可执行文件旁边,直接运行即可。

配置/构建(从 PowerShell 或任意 shell,确保 w64devkit 的 `bin` 在 PATH 中):

```powershell
$env:PATH = "D:\Program Files\MinGW\w64devkit\bin;" + $env:PATH
cmake -S . -B build -G "MinGW Makefiles" `
      -DCMAKE_TOOLCHAIN_FILE=w64devkit-toolchain.cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

> `w64devkit-toolchain.cmake` 指定 w64devkit 的 `gcc/g++/mingw32-make`,并**跳过 CMake 的编译器探测**。
> 这是因为本机 Windows Defender 会在编译器刚生成 `a.exe` 时短暂占用文件,导致 CMake 读取该探测
> 可执行文件失败(`file failed to open for reading (Invalid argument)`);预置编译器标识即可绕过。
> 此外该工程也在 `external/eigen3` 存在时优先使用自带 Eigen,避免被用户级 CMake 包注册表里指向其他
> 工程构建目录的陈旧 `Eigen3_DIR` 劫持。

### macOS + Apple Clang 的 OpenMP 注意

> Apple 自带的 `clang` 需要 `libomp` 才能用
> OpenMP。若 `find_package(OpenMP)` 没找到,装配会以单线程运行(代码已优雅处理),
> 不影响结果正确性。如需开启,可安装 `libomp` 后给 CMake 传入对应的
> `-DOpenMP_CXX_FLAGS`/`-DOpenMP_CXX_LIB_NAMES` 提示,或改用 `brew install llvm` 的编译器。

---

## 3. 构建

在 `Poisson-cpp/` 目录下:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

如果 Eigen / SuiteSparse 装在非标准前缀,补一个前缀提示即可:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/homebrew
```

> **Eigen 版本**:`CMakeLists.txt` 会优先查找 3.x 系列的 Eigen,找不到再退而接受任意版本
> (含 Homebrew 现在默认的 **5.x**),并统一通过导入目标 `Eigen3::Eigen` 传递头文件路径。
> 代码在 Eigen 3.4 与 5.0 下都已验证可编译运行,通常无需额外操作。若想固定用 3.x,可装
> `eigen@3` 并把它放在前缀首位:
> ```bash
> brew install eigen@3
> cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
>       -DCMAKE_PREFIX_PATH="$(brew --prefix eigen@3);/opt/homebrew"
> ```
> Windows + w64devkit 下另有一套自带依赖(`external/` + 工具链文件),见本节后文。

构建产物(八个可执行文件 + 四个静态库)位于 `build/`:

```
build/poisson_pk
build/poisson_ipdg
build/biharmonic_ipcg
build/biharmonic_argyris
build/poisson_mixed_hdg
build/poisson_hdg_hybrid
build/cahn_hilliard
build/cahn_hilliard_convergence
```

> Release 模式默认带 `-O3 -march=native`(见 `CMakeLists.txt`)。`-march=native` 会针对
> 本机 CPU 优化,**不要把编译产物拷到不同架构的机器上运行**。

---

## 4. 运行

可执行文件**不接收命令行参数**,直接运行即可:

```bash
./build/poisson_pk
./build/poisson_ipdg
./build/biharmonic_ipcg
./build/poisson_mixed_hdg
./build/poisson_hdg_hybrid           # 可选参数: [k] [tau] [Nref] [solver]
./build/biharmonic_argyris           # 可选参数: [Nref] [solver]
# biharmonic_ipcg 也支持可选参数: [ord] [bc_type] [Nref] [ip_type]
```

### Argyris 元 vs. C⁰-IPCG 对比文档

`biharmonic_argyris`(Argyris C¹ 协调元)与 `biharmonic_ipcg`(C⁰ 内罚)
求解同一双调和问题的优劣对比(自由度规模、条件数、收敛阶、实现复杂度、
适用场景,含实测数据)见 [`docs/Argyris_vs_C0IPCG.md`](docs/Argyris_vs_C0IPCG.md)。

### 输出示例(`poisson_pk`)

```
Poisson conforming Pk FEM (C++)
ord = 4
Solver: CholmodSupernodalLLT
...
Domain: Polygon (regular hexagon)
Lv 0 (h=0.5000):
  Init: 0.001s
  Assemble: 0.004s
  Solve: 0.002s
  Error: 0.001s
  nDof=... L2=... H1=...
...
Rates (Polygon (regular hexagon)):
Lv 0->1: L2 rate=5.00 H1 rate=4.00
...
```

- **L2 / H1**:制造解与数值解的 L² 与 H¹ 误差。
- **rate**:相邻两层误差比的以 2 为底对数,即网格加密一倍时的收敛阶。
  对 k 阶元,Poisson/IPDG 理论上 L² 阶为 k+1、H¹ 阶为 k。
- `poisson_mixed_hdg` 额外输出 `σ`、`u`、`u*` 三者的 L² 误差与阶,其中 `u*` 应比 `u` 高约一阶(超收敛)。

---

## 5. 配置参数

所有参数都**硬编码在各 `*_main.cpp` 顶部**,改完需要**重新编译**(`cmake --build build -j`)。
常用开关如下:

| 参数 | 含义 | 出现位置 |
|---|---|---|
| `ord` | 多项式次数 k(双调和要求 `ord ≥ 2`) | 全部 main |
| `Nref` | 网格加密层数 | 全部 main |
| `h0` | 初始网格尺寸 | 全部 main |
| `sigma` / `alpha` | 内罚 / 杂交罚系数 | IPDG、biharmonic、mixed |
| `beta` 或 `ip_type` | 内罚变体:`SIPG`/`NIPG`/`IIPG` | IPDG、biharmonic |
| `bc_type` | 双调和边界条件:`CLAMPED`(固支)/`SIMPLY_SUPPORTED`(简支,默认) | biharmonic |
| `solver_type` | 线性求解器选择(见下) | 全部 main |
| `ExactSolution sol(0.3)` | 制造解中心 c | 全部 main |
| `cases` 向量 / `regularPolygon(n)` | 计算区域 | 全部 main |

### 线性求解器(`solver_type`)

| 值 | 求解器 | 适用 |
|---|---|---|
| 0 | `SimplicialLDLT` | 对称矩阵 |
| 1 | `SimplicialLLT`(Eigen 自带 Cholesky) | 对称正定 |
| 2 | `ConjugateGradient` | 对称正定,迭代 |
| 3 | `CholmodSupernodalLLT` | 对称正定,需 CHOLMOD;未找到时自动回退到 `SimplicialLLT` |
| 4 | `SparseLU` | **仅 `poisson_mixed_hdg`**,鞍点系统不定,默认用它 |

`poisson_pk` / `poisson_ipdg` / `biharmonic_ipcg` 默认 `solver_type=3`(CHOLMOD);
`poisson_mixed_hdg` 默认 `solver_type=4`(LU,因为混合系统是对称不定的鞍点问题,
不能用 Cholesky)。

### 切换计算区域

main 里用一个 `cases` 向量描述区域,例如 `poisson_pk` / `poisson_ipdg`:

```cpp
vector<CaseCfg> cases;
cases.push_back({"Polygon (regular hexagon)", true, regularPolygon(6), 0.5});
// 改成单位正方形:
// cases.push_back({"Unit square", false, MatrixXd(), 0.5});
```

- `isPolygon == false` → `Mesh::getMesh(h)` 生成 `[0,1]²` 单位正方形三角网格。
- `isPolygon == true`  → `Mesh::getPolygonMesh(vertices, h)`,`regularPolygon(n)` 生成正 n 边形;
  也可以传任意 `n×2` 顶点矩阵。

`biharmonic_ipcg` 默认跑单位正方形,正六边形那行被注释掉,按需切换即可。

---

## 6. Cahn-Hilliard 演化算例

`cahn_hilliard` 与上面四个收敛性实验不同,它是一个**时间演化**算例:在单位正方形 $[0,1]^2$ 上
求解 Cahn-Hilliard 相场方程(无通量边界),空间用已有的 SIPG 间断 Galerkin 离散 $-\Delta$、
时间用稳定化半隐(IMEX / 凸分裂)格式——`time_order=1`(稳定化后向欧拉,无条件能量稳定)或
`time_order=2`(SBDF2,二阶),逐帧把相场 $c$ 写成 PPM 图片并用 ffmpeg 合成 MP4,直观展示
**spinodal decomposition(相分离)与粗化(coarsening)**。求解保持**质量精确守恒**与**自由能耗散**
(运行时可见 `mass_drift` $\sim10^{-14}$、一阶下能量单调下降)。

默认开启**自适应时间步**(`adaptive=true`):早期相分离快时 $\tau$ 自动变小、后期粗化慢时自动变大,
于是早期出帧密、播放慢,后期快进——避免"前几帧一闪而过"。$\tau$ snap 到离散档、每档矩阵只分解一次。

**计算区域**可选(`domain`):`"square"`(单位正方形)、`"hexagon"`(正六边形)、`"disk"`(圆盘)。
无通量自然边界对任意形状都成立,质量守恒与能量耗散照旧。圆盘网格由 `Mesh::getDiskMesh` 按 $h$ 自动生成
(六边形等边网格光滑映射到圆,高质量、边界精确在圆上)。

> 📄 **方程、边界条件、完整时空离散方案(含二阶 SBDF2)、全部可调参数及含义,见专门文档
> [`docs/cahn-hilliard.md`](docs/cahn-hilliard.md)。**

参数全部写在 **JSON 配置文件**，参考 `examples/ch_config.json`。

```bash
# 演示:自动读取 ch_config.json,帧写入 ./ch_frames/,再合成 MP4
./build/cahn_hilliard examples/ch_config.json
ffmpeg -y -framerate 25 -i ch_frames/frame_%05d.ppm \
       -c:v libx264 -pix_fmt yuv420p -crf 18 cahn_hilliard.mp4

# 时空收敛阶测试(MMS,不产生视频):验证空间 L2 阶 k+1 / H1 阶 k、时间 1 与 2 阶
./build/cahn_hilliard_convergence
```

视频与帧图为生成产物,已在 `.gitignore` 中忽略;ffmpeg 仅用于合成视频,求解/出帧本身只依赖 Eigen。

---

## 7. 不可压 Navier–Stokes 圆柱绕流算例

`navier_stokes_cylinder` 是一个**时间演化**算例:在矩形挖去一个圆盘(圆柱截面)的区域上求解二维
**不可压 Navier–Stokes 方程**,用经典的**圆柱绕流**展示 **von Kármán 涡街**——逐帧把**涡量**
栅格化为 PPM 并用 ffmpeg 合成 MP4。

- **网格自动生成**:只给一个目标尺度 $h$,程序用从零实现的 **DistMesh**(Persson–Strang)式生成器
  (符号距离函数 + 力平衡 + 自带的 **Bowyer–Watson** Delaunay)生成**高质量、按到圆柱距离渐变加密**的
  非结构三角网格(典型最小内角 $\sim30^\circ$、圆周精确分辨)。
- **边界条件**(用户设定):左**入流** $u=(U_\infty,0)$、右**出流**(do-nothing)、上下**可滑移**
  侧壁($v=0$、切向无应力,即对称面)、圆柱面**无滑移**。
- **空间**:任意阶 $\mathbb{dP}_k$ 间断 Galerkin,黏性/压力用 SIPG $-\Delta$(Dirichlet 用 Nitsche 弱加),
  对流用显式 **Lax–Friedrichs** 通量,散度/梯度用一对自洽的弱导算子。
- **时间**:**二阶 IMEX/PPE 格式**:时间导数 BDF2、对流外插 EX2,压力默认用直接 PPE
  (`pressure_mode="direct_ppe"`)和 KIO/Gresho-Sani 高阶 Neumann 边界;默认用 `ppe_div_damping`
  在压力方程中阻尼散度,保持逐分量 Helmholtz 快速回代。可选 `grad_div>0` 会启用 u-v 耦合稳定化。
- **后处理**:在圆柱面积分应力得阻力/升力系数 $C_D,C_L$,由升力振荡估计 **Strouhal 数**;力的时间序列
  写入 `ns_forces.csv`。$\mathrm{Re}=100$ 默认参数下给出 $C_D\approx1.5$、$\mathrm{St}\approx0.18$。

> 📄 **方程、计算域与边界、网格生成算法、完整时空离散方案(含高阶压力边界条件)、全部可调参数及
> Taylor–Green 收敛阶验证,见专门文档 [`docs/navier-stokes.md`](docs/navier-stokes.md)。**

参数全写在 **JSON 配置文件**里，参考 `examples/ns_config.json`。

```bash
# 圆柱绕流影片:自动读取 ns_config.json,帧写入 ./ns_frames/,再合成 MP4
./build/navier_stokes_cylinder examples/ns_config.json
ffmpeg -y -framerate 25 -i ns_frames/frame_%05d.ppm \
       -c:v libx264 -pix_fmt yuv420p -crf 18 cylinder_vortex.mp4

# Taylor–Green 时空收敛阶测试(无视频):验证速度 L2 空间阶 k+1、时间阶 2
./build/navier_stokes_convergence

# 圆柱+小蝌蚪:刚性尾巴(36c0dba 保留版)
./build/navier_stokes_tadpole examples/ns_tadpole_config.json

# 圆柱+小蝌蚪:弹性尾巴(Cosserat 杆 + immersed-boundary 反作用力)
./build/navier_stokes_tadpole_elastic examples/ns_tadpole_elastic_config.json
```

视频与帧图为生成产物,已在 `.gitignore` 中忽略;ffmpeg 仅用于合成视频。

---

## 8. 可压缩 Euler 双马赫反射算例

二维**可压缩 Euler 方程** $U_t+\nabla\cdot F(U)=0$ 的间断 Galerkin($\mathbb{dP}_k$)求解器,时间用
**二阶 IMEX-RK(ARS(2,2,2))**:无粘通量(默认 **HLLC**,可选 Rusanov)**显式**、激波捕捉的
**人工粘性**(Persson–Peraire,$\{\varepsilon\}$-加权 SIPG)**隐式**——刚性扩散不再限制时间步。

- **空间**:守恒变量 $(\rho,\rho u,\rho v,E)$ 各在破碎 $P_k$ 空间;体积分 $\int F:\nabla\phi$ +
  **HLLC** 界面通量(解析接触波→滑移线锐利、涡街清晰);块对角质量阵(仿射,一次参考逆质量阵即可)。
- **时间**:ARS(2,2,2),$\gamma=1-\tfrac{\sqrt2}{2}$,两个隐式级共用 $M+\gamma\,\mathrm dt\,A$(每
  刷新一次 Cholesky 分解);人工粘性为零时退化为稳定的二阶显式 RK。
- **激波捕捉**:Persson–Peraire 人工粘性(**压力传感器**,只作用激波、保滑移线锐利;per-side $\varepsilon$
  加权 SIPG,严格半正定)+ Zhang–Shu **保正限制器**(密度/压力复合压缩,精确守恒,光滑区不降阶)。
- **边界**:统一用"鬼状态"经数值通量弱施加(精确解 / Dirichlet / 零梯度出流 / 反射壁)。
- **并行/高效**:残差/限制器/质量 apply 按单元 `std::thread` 并行;**分区隐式解**利用人工粘性的局部性
  (仅激波区 ~1% 自由度"活跃",$K$ 在活跃/解耦间严格块对角→其余退化为并行块逆质量、只解极小活跃块,
  与整体分解解一致到 ~1e-16);纹影整数幂查表。18 核上 `ny=120` 全程约 **14–15 分钟**(纯效率优化、
  输出不变;优化前 ~35 分钟)。

**收敛阶**(`euler_convergence`,等熵涡精确解):**时间阶严格为 2**(2.00/2.00/2.00);空间 $L^2$ 阶
**P1≈2.0、P2≈3.0、P3≈4.0**(HLLC 低耗散,各阶干净达到 $k+1$)。

**双马赫反射**(`euler_dmr`):马赫 10 激波、域 $[0,4]\times[0,1]$、$T=0.2$、$\mathbb{dP}_2$、$n_y=120$。
输出**三段影片**(全域密度 / 滑移线**局部放大·密度** / 同窗口**局部放大·数值纹影**)+ 最终静帧;
密度与纹影均**逐像素**用 $\mathbb{dP}_k$ 多项式(纹影用真实梯度 $|\nabla\rho|$),清晰呈现三波点与尾部
**涡街**(Kelvin–Helmholtz 卷起)的丰富细节。

**自适应网格加密**(`euler_dmr_amr`):同一套数值,网格改为 **h-AMR**——用**最新顶点二分(NVB)**保持
**协调**(无悬挂结点,求解器一字不改),由**密度梯度指示子**(同时抓住激波**与**滑移线)只在激波/滑移线
附近加密、光滑区粗化(迟滞 + 缓冲),解的传递加密**精确**、粗化**守恒**(实测漂移 $\sim10^{-14}$),并在
每次重剖后按最细单元**重算 $\mathrm dt$ 控制 CFL**。base $n_y{=}50$ + `max_gen=3`(激波处 $\approx$ 均匀
$n_y{=}150$):约 **2.6–4 万**三角形即达均匀 **20 万**三角形的同量级效果——**单元数约 $1/8$**。额外输出一段
**网格叠加影片**展示自适应过程。

> 📄 方程、完整时空格式(含 ARS(2,2,2) 逐级公式)、人工粘性与保正限制器细节、双马赫反射设置、
> **自适应网格加密(NVB / 守恒传递 / CFL 控制)**与全部参数,见 [`docs/euler.md`](docs/euler.md)(§9)。

```bash
# 等熵涡时空收敛阶(约 2 分钟)
./build/euler_convergence

# 自适应网格双马赫反射 + 网格叠加影片
./build/euler_dmr_amr examples/dmr_amr_config.json
ffmpeg -y -framerate 25 -i dmr_amr_frames/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 dmr_amr.mp4

# 均匀网格双马赫反射影片 + 放大 + 纹影 (耗时较久，评测时建议跳过)
./build/euler_dmr examples/dmr_config.json
ffmpeg -y -framerate 25 -i dmr_frames/frame_%05d.ppm -c:v libx264 -pix_fmt yuv420p -crf 16 dmr.mp4
```

### 8.1 更多可压缩 Euler 算例（共用 AMR 场景驱动器）

在双马赫反射之外，新增 4 个经典可压缩 Euler 算例，复用同一套 **h-AMR 场景驱动器**
`src/euler/euler_amr_scene.h`（NVB 自适应 + 守恒传递 + 密度指示子 + CFL 控制 + HLLC + 人工粘性 +
Zhang–Shu 保正，与双马赫反射一字不改），各算例只给出**初值、鬼状态边界、边界标记与基础网格**：

| 可执行 | 物理 | 看点 |
|---|---|---|
| `euler_riemann` | **二维黎曼问题**（Schulz-Rinne/Lax-Liu，`config` 选 3/6/12） | 中心蘑菇喷流 / 四臂剪切涡片 / 螺旋接触面 |
| `euler_shock_bubble` | **激波-气泡相互作用**（Haas–Sturtevant，单 $\gamma$ 代理） | 重气泡（默认 $\rho_b/\rho{=}3$）汇聚透镜→激波聚焦 + 空气射流穿刺 + 界面 KH 涡卷 |
| `euler_rmi_amr` | **Richtmyer–Meshkov 不稳定** | 激波冲过正弦扰动的轻/重气界面，斜压涡量→蘑菇钉 + KH 卷起（脉冲驱动、无需重力） |
| `euler_corner_diffraction` | **激波绕 90° 凸角衍射**（Zhang–Shu） | 马赫 5.09 在 L 形域绕凸角：弯曲衍射激波 + 角点大涡 + **角点近真空**（保正招牌算例） |

**传感器选择（实测关键）**：corner 衍射是**激波主导**→用**压力传感器**（`av_indicator=1`）保激波/剪切锐利、
角点近真空交给 Zhang–Shu 保正。其余三个（二维黎曼、激波-气泡、RMI）是**接触/剪切主导**——保正只保正、
不抑制（保持正值的）高阶振荡，且求解器**无斜率限制器**（`tvbLimit` 为空操作），欠分辨剪切层会无界增长
（$\rho\to10^5$–$10^6$）。对策用**密度传感器**：二维黎曼/RMI 再配低 Persson 阈值 `av_s0=-3`（及早起粘）；
激波-气泡用**温和**密度粘性（过激的低 `av_s0` 会破坏 $\{\varepsilon\}$-SIPG 强制性→反扩散发散）。

**出流边界（长时间跑的关键）**：
- **二维黎曼**：四象限交点放在**离心 (0.8,0.8)**（Schulz-Rinne），主结构往左下大区域发展 + **弱远场鬼态**
  （恒取象限常状态、经 HLLC 逐特征处理）→ 固定 $[0,1]^2$ 干净跑到 **t=0.8**。
- **激波-气泡**：右出流用**特征型非反射出流**（出射不变量取内部、入射 $R^-$ 取环境远场）——裸零梯度会在激波
  抵壁后反射出向左生长的稀疏波黑楔→NaN；**默认**还把气泡整体**左移**（$c_x{=}0.4$、激波 $0.15$）给慢碎片腾空间，
  固定域 $[0,2.5]$ 干净跑到 **t≈4.8**（居中经典布局 $c_x{=}0.9$、$t{=}1.2$ 见 config 注释）。
- **RMI**：界面/激波左移（$x_0{=}0.7/0.3$）+ 密度传感器，固定 $[0,4]$ 干净到 **t≈3.5**；之后透射激波抵右壁
  （亚声速出流）长稀疏波楔、网格爆涨——RMI 是扫掠结构，左移到顶也无法更久，要更长需**加宽域**。
- 另有 `cfl_rho_floor`（仅给 dt 的波速估计设密度地板，防近真空节点压崩 dt）；L 形域由 `makeMaskedRectMesh`
  整格删孔生成，`base_cell` 自动吸附到 $1/k$ 保证孔角对齐网格线。

```bash
./build/euler_riemann examples/riemann_config.json            # 二维黎曼 config 3（config 6/12 见配置）
./build/euler_shock_bubble examples/shock_bubble_config.json  # 激波-气泡（改 bubble_ratio=0.138 得轻气泡涡环）
./build/euler_rmi_amr examples/rmi_config.json                # Richtmyer–Meshkov（max_gen=6 → 界面 ~1/280）
./build/euler_corner_diffraction examples/corner_config.json  # 激波绕 90° 凸角衍射（L 形域）
# 每个跑完会打印对应的 ffmpeg 合成命令
```

## 9. 目录结构

```
Poisson-cpp/
├── CMakeLists.txt
├── README.md                ← 本文件
└── src/
    ├── common/              # 公共库 poisson_common
    │   ├── Mesh.{h,cpp}        网格生成(正方形 / 多边形)、加密、边-单元关系
    │   ├── FEM.{h,cpp}         标量 Pk 元、连续/间断自由度编号、误差计算
    │   ├── Quadrature.{h,cpp}  三角形数值积分
    │   ├── ExactSolution.h     制造解 u、∇u、Δu、Δ²u
    │   └── utils.h
    ├── ipdg/                # 库 dg_assembly + 两个可执行
    │   ├── DGAssembly.{h,cpp}  刚度阵、内罚项、载荷、Nitsche、双调和装配
    │   ├── ipdg_main.cpp       → poisson_ipdg
    │   └── biharmonic_ipcg_main.cpp → biharmonic_ipcg
    ├── cg/
    │   └── pk_main.cpp         → poisson_pk
    ├── mixed/               # 库 hdg + 两个可执行
    │   ├── VecFEM.{h,cpp}      向量值间断元(σ / q)
    │   ├── AssemblyHDG.{h,cpp} 质量阵、混合算子、杂交内罚、弱边界、L2 误差
    │   ├── PostProcess.{h,cpp} 局部 Poisson 求解(超收敛后处理 u*)
    │   ├── HybridHDG.{h,cpp}   杂交化:1D trace 基、单元静态凝聚、SPD Schur 补、CHOLMOD 求解、回代
    │   ├── mixed_main.cpp      → poisson_mixed_hdg(直接解鞍点系统)
    │   └── hybrid_main.cpp     → poisson_hdg_hybrid(杂交化 + Schur 补)
    ├── argyris/             # 库 argyris + 一个可执行
    │   ├── Argyris.{h,cpp}     Argyris C¹ 协调元(P5,21 自由度):逐单元节点基、
    │   │                       全局自由度编号、双调和能量装配、强边界、L2/H1/H2 误差
    │   └── argyris_main.cpp    → biharmonic_argyris
    ├── cahn_hilliard/       # 复用 dg_assembly,两个可执行(演化+视频 / 收敛测试)
    │   ├── CahnHilliard.{h,cpp}   质量阵、非线性载荷、能量/质量、PPM 栅格化、CHIntegrator(一阶/二阶)
    │   ├── Json.h                 无依赖 JSON 配置读取器(支持注释)
    │   ├── ExactSolutionCH.h      收敛测试的制造解 + 源项
    │   ├── ch_main.cpp            → cahn_hilliard(分块 IMEX 时间推进 + 逐帧出图)
    │   └── ch_convergence_main.cpp → cahn_hilliard_convergence(MMS 时空收敛阶测试,无视频)
    └── navier_stokes/       # 库 navier_stokes + 两个可执行(圆柱绕流影片 / 收敛测试)
        ├── MeshGen.{h,cpp}        DistMesh 式网格生成 + Bowyer–Watson Delaunay + 边界分类
        ├── NavierStokes.{h,cpp}   DG 质量阵、Nitsche 边界、弱导算子、LF 对流、BDF2/EX2 分裂积分器、受力、栅格化
        ├── ns_main.cpp            → navier_stokes_cylinder(圆柱绕流 von Kármán 涡街影片 + 力/Strouhal)
        └── ns_convergence_main.cpp → navier_stokes_convergence(Taylor–Green 时空收敛阶测试,无视频)
```

CMake 目标关系:`poisson_common`(基础)← `dg_assembly`(IPDG 装配)← `hdg`(混合元);
`argyris`(C¹ 元)直接依赖 `poisson_common`;`cahn_hilliard*` 与 `navier_stokes*` 可执行复用
`dg_assembly`(借其 SIPG 刚度阵与内罚装配)。十个可执行各自链接所需的库。
