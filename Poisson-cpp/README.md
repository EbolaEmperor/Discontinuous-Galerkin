# Poisson-cpp

二维(间断)Galerkin 有限元求解器的 **C++ / Eigen** 实现,是上层 MATLAB 仓库
[`Discontinuous-Galerkin`](../README.md) 中 2D 求解器的移植版。

每个可执行文件都是一个**收敛性实验**:在一系列逐次加密的三角网格上,用制造解
(manufactured solution)装配并求解,输出每一层的自由度数、L²/H¹ 误差,并拟合收敛阶。

---

## 1. 包含哪些求解器

均在区域 $\Omega \subset \mathbb{R}^2$ 上求解,$f$ 由制造解反推。

| 可执行文件 | 方程 | 边界条件 | 方法 | 默认区域 |
|---|---|---|---|---|
| `poisson_pk` | $-\Delta u = f$ | $u = g \ \text{on}\ \partial\Omega$(强/本质 Dirichlet) | 协调 Lagrange $P_k$ | 正六边形 |
| `poisson_ipdg` | $-\Delta u = f$ | $u = g \ \text{on}\ \partial\Omega$(强 Dirichlet;内罚仅作用于内部边) | 内罚间断 Galerkin(IPDG) | 正六边形 |
| `biharmonic_ipcg` | $\Delta^2 u = f \quad (k \ge 2)$ | $u = g_1$ 且 $\partial_{\mathbf{n}} u = g_2 \ \text{on}\ \partial\Omega$($u$ 强 Dirichlet,$\partial_{\mathbf{n}} u$ 用 Nitsche 弱施加) | C⁰ 内罚 Galerkin(C⁰-IPCG) | 单位正方形 |
| `poisson_mixed_hdg` | $\boldsymbol{\sigma} = \nabla u,\quad -\nabla\!\cdot\!\boldsymbol{\sigma} = f$ | $u = g \ \text{on}\ \partial\Omega$(弱 Dirichlet) | HDG 混合元 $(\mathrm{vec}DP_k,\, DP_k)$ | 单位正方形 |

- **边界数据全部取自制造解**:$g = u_{\text{exact}}\big|_{\partial\Omega}$,双调和的法向导数 $g_2 = \nabla u_{\text{exact}}\cdot\mathbf{n}$。
  因此四个程序本质上都是**非齐次 Dirichlet** 类问题的收敛性验证;目前 C++ 版未实现 Neumann 边界
  (MATLAB 版的 1D 混合求解器才支持逐端选 D/N)。
- **IPDG / C⁰-IPCG** 支持三种内罚变体:**SIPG**($\beta=1$)、**NIPG**($\beta=-1$)、**IIPG**($\beta=0$)。
- **mixed HDG** 在解出 $(\boldsymbol{\sigma}_h, u_h)$ 后,逐单元求解局部 Poisson 问题把 $u_h$ 投影到
  $DP_{k+1}$ 空间得到 $u_h^{\ast}$,用于展示超收敛(误差阶比 $u_h$ 高一阶)。
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

> **macOS + Apple Clang 的 OpenMP 注意**:Apple 自带的 `clang` 需要 `libomp` 才能用
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

构建产物(四个可执行文件 + 三个静态库)位于 `build/`:

```
build/poisson_pk
build/poisson_ipdg
build/biharmonic_ipcg
build/poisson_mixed_hdg
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
```

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

## 6. 目录结构

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
    └── mixed/               # 库 hdg + 一个可执行
        ├── VecFEM.{h,cpp}      向量值间断元(σ)
        ├── AssemblyHDG.{h,cpp} 质量阵、混合算子、杂交内罚、弱边界
        ├── PostProcess.{h,cpp} 局部 Poisson 求解(超收敛后处理)
        └── mixed_main.cpp      → poisson_mixed_hdg
```

CMake 目标关系:`poisson_common`(基础)← `dg_assembly`(IPDG 装配)← `hdg`(混合元);
四个可执行各自链接所需的库。

---

## 7. 与 MATLAB 版的对应关系

| C++ 可执行 | 对应 MATLAB 主程序 |
|---|---|
| `poisson_pk` | `main_Poisson_Pk2D.m` |
| `poisson_ipdg` | `main_Poisson_DG2D.m` |
| `biharmonic_ipcg` | `main_Biharmonic_C0G2D.m` |
| `poisson_mixed_hdg` | `main_PoissonMixed_2D.m` |

C++ 版只覆盖 **2D** 求解器,不含 MATLAB 仓库里的 1D 程序、Hermite 单元和 Bratu
非线性求解器。MATLAB 版还提供求解结果与收敛曲线的可视化;C++ 版只在终端输出
误差与收敛阶表格。
