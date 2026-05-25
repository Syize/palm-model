# PALM 中建筑物壁面对空气动力与热力影响的代码分析报告

## 1. 报告目的

本文基于 PALM 源码中城市表面模型（USM, Urban Surface Model）相关实现，对建筑物壁面与邻近空气之间的动力和热力耦合进行梳理，重点关注以下问题：

1. 太阳辐射如何加热建筑物壁面；
2. 加热后的壁面如何通过感热通量影响邻近空气；
3. 这些热力变化如何进一步影响近壁湍流与边界层；
4. 相关变量、子程序与公式在代码中的实现位置。

本文主要依据如下文件：

- [src/urban_surface_mod.f90](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90)
- [src/surface_layer_fluxes_mod.f90](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_layer_fluxes_mod.f90)
- [src/diffusion_s.f90](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/diffusion_s.f90)
- [src/diffusion_u.f90](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/diffusion_u.f90)
- [src/diffusion_v.f90](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/diffusion_v.f90)
- [src/turbulence_closure_mod.f90](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/turbulence_closure_mod.f90)

---

## 2. 总体结论

从代码实现上看，PALM 中建筑壁面对空气的影响可以分成两部分：

1. **动力影响**  
   建筑壁面对流动的直接机械作用，主要通过壁面动量通量和壁面边界条件体现；这一部分主要不在 `urban_surface_mod` 内完成，而是在 `surface_layer_fluxes_mod`、`diffusion_u`、`diffusion_v` 等模块中完成。

2. **热力影响**  
   建筑壁面对空气的热作用，主要通过：
   - 辐射收支决定表皮温度；
   - 表皮温度与空气温度差决定感热交换；
   - 感热通量进入标量扩散方程；
   - 局地稳定度和浮力状态改变后，再反馈到近壁湍流和边界层结构。

对你关心的“**太阳辐射加热墙面后导致的额外湍流影响**”，代码的处理方式是：

- **有显式的墙面受热和墙-空气热交换**；
- **有热通量进入空气标量方程**；
- **有通过稳定度参数和浮力间接影响湍流的机制**；
- 但**没有看到一个专门针对竖直受热墙面自然对流羽流的独立湍流增强参数化**。

换句话说，PALM/USM 对这一过程的描述是：

> 太阳辐射先改变壁面温度，再通过壁面感热通量改变邻近空气温度和稳定度，随后通过已有的流动与湍流方程间接体现其对湍流的增强或抑制。

---

## 3. 相关数据结构与核心变量

### 3.1 建筑表面分类

在 USM 中，一个城市表面单元通常由 3 类子表面 fraction 构成：

- `ind_veg_wall = 0`：墙面 fraction（对 USM 来说这里的命名沿用了通用 surface 结构）
- `ind_pav_green = 1`：绿化墙/绿化表面 fraction
- `ind_wat_win = 2`：窗面 fraction

相关索引和参数定义见：

- [src/urban_surface_mod.f90:285](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:285)
- [src/urban_surface_mod.f90:323](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:323)

### 3.2 墙体层结构

墙体内部采用固定 4 层结构：

- `nzb_wall = 0`
- `nzt_wall = 3`
- `nzw = 4`

定义见：

- [src/urban_surface_mod.f90:279](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:279)

墙体、窗体和绿化层的主要热学变量包括：

- `t_surf_wall`, `t_surf_window`, `t_surf_green`：表皮温度
- `t_wall`, `t_window`, `t_green`：内部层温度
- `lambda_h`, `lambda_h_window`, `lambda_h_green`：热导率
- `rho_c_wall`, `rho_c_window`, `rho_c_green`：体积热容
- `dz_wall`, `dz_window`, `dz_green`：层厚
- `lambda_surf`, `lambda_surf_window`, `lambda_surf_green`：表皮与第一层之间的导热系数

变量所在数据结构见：

- [src/surface_mod.f90:388](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_mod.f90:388)
- [src/surface_mod.f90:405](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_mod.f90:405)
- [src/surface_mod.f90:432](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_mod.f90:432)

### 3.3 与空气交换相关的关键变量

- `rad_sw_in`, `rad_sw_out`：短波入射、出射
- `rad_lw_in`, `rad_lw_out`：长波入射、出射
- `rad_net_l`：净辐射（代码中汇总 SW/LW）
- `r_a`：空气动力阻力
- `wghf_eb`：表皮到墙体内部的导热通量
- `wshf_eb`：壁面法向感热通量
- `shf`：供大气模块使用的感热通量
- `pt_surface`：表面位温
- `vpt_surface`：表面虚位温
- `rib`：体 Richardson 数
- `ol`：Obukhov 长度

---

## 4. 变量/子程序调用链

这一部分按“太阳辐射加热墙面 -> 形成热通量 -> 影响空气 -> 影响湍流”来梳理。

### 4.1 总体调用链

```text
辐射模块
  -> rad_sw_in, rad_sw_out, rad_lw_in, rad_lw_out
  -> urban_surface_mod: usm_surface_energy_balance
     -> 更新 t_surf_wall / t_surf_window / t_surf_green
     -> 计算 wghf_eb, wshf_eb, shf, pt_surface
  -> urban_surface_mod: usm_wall_heat_model
     -> 更新 t_wall / t_window 内部层温度
  -> diffusion_s
     -> 将 surf_usm%shf 作为城市表面热通量写入标量扩散
  -> surface_layer_fluxes_mod
     -> 用 pt_surface 与近壁空气状态计算 rib, ol, 摩擦速度与壁面通量
  -> diffusion_u / diffusion_v / turbulence_closure_mod
     -> 动量、TKE、近壁湍流响应
```

### 4.2 建立墙体导热离散：`usm_init_wall_heat_model`

子程序：

- [src/urban_surface_mod.f90:1892](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:1892)

功能：

1. 根据墙体、窗体、绿化层厚度建立网格；
2. 计算层中心距 `dz_*_center`；
3. 计算层间等效热导率 `lambda_h_layer`；
4. 初始化绿化层参数。

核心关系：

```text
zw(k) -> dz(k) -> dz_center(k) -> lambda_h_layer(k)
```

对应代码位置：

- 网格与层厚：[src/urban_surface_mod.f90:1901](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:1901)
- 层间导热率：[src/urban_surface_mod.f90:2037](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:2037)

### 4.3 表皮能量平衡：`usm_surface_energy_balance`

子程序：

- [src/urban_surface_mod.f90:4334](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4334)

这是最关键的热力耦合子程序，完成：

1. 计算空气动力阻力 `r_a`；
2. 汇总净辐射；
3. 隐式更新墙/窗/绿化表皮温度；
4. 计算感热通量、潜热通量、向墙体内部的导热通量；
5. 更新 `pt_surface`，供近壁稳定度和湍流模块使用。

关键步骤如下。

#### 步骤 A：计算空气动力阻力 `r_a`

位置：

- [src/urban_surface_mod.f90:4438](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4438)

对上向水平面，`r_a` 来自类似 LSM/MOST 的形式：

```math
r_a \approx \frac{\theta_1 - \theta_s}{\theta_* u_*}
```

代码对应：

- [src/urban_surface_mod.f90:4441](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4441)

对竖直和下向表面，`r_a` 使用 TUF3D 风格强迫对流换热公式：

```math
H = h_{ttc}(T_s - T_a)
```

```math
h_{ttc} = r_w (11.8 + 4.2 U_{eff}) - 4.0
```

```math
r_a = \frac{\rho c_p}{h_{ttc}}
```

对应代码：

- [src/urban_surface_mod.f90:4454](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4454)
- [src/urban_surface_mod.f90:4475](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4475)
- [src/urban_surface_mod.f90:4482](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4482)

这里有一个对你的研究很重要的点：

- 注释里提到了 `wstar` 和对流速度尺度；
- 但**实际公式实现只显式依赖 `Ueff` 和粗糙度 `z0`**；
- 也就是说，**竖直受热墙面引起的自然对流羽流并没有直接写进这个 `r_a` 参数化里**。

#### 步骤 B：净辐射进入表皮能量平衡

位置：

- [src/urban_surface_mod.f90:4592](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4592)

净辐射为：

```math
R_n = SW_{in} - SW_{out} + LW_{in} - LW_{out}
```

代码中存入：

- `surf % rad_net_l`

#### 步骤 C：更新表皮温度

位置：

- 墙面：[src/urban_surface_mod.f90:4652](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4652)
- 窗面：[src/urban_surface_mod.f90:4607](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4607)
- 绿化面：[src/urban_surface_mod.f90:4624](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4624)

墙面表皮的离散方程可概括为：

```math
C_s \frac{T_s^{n+1} - T_s^n}{\Delta t}
= R_n - \epsilon \sigma T_s^4 - H - G
```

其中：

- `C_s`：表皮热容，对应 `c_surface`
- `R_n`：净辐射
- `\epsilon \sigma T_s^4`：向外长波辐射
- `H`：感热交换
- `G`：向墙体内部导热

代码把这一非线性式线性化成：

```math
T_s^{n+1} = \frac{coef_1 \Delta t + C_s T_s^n}{C_s + coef_2 \Delta t}
```

对应代码：

- [src/urban_surface_mod.f90:4660](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4660)

`coef_1` 和 `coef_2` 分别包含：

- 净辐射
- 长波线性化项
- 感热项
- 潜热项（仅绿化面）
- 与第一层材料之间的导热项

#### 步骤 D：表皮向墙体内部导热

位置：

- [src/urban_surface_mod.f90:4750](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4750)

导热通量定义为：

```math
G = \lambda_{surf}(T_s - T_{wall,1})
```

对应变量：

- 墙面：`wghf_eb`
- 绿化面：`wghf_eb_green`
- 窗面：`wghf_eb_window`

#### 步骤 E：表皮向空气的感热通量

位置：

- [src/urban_surface_mod.f90:4756](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4756)

代码形式可写成：

```math
H = \frac{\rho c_p}{r_a} (T_s - T_a)
```

在代码中，墙、窗、绿化各自算一部分，再按面积 fraction 汇总为：

```math
H_{tot} = f_{wall} H_{wall} + f_{win} H_{win} + f_{green} H_{green}
```

然后：

- `wshf_eb`：总感热通量
- `shf = wshf_eb / c_p`

见：

- [src/urban_surface_mod.f90:4757](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4757)
- [src/urban_surface_mod.f90:4765](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4765)

竖直面还有密度修正：

- [src/urban_surface_mod.f90:4770](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4770)

这一步的意义是：  
USM 内部算出来的壁面热通量，要转换成与大气扩散方程一致的形式，否则竖直壁面热通量会偏大。

#### 步骤 F：更新表面位温 `pt_surface`

位置：

- [src/urban_surface_mod.f90:4680](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4680)

```math
\theta_{surface}
= \frac{f_{wall} T_{surf,wall} + f_{win} T_{surf,win} + f_{green} T_{surf,green}}{\Pi}
```

这里 `\Pi` 为 Exner 函数，对应 `exner(k)`。

`pt_surface` 是连接 USM 和近壁湍流/稳定度计算的核心变量。

### 4.4 更新墙体内部温度：`usm_wall_heat_model`

子程序：

- [src/urban_surface_mod.f90:3430](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:3430)

功能：

1. 用前面得到的 `wghf_eb` 驱动墙体表层升温或降温；
2. 在墙体内部求解 1D 热传导；
3. 更新 `t_wall`、`t_window`；
4. 把太阳透过窗的短波吸收分布到窗层内部。

#### 墙体内部导热方程

代码本质上是 1D 热传导方程：

```math
\rho c \frac{\partial T}{\partial t}
= \frac{\partial}{\partial z}
\left( \lambda \frac{\partial T}{\partial z} \right)
```

表面第一层离散时，外边界项由 `wghf_eb` 提供：

- [src/urban_surface_mod.f90:3490](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:3490)

内部层：

- [src/urban_surface_mod.f90:3505](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:3505)

内边界：

- [src/urban_surface_mod.f90:3513](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:3513)

#### 窗体内部短波吸收

位置：

- [src/urban_surface_mod.f90:3532](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:3532)

代码中先由透射率和反照率构造单侧非反射比例，再用指数衰减给出窗体内部的层吸收：

```math
Q_{sw}(z) \propto e^{-k z}
```

这说明：

- 窗面短波吸收是**体吸收**；
- 墙面短波吸收主要先进入**表皮净辐射**再传入墙体。

### 4.5 热通量如何进入大气：`diffusion_s`

子程序：

- [src/diffusion_s.f90:86](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/diffusion_s.f90:86)

这是墙面热作用真正进入空气标量方程的地方。

该程序先构造 6 个方向的扩散通量：

- `flux_r`, `flux_l`
- `flux_n`, `flux_s`
- `flux_t`, `flux_d`

然后如果使用表面通量，就用 surface 通量覆盖相应方向上的扩散通量。对 USM：

- [src/diffusion_s.f90:184](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/diffusion_s.f90:184)
- [src/diffusion_s.f90:190](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/diffusion_s.f90:190)

代码显示 `surf_usm` 的标量通量对所有表面方向都可生效：

```math
\frac{\partial s}{\partial t}
= - \nabla \cdot \mathbf{F}_s
```

对应离散更新：

- [src/diffusion_s.f90:215](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/diffusion_s.f90:215)

这意味着：

- 对屋顶，热通量以垂直通量进入空气；
- 对竖直墙面，热通量以侧向通量进入邻近空气单元。

这一步是“墙面热羽流”物理机制的第一层数值体现。

### 4.6 表面稳定度与近壁湍流反馈：`surface_layer_fluxes_mod`

相关位置：

- `calc_rib`：[src/surface_layer_fluxes_mod.f90:1116](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_layer_fluxes_mod.f90:1116)
- `calc_usws/calc_vsws`：[src/surface_layer_fluxes_mod.f90:1683](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_layer_fluxes_mod.f90:1683)
- 主调部分：[src/surface_layer_fluxes_mod.f90:384](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_layer_fluxes_mod.f90:384)

#### 体 Richardson 数

稳定度诊断公式：

```math
Ri_b = \frac{g z_{mo} (\theta_1 - \theta_{surface})}
             {U^2 \theta_1 + \varepsilon}
```

对应代码：

- [src/surface_layer_fluxes_mod.f90:1137](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_layer_fluxes_mod.f90:1137)

意义：

- 若墙面被太阳加热，`pt_surface` 升高；
- 则 `pt1 - pt_surface` 变小，甚至为负；
- `Ri_b` 变得更不稳定；
- 进而影响 `ol`、摩擦速度、壁面通量和湍流交换。

#### 自由对流速度尺度

程序 `calc_uvw_abs_s` 中有一个局地自由对流速度尺度：

```math
w_{lfc} \sim \left( \frac{g}{\theta} z_{mo} shf \right)^{1/3}
```

对应代码：

- [src/surface_layer_fluxes_mod.f90:580](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_layer_fluxes_mod.f90:580)
- [src/surface_layer_fluxes_mod.f90:589](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_layer_fluxes_mod.f90:589)

但注释明确指出：

- 该自由对流速度尺度**只用于水平面**；
- 对竖直墙面并不直接引入这个 free-convection enhancement。

这也是为什么说：  
PALM 在这里对“受热竖直墙面羽流”的处理更偏**间接反馈**而不是**直接参数化**。

### 4.7 墙面对动量方程的直接影响

动量通量不是由 `urban_surface_mod` 直接计算的，而是由 `surface_layer_fluxes_mod` 给出，再进入 `diffusion_u` 与 `diffusion_v`。

#### 水平面相关动量通量

- `usws`, `vsws`

相关公式在：

- [src/surface_layer_fluxes_mod.f90:1683](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_layer_fluxes_mod.f90:1683)
- [src/surface_layer_fluxes_mod.f90:1744](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_layer_fluxes_mod.f90:1744)

进入动量扩散：

- [src/diffusion_u.f90:224](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/diffusion_u.f90:224)
- [src/diffusion_v.f90:222](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/diffusion_v.f90:222)

#### 竖直面相关动量通量

主调位置说明：

- [src/surface_layer_fluxes_mod.f90:421](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_layer_fluxes_mod.f90:421)

也就是说：

- 建筑壁面的**机械阻挡和壁面剪切**确实存在；
- 但它属于 surface-layer / diffusion 模块的职责；
- `urban_surface_mod` 更偏向**热学闭合**。

---

## 5. 相关科学公式与变量解释

下面把和“太阳加热墙面 -> 影响空气 -> 影响湍流”最相关的公式单独整理。

### 5.1 表面净辐射

```math
R_n = SW_{in} - SW_{out} + LW_{in} - LW_{out}
```

变量：

- `rad_sw_in`：入射短波辐射
- `rad_sw_out`：反射短波辐射
- `rad_lw_in`：入射长波辐射
- `rad_lw_out`：出射长波辐射
- `rad_net_l`：代码中汇总后的净辐射

代码：

- [src/urban_surface_mod.f90:4594](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4594)

### 5.2 墙面表皮能量平衡

概念式：

```math
C_s \frac{\partial T_s}{\partial t}
= R_n - \epsilon \sigma T_s^4 - H - G
```

变量解释：

- `C_s`：表皮热容，`c_surface`
- `T_s`：墙面表皮温度，`t_surf_wall`
- `\epsilon`：表面发射率，`emissivity`
- `\sigma`：Stefan-Boltzmann 常数，`sigma_sb`
- `H`：墙面到空气的感热通量
- `G`：表皮到墙体第一层的导热通量

对应代码：

- [src/urban_surface_mod.f90:4652](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4652)

### 5.3 表皮到墙体内部的导热

```math
G = \lambda_{surf}(T_s - T_{wall,1})
```

变量：

- `lambda_surf`：表皮与墙体首层之间的导热系数
- `t_wall(nzb_wall,m)`：墙体第一层温度
- `wghf_eb`：导热通量

代码：

- [src/urban_surface_mod.f90:4750](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4750)

### 5.4 墙体内部热传导

```math
\rho c \frac{\partial T}{\partial t}
= \frac{\partial}{\partial z}
\left( \lambda \frac{\partial T}{\partial z} \right)
```

变量：

- `rho_c_wall`：体积热容
- `lambda_h_layer`：层间导热率
- `t_wall`：墙体内部各层温度

代码：

- [src/urban_surface_mod.f90:3490](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:3490)
- [src/urban_surface_mod.f90:3505](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:3505)

### 5.5 墙面对空气的感热通量

近似写成：

```math
H = \frac{\rho c_p}{r_a}(T_s - T_a)
```

对应代码中先定义：

```math
f_{shf} = \frac{\rho c_p}{r_a}
```

再计算：

```math
H = f_{shf}(T_s - T_a)
```

代码位置：

- `f_shf`：[src/urban_surface_mod.f90:4495](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4495)
- `wshf_eb`：[src/urban_surface_mod.f90:4757](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4757)

### 5.6 空气动力阻力

对竖直墙面：

```math
r_a = \frac{\rho c_p}{r_w (11.8 + 4.2 U_{eff}) - 4.0}
```

变量：

- `r_a`：空气动力阻力
- `Ueff`：近壁有效风速
- `r_w`：相对粗糙度，代码中由 `z0` 和参考混凝土粗糙度隐式表达
- `z0`：动量粗糙度长度

代码：

- [src/urban_surface_mod.f90:4475](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4475)
- [src/urban_surface_mod.f90:4482](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4482)

### 5.7 体 Richardson 数

```math
Ri_b = \frac{g z_{mo} (\theta_1 - \theta_s)}
             {U^2 \theta_1 + \varepsilon}
```

变量：

- `g`：重力加速度
- `z_mo`：近壁层参考高度
- `theta_1`：第一层空气位温，`pt1`
- `theta_s`：表面位温，`pt_surface`
- `U`：近壁平行风速，`uvw_abs`

物理意义：

- `Ri_b < 0`：不稳定，有利于热对流和湍流增强；
- `Ri_b > 0`：稳定，抑制湍流交换。

代码：

- [src/surface_layer_fluxes_mod.f90:1137](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_layer_fluxes_mod.f90:1137)

### 5.8 自由对流速度尺度

```math
w_* \sim \left( \frac{g}{\theta} z_{mo} shf \right)^{1/3}
```

代码中类似量：

- `w_lfc`

代码：

- [src/surface_layer_fluxes_mod.f90:589](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/surface_layer_fluxes_mod.f90:589)

但注意：  
这个自由对流增强只对**水平面**显式生效。

### 5.9 标量扩散方程中的壁面热通量

概念上，大气标量方程使用：

```math
\frac{\partial s}{\partial t} = - \nabla \cdot \mathbf{F}_s
```

其中 `s` 可以是位温、湿度等标量。

代码中，USM 壁面热通量会覆盖对应方向上的扩散通量：

- [src/diffusion_s.f90:190](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/diffusion_s.f90:190)
- [src/diffusion_s.f90:215](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/diffusion_s.f90:215)

这一步是墙面热量进入空气控制方程的关键接口。

---

## 6. 对“太阳加热墙面导致额外湍流”的物理解释

基于上述代码，可以把物理链条总结为：

### 6.1 正向过程

1. 墙面接受太阳短波和环境长波；
2. 表皮温度 `t_surf_wall` 升高；
3. 墙-气温差增大；
4. 壁面对空气的感热通量 `shf` 增强；
5. 邻近空气升温，近壁稳定度减弱甚至转为不稳定；
6. 局地湍流交换增强，近壁边界层结构被改变。

### 6.2 代码中这种影响主要体现在哪些环节

主要体现在三个层次：

1. **USM 层面**  
   辐射和导热决定壁面热通量。

2. **标量扩散层面**  
   壁面热通量进入空气温度方程，尤其对竖直面是“侧向热通量”。

3. **近壁稳定度/湍流层面**  
   `pt_surface` 通过 `Ri_b` 和 `ol` 影响近壁交换参数。

### 6.3 需要特别注意的限制

从实现看，PALM 对这一过程的处理仍有几条重要限制：

1. **竖直受热墙面的自然对流羽流没有独立参数化**  
   竖直壁面 `r_a` 主要依赖 `Ueff`，不是显式依赖壁面浮力羽流强度。

2. **自由对流速度尺度只用于水平面**

3. **竖直绿墙是 workaround**
   当前实现里竖直绿墙温度直接取墙温：
   - [src/urban_surface_mod.f90:3953](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:3953)

因此，如果你研究的是：

- 屋顶太阳加热后的热湍流增强；
- 街谷中整体热力不稳定度增强；

那么当前代码框架是比较合理的。

但如果你研究的是：

- 竖直受热墙面产生的贴壁热羽流；
- 该羽流对街谷湍流的局地增强；
- 朝阳面/背阳面之间由壁面热羽流引起的非对称近壁湍流结构；

那么需要意识到：  
**PALM 现有实现更像是“通过热通量和稳定度的间接影响”来表现它，而不是显式解析该壁面热羽流过程。**

---

## 7. 可直接引用的总结性表述

可以将本报告浓缩为如下表述：

> 在 PALM 的城市表面模型中，建筑物壁面的太阳辐射加热主要通过表皮能量平衡和墙体导热模型来描述。短波和长波辐射首先决定墙面表皮温度，随后墙面与邻近空气之间的温差通过空气动力阻力参数化转化为感热通量，该通量再作为边界通量进入大气标量扩散方程。由此产生的近壁增温会改变表面位温、局地 Richardson 数和 Obukhov 长度，并进一步通过近壁稳定度与浮力反馈影响湍流交换。需要注意的是，现有实现对竖直受热墙面的影响主要表现为热通量和稳定度的间接作用，而不是通过显式的壁面自然对流羽流参数化来表示。

---

## 8. 后续建议

如果你准备继续深入分析，我建议下一步优先做两件事：

1. 沿着 `surf_usm%shf -> 位温方程 -> 浮力项 -> TKE/速度方程` 再追一次完整调用链；
2. 对比太阳照射强/弱时街谷内靠墙网格单元的 `pt`, `shf`, `rib`, `km`, `e`（若使用 TKE 闭合）输出，验证“额外湍流影响”到底主要表现在哪个量上。

如果后续需要，我可以继续在这份报告基础上补一章：

- “`shf` 进入位温方程和浮力项的更深层调用链”
- 或者补一张“变量依赖关系图/流程图”。
