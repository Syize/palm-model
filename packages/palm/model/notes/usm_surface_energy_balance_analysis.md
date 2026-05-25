# `usm_surface_energy_balance` 实现细节解析

## 1. 目的与位置

目标子程序位于：

- [src/urban_surface_mod.f90:4334](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4334)

调用位置位于：

- [src/urban_surface_mod.f90:4316](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4316)

调用顺序是：

```text
usm_energy_balance
  -> usm_surface_energy_balance
  -> usm_green_heat_model
  -> usm_wall_heat_model
```

这说明 `usm_surface_energy_balance` 的职责不是更新墙体内部温度，而是更靠前一步：

1. 读取当前辐射与近壁空气状态；
2. 计算表皮层（wall/window/green）的能量平衡；
3. 得到新的表皮温度；
4. 由表皮温度反推出感热通量、潜热通量和向内部材料的导热通量；
5. 把这些通量和表面状态写回 `surf_usm`，供后续墙体导热、标量扩散和湍流模块使用。

可以把它理解为：  
**这是 USM 中“建筑表皮与大气交换”的核心闭合子程序。**

---

## 2. 总体结构

按实现顺序，这个 subroutine 可以拆成 10 个逻辑块：

1. 初始化本地变量和控制量
2. 设置 surface fraction（wall / window / green）
3. 预计算热力常数与近壁湿度
4. 计算空气动力阻力 `r_a`
5. 为绿化表面计算蒸散相关参数
6. 汇总净辐射
7. 构造 wall / window / green 的表皮能量平衡方程
8. 更新表皮温度与 Runge-Kutta 倾向项
9. 由新表皮温度反算各种通量与表面诊断量
10. 判断是否需要强制再调用辐射模块

下面按源码顺序详细解释。

---

## 3. 逻辑块 1：初始化与基本控制量

对应代码：

- [src/urban_surface_mod.f90:4338](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4338)
- [src/urban_surface_mod.f90:4396](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4396)

### 3.1 变量类别

局部变量主要分为几类：

- 控制量：
  - `runge_l`：是否使用 Runge-Kutta 时间推进
  - `horizontal`：surface 是否为水平面
  - `force_radiation_call_l_v`：是否触发额外辐射调用

- 几何/索引量：
  - `i, j, k`
  - `i_off, j_off, k_off`
  - `m`：surface element 编号

- 热力学与交换系数：
  - `rho_cp`
  - `rho_lv`
  - `f_shf`, `f_shf_window`, `f_shf_green`
  - `coef_1`, `coef_2`
  - `coef_window_1`, `coef_window_2`
  - `coef_green_1`, `coef_green_2`

- 水分与植被相关：
  - `qv1`, `q_s`
  - `f1`, `f2`, `f3`
  - `r_canopy`
  - `f_qsws`, `f_qsws_veg`, `f_qsws_liq`

### 3.2 指针别名

```fortran
surf => surf_usm
```

这里把 `surf_usm` 绑定到局部指针 `surf`，后续代码就统一通过 `surf % ...` 访问 USM 表面数据结构。这是典型的简化写法，不改变物理含义。

### 3.3 时间推进控制

```fortran
runge_l = (timestep_scheme(1:5) == 'runge')
```

后续如果是 Runge-Kutta，会维护表皮温度和液态水储量的历史倾向项。

---

## 4. 逻辑块 2：设置 wall / window / green fraction

对应代码：

- [src/urban_surface_mod.f90:4400](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4400)

这里先做两件事：

1. 判断 surface 是否是水平面：

```fortran
horizontal = upward .or. downward
```

用途是后面在竖直面上做密度修正。

2. 设置每个 surface element 的三种 tile fraction：

- `frac_wall`
- `frac_win`
- `frac_green`

### 4.1 spinup 特殊处理

如果 `spinup_phase = .true.`，代码强制：

- `frac_wall = 1`
- `frac_win = 0`
- `frac_green = 0`

也就是在 spinup 阶段把所有建筑表面都当成“纯墙面”处理。

这有两个直接效果：

1. 窗面和绿化面的热容量、潜热过程、辐射响应都先不参与；
2. 可降低模型刚性，让初始旋转过程更稳定。

这是一个数值上的简化，而不是物理上的真实建筑表面构成。

---

## 5. 逻辑块 3：预计算热力学常数与空气湿度

对应代码：

- [src/urban_surface_mod.f90:4416](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4416)

### 5.1 `rho_cp`

```fortran
rho_cp(m) = c_p * hyp(k) / (r_d * surf % pt1(m) * exner(k))
```

这实际是在计算：

```math
\rho c_p
```

其中空气密度通过状态方程重构：

```math
\rho = \frac{p}{R_d T}
```

这里：

- `hyp(k)`：与气压相关
- `pt1(m)`：第一层空气位温
- `exner(k)`：Exner 函数

这个量后面会用于：

- 空气动力阻力换算
- 感热通量系数 `f_shf`

### 5.2 `rho_lv`

如果有绿化面，则计算：

```fortran
rho_lv = rho_cp / c_p * l_v
```

也就是：

```math
\rho l_v
```

用于把蒸发质量通量换成潜热通量。

### 5.3 `qv1`

```fortran
qv1(m) = q(k,j,i)
```

这是相邻第一层大气网格中的比湿，用于：

- 绿化面的蒸散
- 饱和亏缺计算
- `qsws`（潜热/水汽通量）计算

如果 `humidity = .false.`，则统一设为 0。

---

## 6. 逻辑块 4：计算空气动力阻力 `r_a`

对应代码：

- [src/urban_surface_mod.f90:4438](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4438)

这是整个 subroutine 里最关键的交换参数之一。它决定墙面和空气之间的感热、潜热交换强度。

### 6.1 上向水平面：沿用 LSM 风格

如果 surface 是 `upward`，则：

```fortran
surf % r_a = (pt1 - pt_surface) / (ts * us)
```

可以理解为：

```math
r_a \approx \frac{\theta_1 - \theta_s}{\theta_* u_*}
```

其中：

- `pt1`：第一层空气位温
- `pt_surface`：表面位温
- `ts`：温度尺度
- `us`：摩擦速度

特点：

- 更接近 MOST/LSM 的近地层做法；
- 使用上一时刻已有的 `ts`、`us`，因为当前预报步尚未完成所有更新。

### 6.2 竖直面和下向面：强迫对流换热公式

如果不是 `upward`，代码使用 TUF3D 风格的换热系数：

```math
H = h_{ttc}(T_s - T_a)
```

```math
h_{ttc} = r_w (11.8 + 4.2 U_{eff}) - 4.0
```

再换成阻力形式：

```math
r_a = \frac{\rho c_p}{h_{ttc}}
```

代码实现：

```fortran
ueff = max( velocity_magnitude , lower_limit )
surf % r_a = rho_cp / (z0 * d_roughness_concrete * (11.8 + 4.2 * ueff) - 4.0)
```

### 6.3 `ueff` 的物理含义

`ueff` 用的是邻近标量点上的三维速度模长：

- `u`
- `v`
- `w`

这代表竖直墙附近由大气 resolved flow 提供的强迫通风能力。

### 6.4 重要实现特征

虽然注释提到了 `wstar` 和对流速度尺度，但**实际公式中没有显式把墙面热浮力引起的自然对流速度加进去**。  
因此：

- 竖直受热墙面的换热增强，主要通过 `Tsfc - Tair` 进入感热通量；
- 但没有单独的“热壁羽流换热增强项”进入 `r_a`。

### 6.5 数值限制

代码把 `r_a` 限制在：

- 最小 `1 s/m`
- 最大 `300 s/m`

目的是防止：

- 风速过小时换热过强；
- 阻力过大导致数值不稳定。

### 6.6 导出给 wall/window/green

后面直接令：

```fortran
r_a_window = r_a
r_a_green  = r_a
```

也就是说，当前实现中 wall、window、green 在空气动力阻力上默认共用同一交换尺度，差异主要来自辐射、热容和水分过程。

---

## 7. 逻辑块 5：绿化面蒸散参数

对应代码：

- [src/urban_surface_mod.f90:4502](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4502)
- [src/urban_surface_mod.f90:4530](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4530)

这一段只在 `frac_green > 0` 时有意义。

### 7.1 `f2`：土壤水分胁迫因子

对上向绿化面，先用根区含水量构造 `m_total`，再根据：

- `wilt`
- `fc`

得到 `f2`：

- 太干时接近 0
- 足够湿时为 1

这体现了植物缺水时蒸腾减弱。

对竖直绿化面，代码直接设：

```fortran
f2 = 1
```

即不考虑竖直绿墙根区水分限制。这是一个简化。

### 7.2 `f1`：短波辐射控制因子

```fortran
f1 = min(1, radiation_function(rad_sw_in))
```

物理意义：

- 白天短波强，气孔更容易打开；
- 夜间短波弱，蒸腾受到抑制。

### 7.3 `f3`：饱和亏缺因子

先根据表皮温度求饱和水汽压 `e_s`，再与空气水汽压 `e` 比较：

```fortran
f3 = exp(-g_d * (e_s - e))
```

物理意义：

- 空气越干，饱和亏缺越大；
- 植物蒸腾越受抑制。

### 7.4 `r_canopy`

```fortran
r_canopy = r_canopy_min / (lai * f1 * f2 * f3)
```

这是典型的植被气孔阻力构造方式。

### 7.5 `q_s`, `dq_s_dt`

代码计算：

- 表面饱和比湿 `q_s`
- 其对温度的导数 `dq_s_dt`

后者用于表皮能量平衡中的线性化。

### 7.6 `f_qsws`

这一组量：

- `f_qsws_veg`
- `f_qsws_liq`
- `f_qsws`

本质上是把潜热/蒸发交换写成线性系数：

```math
LE \sim f_{qsws} (q_{air} - q_s + correction)
```

对上向绿化面：

- 同时考虑植被蒸腾和表面液态水蒸发；
- 用 `c_liq` 把两部分分开。

对竖直绿化面：

- 不考虑液态水储层蒸发；
- 只保留植被蒸腾那一部分。

---

## 8. 逻辑块 6：汇总净辐射

对应代码：

- [src/urban_surface_mod.f90:4592](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4592)

这里定义：

```fortran
rad_net_l = rad_sw_in - rad_sw_out + rad_lw_in - rad_lw_out
```

物理上就是：

```math
R_n = SW_{in} - SW_{out} + LW_{in} - LW_{out}
```

这是 wall/window/green 三类表皮能量平衡的共同辐射驱动项。

注意注释写的是 “Add LW up so that it can be removed in prognostic equation”，说明作者后面用线性化方式处理 `\sigma T^4`，因此在构造离散方程时会额外把某些长波项拆开重组。

---

## 9. 逻辑块 7：构造 wall / window / green 的表皮能量平衡方程

对应代码：

- 窗面：[src/urban_surface_mod.f90:4598](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4598)
- 绿化面：[src/urban_surface_mod.f90:4618](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4618)
- 墙面：[src/urban_surface_mod.f90:4644](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4644)

这一段是整个 subroutine 的数学核心。

### 9.1 方程的共同形式

每一种表皮都在求解近似形式：

```math
C_s \frac{T_s^{n+1} - T_s^n}{\Delta t}
= R_n - \epsilon \sigma (T_s^{n+1})^4 - H(T_s^{n+1}) - LE(T_s^{n+1}) - G(T_s^{n+1})
```

由于：

- 长波辐射是四次非线性；
- 饱和比湿随温度变化也是非线性；

代码把它线性化为：

```math
T_s^{n+1} =
\frac{coef_1 \Delta t + C_s T_s^n}
     {C_s + coef_2 \Delta t}
```

其中：

- `coef_1`：把已知项和线性化后的常数项合并
- `coef_2`：把线性化后的温度系数合并

### 9.2 窗面方程

窗面只有：

- 净辐射
- 长波线性化
- 感热
- 向窗体内部导热

没有蒸散项。

对应：

```fortran
coef_window_1 = rad_net_l + 4 eps sigma T^4 + f_shf_window * pt1 + lambda_surf_window * t_window(1)
coef_window_2 = 4 eps sigma T^3 + lambda_surf_window + f_shf_window / exner
```

### 9.3 绿化面方程

绿化面多了一项潜热交换：

```fortran
+ f_qsws * (qv1 - q_s + dq_s_dt * t_surf_green)
```

并在 `coef_green_2` 中多出：

```fortran
+ f_qsws * dq_s_dt
```

这说明作者采用的是“对 `q_s(T)` 做一阶 Taylor 展开”的线性化方式。

### 9.4 墙面方程

墙面没有水分过程，因此最简单，只包含：

- 净辐射
- 长波
- 感热
- 向墙体首层导热

写成概念式：

```math
C_{wall,surf}\frac{dT_{wall,surf}}{dt}
= R_n - LW_{up} - H - G
```

### 9.5 一个实现上的细节

墙面注释写道：

> `Todo: Adjust to tile approach. So far, emissivity for wall (element 0) is used`

说明当前墙面能量平衡在 emissivity 处理上仍有简化，tile-level 的更细致辐射参数化还没有完全展开。

---

## 10. 逻辑块 8：更新表皮温度与 Runge-Kutta 倾向项

对应代码：

- [src/urban_surface_mod.f90:4660](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4660)

### 10.1 先更新表皮温度

分别得到：

- `t_surf_wall_p`
- `t_surf_window_p`
- `t_surf_green_p`

其中 `_p` 表示新的时间层。

### 10.2 再叠加 RK 历史项

```fortran
t_surf_*_p = t_surf_*_p + dt_3d * tsc(3) * tt_surface_*_m
```

所以这个子程序实际上不是单纯的一步显式更新，而是与 PALM 整体时间推进一致地使用多步/多阶段结构。

### 10.3 生成新的 `pt_surface`

```fortran
pt_surface =
    ( frac_wall  * t_surf_wall_p
    + frac_win   * t_surf_window_p
    + frac_green * t_surf_green_p ) / exner
```

这一步非常关键，因为：

- `pt_surface` 是后续稳定度计算直接使用的量；
- 它把三类 tile 的温度按面积分数混合成一个给近壁大气使用的“等效表面位温”。

### 10.4 `vpt_surface`

当前代码直接令：

```fortran
vpt_surface = pt_surface
```

并在注释中承认：

- 这并不完全正确；
- 尤其对 walls/windows，并没有显式构造真实的 `q_surface`。

也就是说，这里对表面虚位温的表示是近似的。

### 10.5 更新 Runge-Kutta 倾向项

通过：

```fortran
stend_* = (new - old - RK_history_part) / (dt * tsc(2))
```

得到真实倾向后，更新：

- `tt_surface_wall_m`
- `tt_surface_window_m`
- `tt_surface_green_m`

供下一 RK 子步使用。

---

## 11. 逻辑块 9：检查是否需要强制重新调用辐射模块

对应代码：

- [src/urban_surface_mod.f90:4717](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4717)

判断条件是：

- 若任一表皮温度在一个子步内变化超过 `1 K`
- 且 `unscheduled_radiation_calls = .true.`

则设置：

```fortran
force_radiation_call_l_v(m) = .true.
```

其目的很明确：

- 表皮温度变化太快时，原来的辐射收支已经过时；
- 如果不更新辐射，能量平衡会失真，甚至导致数值不稳定。

这个机制本质上是一个 **subcycling / adaptive recoupling** 思想：  
当表皮温度快变时，动态增加辐射调用频率。

---

## 12. 逻辑块 10：由新表皮温度反算通量与诊断量

对应代码：

- [src/urban_surface_mod.f90:4730](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4730)

这段代码的作用不是再求温度，而是把前面得到的 `t_surf_*_p` 翻译成各种物理通量。

### 12.1 修正后的净辐射

```fortran
rad_net_l = rad_net_l + sigma * emissivity * (T_new^4 - T_old^4)
```

虽然注释写着 “rad_net_l is never used!”，但这一步是在用新的表皮温度修正净长波部分，使诊断量在当前子步末更一致。

### 12.2 向内部材料的导热通量 `wghf_eb`

```fortran
wghf_eb = lambda_surf * (t_surf_wall_p - t_wall(first_layer))
```

同理还有：

- `wghf_eb_green`
- `wghf_eb_window`

这些量是下一步 `usm_wall_heat_model` 和 `usm_green_heat_model` 的边界驱动。

### 12.3 感热通量 `wshf_eb`

代码写法：

```fortran
wshf_eb =
  - f_shf        * (pt1 - t_surf_wall_p   / exner) * frac_wall
  - f_shf_window * (pt1 - t_surf_window_p / exner) * frac_win
  - f_shf_green  * (pt1 - t_surf_green_p  / exner) * frac_green
```

这等价于按 tile 分别算：

```math
H_i = \frac{\rho c_p}{r_{a,i}} (T_{s,i} - T_a)
```

再按面积分数组合。

### 12.4 `shf`

```fortran
shf = wshf_eb / c_p
```

这是传给扩散和 surface layer 模块使用的运动学感热通量形式。

若开启室内模型，还会加上：

```fortran
waste_heat / c_p
```

即建筑废热的附加感热输入。

### 12.5 竖直面的密度修正

如果不是水平面：

```fortran
shf = shf * (r_d * pt1 * exner) / hyp
```

原因在注释里说得很清楚：

- 水平面进入 `diffusion_s` 时，密度会自然抵消；
- 竖直面不会自动抵消；
- 不修正的话，壁面热通量会高估约 15%-20%。

这说明作者对不同方向 surface flux 在标量扩散方程中的单位处理是有意识区分的。

---

## 13. 绿化面潜热与液态水储层的后处理

对应代码：

- [src/urban_surface_mod.f90:4780](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4780)

这一大段只对 `humidity = .true.` 且 `frac_green > 0` 的 surface 生效。

### 13.1 计算 `qsws`

形式上：

```fortran
qsws = - f_qsws * (qv1 - q_s + dq_s_dt * T_old - dq_s_dt * T_new)
qsws = qsws / l_v
```

这里：

- 若空气比湿低于表面饱和比湿，则蒸发为正；
- 若空气过饱和，则会出现结露趋势。

### 13.2 分解为 `qsws_veg` 与 `qsws_liq`

对上向绿化面：

- `qsws_veg`：植被蒸腾
- `qsws_liq`：液态水储层蒸发/结露

这有利于后续分别更新液态水库和土壤水分。

### 13.3 `r_s`

代码还反推出总表面阻力：

```fortran
r_s = inferred_total_resistance - r_a_green
```

这是一个诊断量，用于表征绿化面从空气到叶面/表面的总水汽交换阻力。

### 13.4 降水、结露与液态水库

如果 `precipitation = .true.`：

- 雨水可以补给液态水库 `m_liq_usm`
- 如果空气饱和导致结露，结露水也会进入液态水库

### 13.5 液态水库时间推进

通过：

```fortran
tend = -qsws_liq / (rho_l * l_v)
```

更新：

- `m_liq_usm_p`

再限制在：

- `0 <= m_liq_usm_p <= m_liq_max`

这本身并不严格守恒，代码注释也明确承认超界裁剪会破坏守恒。

### 13.6 竖直绿化面的简化

对竖直或下向绿化面：

- 仍可算 `qsws` 与 `qsws_veg`
- 但不考虑液态水库蒸发
- 竖直面 `qsws` 同样要做密度修正

这再次说明：  
竖直绿墙在现有实现中比水平绿化面简单得多。

### 13.7 没有绿化面时

如果 `frac_green = 0`，则直接：

```fortran
r_s = 1.0e10
```

表示几乎无水汽交换。

---

## 14. 最后一段：跨进程汇总 `force_radiation_call`

对应代码：

- [src/urban_surface_mod.f90:4881](/home/syize/Documents/Python/PALM/model-code/palm_model_system-v25.04/packages/palm/model/src/urban_surface_mod.f90:4881)

这里把本进程上的 `force_radiation_call_l_v` 汇总成单个逻辑量：

```fortran
force_radiation_call_l = any(force_radiation_call_l_v)
```

如果并行运行，则再用：

```fortran
MPI_ALLREDUCE(..., MPI_LOR, ...)
```

把所有进程上的需求做逻辑或。

效果是：

- 只要有任何一个 surface element 温度变化过快；
- 整个并行域下一子步都可以触发额外辐射调用。

这是保证并行一致性的必要步骤。

---

## 15. 从物理与数值角度看，这个 subroutine 在做什么

把全部实现压缩成一句话：

> `usm_surface_energy_balance` 先用辐射、近壁空气状态和材料参数构造 wall/window/green 三类建筑表皮的线性化能量平衡方程，求出新的表皮温度，再据此反算感热、潜热和向内部材料的导热通量，并把这些通量写回城市表面数据结构，供后续墙体导热、标量扩散和近壁湍流参数化使用。

如果按“输入-计算-输出”来总结：

### 输入

- 辐射：
  - `rad_sw_in/out`
  - `rad_lw_in/out`

- 近壁空气：
  - `pt1`
  - `q`
  - `u, v, w`

- 表面状态：
  - `t_surf_wall`, `t_surf_window`, `t_surf_green`
  - `t_wall`, `t_window`, `t_green`
  - `m_liq_usm`
  - `swc`, `fc`, `wilt`

- 参数：
  - `emissivity`
  - `z0`
  - `lambda_surf`
  - `c_surface`
  - `lai`
  - `r_canopy_min`

### 计算

1. `r_a`
2. 绿化蒸散阻力与线性化系数
3. `rad_net_l`
4. `t_surf_*_p`
5. `pt_surface`
6. `wghf_eb`, `wshf_eb`, `shf`
7. `qsws`, `qsws_veg`, `qsws_liq`
8. `m_liq_usm_p`

### 输出

写回 `surf_usm` 的关键量包括：

- `pt_surface`
- `vpt_surface`
- `r_a`, `r_a_window`, `r_a_green`
- `wghf_eb`, `wghf_eb_green`, `wghf_eb_window`
- `wshf_eb`
- `shf`
- `qsws`, `qsws_veg`, `qsws_liq`
- `r_s`
- `tt_surface_*_m`
- `m_liq_usm_p`
- `force_radiation_call_l`

---

## 16. 对你研究问题最关键的几点理解

如果你的问题是“太阳辐射加热建筑壁面后，代码里究竟发生了什么”，那么 `usm_surface_energy_balance` 给出的答案是：

1. 太阳辐射先进入 `rad_net_l`；
2. `rad_net_l` 通过表皮能量平衡抬升 `t_surf_wall`；
3. 更高的 `t_surf_wall` 会增大壁面到空气的感热通量 `shf`；
4. `shf` 随后进入标量扩散与近壁稳定度计算；
5. 于是热力影响再间接反馈到湍流和边界层结构。

但还要加一句：

> 这个 subroutine 对竖直受热墙面的处理核心是“更新表皮温度和感热通量”，而不是显式求解一个贴壁自然对流羽流模型。

也就是说，它对“热墙增强湍流”的表达是：

- **显式处理热通量**
- **间接影响稳定度与湍流**

而不是：

- **显式参数化竖直热羽流本身**

---

## 17. 一句话总结

`usm_surface_energy_balance` 是 PALM 城市表面模型中把**辐射、表皮储热、墙-气感热、绿化潜热、表皮-内部导热**统一闭合到同一个时间步里的核心子程序；它的输出决定了建筑表面对空气的即时热力强迫强度。

