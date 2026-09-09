# 演示 II 型超新星光变曲线的简化输运模型

本仓库给出一个教学用简化模型的完整数值实现，用于演示 II 型超新星光变曲线平台阶段的形成。模型把辐射输运简化为直杆内两束反向飞行的粒子，粒子以固定速度自由飞行，发生碰撞后等概率地改变方向。单位长度上的碰撞概率随粒子密度变化，在某一阈值附近陡然升高，用来模仿恒星抛射物中氢复合前后不透明度的巨大反差。该模型的演化规律在数学上是一对非线性电报方程，采用算子分裂格式求解，对流步使用显式迎风差分，碰撞步有精确解析解，普通个人电脑数秒内即可完成全部计算。

数值结果表明，杆端的出射通量曲线会自然出现早期弛豫、近似恒定的平台、快速下降三个阶段，与 II-P 型超新星光变曲线的形态定性相符。改变总权重（对应总辐射能）与杆长，还可以演示平台长短的变化以及向 II-L 型形态的过渡。

## 文件说明

- `telegraph_solver.py`：模型求解器。运行后依次完成三组实验（基准算例、改变总权重、固定密度改变杆长），数据以 npz 格式存入 `data/`，并在终端打印守恒检验等核验信息。
- `make_figures.py`：绘图脚本。读取 `data/` 中的数据，在 `figures/` 下生成五张图。

## 运行环境

- Python 3（建议 3.8 及以上），依赖 numpy 与 matplotlib。
- 图的坐标轴标签为中文，系统需装有中文字体（如 SimHei 或 Droid Sans Fallback），否则标签可能显示为方框，但不影响计算与数据。

## 使用方法

```bash
python3 telegraph_solver.py
python3 make_figures.py
```

第一条命令完成全部三组实验的计算，第二条命令生成全部插图。模型参数（杆长、阈值密度、网格数、终止时间等）集中在 `telegraph_solver.py` 中 `solve()` 函数的参数表里，可直接修改后重新运行。

## 许可证

本仓库代码以 MIT 许可证发布，详见 [LICENSE](LICENSE)。

---

# A Simplified Transport Model Demonstrating Type II Supernova Light Curves

This repository provides a complete numerical implementation of a simple heuristic model that demonstrates how the plateau phase of Type II supernova light curves arises. Radiation transport is reduced to two counter-propagating beams of particles inside a one-dimensional rod. Particles stream freely at a fixed speed and, upon collision, reverse direction with equal probability. The collision probability per unit length depends on the local particle density and rises sharply around a threshold value, mimicking the strong opacity contrast of a stellar ejecta before and after hydrogen recombination. The model is governed by a pair of nonlinear telegraph equations and is solved with an operator-splitting scheme, using explicit upwind differencing for the streaming step and an exact analytic solution for the collision step. A full run takes only seconds on an ordinary computer.

The flux escaping from the rod ends naturally exhibits three successive stages, namely an early relaxation, a nearly constant plateau, and a rapid decline, in qualitative agreement with the shape of Type II-P light curves. Varying the total weight (which plays the role of the total radiated energy) and the rod length further demonstrates how the plateau length changes and how the morphology transitions toward the Type II-L case.

## Files

- `telegraph_solver.py`: the model solver. Running it performs three sets of experiments in sequence (a baseline case, a total-weight scan, and a rod-length scan at fixed density), saves the data as npz files into `data/`, and prints conservation checks and other diagnostics to the terminal.
- `make_figures.py`: the plotting script. It reads the data in `data/` and produces five figures under `figures/`.

## Requirements

- Python 3 (3.8 or later recommended), with numpy and matplotlib installed.
- The axis labels in the figures are in Chinese, so a Chinese font (e.g. SimHei or Droid Sans Fallback) is needed for proper rendering; missing fonts only affect the labels, not the computation or the data.

## Usage

```bash
python3 telegraph_solver.py
python3 make_figures.py
```

The first command runs all three experiments, and the second produces all the figures. Model parameters (rod length, threshold density, grid resolution, final time, and so on) are collected in the argument list of `solve()` in `telegraph_solver.py` and can be modified directly before rerunning.

## License

The code in this repository is released under the MIT License; see [LICENSE](LICENSE).
