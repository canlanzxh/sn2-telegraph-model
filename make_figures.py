#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
由 data/ 中的计算结果生成全部插图（figures/ 下 5 张 PNG）。
先运行 telegraph_solver.py 生成数据，再运行本脚本。

图件清单：
  fig1_chi.png         碰撞概率随粒子权重密度的变化
  fig2a_density.png    基准算例：密度分布的演化
  fig2b_flux.png       基准算例：杆端总出射通量
  fig3_flux_L.png      固定密度、不同杆长 L 下的杆端总出射通量
  fig4_flux_wtot.png   不同总权重 W_tot 下的杆端总出射通量

作图约定：
- 通量一律画两端出射之和（两侧对称，为单侧的两倍），
  通量图一律采用线性纵轴并收窄显示范围（突出平台期）；
- 曲线一律用灰度 + 线型区分，黑白打印也可辨认；
- 坐标轴标签含中文，需要系统装有中文字体（如 SimHei 或 Droid Sans Fallback）。
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["SimHei", "Droid Sans Fallback"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "stix"   # 新罗马斜体风格的数学符号
plt.rcParams["font.size"] = 7.5        # 图内所有文字统一为六号（7.5 pt）
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["xtick.major.width"] = 0.8    # 刻度线宽与边框一致
plt.rcParams["ytick.major.width"] = 0.8

CM = 1.0 / 2.54                        # 厘米 → 英寸
DPI = 400

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(HERE, "data"))
FIG_DIR = os.path.normpath(os.path.join(HERE, "figures"))

# 灰度下也可区分的线型与配色（5 条曲线用）
STYLES = [
    dict(color="0.10", ls="-"),
    dict(color="0.35", ls="--"),
    dict(color="0.50", ls="-."),
    dict(color="0.65", ls=":"),
    dict(color="0.25", ls=(0, (3, 1, 1, 1, 1, 1))),
]

# 8 条曲线用的互不相同组合（深浅 × 线型）
STYLES8 = [
    dict(color="0.05", ls="-"),
    dict(color="0.42", ls="-"),
    dict(color="0.05", ls="--"),
    dict(color="0.42", ls="--"),
    dict(color="0.05", ls="-."),
    dict(color="0.42", ls="-."),
    dict(color="0.05", ls=":"),
    dict(color="0.42", ls=":"),
]


def save(fig, name):
    # bbox_inches="tight" 裁掉图四周的白边
    fig.savefig(os.path.join(FIG_DIR, name), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def fig1_chi():
    """图1：碰撞概率随粒子权重密度的变化。"""
    fig, ax = plt.subplots(figsize=(6.0 * CM, 4.5 * CM))

    w = np.linspace(0, 1.6, 600)
    w_c, chi0, delta = 0.5, 1.0, 0.05
    chi = 0.5 * chi0 * (1 + np.tanh((w - w_c) / delta))
    ax.plot(w, chi / chi0, color="0.1", lw=1.8, label=rf"$\delta={delta}$")
    ax.plot(w, (w > w_c).astype(float), color="0.45", lw=1.3, ls="--",
            label=r"$\delta\to0$")
    ax.axvline(w_c, color="0.6", lw=0.9, ls=":")
    ax.text(w_c + 0.03, 0.10, "$w_c$")
    ax.set_xlabel(r"粒子权重密度 $w$")
    ax.set_ylabel(r"$\chi(w)/\chi_0$")
    ax.set_xlim(0, 1.6)
    ax.set_ylim(-0.05, 1.1)
    ax.tick_params(labelsize=7.5)
    ax.legend(loc="lower right", fontsize=7.5, frameon=False)

    fig.tight_layout()
    save(fig, "fig1_chi.png")


def fig2a_density():
    """图2a：基准算例，密度分布的演化。"""
    d = np.load(os.path.join(DATA_DIR, "baseline.npz"))
    x, snap_times, snaps = d["x"], d["snap_times"], d["snapshots"]
    w_c = float(d["par_w_c"])

    fig, ax = plt.subplots(figsize=(6.0 * CM, 4.5 * CM))

    show = [0.0, 2.0, 5.0, 8.0, 10.0, 12.0, 15.0, 20.0]
    for i, ts in enumerate(show):
        j = int(np.argmin(np.abs(snap_times - ts)))
        st = STYLES8[i % len(STYLES8)]
        ax.plot(x, snaps[j], lw=1.3, label=f"$t={ts:g}$", **st)
    ax.axhline(w_c, color="0.2", lw=1.1, ls=(0, (5, 2)))
    ax.set_xlabel("位置 $x$")
    ax.set_ylabel(r"粒子权重密度 $w(x,t)$")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 1.35)
    ax.tick_params(labelsize=7.5)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4,
              fontsize=7.5, frameon=False, columnspacing=0.6,
              handlelength=1.3, labelspacing=0.25, borderaxespad=0.2)

    fig.tight_layout()
    save(fig, "fig2a_density.png")


def fig2b_flux():
    """图2b：基准算例，杆端总出射通量（灰色底纹为平台期，按通量曲线表观形态划分）。"""
    d = np.load(os.path.join(DATA_DIR, "baseline.npz"))
    t, flux = d["time"], d["flux_total"]

    fig, ax = plt.subplots(figsize=(6.0 * CM, 4.5 * CM))

    ax.axvspan(2.0, 12.0, color="0.94", zorder=0)
    ax.plot(t, flux, color="0.1", lw=1.5)
    ax.set_xlabel("时间 $t$")
    ax.set_ylabel(r"总出射通量 $F(t)$")
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 1.2)
    ax.tick_params(labelsize=7.5)

    fig.tight_layout()
    save(fig, "fig2b_flux.png")


def fig3_flux_L():
    """图3：固定密度、不同杆长 L 下的杆端总出射通量。"""
    d = np.load(os.path.join(DATA_DIR, "L_scan.npz"))
    L_vals = d["L_values"]

    fig, ax = plt.subplots(figsize=(6.0 * CM, 4.5 * CM))
    for i, Lv in enumerate(L_vals):
        st = STYLES[i % len(STYLES)]
        ax.plot(d[f"time_{i}"], d[f"flux_{i}"], lw=1.4,
                label=f"$L={Lv:g}$", **st)
    ax.set_xlabel("时间 $t$")
    ax.set_ylabel(r"总出射通量 $F(t)$")
    ax.set_xlim(0, 40)
    ax.set_ylim(0, 1.6)
    ax.tick_params(labelsize=7.5)
    ax.legend(loc="upper right", fontsize=7.5, frameon=False)

    fig.tight_layout()
    save(fig, "fig3_flux_L.png")


def fig4_flux_wtot():
    """图4：不同总权重 W_tot 下的杆端总出射通量。"""
    d = np.load(os.path.join(DATA_DIR, "wtot_scan.npz"))
    W_vals = d["W_tot_values"]

    fig, ax = plt.subplots(figsize=(6.0 * CM, 4.5 * CM))
    for i, W in enumerate(W_vals):
        st = STYLES[i % len(STYLES)]
        ax.plot(d[f"time_{i}"], d[f"flux_{i}"], lw=1.4,
                label=f"$W_{{\\mathrm{{tot}}}}={W:g}$", **st)
    ax.set_xlabel("时间 $t$")
    ax.set_ylabel(r"总出射通量 $F(t)$")
    ax.set_xlim(0, 45)
    ax.set_ylim(0, 3.8)
    ax.tick_params(labelsize=7.5)
    ax.legend(loc="upper right", fontsize=7.5, frameon=False)

    fig.tight_layout()
    save(fig, "fig4_flux_wtot.png")


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    fig1_chi()
    fig2a_density()
    fig2b_flux()
    fig3_flux_L()
    fig4_flux_wtot()
    print("已生成 5 张图，保存在", FIG_DIR)


if __name__ == "__main__":
    main()
