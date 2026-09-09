#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
一维杆中非线性双束输运模型（非线性电报方程）的数值求解器

模型方程：
    ∂w+/∂t + v ∂w+/∂x = (v χ(w) / 2) (w- - w+)
    ∂w-/∂t - v ∂w-/∂x = (v χ(w) / 2) (w+ - w-)
其中 w = w+ + w- 为总权重密度，碰撞概率随密度变化：
    χ(w) = (χ0 / 2) [1 + tanh((w - w_c) / δ)],  χ0 = 1 / l0
边界为吸收边界：w-(0,t) = 0, w+(L,t) = 0。
初始条件为均匀分布：w(x,0) = w0 = W_tot / L，各向同性 w+ = w- = w0/2。

数值方法：Lie–Trotter 算子分裂。
  对流步：显式迎风格式，CFL 数取 0.4；
  碰撞步：碰撞不改变总密度，步内 χ 为常数，退化为常系数常微分方程，精确解析求解。

运行本文件即完成全部三组实验（基准、改变总权重、固定密度改变杆长），
结果存入 data/，并打印守恒检验与相变时刻等核验数据。
"""

import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(HERE, "data"))


def chi_of(w, chi0, w_c, delta):
    """密度依赖的碰撞概率（单位长度碰撞率）。"""
    return 0.5 * chi0 * (1.0 + np.tanh((w - w_c) / delta))


def solve(L=10.0, v=1.0, l0=1.0, w_c=0.5, delta=0.05,
          N=400, cfl=0.4, T_final=20.0, w0=1.25,
          snapshot_times=None):
    """
    求解非线性双束输运方程（均匀初始条件）。

    参数
    ----
    L        杆长
    v        粒子飞行速度
    l0       正常平均自由程（高密度区）
    w_c      阈值权重密度
    delta    tanh 过渡宽度
    N        空间网格数
    cfl      CFL 数
    T_final  终止时间
    w0       初始均匀密度，w0 = W_tot / L，需 w0 > w_c

    返回字典：x、time、flux_single（单侧边界通量）、flux_total（两侧之和）、
    weight（杆内总权重随时间变化）、snapshots（密度快照）、t_transition
    （全面弹道化时刻）、以及守恒检验数据。
    """
    chi0 = 1.0 / l0
    dx = L / N
    dt = cfl * dx / v
    n_steps = int(round(T_final / dt))
    x = np.linspace(dx / 2, L - dx / 2, N)

    # 均匀初始条件，各向同性
    w_plus = np.full(N, 0.5 * w0)
    w_minus = np.full(N, 0.5 * w0)
    assert w0 > w_c, "初始密度需高于阈值"

    if snapshot_times is None:
        snapshot_times = [0.0, 2.0, 5.0, 8.0, 10.0, 12.0, 15.0, 20.0]
    snapshot_times = sorted(t for t in snapshot_times if t <= T_final)
    snapshots = {0.0: w_plus + w_minus}
    snap_idx = 1

    time_hist, flux_single_hist, flux_total_hist, weight_hist = [], [], [], []
    t_now = 0.0
    t_transition = None  # 首次全部网格 w < w_c 的时刻

    for _ in range(n_steps):
        w_before = np.sum(w_plus + w_minus) * dx

        # ---- 对流步：显式迎风 ----
        wp = w_plus.copy()
        wm = w_minus.copy()
        wp[0] -= (v * dt / dx) * w_plus[0]
        wp[1:] -= (v * dt / dx) * (w_plus[1:] - w_plus[:-1])
        wm[:-1] += (v * dt / dx) * (w_minus[1:] - w_minus[:-1])
        wm[-1] -= (v * dt / dx) * w_minus[-1]
        # 吸收边界：飞出边界的粒子消失
        wp[-1] = 0.0
        wm[0] = 0.0

        # 边界通量 = 本步对流损失 / 时间步长（两端对称，单侧取半）
        loss = w_before - np.sum(wp + wm) * dx
        time_hist.append(t_now + dt / 2.0)
        flux_total_hist.append(loss / dt)
        flux_single_hist.append(loss / (2.0 * dt))
        weight_hist.append(np.sum(wp + wm) * dx)

        # ---- 碰撞步：精确解析求解 ----
        w_star = wp + wm
        dw_star = wp - wm
        lam = v * chi_of(w_star, chi0, w_c, delta)  # 步内 w 不变，故为常系数
        e = np.exp(-lam * dt)
        w_plus = 0.5 * w_star + 0.5 * dw_star * e
        w_minus = 0.5 * w_star - 0.5 * dw_star * e
        # 边界处碰撞不应产生向内的粒子，再次强制吸收边界
        w_plus[-1] = 0.0
        w_minus[0] = 0.0

        t_now += dt

        if t_transition is None:
            if np.all((w_plus + w_minus) < w_c):
                t_transition = t_now

        if snap_idx < len(snapshot_times) and abs(t_now - snapshot_times[snap_idx]) < dt / 2.0:
            snapshots[snapshot_times[snap_idx]] = (w_plus + w_minus).copy()
            snap_idx += 1

    time_arr = np.array(time_hist)
    flux_single = np.array(flux_single_hist)
    flux_total = np.array(flux_total_hist)
    weight = np.array(weight_hist)

    # 守恒检验
    W_init = w0 * L
    W_escaped = np.sum(flux_total) * dt
    W_final = weight[-1] if len(weight) else W_init
    conserv_err = abs(W_escaped + W_final - W_init) / W_init

    return {
        "x": x,
        "time": time_arr,
        "flux_single": flux_single,
        "flux_total": flux_total,
        "weight": weight,
        "snapshots": snapshots,
        "snapshot_times": sorted(snapshots.keys()),
        "t_transition": t_transition,
        "W_init": W_init,
        "W_escaped": W_escaped,
        "W_final": W_final,
        "conserv_err": conserv_err,
        "params": dict(L=L, v=v, l0=l0, chi0=chi0, w_c=w_c, delta=delta,
                       N=N, cfl=cfl, T_final=T_final, w0=w0, W_tot=w0 * L),
    }


def print_check(tag, res):
    p = res["params"]
    print(f"\n[{tag}] L={p['L']}, w0={p['w0']}, W_tot={p['W_tot']:.2f}, "
          f"L/l0={p['L']/p['l0']:.0f}, T_final={p['T_final']}")
    print(f"  相变时刻 t* = {res['t_transition']}")
    print(f"  守恒: 初始 {res['W_init']:.4f} = 逃逸 {res['W_escaped']:.4f} "
          f"+ 剩余 {res['W_final']:.4f}, 相对误差 {res['conserv_err']:.2e}")
    # 平台期平均通量（t 在 [2, 8] 内，仅供核验）
    m = (res["time"] >= 2.0) & (res["time"] <= 8.0)
    if np.any(m):
        print(f"  平台期平均单侧通量(t∈[2,8]): {np.mean(res['flux_single'][m]):.4f}")
    print(f"  初始尖峰(单侧): {np.max(res['flux_single']):.4f}")


def main():
    os.makedirs(DATA_DIR, exist_ok=True)

    # ---------- 1. 基准算例 ----------
    snap = [0.0, 2.0, 5.0, 8.0, 9.5, 10.0, 10.5, 11.0, 12.0, 15.0, 20.0]
    base = solve(T_final=20.0, w0=1.25, snapshot_times=snap)
    print_check("基准", base)
    snap_times = np.array(base["snapshot_times"])
    snap_mat = np.array([base["snapshots"][t] for t in snap_times])
    np.savez(os.path.join(DATA_DIR, "baseline.npz"),
             x=base["x"], time=base["time"],
             flux_single=base["flux_single"], flux_total=base["flux_total"],
             weight=base["weight"], snap_times=snap_times, snapshots=snap_mat,
             t_transition=base["t_transition"],
             **{f"par_{k}": v for k, v in base["params"].items()})

    # ---------- 2. 改变总权重 W_tot ----------
    w0_list = [0.625, 1.25, 2.5, 5.0, 6.25]
    Tf_list = [20.0, 20.0, 30.0, 50.0, 60.0]
    wtot_data = {}
    t_trans_list = []
    for i, (w0, Tf) in enumerate(zip(w0_list, Tf_list)):
        res = solve(T_final=Tf, w0=w0)
        print_check(f"W_tot 扫描 {i}", res)
        wtot_data[f"time_{i}"] = res["time"]
        wtot_data[f"flux_{i}"] = res["flux_total"]
        t_trans_list.append(res["t_transition"])
    np.savez(os.path.join(DATA_DIR, "wtot_scan.npz"),
             w0_values=np.array(w0_list),
             W_tot_values=np.array(w0_list) * 10.0,
             t_transition=np.array(t_trans_list),
             **wtot_data)

    # ---------- 3. 固定密度改变杆长 L ----------
    L_list = [5.0, 10.0, 15.0, 20.0]
    Tf_list = [20.0, 20.0, 40.0, 60.0]
    L_data = {}
    t_trans_list = []
    for i, (L, Tf) in enumerate(zip(L_list, Tf_list)):
        res = solve(L=L, T_final=Tf, w0=1.25)
        print_check(f"L 扫描 {i}", res)
        L_data[f"time_{i}"] = res["time"]
        L_data[f"flux_{i}"] = res["flux_total"]
        t_trans_list.append(res["t_transition"])
    np.savez(os.path.join(DATA_DIR, "L_scan.npz"),
             L_values=np.array(L_list),
             W_tot_values=1.25 * np.array(L_list),
             t_transition=np.array(t_trans_list),
             **L_data)

    print("\n全部实验完成，数据已保存至", DATA_DIR)


if __name__ == "__main__":
    main()
