# -*- coding: utf-8 -*-
"""Lecture 4 插图：采样手段对比、拉丁方设计、单样本 t 检验。

用法：python make_figures.py  （输出到 figures/）
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from scipy.stats import qmc

sns.set_style("whitegrid")
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

OUT = "figures"
import os
os.makedirs(OUT, exist_ok=True)

C_RAND, C_GRID, C_LHS = "#d62728", "#7f7f7f", "#1f77b4"


# ---------------------------------------------------------------- fig 1: 2D
def fig_sampling_2d():
    n = 9
    rng = np.random.default_rng(42)
    pts_rand = rng.random((n, 2))
    g = np.linspace(1 / 6, 5 / 6, 3)
    pts_grid = np.array([(x, y) for x in g for y in g])
    pts_lhs = qmc.LatinHypercube(d=2, seed=42).random(n)
    # Latin square 型：与 LHS 相同的 bin 占用结构（每行/列 bin 恰一点），
    # 但层号配对与落点位置是确定性的（bin 中心、固定错位配对）
    xs_ls = (np.arange(n) + 0.5) / n
    ys_ls = (((np.arange(n) * 4 + 2) % n) + 0.5) / n
    pts_ls = np.column_stack([xs_ls, ys_ls])

    fig, axes = plt.subplots(1, 4, figsize=(17.5, 4.2), sharex=True, sharey=True)
    for ax, pts, title, c, bins in zip(
        axes,
        [pts_rand, pts_grid, pts_lhs, pts_ls],
        ["Random sampling（随机采样）", "Grid sampling（网格采样 = 全组合）",
         "Latin hypercube（拉丁超立方采样）", "Latin square 型（确定性 bin 配对）"],
        [C_RAND, C_GRID, C_LHS, "#9467bd"],
        [False, False, True, True],
    ):
        if bins:
            for k in range(1, n):
                ax.axhline(k / n, ls="--", color="gray", lw=0.6, zorder=1)
                ax.axvline(k / n, ls="--", color="gray", lw=0.6, zorder=1)
        ax.scatter(pts[:, 0], pts[:, 1], s=90, c=c, edgecolors="white", zorder=3)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xticks(np.linspace(0, 1, 5)); ax.set_yticks(np.linspace(0, 1, 5))
        ax.set_title(f"{title}\n(n = {n})", fontsize=11)
        ax.set_xlabel("factor X"); ax.set_ylabel("factor Y")
        ax.set_aspect("equal")
    fig.suptitle("同样 9 个样本点：对因子空间的覆盖能力对比\n"
                 "（后两幅 bin 占用结构相同——每行/列 bin 恰一点；LHS 随机配对，Latin square 型为确定性配对）",
                 fontsize=12.5, y=1.06)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_sampling_2d.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------- fig 2: 3D 投影覆盖
def fig_sampling_3d():
    """3 个因子，对比四种采样：
    - grid 2^3=8 点（立方体顶点）：投影只有 4 个不同位置；
    - LHS 9 点：一维分层均匀，但二维投影只是 9 个散点，无组合均衡；
    - Latin square（3 阶，X=行, Y=列, Z=符号）9 点：任意二维投影恰为均匀 3×3 网格；
    - orthogonal array L9(3^3) 9 点：同样任意二维投影恰为均匀 3×3 网格。
    OA(9,3,3,2) 读作 (行, 列, 符号) 三元组时就是一个 3 阶 Latin square
    （c3 = c1 + c2 mod 3 正是循环拉丁方）——第三、四行点集完全相同，殊途同归。
    """
    pts_grid = np.array([(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)],
                        dtype=float)
    pts_lhs = qmc.LatinHypercube(d=3, seed=7).random(9)
    # 手工构造 OA(9, 3, 3, 2)：c3 = c1 + c2 (mod 3)，任意两列含全部 9 种水平组合
    c1 = np.repeat([0, 1, 2], 3)
    c2 = np.tile([0, 1, 2], 3)
    c3 = (c1 + c2) % 3
    pts_ls3d = (np.column_stack([c1, c2, c3]) + 0.5) / 3.0
    pts_oa = pts_ls3d.copy()

    fig = plt.figure(figsize=(14.5, 16.2))
    titles = ["Grid sampling：$2^3 = 8$ 点",
              "Latin hypercube：$n = 9$ 点（仅保证一维分层）",
              "Latin square（3 阶：X=行, Y=列, Z=符号）：$n = 9$ 点",
              "Orthogonal array $L_9(3^3)$：$n = 9$ 点\n"
              "（与 Latin square 行殊途同归——同一点集）"]
    data = [pts_grid, pts_lhs, pts_ls3d, pts_oa]
    colors = [C_GRID, C_LHS, "#9467bd", "#2ca02c"]
    proj_note = ["{} 个不同位置", "{} 个不同位置（散布，无组合均衡）",
                 "{} 个不同位置（均匀 3×3）", "{} 个不同位置（均匀 3×3）"]
    for row in range(4):
        pts = data[row]
        ax3d = fig.add_subplot(4, 3, row * 3 + 1, projection="3d")
        ax3d.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=70, c=colors[row],
                     edgecolors="white", depthshade=False)
        ax3d.set_title(titles[row], fontsize=11)
        ax3d.set_xlabel("X"); ax3d.set_ylabel("Y"); ax3d.set_zlabel("Z")
        for col, (i, j, name) in enumerate([(0, 2, "XZ"), (1, 2, "YZ")], start=1):
            ax = fig.add_subplot(4, 3, row * 3 + col + 1)
            ax.scatter(pts[:, i], pts[:, j], s=70, c=colors[row],
                       edgecolors="white", zorder=3)
            ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05)
            n_pos = len(np.unique(np.round(pts[:, [i, j]], 6), axis=0))
            ax.set_title(f"投影到 {name} 平面：" + proj_note[row].format(n_pos),
                         fontsize=11)
            ax.set_xlabel(name[0]); ax.set_ylabel(name[1])
            ax.set_aspect("equal")
    fig.suptitle("三因子空间中的采样对比：grid 投影重叠；LHS 只保一维分层；"
                 "Latin square 与 orthogonal array 殊途同归——二维投影都是均匀 3×3",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{OUT}/fig_sampling_3d.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------ fig 3: 拉丁方
def fig_latin_square():
    latin = np.array([[0, 1, 2, 3],
                      [1, 2, 3, 0],
                      [2, 3, 0, 1],
                      [3, 0, 1, 2]])
    labels = np.array([["A", "B", "C", "D"][(v)] for v in latin.flat]).reshape(4, 4)
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    sns.heatmap(latin, annot=labels, fmt="", cmap="Set2", cbar=False,
                linewidths=1.5, linecolor="white",
                annot_kws={"fontsize": 20, "weight": "bold"}, ax=ax,
                xticklabels=["位置 1", "位置 2", "位置 3", "位置 4"],
                yticklabels=["11 时", "12 时", "13 时", "14 时"])
    ax.set_xlabel("猪圈位置（confounding factor 2）")
    ax.set_ylabel("喂食顺序（confounding factor 1）")
    ax.set_title("Latin square design：4 种饲料（treatment factor 的 4 个 level）\n"
                 "每行、每列都恰好各含 A/B/C/D 一次 → 同时均衡两个干扰因素",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_latin_square.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------- fig 3b: Latin square vs LHS 的区分
def fig_latin_vs_lhs():
    """同取 n=5，直观对比：
    - Latin square：离散 treatment level 的确定性组合设计；
    - LHS：连续空间上的随机分层采样（bin 占用模式与拉丁方同构，故名）。
    """
    n = 5
    latin = np.array([[(i + j) % n for j in range(n)] for i in range(n)])
    letters = np.array([[chr(ord("A") + v) for v in row] for row in latin])

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.8))

    # 左：Latin square（组合设计，离散水平）
    ax = axes[0]
    sns.heatmap(latin, annot=letters, fmt="", cmap="Set3", cbar=False,
                linewidths=1.5, linecolor="white",
                annot_kws={"fontsize": 17, "weight": "bold"}, ax=ax,
                xticklabels=[f"列 {j + 1}" for j in range(n)],
                yticklabels=[f"行 {i + 1}" for i in range(n)])
    ax.set_title("Latin square（组合设计）\n"
                 "格内是 treatment factor 的离散 level（A–E）；\n"
                 "确定性构造：每行、每列各含每个 level 一次", fontsize=11)

    # 右：LHS（随机分层采样，连续坐标）
    ax = axes[1]
    rng = np.random.default_rng(3)
    px, py = rng.permutation(n), rng.permutation(n)
    xs = (px + rng.random(n)) / n
    ys = (py + rng.random(n)) / n
    for k in range(1, n):
        ax.axhline(k / n, ls="--", color="gray", lw=0.8, zorder=1)
        ax.axvline(k / n, ls="--", color="gray", lw=0.8, zorder=1)
    ax.scatter(xs, ys, s=110, c="#1f77b4", edgecolors="white", zorder=3)
    for x, y in zip(xs, ys):
        ax.annotate("", (x, y))  # no-op, keeps loop explicit
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks((np.arange(n) + 0.5) / n)
    ax.set_xticklabels([f"bin {j + 1}" for j in range(n)])
    ax.set_yticks((np.arange(n) + 0.5) / n)
    ax.set_yticklabels([f"bin {i + 1}" for i in range(n)])
    ax.set_aspect("equal")
    ax.set_title("Latin hypercube sampling（随机分层采样）\n"
                 "点是连续坐标；每维等分 n 层，层号随机配对：\n"
                 "每个行 bin、列 bin 恰含一个点（二维时即拉丁方图案）", fontsize=11)

    fig.suptitle("同名不同物：Latin square 是离散水平的确定性设计，"
                 "LHS 是连续空间的随机采样——二维 bin 占用模式同构是命名的由来",
                 fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(f"{OUT}/fig_latin_vs_lhs.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------- fig 4: 单样本 t 检验
def fig_t_test():
    df, t_calc, t_crit = 24, 1.692, stats.t.ppf(0.95, 24)
    p_val = 1 - stats.t.cdf(t_calc, df)
    x = np.linspace(-4, 4.5, 1200)
    y = stats.t.pdf(x, df)

    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    ax.plot(x, y, color="#1f77b4", lw=2)
    # α = 0.05 的拒绝域（临界值右侧）
    ax.fill_between(x, y, where=(x >= t_crit), color="#d62728", alpha=0.45,
                    label=f"拒绝域：$t \\geq t_{{0.05,24}} = {t_crit:.3f}$，面积 $= \\alpha = 0.05$")
    # 计算值右侧的尾部面积 = P 值
    ax.fill_between(x, y, where=(x >= t_calc), color="#ff7f0e", alpha=0.30, hatch="//",
                    edgecolor="#ff7f0e", linewidth=0,
                    label=f"$P$ 值：$t \\geq 1.692$ 的面积 $= {p_val:.3f}$")
    ax.axvline(t_calc, color="#ff7f0e", ls="--", lw=1.5)
    ax.axvline(t_crit, color="#d62728", ls="--", lw=1.5)
    ax.set_xlabel("$t$"); ax.set_ylabel("density")
    ax.set_title("Example 6.3：$\\nu = 24$ 的 $t$ 分布，one-sided $\\alpha = 0.05$\n"
                 "$t_{\\mathrm{calc}} = 1.692 < 1.711 = t_{0.05,24}$ → $P = "
                 f"{p_val:.3f} > 0.05$ → 不拒绝 $H_0$", fontsize=12)
    ax.legend(fontsize=10, loc="upper left")

    # 尾部放大 inset：1.692 与 1.711 在主图上几乎重合，放大后可见 P 值区域略宽于 α
    axins = ax.inset_axes([0.55, 0.42, 0.42, 0.5])
    xs = np.linspace(1.5, 2.4, 600)
    ys = stats.t.pdf(xs, df)
    axins.plot(xs, ys, color="#1f77b4", lw=2)
    axins.fill_between(xs, ys, where=(xs >= t_crit), color="#d62728", alpha=0.45)
    axins.fill_between(xs, ys, where=(xs >= t_calc), color="#ff7f0e", alpha=0.30,
                       hatch="//", edgecolor="#ff7f0e", linewidth=0)
    axins.axvline(t_calc, color="#ff7f0e", ls="--", lw=1.5)
    axins.axvline(t_crit, color="#d62728", ls="--", lw=1.5)
    axins.annotate("$t_{\\mathrm{calc}}$\n$=1.692$", (t_calc, 0.008), xytext=(1.53, 0.030),
                   fontsize=10, color="#ff7f0e",
                   arrowprops=dict(arrowstyle="->", color="#ff7f0e"))
    axins.annotate("$t_{0.05,24}$\n$=1.711$", (t_crit, 0.008), xytext=(1.85, 0.045),
                   fontsize=10, color="#d62728",
                   arrowprops=dict(arrowstyle="->", color="#d62728"))
    axins.set_title("尾部放大：$P$ 值区域（橙）比 $\\alpha$ 区域（红）宽出一小段",
                    fontsize=10)
    axins.set_xlim(1.5, 2.4)
    ax.indicate_inset_zoom(axins, edgecolor="gray", alpha=0.5)

    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_t_test.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_sampling_2d()
    fig_sampling_3d()
    fig_latin_square()
    fig_latin_vs_lhs()
    fig_t_test()
    print("done ->", OUT)
