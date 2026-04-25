# ===============================================================
# 抽取自仓库 [current]: examples/plot_validation_curves.py
# 行号区段：91..167
# 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
# ===============================================================
def make_combined_plot(
    *,
    it11_theta: list[float],
    it11_exp: list[float],
    it11_model: list[float],
    ds_theta: list[float],
    ds_exp: list[float],
    ds_model: list[float],
    out_file: Path,
    annotation: Dict[str, Any],
    dsigma_y_label: str,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4), dpi=140)
    fig.suptitle("Experiment vs Faddeev Simulation (solver Tlab benchmark)", fontsize=15, fontweight="bold")

    best_tlab = annotation.get("best_tlab", math.nan)
    target_tlab = annotation.get("target_tlab", math.nan)
    fig.text(
        0.5,
        0.945,
        f"Target Tlab = {target_tlab:.3f} MeV, Best solver Tlab = {best_tlab:.3f} MeV",
        ha="center",
        va="center",
        fontsize=10,
    )

    ax0 = axes[0]
    ax0.scatter(it11_theta, it11_exp, s=26, color="#111827", label="Experiment data", zorder=3)
    ax0.plot(it11_theta, it11_model, color="#dc2626", linewidth=2.2, label="Faddeev simulation", zorder=2)
    ax0.set_title("iT11")
    ax0.set_xlabel("theta_cm (deg)")
    ax0.set_ylabel("iT11")
    ax0.grid(True, alpha=0.25, linestyle="--", linewidth=0.6)
    ax0.legend(loc="best")
    ax0.text(
        0.03,
        0.97,
        (
            f"RMSE = {annotation.get('it11_rmse', math.nan):.5f}\n"
            f"MAE = {annotation.get('it11_mae', math.nan):.5f}\n"
            f"Max |err| = {annotation.get('it11_max_abs_error', math.nan):.5f}"
        ),
        transform=ax0.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#f3f4f6", edgecolor="#9ca3af", alpha=0.95),
    )

    ax1 = axes[1]
    ax1.scatter(ds_theta, ds_exp, s=26, color="#111827", label="Experiment data", zorder=3)
    ax1.plot(ds_theta, ds_model, color="#dc2626", linewidth=2.2, label="Faddeev simulation", zorder=2)
    ax1.set_yscale("log")
    ax1.set_title("dSigma/dOmega")
    ax1.set_xlabel("theta_cm (deg)")
    ax1.set_ylabel(dsigma_y_label)
    ax1.grid(True, alpha=0.25, linestyle="--", linewidth=0.6)
    ax1.legend(loc="best")
    ax1.text(
        0.03,
        0.97,
        (
            f"RMSE = {annotation.get('dsigma_rmse', math.nan):.5f}\n"
            f"Rel. RMSE = {annotation.get('dsigma_rel_rmse', math.nan):.5f}\n"
            f"MAE = {annotation.get('dsigma_mae', math.nan):.5f}"
        ),
        transform=ax1.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#f3f4f6", edgecolor="#9ca3af", alpha=0.95),
    )

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.92))
    fig.savefig(out_file)
    plt.close(fig)

