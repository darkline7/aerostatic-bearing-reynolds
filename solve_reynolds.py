"""Solve the 1D modified Reynolds equation of an aerostatic bearing (Miettinen et al., 2024).

    (1/r) d/dr( p r h^3 dp/dr ) = (6 k / h_p) (p^2 - p_s^2)

BCs: dp/dr = 0 at r = 0 (applied at r_min = 1 um to avoid the 1/r singularity),
     p = p_a at r = R.
First-order system with y = [p, dp/dr]:
    dy0/dr = y1
    dy1/dr = (6 k / (h^3 h_p y0)) (y0^2 - p_s^2) - (1/r + y1/y0) y1
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_bvp

H = 10e-6
R = 18.415e-3
HP = 4.5e-3
PA = 0.1e6
RHO = 1.225
MU = 1.78e-5
EPS = 0.2
K_DEFAULT = 1.44e-15
R_MIN = 1e-6
PS_LIST = [0.6e6, 0.4e6, 0.2e6]
MESH_SIZES = [10, 100, 1000, 10000]
TOL = 1e-4


def reynolds_ode(r, y, ps, k):
    p, dp = y
    d2p = (6.0 * k / (H**3 * HP * p)) * (p**2 - ps**2) - (1.0 / r + dp / p) * dp
    return np.vstack((dp, d2p))


def boundary_conditions(ya, yb, ps):
    return np.array([ya[1] * R / ps, (yb[0] - PA) / ps])


def initial_guess(r, ps):
    p = PA + (ps - PA) * (1.0 - (r / R) ** 2)
    dp = -2.0 * (ps - PA) * r / R**2
    return np.vstack((p, dp))


def solve_pressure(ps, n, k=K_DEFAULT, tol=TOL, adaptive=True):
    r = np.linspace(R_MIN, R, n)
    max_nodes = max(10000, 10 * n) if adaptive else n
    sol = solve_bvp(
        lambda r, y: reynolds_ode(r, y, ps, k),
        lambda ya, yb: boundary_conditions(ya, yb, ps),
        r,
        initial_guess(r, ps),
        tol=tol,
        max_nodes=max_nodes,
    )
    return sol


def load_capacity(sol, ps):
    r = np.linspace(R_MIN, R, 20001)
    p = sol.sol(r)[0]
    return 2.0 * np.pi * np.trapezoid((p - PA) * r, r)


def plot_pressure_distribution(results, outdir, paper_csv=None):
    fig, ax = plt.subplots(figsize=(8, 5))
    for ps, sols in results.items():
        sol = sols[max(MESH_SIZES)]
        r = np.linspace(R_MIN, R, 2000)
        ax.plot(r * 1e3, sol.sol(r)[0] / 1e6, lw=2, label=f"Simulation $P_s$ = {ps/1e6:g} MPa")
    if paper_csv and Path(paper_csv).exists():
        data = np.genfromtxt(paper_csv, delimiter=",", names=True)
        for ps_val in np.unique(data["ps_MPa"]):
            d = data[data["ps_MPa"] == ps_val]
            ax.plot(d["r_mm"], d["p_MPa"], "^", ms=5,
                    label=f"Paper data $P_s$ = {ps_val:g} MPa")
    ax.set_xlabel("Radius r (mm)")
    ax.set_ylabel("Pressure p (MPa)")
    ax.set_title("Pressure distribution along the radius (N = 10000)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "fig1_pressure_distribution.png", dpi=200)
    plt.close(fig)


def plot_mesh_convergence(results, outdir):
    fig, axes = plt.subplots(1, len(PS_LIST), figsize=(16, 4.5), sharey=True)
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(MESH_SIZES)))
    for ax, ps in zip(axes, PS_LIST):
        for n, color in zip(MESH_SIZES, colors):
            sol = results[ps][n]
            r = np.linspace(R_MIN, R, 2000)
            ax.plot(r * 1e3, sol.sol(r)[0] / 1e6, color=color, lw=1.5,
                    label=f"N = {n}")
        ax.set_title(f"$P_s$ = {ps/1e6:g} MPa")
        ax.set_xlabel("Radius r (mm)")
        ax.grid(True, alpha=0.3)
        ax.legend()
    axes[0].set_ylabel("Pressure p (MPa)")
    fig.suptitle("Mesh convergence (N = 10, 100, 1000, 10000)")
    fig.tight_layout()
    fig.savefig(outdir / "fig2_mesh_convergence.png", dpi=200)
    plt.close(fig)


def plot_error_vs_mesh(results, outdir):
    r_fine = np.linspace(R_MIN, R, 5001)
    ref = {ps: results[ps][max(MESH_SIZES)].sol(r_fine)[0] for ps in PS_LIST}
    fig, ax = plt.subplots(figsize=(7, 5))
    for ps in PS_LIST:
        errs = []
        for n in MESH_SIZES:
            p_n = results[ps][n].sol(r_fine)[0]
            errs.append(np.max(np.abs(p_n - ref[ps])) / 1e6)
        ax.semilogy(MESH_SIZES, errs, "o-", label=f"$P_s$ = {ps/1e6:g} MPa")
    ax.set_xlabel("Number of grid points N")
    ax.set_ylabel("Max |p_N - p_10000| (MPa)")
    ax.set_title("Error vs. mesh resolution (reference N = 10000)")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "fig3_error_vs_mesh.png", dpi=200)
    plt.close(fig)
    return {ps: {n: float(np.max(np.abs(results[ps][n].sol(r_fine)[0] - ref[ps])) / 1e6)
                for n in MESH_SIZES} for ps in PS_LIST}


def save_csv(results, outdir):
    n_ref = max(MESH_SIZES)
    r = np.linspace(R_MIN, R, 2000)
    for ps, sols in results.items():
        p = sols[n_ref].sol(r)[0]
        np.savetxt(
            outdir / f"pressure_Ps_{ps/1e6:g}MPa_N{n_ref}.csv",
            np.column_stack([r * 1e3, p / 1e6]),
            delimiter=",",
            header="r_mm,p_MPa",
            comments="",
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=float, default=K_DEFAULT,
                        help="permeability k in m^2 (default 1.44e-15; tuned value 1.20e-15)")
    parser.add_argument("--outdir", type=Path, default=Path(__file__).parent / "outputs")
    parser.add_argument("--paper-csv", type=str, default=None,
                        help="optional digitized paper data: r_mm,p_MPa,ps_MPa")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    results = {}
    print(f"{'Ps (MPa)':>9} {'N':>6} {'p(0) (MPa)':>12} {'W (N)':>10}")
    for ps in PS_LIST:
        results[ps] = {}
        for n in MESH_SIZES:
            sol = solve_pressure(ps, n, k=args.k, adaptive=False)
            results[ps][n] = sol
            p0 = sol.sol(R_MIN)[0] / 1e6
            w = load_capacity(sol, ps)
            print(f"{ps/1e6:>9g} {n:>6} {p0:>12.6f} {w:>10.3f}")

    errs = plot_error_vs_mesh(results, args.outdir)
    for ps in PS_LIST:
        for n in MESH_SIZES:
            if n != max(MESH_SIZES):
                print(f"Ps={ps/1e6:g} MPa, N={n:>5}: max err = {errs[ps][n]:.3e} MPa")

    plot_pressure_distribution(results, args.outdir, args.paper_csv)
    plot_mesh_convergence(results, args.outdir)
    save_csv(results, args.outdir)
    print(f"Figures and CSV saved to {args.outdir.resolve()}")


if __name__ == "__main__":
    main()
