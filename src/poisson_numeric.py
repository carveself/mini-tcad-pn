"""Finite-difference Poisson solve with prescribed depletion charge.

This is a linear discretization check, not a self-consistent carrier solver.
Only charge and contact potentials use the junction's analytical parameters.
Interior numerical potentials are obtained by solving a tridiagonal system.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pn_junction import (EPS_SI, Q, Junction, analytical, argument_parser,
                         plt, save_figure)


def solve_poisson(x: np.ndarray, rho: np.ndarray, left: float, right: float,
                  permittivity: float = EPS_SI) -> np.ndarray:
    """Solve -psi''=rho/eps on a uniform grid with Dirichlet contact values.

    Central differences: -psi[i-1]+2psi[i]-psi[i+1]=h^2 rho[i]/eps.
    Thomas elimination takes O(N) time and memory (no dense matrix/inverse).
    """
    x, rho = np.asarray(x, dtype=float), np.asarray(rho, dtype=float)
    if x.ndim != 1 or x.size < 5 or rho.shape != x.shape:
        raise ValueError("x and rho must be 1D, same length, with at least 5 nodes")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(rho)):
        raise ValueError("Grid and charge must be finite")
    if not all(np.isfinite(v) for v in (left, right, permittivity)) or permittivity <= 0:
        raise ValueError("Boundary values must be finite; permittivity must be positive")
    steps = np.diff(x)
    if steps[0] <= 0 or not np.allclose(steps, steps[0], rtol=1e-9, atol=0):
        raise ValueError("Use a strictly increasing, uniformly spaced grid")
    h = steps[0]
    diagonal = np.full(x.size - 2, 2.0)
    rhs = h**2 * rho[1:-1] / permittivity
    rhs[0] += left
    rhs[-1] += right
    for i in range(1, len(rhs)):
        factor = -1.0 / diagonal[i - 1]
        diagonal[i] += factor
        rhs[i] -= factor * rhs[i - 1]
    psi = np.empty_like(x)
    psi[0], psi[-1] = left, right
    psi[-2] = rhs[-1] / diagonal[-1]
    for i in range(len(rhs) - 2, -1, -1):
        psi[i + 1] = (rhs[i] + psi[i + 2]) / diagonal[i]
    return psi


def depletion_charge(j: Junction, x: np.ndarray) -> np.ndarray:
    """Average piecewise-constant charge over each node's dual cell.

    Exact overlap lengths avoid arbitrary sampling of a charge discontinuity.
    This uses depletion edges, not the analytical interior potential.
    """
    h = x[1] - x[0]
    lo, hi = x - h / 2, x + h / 2
    p_length = np.maximum(0, np.minimum(hi, 0) - np.maximum(lo, -j.xp))
    n_length = np.maximum(0, np.minimum(hi, j.xn) - np.maximum(lo, 0))
    return Q * (-j.na * p_length + j.nd * n_length) / h


def compare(j: Junction, nodes: int) -> dict:
    if not isinstance(nodes, int) or nodes < 21:
        raise ValueError("Use at least 21 integer grid nodes")
    x = np.linspace(-j.xp - .25 * j.width, j.xn + .25 * j.width, nodes)
    h = x[1] - x[0]
    rho = depletion_charge(j, x)
    psi = solve_poisson(x, rho, 0.0, j.vbi)
    field = -np.gradient(psi, h, edge_order=2)
    exact = analytical(j, x)
    source = h**2 * rho[1:-1] / EPS_SI
    residual = 2 * psi[1:-1] - psi[:-2] - psi[2:] - source
    return {"x": x, "rho": rho, "psi": psi, "field": field, "exact": exact,
            "metrics": {"nodes": nodes, "h_m": h,
                        "max_potential_error_V": float(np.max(np.abs(psi - exact["psi"]))),
                        "rms_potential_error_V": float(np.sqrt(np.mean((psi - exact["psi"])**2))),
                        "max_field_error_V_m": float(np.max(np.abs(field - exact["field"]))),
                        "scaled_discrete_residual": float(np.max(np.abs(residual)) / np.max(np.abs(source))),
                        "net_sheet_charge_C_m2": float(np.sum(rho) * h)}}


def plot_comparison(j: Junction, out: Path, nodes: int = 1601) -> None:
    result = compare(j, nodes)
    x = result["x"] * 1e6
    exact = result["exact"]
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True, constrained_layout=True)
    for ax, key, label in zip(axes[:2], ("psi", "field"), ("Potential (V)", "Electric field (V/m)")):
        ax.plot(x, exact[key], label="Analytical depletion model", lw=2)
        ax.plot(x, result[key], "--", label="Finite difference", lw=1.6)
        ax.set_ylabel(label)
        ax.legend()
    axes[2].plot(x, (result["psi"] - exact["psi"]) * 1e6)
    axes[2].set(ylabel="Potential error (uV)", xlabel=r"Position x ($\mu$m)")
    axes[0].set_title(f"Poisson comparison: prescribed depletion charge, {nodes} nodes")
    for ax in axes:
        ax.grid(alpha=.25)
    save_figure(fig, out, "poisson_comparison.png")

    studies = [compare(j, n)["metrics"] for n in (201, 401, 801, 1601, 3201)]
    hs = np.array([r["h_m"] for r in studies])
    errors = np.array([r["max_potential_error_V"] for r in studies])
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    ax.loglog(hs * 1e9, errors, "o-", label="Max potential error")
    ax.loglog(hs * 1e9, errors[0] * (hs / hs[0])**2, "--", label="O(h^2) guide")
    ax.set(xlabel="Grid spacing h (nm)", ylabel="Max potential error (V)",
           title="Mesh refinement (charge-cell alignment can change the error)")
    ax.grid(True, which="both", alpha=.25)
    ax.legend()
    save_figure(fig, out, "poisson_convergence.png")
    (out / "poisson_metrics.json").write_text(json.dumps(
        {"model": "linear Poisson with prescribed depletion charge",
         "parameters": j.summary(), "comparison": result["metrics"], "refinement": studies}, indent=2) + "\n")
    np.savetxt(out / "poisson_profiles.csv", np.column_stack(
        (result["x"], result["rho"], exact["psi"], result["psi"], exact["field"], result["field"])),
        delimiter=",", header="x_m,rho_cell_C_m3,psi_analytic_V,psi_numeric_V,E_analytic_V_m,E_numeric_V_m", comments="")
    print(json.dumps(result["metrics"], indent=2))


def main() -> None:
    parser = argument_parser("Finite-difference Poisson solution and analytical comparison")
    parser.add_argument("--nodes", type=int, default=1601, help="Uniform grid nodes, >=21")
    args = parser.parse_args()
    try:
        plot_comparison(Junction(args.na, args.nd, args.temperature), args.output, args.nodes)
    except ValueError as exc:
        parser.error(str(exc))
    print(f"Saved Poisson comparison and convergence to {args.output.resolve()}")


if __name__ == "__main__":
    main()
