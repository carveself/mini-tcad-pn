"""Analytical equilibrium silicon PN junction and ideal diffusion diode.

CLI doping inputs are in cm^-3. All electrostatics use SI units internally;
energies are in eV. Run from any directory; figures default to project/figures.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # Works on servers and computers without a display.
import matplotlib.pyplot as plt
import numpy as np

Q = 1.602176634e-19  # Elementary charge, C
KB = 1.380649e-23  # Boltzmann constant, J/K
EPS_SI = 11.7 * 8.8541878128e-12  # Silicon permittivity, F/m
ROOT = Path(__file__).resolve().parents[1]


def band_gap(temperature: float) -> float:
    """Varshni silicon band gap, eV (teaching parametrization)."""
    return 1.17 - 4.73e-4 * temperature**2 / (temperature + 636.0)


def intrinsic_density(temperature: float) -> float:
    """Approximate n_i(T), m^-3, calibrated to 1e10 cm^-3 at 300 K."""
    kb_ev = KB / Q
    return float(1e16 * (temperature / 300.0)**1.5 * np.exp(
        band_gap(300.0) / (2 * kb_ev * 300.0)
        - band_gap(temperature) / (2 * kb_ev * temperature)))


@dataclass(frozen=True)
class Junction:
    na_cm3: float = 1e16
    nd_cm3: float = 1e17
    temperature: float = 300.0

    def __post_init__(self) -> None:
        if not np.isfinite(self.temperature) or not 250 <= self.temperature <= 400:
            raise ValueError("Teaching model requires 250 <= T <= 400 K")
        for name, value in (("NA", self.na_cm3), ("ND", self.nd_cm3)):
            if not np.isfinite(value) or not 1e14 <= value <= 1e18:
                raise ValueError(f"{name} must be between 1e14 and 1e18 cm^-3")
            if value * 1e6 < 100 * self.ni:
                raise ValueError(f"{name} must be at least 100 times n_i(T)")

    @property
    def na(self) -> float:
        return self.na_cm3 * 1e6

    @property
    def nd(self) -> float:
        return self.nd_cm3 * 1e6

    @property
    def vt(self) -> float:
        return KB * self.temperature / Q

    @property
    def ni(self) -> float:
        return intrinsic_density(self.temperature)

    @property
    def vbi(self) -> float:
        return float(self.vt * np.log(self.na * self.nd / self.ni**2))

    @property
    def width(self) -> float:
        return float(np.sqrt(2 * EPS_SI * self.vbi / Q * (1 / self.na + 1 / self.nd)))

    @property
    def xp(self) -> float:
        return self.width * self.nd / (self.na + self.nd)

    @property
    def xn(self) -> float:
        return self.width * self.na / (self.na + self.nd)

    def summary(self) -> dict:
        return {"NA_cm-3": self.na_cm3, "ND_cm-3": self.nd_cm3,
                "T_K": self.temperature, "ni_cm-3": self.ni / 1e6,
                "Eg_eV": band_gap(self.temperature), "VT_V": self.vt,
                "Vbi_V": self.vbi, "W_um": self.width * 1e6,
                "xp_um": self.xp * 1e6, "xn_um": self.xn * 1e6,
                "peak_field_magnitude_V_per_m": Q * self.na * self.xp / EPS_SI}


def analytical(junction: Junction, x: np.ndarray) -> dict[str, np.ndarray]:
    """P side is x<0. psi(-xp)=0, psi(xn)=Vbi, E=-d(psi)/dx."""
    j = junction
    x = np.asarray(x, dtype=float)
    p = (x >= -j.xp) & (x < 0)
    n = (x >= 0) & (x <= j.xn)
    rho = np.zeros_like(x)
    field = np.zeros_like(x)
    psi = np.where(x >= j.xn, j.vbi, 0.0)
    rho[p], rho[n] = -Q * j.na, Q * j.nd
    field[p] = -Q * j.na * (x[p] + j.xp) / EPS_SI
    field[n] = Q * j.nd * (x[n] - j.xn) / EPS_SI
    psi[p] = Q * j.na * (x[p] + j.xp)**2 / (2 * EPS_SI)
    psi[n] = j.vbi - Q * j.nd * (j.xn - x[n])**2 / (2 * EPS_SI)
    # Choosing EF=0 fixes an arbitrary energy reference. Midgap Ei is approximate.
    ei = j.vt * np.log(j.na / j.ni) - psi
    return {"rho": rho, "field": field, "psi": psi, "Ei": ei,
            "Ec": ei + band_gap(j.temperature) / 2,
            "Ev": ei - band_gap(j.temperature) / 2, "Ef": np.zeros_like(x)}


def ideal_iv(junction: Junction, voltage: np.ndarray, area_m2: float = 1e-8,
             tau_n_s: float = 1e-6, tau_p_s: float = 1e-6) -> tuple[np.ndarray, float]:
    """Long-base diffusion diode, ideality factor=1; returns I [A], Is [A].

    Fixed teaching mobilities: electrons=1350, holes=480 cm^2/(V s).
    These do not model doping or temperature dependent mobility.
    """
    voltage = np.asarray(voltage, dtype=float)
    if not all(np.isfinite(v) and v > 0 for v in (area_m2, tau_n_s, tau_p_s)):
        raise ValueError("Area and lifetimes must be finite and positive")
    if not np.all(np.isfinite(voltage)) or np.any(np.abs(voltage / junction.vt) > 600):
        raise ValueError("Voltage is nonfinite or outside safe exponential range")
    dn, dp = 0.135 * junction.vt, 0.048 * junction.vt
    ln, lp = np.sqrt(dn * tau_n_s), np.sqrt(dp * tau_p_s)
    isat = Q * area_m2 * junction.ni**2 * (
        dn / (ln * junction.na) + dp / (lp * junction.nd))
    return isat * np.expm1(voltage / junction.vt), float(isat)


def save_figure(fig: plt.Figure, directory: Path, name: str) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(directory / name, dpi=170, bbox_inches="tight")
    plt.close(fig)


def plot_junction(j: Junction, out: Path) -> None:
    x = np.linspace(-j.xp - .3 * j.width, j.xn + .3 * j.width, 2001)
    data = analytical(j, x)
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True, constrained_layout=True)
    for ax, key, label in zip(axes, ("rho", "field", "psi"),
                              (r"Charge density (C/m$^3$)", "Electric field (V/m)", "Potential (V)")):
        ax.plot(x * 1e6, data[key], lw=2)
        ax.axvspan(-j.xp * 1e6, j.xn * 1e6, color="orange", alpha=.12)
        ax.axvline(0, color="gray", ls=":")
        ax.set_ylabel(label)
        ax.grid(alpha=.25)
    axes[-1].set_xlabel(r"Position x ($\mu$m), P side < 0 < N side")
    axes[0].set_title(f"Equilibrium silicon PN junction | T={j.temperature:g} K")
    save_figure(fig, out, "electrostatics.png")

    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    for key, label, style in (("Ec", "Ec: conduction band edge", "-"),
                              ("Ei", "Ei: intrinsic level (midgap approximation)", "--"),
                              ("Ev", "Ev: valence band edge", "-"),
                              ("Ef", "EF: equilibrium Fermi level", ":")):
        ax.plot(x * 1e6, data[key], style, lw=2, label=label)
    ax.axvspan(-j.xp * 1e6, j.xn * 1e6, color="orange", alpha=.12)
    ax.set(xlabel=r"Position x ($\mu$m)", ylabel="Energy relative to EF (eV)",
           title="Energy bands: increasing potential lowers electron energy")
    ax.legend(fontsize=9)
    ax.grid(alpha=.25)
    save_figure(fig, out, "energy_bands.png")

    voltage = np.linspace(-.5, .6, 1101)
    current, isat = ideal_iv(j, voltage)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    axes[0].plot(voltage, current * 1e3)
    axes[0].set(ylabel="Current (mA)", title="Ideal diffusion diode (A=100 x 100 um^2)")
    axes[1].plot(voltage, current)
    axes[1].set_yscale("symlog", linthresh=isat)
    axes[1].set(ylabel="Current (A), symmetric log", title=f"Reverse saturation Is={isat:.3e} A")
    for ax in axes:
        ax.set_xlabel("Applied voltage (V), forward bias positive")
        ax.grid(alpha=.25)
    save_figure(fig, out, "ideal_iv.png")
    np.savetxt(out / "ideal_iv.csv", np.column_stack((voltage, current)), delimiter=",",
               header="voltage_V,current_A", comments="")
    np.savetxt(out / "equilibrium_profiles.csv",
               np.column_stack((x, *(data[k] for k in ("rho", "field", "psi", "Ec", "Ei", "Ev", "Ef")))),
               delimiter=",", header="x_m,rho_C_m3,E_V_m,psi_V,Ec_eV,Ei_eV,Ev_eV,Ef_eV", comments="")
    (out / "junction_summary.json").write_text(
        json.dumps({**j.summary(), "Is_A": isat, "diode_area_m2": 1e-8,
                    "electron_lifetime_s": 1e-6, "hole_lifetime_s": 1e-6}, indent=2) + "\n")


def argument_parser(description: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--na", type=float, default=1e16, help="Acceptor density [cm^-3]")
    parser.add_argument("--nd", type=float, default=1e17, help="Donor density [cm^-3]")
    parser.add_argument("--temperature", type=float, default=300, help="Temperature [K]")
    parser.add_argument("--output", type=Path, default=ROOT / "figures")
    return parser


def main() -> None:
    parser = argument_parser("Analytical PN junction electrostatics, bands, and ideal I-V")
    args = parser.parse_args()
    try:
        j = Junction(args.na, args.nd, args.temperature)
    except ValueError as exc:
        parser.error(str(exc))
    plot_junction(j, args.output)
    print(json.dumps(j.summary(), indent=2))
    print(f"Saved electrostatics, bands, and ideal I-V to {args.output.resolve()}")


if __name__ == "__main__":
    main()
