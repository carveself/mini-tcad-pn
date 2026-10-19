"""Separate NA, ND, and temperature scans with reproducible figures and CSV."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from pn_junction import (Junction, analytical, argument_parser, ideal_iv,
                         plt, save_figure)


def run_sweeps(base: Junction, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    doping = np.geomspace(1e15, 1e18, 41)
    temperatures = np.linspace(250, 400, 41)
    groups = [
        ("NA", doping, [Junction(n, base.nd_cm3, base.temperature) for n in doping]),
        ("ND", doping, [Junction(base.na_cm3, n, base.temperature) for n in doping]),
        ("T", temperatures, [Junction(base.na_cm3, base.nd_cm3, t) for t in temperatures])]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    rows = []
    for column, (name, values, junctions) in enumerate(groups):
        axes[0, column].plot(values, [j.vbi for j in junctions])
        axes[1, column].plot(values, [j.width * 1e6 for j in junctions])
        for ax in axes[:, column]:
            if name != "T":
                ax.set_xscale("log")
            ax.set_xlabel(f"{name} ({'K' if name == 'T' else 'cm^-3'})")
            ax.grid(alpha=.25)
        axes[0, column].set_ylabel("Built-in potential (V)")
        axes[1, column].set_ylabel("Depletion width (um)")
        for value, j in zip(values, junctions):
            rows.append({"sweep": name, "sweep_value": value, **j.summary()})
    fig.suptitle("One parameter changes at a time; other values follow the CLI inputs")
    save_figure(fig, out, "parameter_sweep.png")
    with (out / "parameter_sweep.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    # Common physical x axis allows depletion width and band bending comparison.
    temperatures_for_bands = (250, 300, 400)
    band_junctions = [Junction(base.na_cm3, base.nd_cm3, t) for t in temperatures_for_bands]
    extent = max(j.width for j in band_junctions)
    x = np.linspace(-1.3 * extent, 1.3 * extent, 2001)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    voltage = np.linspace(-.3, .55, 700)
    for j in band_junctions:
        data = analytical(j, x)
        line, = axes[0].plot(x * 1e6, data["Ec"], label=f"Ec, {j.temperature:g} K")
        axes[0].plot(x * 1e6, data["Ev"], "--", color=line.get_color(), alpha=.75)
        current, _ = ideal_iv(j, voltage)
        axes[1].plot(voltage, current, label=f"{j.temperature:g} K")
    axes[0].axhline(0, color="gray", ls=":", label="EF = 0 eV; Ev dashed")
    axes[0].set(xlabel="Position x (um)", ylabel="Energy relative to EF (eV)",
                title="Band edges at fixed doping")
    axes[1].set_yscale("symlog", linthresh=1e-12)
    axes[1].set(xlabel="Applied voltage (V)", ylabel="Ideal current (A)",
                title="Temperature scan; fixed mobility and lifetime")
    for ax in axes:
        ax.legend(fontsize=9)
        ax.grid(alpha=.25)
    save_figure(fig, out, "temperature_bands_iv.png")
    print(f"Saved {len(rows)} scan cases and temperature bands/I-V to {out.resolve()}")


def main() -> None:
    parser = argument_parser("Scan NA, ND, T; plot built-in voltage, width, bands and I-V")
    args = parser.parse_args()
    try:
        run_sweeps(Junction(args.na, args.nd, args.temperature), args.output)
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
