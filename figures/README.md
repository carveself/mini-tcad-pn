# Generated figures and data

Run `python run_all.py` from the repository root to regenerate the seven PNG
figures, four CSV tables, and two JSON summaries in this directory.

All checked-in results use NA=1e16 cm^-3, ND=1e17 cm^-3 and T=300 K as the
baseline. The sweep script varies one parameter at a time as described in the
main README. The default Poisson comparison uses 1601 grid nodes.

CSV headers specify SI units for electrostatics and eV for energies. The
Poisson solution uses prescribed depletion charge; it is not self-consistent
carrier transport. Regenerating with custom CLI inputs overwrites the outputs
unless you specify a different `--output` directory.
