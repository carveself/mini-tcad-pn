# Source code

Run these files with Python from any working directory. Default outputs are
resolved relative to the project root. CLI doping units are cm^-3; internal
electrostatics use SI units, and band energies use eV.

- `pn_junction.py`: material parameters, analytical depletion model, bands,
  ideal diffusion diode, common CLI and plotting helpers.
- `parameter_sweep.py`: independent NA, ND and temperature sweeps.
- `poisson_numeric.py`: tridiagonal finite-difference solver, prescribed
  charge, analytical comparison and mesh-refinement study.

The numerical interior potential is obtained by linear-system elimination.
Analytical depletion edges and boundary potentials define the comparison
problem; this is not self-consistent carrier transport.
