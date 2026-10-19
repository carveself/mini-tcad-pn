# Physics validation

From the repository root, run:

```bash
python -m unittest discover -s tests -v
```

Seven tests cover charge neutrality, equilibrium boundary and band values,
ideal-diode limits, parameter trends, and symmetric/asymmetric junctions.
Independent quadratic and sinusoidal manufactured solutions check the
Poisson solver and its second-order convergence on smooth charge.

Only Python's standard unittest library and the project dependencies are
needed. The tests do not change any generated figures or data files.
