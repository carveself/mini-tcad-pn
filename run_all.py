"""Run all three demonstrations, propagating any failure to the caller."""
from pathlib import Path
import subprocess
import sys


def main() -> None:
    root = Path(__file__).resolve().parent
    for script in ("pn_junction.py", "parameter_sweep.py", "poisson_numeric.py"):
        subprocess.run([sys.executable, str(root / "src" / script), *sys.argv[1:]], check=True)
    print("All demonstrations completed successfully.")


if __name__ == "__main__":
    main()
