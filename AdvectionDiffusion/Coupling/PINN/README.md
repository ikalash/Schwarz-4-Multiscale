
# Overview

This package supplies utilities for training physics-informed neural networks (PINNs) coupled by the Schwarz alternating method for the 1D advection-diffusion equation. The code was originally developed by Will Snyder of Virginia Tech for research conducted at Sandia National Laboratories under the direction of Irina Tezaur and Chris Wentland. Refer to the CSRI proceedings paper draft at `docs/Will_Snyder_CSRI_Summer_Proceedings_2023_Paper.pdf` for theory and implementation details.

# Prerequisites and testing

This package requires `numpy`, `matplotlib`, `pandas`, and `tensorflow`. You can install these yourself with Conda, or allow the `pip` installation below to download them (include the `--user` flag if using the system Python installation, which is not recommended). The package can be downloaded and installed with the following commands:

```
git clone git@github.com:ikalash/Schwarz-4-Multiscale.git
cd Schwarz-4-Multiscale/AdvectionDiffusion/Coupling/PINN/
pip install -e .
```
# Add yaml input file

This package thus also requires yaml.

After installation has completed, test the installation by executing the following command:

```
python pinnschwarz/driver.py cases/example/input.csv cases/example/
```
