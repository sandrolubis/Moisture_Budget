# Moisture Budget

**Sandro W. Lubis, Ph.D.**
Pacific Northwest National Laboratory (PNNL)

This repository provides **NCL and Python scripts for calculating atmospheric moisture-budget diagnostics**, including horizontal moisture advection, moisture convergence, moisture flux convergence (MFC), vertical moisture transport, and vertically integrated moisture-budget terms.

The scripts are designed for applications involving tropical variability, large-scale circulation, precipitation, atmospheric rivers, and other hydroclimate processes in which changes in atmospheric moisture transport are important.

<p align="center">
  <img src="input/moisture_budget.png" width="800">
</p>

---

## Overview

The atmospheric moisture budget describes how water vapor changes through transport, convergence, and non-advective moisture sources and sinks.

In pressure coordinates, specific humidity $q$ evolves according to contributions from:

* horizontal moisture advection,
* horizontal moisture convergence,
* vertical moisture transport,
* local moisture storage, and
* non-advective moisture sources and sinks associated with processes such as condensation, evaporation of condensate, and surface moisture exchange.

The scripts in this repository diagnose these individual terms at pressure levels and after vertical integration through the atmospheric column.

---

# 1. Moisture Flux Convergence

The horizontal moisture flux is

```math
\mathbf{F}_q = q\mathbf{v},
```

where $q$ is specific humidity and $\mathbf{v}=(u,v)$ is the horizontal wind.

The horizontal **moisture flux convergence (MFC)** is defined as

```math
\mathrm{MFC}
=
-\nabla_h \cdot (q\mathbf{v}).
```

Using the product rule, MFC can be decomposed into moisture-advection and wind-convergence components:

```math
\mathrm{MFC}
=
-\mathbf{v}\cdot\nabla_h q
-
q\nabla_h\cdot\mathbf{v}.
```

Therefore,

```math
\mathrm{MFC}
=
\mathrm{ADV}_q
+
\mathrm{CONV}_q,
```

where

```math
\mathrm{ADV}_q
=
-u\frac{\partial q}{\partial x}
-
v\frac{\partial q}{\partial y}
```

is the **horizontal moisture-advection term**, and

```math
\mathrm{CONV}_q
=
-q
\left(
\frac{\partial u}{\partial x}
+
\frac{\partial v}{\partial y}
\right)
```

is the contribution from **horizontal wind convergence**.

This decomposition is useful for distinguishing whether anomalous moisture convergence arises primarily from the transport of moisture by the circulation or from convergence of the flow itself.

---

# 2. Vertical Moisture Transport

In pressure coordinates, the vertical moisture-advection term is

```math
\mathrm{ADV}_{q,p}
=
-\omega\frac{\partial q}{\partial p},
```

where $\omega = dp/dt$ is pressure vertical velocity.

The corresponding vertical convergence contribution is

```math
\mathrm{CONV}_{q,p}
=
-q\frac{\partial\omega}{\partial p}.
```

Together,

```math
\mathrm{MFC}_{p}
=
-\omega\frac{\partial q}{\partial p}
-
q\frac{\partial\omega}{\partial p}
=
-\frac{\partial(q\omega)}{\partial p}.
```

These terms describe vertical redistribution of atmospheric moisture.

---

# 3. Local Moisture Budget

The local moisture tendency is

```math
\frac{\partial q}{\partial t}.
```

The scripts additionally diagnose the residual non-advective moisture source/sink as

```math
S_q
=
\frac{\partial q}{\partial t}
+
u\frac{\partial q}{\partial x}
+
v\frac{\partial q}{\partial y}
+
\omega\frac{\partial q}{\partial p}.
```

This term represents changes in specific humidity that are not explained by resolved three-dimensional advection and may include effects associated with condensation, evaporation, and other moist physical processes.

The sign convention should always be checked when comparing this residual with precipitation, evaporation, or model-physics tendencies.

---

# 4. Column-Integrated Moisture Budget

Column water vapor is defined as

```math
W
=
\frac{1}{g}
\int_{p_t}^{p_s}
q\,dp,
```

where $g$ is gravitational acceleration, $p_s$ is the lower pressure boundary, and $p_t$ is the upper pressure boundary.

The vertically integrated horizontal moisture transport is

```math
\mathbf{Q}
=
\frac{1}{g}
\int_{p_t}^{p_s}
q\mathbf{v}\,dp.
```

The column moisture budget can then be written as

```math
\frac{\partial W}{\partial t}
=
-\nabla_h\cdot\mathbf{Q}
+
E-P,
```

or equivalently,

```math
P-E
=
\mathrm{MFC}_{\mathrm{column}}
-
\frac{\partial W}{\partial t},
```

where

```math
\mathrm{MFC}_{\mathrm{column}}
=
-\nabla_h\cdot\mathbf{Q}.
```

Here,

* $P$ is precipitation,
* $E$ is surface evaporation,
* $W$ is column water vapor, and
* $\mathrm{MFC}_{\mathrm{column}}$ is vertically integrated moisture flux convergence.

Under approximately steady conditions,

```math
\frac{\partial W}{\partial t} \approx 0,
```

and therefore

```math
P-E
\approx
\mathrm{MFC}_{\mathrm{column}}.
```

Positive vertically integrated MFC therefore indicates net convergence of atmospheric moisture and generally favors positive $P-E$, whereas negative MFC indicates moisture divergence and drying.

For transient systems, however, the moisture-storage term should not be neglected.

---

# Repository Contents

```text
Moisture_Budget/
│
├── cal_mfc_850.ncl
├── cal_mfc_850.py
├── cal_moist_budget_multi_levels.ncl
├── cal_moist_budget_single_level.ncl
├── cal_vint_moist_budget.ncl
│
├── cross_products/
│   ├── advection_term.ncl
│   ├── cal_budget_decomp.ncl
│   ├── cal_filter_eddy.ncl
│   ├── divergence_term.ncl
│   ├── prep_data.sh
│   └── vort_flux_divergence.ncl
│
├── input/
│   ├── input_2001.nc
│   └── moisture_budget.png
│
├── README.md
└── LICENSE
```

---

# Main Scripts

## `cal_moist_budget_multi_levels.ncl`

Calculates the moisture budget at multiple atmospheric pressure levels.

The script diagnoses:

```text
dq_dt      local moisture tendency
adv_q      horizontal moisture advection
conv_q     horizontal moisture convergence
mfc        horizontal moisture flux convergence
mfc_ver    vertical moisture flux convergence
q2         diagnosed moisture source/sink
```

The horizontal MFC satisfies

```math
\mathrm{MFC}
=
\mathrm{ADV}_q
+
\mathrm{CONV}_q.
```

The script also explicitly calculates vertical moisture-advection and convergence terms.

---

## `cal_moist_budget_single_level.ncl`

Calculates the same moisture-budget decomposition but saves the diagnostics at a selected pressure level.

The current script contains a pressure-level selection that can be modified according to the desired analysis level.

Typical output variables include

```text
dq_dt
adv_q
conv_q
mfc
mfc_ver
q2
q
conv_q_zonal
conv_q_meridional
```

This script is useful when the analysis focuses on a particular part of the troposphere rather than the entire atmospheric column.

---

## `cal_vint_moist_budget.ncl`

Calculates **mass-weighted vertically integrated moisture-budget terms**.

The vertical integration follows the general form

```math
\left\langle A \right\rangle
=
\frac{1}{g}
\int A\,dp.
```

The script calculates vertically integrated:

```text
dq_dt                moisture tendency
adv_q                horizontal moisture advection
conv_q               horizontal moisture convergence
conv_q_zonal         zonal convergence contribution
conv_q_meridional    meridional convergence contribution
mfc                   moisture flux convergence
mfc_ver               vertical moisture flux convergence
adv_q_vertical        vertical moisture advection
conv_q_vertical       vertical convergence contribution
q2                    moisture source/sink
q                     column water vapor
```

The vertically integrated tendency and transport terms are written in

```text
kg m^-2 s^-1
```

while vertically integrated specific humidity gives column water vapor in

```text
kg m^-2
```

which is numerically equivalent to millimeters of liquid water.

---

## `cal_mfc_850.ncl`

A compact NCL script for calculating horizontal moisture flux convergence at **850 hPa**.

It evaluates

```math
\mathrm{MFC}
=
-\nabla_h\cdot(q\mathbf{v})
```

through its advection and convergence components:

```math
\mathrm{MFC}
=
-\mathbf{v}\cdot\nabla_hq
-
q\nabla_h\cdot\mathbf{v}.
```

The default input files are

```text
u850.nc
v850.nc
q850.nc
```

and the default output is

```text
mfc_850.nc
```

---

## `cal_mfc_850.py`

Python implementation of the 850-hPa moisture flux convergence calculation.

The script uses **xarray** and **NumPy** and directly evaluates

```math
\mathrm{MFC}
=
-
\left[
\frac{\partial(qu)}{\partial x}
+
\frac{\partial(qv)}{\partial y}
\right]
```

using centered finite differences on a regular latitude-longitude grid.

Default input:

```text
u850.nc
v850.nc
q850.nc
```

Default output:

```text
mfc_850.nc
```

---

# Input Data

The NCL moisture-budget scripts are written for data with dimensions approximately organized as

```text
[time, pressure, latitude, longitude]
```

The current implementations use ERA5-style variable names:

```text
var131    zonal wind, u
var132    meridional wind, v
var133    specific humidity, q
var135    pressure vertical velocity, omega
```

Typical physical units are:

```text
u, v       m s^-1
q          kg kg^-1
omega      Pa s^-1
pressure   Pa or hPa, depending on the input dataset
```

Before running the scripts, users should verify:

* variable names;
* pressure units;
* pressure-level ordering;
* latitude ordering;
* longitude convention;
* temporal resolution;
* input/output directories; and
* coordinate names.

The scripts can be readily modified for datasets that use different variable names or file structures.

---

# Running the Scripts

## NCL

For the multilevel moisture budget:

```bash
ncl cal_moist_budget_multi_levels.ncl
```

For a selected pressure level:

```bash
ncl cal_moist_budget_single_level.ncl
```

For the vertically integrated moisture budget:

```bash
ncl cal_vint_moist_budget.ncl
```

For 850-hPa MFC:

```bash
ncl cal_mfc_850.ncl
```

Input paths, years, pressure levels, and variable names should be modified directly in the scripts as needed.

---

## Python

The Python MFC calculation requires

```text
xarray
numpy
```

Run with

```bash
python cal_mfc_850.py
```

---

# Sign Convention

In this repository,

```math
\mathrm{MFC}
=
-\nabla\cdot(q\mathbf{v}).
```

Therefore:

```text
MFC > 0    moisture flux convergence
MFC < 0    moisture flux divergence
```

Positive MFC indicates convergence of water vapor into a region, whereas negative MFC indicates export of water vapor.

When comparing MFC with $P-E$, remember that the atmospheric moisture-storage term can become important on daily and other transient timescales.

---

# Example Scientific Application

The moisture-budget framework used in this repository was applied to investigate the propagation of the **Madden–Julian Oscillation (MJO)** across the Maritime Continent and the role of cross-equatorial surges in modifying its moisture distribution and southward propagation.

Please cite the following paper when using this code:

**Lubis, S. W., Hagos, S., Chang, C.-C., Balaguru, K., & Leung, L. R. (2023).**
*Cross-equatorial surges boost MJO's southward detour over the Maritime Continent.*
**Geophysical Research Letters, 50**, e2023GL104770.
doi:10.1029/2023GL104770

---

# Citation

If you use these scripts in your research, please cite:

> Lubis, S. W., Hagos, S., Chang, C.-C., Balaguru, K., & Leung, L. R. (2023). Cross-equatorial surges boost MJO's southward detour over the Maritime Continent. *Geophysical Research Letters, 50*, e2023GL104770. doi:10.1029/2023GL104770

### BibTeX

```bibtex
@article{Lubis2023MJO,
  author  = {Lubis, Sandro W. and Hagos, Samson and Chang, Chuan-Chieh and Balaguru, Karthik and Leung, L. Ruby},
  title   = {Cross-Equatorial Surges Boost MJO's Southward Detour over the Maritime Continent},
  journal = {Geophysical Research Letters},
  volume  = {50},
  pages   = {e2023GL104770},
  year    = {2023},
  doi     = {10.1029/2023GL104770}
}
```

---

# Notes

These scripts are intended primarily as research tools. The numerical derivatives and vertical integrations depend on the structure and resolution of the input data, so users should verify budget closure and sign conventions for their particular dataset.

Particular care should be taken when applying the calculations near:

* the poles,
* topography,
* the lower atmospheric boundary,
* missing pressure levels, or
* irregular grids.

For exact moisture-budget closure, surface pressure, boundary moisture fluxes, precipitation, evaporation, and storage should be treated consistently with the dataset being analyzed.

---

# Author and Contact

**Sandro W. Lubis, Ph.D.**
Pacific Northwest National Laboratory (PNNL)

Contact: **[sandro.lubis@pnnl.gov](mailto:sandro.lubis@pnnl.gov)**

---

# License

This repository is distributed under the **MIT License**.
