# Moisture Budget

**Sandro W. Lubis, Ph.D.**
Pacific Northwest National Laboratory (PNNL)

This repository provides **NCL and Python scripts for calculating atmospheric moisture-budget diagnostics**, including horizontal moisture advection, horizontal moisture convergence, moisture flux convergence (MFC), vertical moisture transport, and vertically integrated moisture-budget terms.

These diagnostics are useful for studying tropical variability, atmospheric convection, precipitation, atmospheric rivers, large-scale circulation, and other hydroclimate processes in which moisture transport plays an important role.

<p align="center">
  <img src="input/moisture_budget.png" width="800">
</p>

---

## Citation
If you use or adapt this code, please cite:

**Lubis, S. W., Hagos, S., Chang, C.-C., Balaguru, K., & Leung, L. R. (2023).**
*Cross-equatorial surges boost MJO's southward detour over the Maritime Continent.*
**Geophysical Research Letters, 50**, e2023GL104770.
https://doi.org/10.1029/2023GL104770


## Overview

The atmospheric moisture budget describes the temporal evolution of water vapor through transport, convergence, storage, and physical moisture sources and sinks.

In pressure coordinates, changes in specific humidity $q$ arise from:

* horizontal moisture advection,
* horizontal wind convergence,
* vertical moisture transport,
* local moisture storage, and
* moisture sources and sinks associated with processes such as condensation, evaporation of condensate, and surface moisture exchange.

The scripts in this repository diagnose these processes at individual pressure levels and after vertical integration through the atmospheric column.

Throughout this README,

```math
\mathbf{V}=(u,v)
```

denotes the **horizontal wind vector**, where $u$ and $v$ are the zonal and meridional wind components, respectively.

---

# 1. Moisture Flux Convergence

The horizontal moisture flux is

```math
\mathbf{F}_q = q\mathbf{V},
```

where $q$ is specific humidity and $\mathbf{V}=(u,v)$ is the horizontal wind vector.

The horizontal **moisture flux convergence (MFC)** is defined as

```math
\mathrm{MFC}
=
-\nabla_h\cdot(q\mathbf{V}).
```

Using the product rule, MFC can be decomposed into moisture-advection and wind-convergence components:

```math
\mathrm{MFC}
=
-\mathbf{V}\cdot\nabla_h q
-
q\nabla_h\cdot\mathbf{V}.
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

This decomposition is useful for determining whether anomalous moisture flux convergence is primarily associated with the transport of moisture by the circulation or with convergence of the horizontal wind field.

---

# 2. Vertical Moisture Transport

In pressure coordinates, vertical moisture transport can similarly be separated into advection and convergence components.

The vertical moisture-advection term is

```math
\mathrm{ADV}_{q,p}
=
-\omega\frac{\partial q}{\partial p},
```

where

```math
\omega=\frac{Dp}{Dt}
```

is pressure vertical velocity.

The vertical moisture-convergence contribution is

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
q\frac{\partial\omega}{\partial p},
```

which can equivalently be written in flux form as

```math
\mathrm{MFC}_{p}
=
-\frac{\partial(q\omega)}{\partial p}.
```

These terms describe the redistribution of atmospheric moisture by vertical motion.

---

# 3. Local Moisture Budget

The local moisture tendency is

```math
\frac{\partial q}{\partial t}.
```

In pressure coordinates, the three-dimensional moisture equation can be written schematically as

```math
\frac{\partial q}{\partial t}
=
-\mathbf{V}\cdot\nabla_h q
-
\omega\frac{\partial q}{\partial p}
+
S_q,
```

where $S_q$ represents non-advective moisture sources and sinks.

Rearranging gives

```math
S_q
=
\frac{\partial q}{\partial t}
+
\mathbf{V}\cdot\nabla_h q
+
\omega\frac{\partial q}{\partial p}.
```

In component form,

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

This residual represents changes in atmospheric moisture not explained by resolved three-dimensional advection and may include processes associated with condensation, evaporation, precipitation formation, and other moist physical processes.

The sign convention should be checked carefully when comparing this residual with precipitation, evaporation, or model-physics tendencies.

---

# 4. Column-Integrated Moisture Budget

The atmospheric column water vapor is defined as

```math
W
=
\frac{1}{g}
\int_{p_t}^{p_s}
q\,dp,
```

where

* $g$ is gravitational acceleration,
* $p_s$ is the lower pressure boundary, and
* $p_t$ is the upper pressure boundary.

The vertically integrated horizontal moisture transport is

```math
\mathbf{Q}
=
\frac{1}{g}
\int_{p_t}^{p_s}
q\mathbf{V}\,dp.
```

The vertically integrated moisture flux convergence is then

```math
\mathrm{MFC}_{\mathrm{column}}
=
-\nabla_h\cdot\mathbf{Q}.
```

The column moisture budget can be written as

```math
\frac{\partial W}{\partial t}
=
\mathrm{MFC}_{\mathrm{column}}
+
E-P,
```

where

* $P$ is precipitation,
* $E$ is surface evaporation,
* $W$ is column water vapor, and
* $\mathrm{MFC}_{\mathrm{column}}$ is vertically integrated moisture flux convergence.

Equivalently,

```math
P-E
=
\mathrm{MFC}_{\mathrm{column}}
-
\frac{\partial W}{\partial t}.
```

Thus, precipitation minus evaporation depends on both moisture flux convergence and changes in atmospheric moisture storage.

---

## Steady-State Approximation

If changes in atmospheric moisture storage are small,

```math
\frac{\partial W}{\partial t}
\approx 0,
```

then the column moisture budget approximately reduces to

```math
P-E
\approx
\mathrm{MFC}_{\mathrm{column}}.
```

Under this approximation:

* **positive MFC** indicates net convergence of atmospheric moisture and favors positive $P-E$;
* **negative MFC** indicates moisture divergence and favors drying or reduced net precipitation.

For transient weather systems and daily or subseasonal variability, however, the storage term can be important and should not automatically be neglected.

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

# Input Data

The moisture-budget scripts are designed for atmospheric data organized approximately as

```text
[time, pressure, latitude, longitude]
```

The current NCL implementations use ERA5-style variable names such as

```text
var131    zonal wind, u
var132    meridional wind, v
var133    specific humidity, q
var135    pressure vertical velocity, omega
```

Typical physical units are

```text
u, v       m s^-1
q          kg kg^-1
omega      Pa s^-1
pressure   Pa or hPa, depending on the dataset
```

Before running the scripts, users should verify:

* input and output directories;
* variable names;
* coordinate names;
* pressure units;
* pressure-level ordering;
* latitude ordering;
* longitude convention;
* temporal resolution; and
* missing-value treatment.

The scripts can be modified for datasets that use different variable names or file structures.

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

For 850-hPa moisture flux convergence:

```bash
ncl cal_mfc_850.ncl
```

Input paths, output paths, pressure levels, years, and variable names can be modified directly in the scripts.

---

## Python

The Python implementation requires

```text
xarray
numpy
```

Run with

```bash
python cal_mfc_850.py
```

---

