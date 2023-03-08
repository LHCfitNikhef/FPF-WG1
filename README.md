# FPF-WG1
Working Group 1 "Neutrino interactions and hadron structure" of the Forward Physics Facility. The following branch is specific to NNPDF-related fits.

## $\nu\rm{FPF}$

[nufpf](https://github.com/juanrojochacon/FPF-WG1/tree/main/nufpf) provides the computation of the differential cross section. It is a wrapper around [Yadism](https://github.com/NNPDF/yadism) and [PineAPPL](https://github.com/NNPDF/pineappl). The workflow is the following: Yadism takes care of computing the coefficient functions and dump them as PineAPPL grids, PineAPPL then perform the convolution with the Parton Distribution Function (PDF).

### Dependencies

The following are the packages required to run $\nu\rm{FPF}$:
- [poetry](https://python-poetry.org/): handles the packaging and the dependency management
- [lhapdf](https://lhapdf.hepforge.org/): evaluates the parton density functions to be convoluted with the PineAPPL grids

### Install the package

To install the package, first clone the repository:
```bash
git clone https://github.com/juanrojochacon/FPF-WG1 --depth 1
```
then just run the following command:
```bash
poetry install
```
The [nufpf](https://github.com/juanrojochacon/FPF-WG1/tree/main/nufpf) package provides various command line interface, to know more about the various features just type:
```bash
poetry run nufpf --help
```

### Generating predictions

[nufpf](https://github.com/juanrojochacon/FPF-WG1/tree/main/nufpf) can both generate as predictions the differential cross sections (as a `.txt` file) or the structure functions (as a `LHAPDF` grid). In order to generate the predictions first we need to generate a Yadism theory card containing the input settings and the observables to be computed. To do so, just input in the following command:
```bash
poetry run nufpf xsecs runcards igrids/xyQ2.csv [-a ${A_VALUE}] [--obs ${OBS_TYPE}]
```
The `${OBS_TYPE}` is an optional argument, its value can only be either `XSEC` (for differential cross sections) or `SF` (for structure functions). If it is not specified, by default the differential cross sections will be computed. In the same way, if `${A_VALUE}` is not specified, the proton `A=1` will be set as the default.

Otherwise specified this will generate a folder called `theory` in the current directory which contains the theory run cards. These theory cards can then be used to generate the PineAPPL grids using the following command:
```bash
poetry run nufpf xsecs grids ${Path_to_Cards}
```
This will dump the grids into the same folder as above which can then be convoluted with a PDF set in order to obtain the final predictions. From this point, we can choose to produce directly the values of the differential cross sections or to have access at the individual structure functions.

* **Computing differential cross sections:**
To compute the differential cross sections and dump the results into a `.txt` just run the following command:
```bash
poetry run nufpf xsecs generate_xsecs_datfile ${Path_to_Grids} igrids/xyQ2.csv ${PDFSET_NAME}
```
This will generate a file called `diff_xsecs_a${A_VALUE}.txt` which otherwise specificied is stored in the current path.

* **Computing structure functions:**
Instead, if one wants to compute structure functions, just type the following command:
```bash
poetry run nufpf xsecs generate_sfs_lhapdf ${Path_to_Grids} igrids/xyQ2.csv ${PDFSET_NAME}
```
The resulting object will be a grid following the LHAPDF format in which the entries are the structure functions:
$$F_2^{\nu N}, F_L^{\nu N}, xF_3^{\nu N}, F_2^{\bar{\nu} N}, F_L^{\bar{\nu} N}, xF_3^{\bar{\nu} N}, \langle F_2^N \rangle, \langle F_L^N \rangle, \langle xF_3^N \rangle$$
The next step is then to multiply the structure functions with the corresponding coefficients in order to compute the cross section:
```bash
poetry run nufpf xsecs multiply_sfs ${LHAPDF_SETNAME}
```
For convenience, pre-computed LHAPDF sets for free-proton, iron, and lead target are available in
the following [location](https://data.nnpdf.science/FPF-DIS/).

* **Computing total cross sections:**

Currently, the total cross section is computed by reading the structure functions from the
LHAPDF set and follows the procedure described in this [paper](https://arxiv.org/pdf/1808.02034.pdf).
For instance, to compute the total cross section for the interaction of a *neutrino* to
a free-proton target with an energy $E_\nu = 5 \cdot 10^3~\mathrm{GeV}$, one can simply run the
following command:
```bash
nufpf extra integrate YADISM_SFS_A1 -a 1 -e 5e3 --type neutrino
```
