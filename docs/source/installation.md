# Installation Quick Start

## Through Pypi

AF-Analysis can install easily through Pypi:

```bash
$ pip install af_analysis
```

or you can install the last version directly from the github repo:

```bash
$ pip install git+https://github.com/samuelmurail/af_analysis.git@main
```

## Get sources from the GithubRepo

The sources for AF-Analysis can be downloaded from the GithubRepo.

You can either clone the public repository:

```bash
$ git clone git@github.com:samuelmurail/af_analysis.git
```

Or download the tarball:

```bash
$ curl -OJL https://github.com/samuelmurail/af_analysis/tarball/master
```

Once you have a copy of the source, switch to the `af_analysis` directory.

```bash
$ cd af_analysis
```

##  Install `af_analysis`

Once you have a copy of the source and have created a conda environment, you can install it with:

```bash
$ pip install .
```

## Test the installation

Use `pytest` to check that the installation was successful:

```bash
$ pip install pytest
$ pytest
============================ test session starts =============================
platform linux -- Python 3.10.13, pytest-7.4.3, pluggy-1.3.0
rootdir: /home/murail/Documents/Code/af_analysis
plugins: anyio-4.0.0
collected 9 items                                                            

src/af_analysis/test/test_analysis.py ..                              [ 22%]
src/af_analysis/test/test_clustering.py .                             [ 33%]
src/af_analysis/test/test_data.py ...                                 [ 66%]
src/af_analysis/test/test_docking.py .                                [ 77%]
src/af_analysis/test/test_format.py ..                                [100%]

============================== warnings summary ==============================
../../../miniforge3/envs/docking/lib/python3.10/site-packages/Bio/Application/__init__.py:40
  /home/murail/miniforge3/envs/docking/lib/python3.10/site-packages/Bio/Application/__init__.py:40: BiopythonDeprecationWarning: The Bio.Application modules and modules relying on it have been deprecated.
  
  Due to the on going maintenance burden of keeping command line application
  wrappers up to date, we have decided to deprecate and eventually remove these
  modules.
  
  We instead now recommend building your command line and invoking it directly
  with the subprocess module.
    warnings.warn(

src/af_analysis/test/test_clustering.py::test_cf_1_5_5_relax
src/af_analysis/test/test_clustering.py::test_cf_1_5_5_relax
src/af_analysis/test/test_clustering.py::test_cf_1_5_5_relax
src/af_analysis/test/test_clustering.py::test_cf_1_5_5_relax
  /home/murail/miniforge3/envs/docking/lib/python3.10/site-packages/MDAnalysis/coordinates/base.py:725: UserWarning: Reader has no dt information, set to 1.0 ps
    return self.ts.dt

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
======================= 9 passed, 5 warnings in 6.14s ========================
```