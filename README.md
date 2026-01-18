# GeiloWinterSchoolUQ
Coding examples for the Geilo Winter School 2026.

## Running on Sigma2/Olivia
See [Sigma2 setup guide](https://md.sigma2.no/s/1wxwWlLjG) for instructions on how to set up your environment on Sigma2. 

Once you have everything set up, JupyterLab session from apps.olivia.sigma2.no, open a new notebook and clone this repository:
```
!git clone https://github.com/kjetilly/GeiloWinterSchoolUQ.git
```
Navigate to the cloned directory and open `notebooks/basic_monte_carlo.ipynb` to get started.

## Setup locally
To run the code examples locally, you need to have Julia and Jupyter installed. To install Julia, follow the instructions on the [official Julia website](https://julialang.org/downloads/). 

To install Jupyter, you can use the following command in Julia's REPL:

```julia
using Pkg
Pkg.add("IJulia")
```
After installing Jupyter, you can clone this repository to your local machine:

```bash
git clone https://github.com/kjetilly/GeiloWinterSchoolUQ.git
```

To run the notebooks, navigate to the cloned directory and start JupyterLab:

```bash
julia -e 'using IJulia; notebook()'
```
