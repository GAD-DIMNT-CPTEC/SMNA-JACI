project_dir: ./src
output_dir: ./html
project: BAM Pre-Processor DOCUMENTATION
include: 
project_github: 
project_website: http://www.cptec.inpe.br/
summary: <img src= "../images/logo-cptec.png" width="5%"> The Brazilian Global Atmospheric Model (BAM) Pre-processor<img src= "../images/logo-inpe.png" width="7%">
author: The BAM's Team
author_description: The Brazilian Atmospheric Model (BAM) was fully developed at CPTEC of the National Institute for Space Research (INPE) and was recently implemented at this center for operational purposes and is used to make weather and climate forecasts
github:
email: atende@inpe.br
fpp_extensions: fpp
predocmark: >
media_dir: ./media
docmark_alt: #
predocmark_alt: <
display: public
         protected
         private
source: false
graph: true
search: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
license: CC-GPL
extra_filetypes: sh #

The BAM model is a global scale model, where the primitive equations are discretized using the spectral method and the hydrostatic approximation. It has two formulations of dynamic cores (Eulerian and Semi-Lagrangian), can be easily configured in any of the modules, it is also possible to work with the option of quadratic and triangular grid, full or reduced regular Gaussian grid. The primitive equations of BAM are divided into two parts: spectral space, where the linear processes of the primitive equations (dynamics) are solved, and physical space, where the nonlinear processes (physical parameterizations) are solved. The two vertical coordinate options used for discretization are the sigma and hybrid coordinates. The system for rotating the model consists of three parts: PRE, MODEL and POST. 

This documentation is relative to the PRE (pre-processing) part, which creates the initial conditions (examples of atmospheric conditions: pressure, wind, temperature and humidity) and the lower boundary conditions (example: topography, TSM sea surface temperature, etc.). 

The MODEL part, runs the global model and the POST part moves from spectral space to grid point, and from sigma / hybrid vertical coordinate to pressure coordinate.


@info
**1. Introduction**<br/>
This page will guide You to install software infrastructure. BAM code is distributed under GLP-3.0 license (https://opensource.org/licenses/GPL-3.0). BAMworks with Linux or Unix operational’s systems. As a first approach, we recommend the Linux UBUNTU flavours distribution: UBUNTU <https://www.ubuntu.com/download/desktop?><br/>
Please, read carefully each of the steps below in order to install and run the model.
---
**2. Install Fortran Compilers** <br/>
PRE has been tested with the compilers: PGI compilers <https://www.pgroup.com/products/community.htm> and GNU Fortran compiler (GPL) <https://gcc.gnu.org/fortran/> or <https://askubuntu.com/questions/358907/how-do-i-install-gfortran>. Follow the instructions of each site to install the compilers.
---
**3. Install MPI Libraries and Software** <br/>
PRE works in parallel mode. Run the model using MPI with the MPIRUN command. We recommend download and install the last version of MPICH stable release. Take care to choose the correct version to your OS.

After download the MPICH in <https://www.mpich.org/downloads/>, please, proceed the installation:

* Uncompress and goto mpich directory:
```bash 
$ tar -zxvf mpich-3.1.4.tar.gz
$ cd mpich-3.1.4
```
* Configure mpichs makefile, make and install:
```bash
$ ./configure -disable-fast CFLAGS=-O2 FFLAGS=-O2 CXXFLAGS=-O2 FCFLAGS=-O2 -prefix=/opt/mpich3 CC=gcc FC=gfortran F77=gfortran
$ make; sudo make install
```

Notes:

prefix directory where mpich will be installed when You make the install command;
CC= C compiler according You install on step (2) above;
FC= Fortran compiler according You install on step (2) above;
F77= Use the same as FC
---
**4. Compile ** <br/>
Check instructions on README.txt file in PRE root directory <br/>
<br/>

@endinfo
    @Bug
If You find any bugs please, send an email to:

* <mailto:atende.cptec@inpe.br>

@endbug

