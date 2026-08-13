#!/bin/bash
#---------------------------------------------------------------------------------------------------#
#                           Compilation Environment Setup Script - Version 1.3                      #
#---------------------------------------------------------------------------------------------------#
#BOP                                                                                                #
# Description:                                                                                      #
#     This script configures the compilation environment based on the specified machine and         #
#     compiler. It loads necessary modules, sets environment variables, and verifies dependencies.  #
#                                                                                                   #
# Usage:                                                                                            #
#     source script.sh [Machine (egeon or xc50)] [Compiler (intel or gnu)]                          #
#                                                                                                   #
# Examples:                                                                                         #
#     source script.sh egeon intel                                                                  #
#     source script.sh xc50 gnu                                                                     #
#                                                                                                   #
# Input Parameters:                                                                                 #
#     * Machine:  "egeon" or "xc50" (required)                                                      #
#     * Compiler: "intel" or "gnu" (required)                                                       #
#                                                                                                   #
# Dependencies:                                                                                     #
#     * The script requires the module system and the specified compiler environment.               #
#     * Required tools: `cmake`, `make`, `git`, `gcc`, `gfortran`, `mpif90`.                        #
#                                                                                                   #
# Notes:                                                                                            #
#     * If paths do not exist, warnings will be displayed, and the user must ensure their presence. #
#     * This script **must be sourced** (`source script.sh`) instead of executed directly.          #
#                                                                                                   #
# Revisions:                                                                                        #
#     * 25-06-2024: de Mattos, J. G. - Initial script version                                       #
#     * 26-06-2024: de Mattos, J. G. - Added environment variable checks                            #
#     * 27-06-2024: de Mattos, J. G. - Improved module loading and compiler validation              #
#     * 21-02-2025: de Mattos, J. G. - Enhanced error handling and logging                          #
#     * 21-02-2025: de Mattos, J. G. - Added missing paths validation and dependency checks         #
#EOP                                                                                                #
#---------------------------------------------------------------------------------------------------#

# Input variables
machine=$1
compiler=$2

# Validate the number of input arguments
if [ $# -ne 2 ]; then
    echo ""
    echo "Usage: source $0 [Machine (egeon or xc50)] [Compiler (intel or gnu)]"
    echo "Example: source $0 egeon intel"
    echo "Example: source $0 xc50 gnu"
    echo ""
    return 1
fi

# Print configuration information
echo "[INFO] Configuring compilation environment with:"
echo "[INFO] Machine: $machine"
echo "[INFO] Compiler: $compiler"

module_swap(){
   local atual=$1
   local novo=$2

   if module -t list 2>&1 | grep -q "$atual"; then
       module unload $atual
   fi
   module load $novo

}
#---------------------------------------------------------------------------------------------------#
# Function: assign                                                                                  #
# Description: Sets an environment variable and verifies its existence if it is a directory.        #
#---------------------------------------------------------------------------------------------------#
assign() {
    export "$1=$2"

    # Check if the path exists if it's a directory
    if [ -d "$2" ]; then
        echo "[OK] Path $2 defined by $1 exists."
    else
        echo "[WARNING] Path $2 defined by $1 does not exist."
    fi
}

#---------------------------------------------------------------------------------------------------#
# Function: check_compilers                                                                         #
# Description: Verifies that essential compilers and MPI tools are available.                       #
#---------------------------------------------------------------------------------------------------#
check_compilers() {
    local tools=("$FC" "$F90" "$CC" "$CXX")
    local missing=()

    for tool in "${tools[@]}"; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            missing+=("$tool")
        fi
    done

    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "[ERROR] The following compilers or MPI tools are missing: ${missing[*]}"
        echo "[ERROR] Ensure that the correct module environment is loaded."
        return 1
    fi
}

#---------------------------------------------------------------------------------------------------#
# Machine-specific Configuration                                                                   #
#---------------------------------------------------------------------------------------------------#
echo "[INFO] Loading module environment..."

if [ "${machine,,}" == "egeon" ]; then
    export LC_ALL="en_US.UTF-8"
    module -q purge
    
    # Configuration for Intel compiler
    if [ "${compiler,,}" == "intel" ]; then
        module load intel impi
        export FC=ifort
        export F90=ifort
        export CC=icc
        export CXX=icx
        
    # Configuration for GNU compiler
    elif [ "${compiler,,}" == "gnu" ]; then
        module load gnu9/9.4.0 mpich
        export FC=gfortran
        export F90=gfortran
        export CC=gcc
        export CXX=g++
        
    # If the compiler is not recognized
    else
        echo "[ERROR] Compiler not recognized, use 'intel' or 'gnu'."
        return 1
    fi

    # Loading additional necessary modules
    module load curl-7.85.0-gcc-9.4.0-qbney7y
    module load cmake/3.21.3
    module load openblas
    module load netcdf
    module load netcdf-fortran

elif [ "${machine,,}" == "jaci" ]; then
    export LC_ALL="en_US.UTF-8"
    #module reset
    
    # 1. Carrega o ambiente de programação (PrgEnv) escolhido
    if [ "${compiler,,}" == "intel" ]; then
	module swap PrgEnv-cray PrgEnv-intel
        module load PrgEnv-intel
        
       # -assume byterecl: garante alinhamento de bytes para registros de E/S
       export FFLAGS="-assume byterecl ${FFLAGS}"
       export FCFLAGS="-assume byterecl ${FCFLAGS}" 
        
    elif [ "${compiler,,}" == "gnu" ]; then
        module load PrgEnv-gnu
        
        # Incompatibilidades conhecidas de convenção de nomes no GCC/Gfortran com GSI
        export FFLAGS="-fno-second-underscore ${FFLAGS}"
        export FCFLAGS="-fno-second-underscore ${FCFLAGS}"
        
    else
        echo "[ERROR] Compiler not recognized, use 'intel' or 'gnu'."
        return 1
    fi

    # 2. Define o alvo de arquitetura da CPU (AMD EPYC 9745)
    module load craype-x86-turin

    # 2. Carrega as bibliotecas do ecossistema Cray e dependências extras
    # Nota: cray-mpich e cray-libsci (BLAS/LAPACK) são vinculados automaticamente via wrappers ftn/cc
    module load cray-netcdf
    module load cray-hdf5
    #module load cmake

    # Mapeie para as variáveis que o CMake do GSI costuma procurar
    export NetCDF_Fortran_DIR="${CRAY_NETCDF_DIR}"
    export NetCDF_C_DIR="${CRAY_NETCDF_DIR}"
    export NETCDF_FORTRAN_ROOT="${CRAY_NETCDF_DIR}"
    export NETCDF_ROOT="${CRAY_NETCDF_DIR}"

    # SOLUÇÃO AQUI: Aponta a variável que o GSI pede para a pasta do NetCDF da Cray -> JACI
    export NETCDF_FORTRAN_DIR="${CRAY_NETCDF_DIR}"

    # 3. Define os wrappers da Cray como compiladores padrão
    export FC=ftn
    export F90=ftn
    export F77=ftn
    export CC=cc
    export CXX=CC

    # 4. Exporta caminhos explícitos das bibliotecas (se exigido pelos CMakeLists/Makefiles do GSI/BAM)
    if [ -n "${CRAY_NETCDF_DIR}" ]; then
        export NETCDF="${CRAY_NETCDF_DIR}"
        export NETCDF_ROOT="${CRAY_NETCDF_DIR}"
    fi
    if [ -n "$HDF5_DIR" ]; then
        export HDF5="$HDF5_DIR"
        export HDF5_ROOT="$HDF5_DIR"
    fi
    
elif [ "${machine,,}" == "jaci2" ]; then
    export LC_ALL="en_US.UTF-8"
    #module -q purge
    
    # Configuration for Intel compiler
    if [ "${compiler,,}" == "intel" ]; then
	module list
	module swap  PrgEnv-cray/8.6.0 PrgEnv-intel/8.6.0
	module load cray-libpals/1.6.1 cray-pals/1.6.1
        module load intel impi
        export FC=ifort
        export F90=ifort
        export CC=icc
        export CXX=icx
        
    # Configuration for GNU compiler
    elif [ "${compiler,,}" == "gnu" ]; then
        module load gnu9/9.4.0 mpich
        export FC=gfortran
        export F90=gfortran
        export CC=gcc
        export CXX=g++
        
    # If the compiler is not recognized
    else
        echo "[ERROR] Compiler not recognized, use 'intel' or 'gnu'."
        return 1
    fi

    # Loading additional necessary modules
    module load curl-7.85.0-gcc-9.4.0-qbney7y
    module load cmake/3.21.3
    module load openblas
    module load netcdf
    module load netcdf-fortran

elif [ "${machine,,}" == "xc50" ]; then
    . /opt/modules/default/etc/modules.sh
    module load pbs
    module load craype-x86-skylake
    module load craype-network-aries
    module_swap PrgEnv-cray PrgEnv-gnu
    module load cray-hdf5
    module load cray-netcdf

    # Setting compiler variables
    export FC=ftn
    export F90=ftn
    export CC=cc
    export CXX=g++

# If the machine is not recognized
else
    echo "[ERROR] Unsupported machine: $machine"
    return 1
fi

# Check if compilers are correctly loaded
check_compilers || return 1

# List loaded modules safely
module avail > /dev/null 2>&1 && module list || echo "[WARNING] Unable to list modules."

SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
RootDir="$(dirname "$SCRIPT_PATH")"

# Assign necessary paths
if [ "${machine,,}" == "jaci" ]; then
    assign DIRGSI "/p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717/SMG/cptec/gsi"
else 
    assign DIRGSI "$(pwd)"
fi

assign DIRLIB "${DIRGSI}/libsrc"
assign install_dir "${DIRGSI}"

# Additional environment paths
declare -A paths=(
   ["NETCDF"]="$NETCDF_DIR"
   ["NetCDF_DIR"]="$NETCDF_DIR"
   ["NetCDF_C_DIR"]="$NETCDF_DIR"
   ["NETCDF_F90"]="$NETCDF_FORTRAN_DIR"
   ["NetCDF_Fortran_DIR"]="$NETCDF_FORTRAN_DIR"
   ["DIRTEST"]="${DIRLIB}/NCEPLIBS-bufr-develop/test"
   ["bacio_DIR"]="${install_dir}/lib64/cmake/bacio"
   ["sigioBAMMod_DIR"]="${install_dir}/lib64/cmake/sigioBAMMod"
   ["sigio_DIR"]="${install_dir}/lib64/cmake/sigio"
   ["sfcio_DIR"]="${install_dir}/lib64/cmake/sfcio"
   ["ncio_DIR"]="${install_dir}/lib64/cmake/ncio"
   ["bufr_DIR"]="${install_dir}/lib64/cmake/bufr"
   ["w3emc_DIR"]="${install_dir}/lib64/cmake/w3emc"
   ["nemsio_DIR"]="${install_dir}/lib64/cmake/nemsio"
   ["sp_DIR"]="${install_dir}/lib64/cmake/sp"
   ["ip_DIR"]="${install_dir}/lib64/cmake/ip"
   ["crtm_DIR"]="${install_dir}/lib64/cmake/crtm"
)

# Assign and verify additional paths
for var_name in "${!paths[@]}"; do
    if [[ -n "${paths[$var_name]}" ]]; then
        assign "$var_name" "${paths[$var_name]}"
    else
        echo "[WARNING] Variable '$var_name' is empty. Skipping export."
    fi
done

echo "[OK] Environment setup completed successfully."
#---------------------------------------------------------------------------------------------------#
# Special Case: CMake for XC50                                                                      #
#---------------------------------------------------------------------------------------------------#
#if [ "${machine,,}" == "xc50" ]; then
#    cmake_path="${DIRGSI}/util/cmake-3.20.5-linux-x86_64/bin/cmake"
#    if [ -x "$cmake_path" ]; then
#        alias cmake="$cmake_path"
#        export PATH="${DIRGSI}/util/cmake-3.20.5-linux-x86_64/bin:$PATH"
#    else
#        echo "[WARNING] CMake not found in expected path: $cmake_path"
#    fi
#fi

export PATH="$PATH:${DIRGSI}/util/ecbuild/bin"

