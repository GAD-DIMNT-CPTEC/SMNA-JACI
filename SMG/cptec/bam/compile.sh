#!/bin/bash
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
#
# !SCRIPT: setup_bam.sh
#
# !DESCRIPTION:
#    This script is responsible for setting up the compilation environment for 
#    the BAM model. It detects the available compiler, identifies the HPC 
#    system, loads the necessary environment modules, and compiles BAM 
#    components.
#
# !CALLING SEQUENCE:
#    ./compile.sh
#
# !REVISION HISTORY: 
#    2025-03-04 - J. G. de Mattos - Initial Version.
#
# !REMARKS:
#    - This script is designed to be run on HPC clusters such as EGEON and Cray.
#    - Ensure that required modules are available before execution.
#
#EOP
#-----------------------------------------------------------------------------#
#BOC

#BOP
# !FUNCTION: detect_compiler
# !DESCRIPTION:
#    Detects the available compiler on the system by checking common compiler 
#    commands. It sets an environment variable "compiler" with the detected name.
#
# !ARGUMENTS:
#    None (relies on system availability of compilers).
#
# !RETURN VALUE:
#    None (sets the environment variable "compiler").
#
#EOP
#!/bin/bash

detect_compiler() {
    if [[ -n "$CRAY_PRGENVINTEL" || -n "$CRAY_PRGENVCRAY" || -n "$CRAY_PRGENVPGI" || -n "$CRAY_PRGENVGNU" ]]; then
        if [[ -n "$CRAY_PRGENVINTEL" ]]; then
            echo "[INFO] Cray Intel Environment Detected"
            export compiler='intel'
        elif [[ -n "$CRAY_PRGENVCRAY" ]]; then
            echo "[INFO] Cray Native Compiler Detected"
            export compiler='cray'
        elif [[ -n "$CRAY_PRGENVPGI" ]]; then
            echo "[INFO] Cray PGI Environment Detected"
            export compiler='pgi'
        elif [[ -n "$CRAY_PRGENVGNU" ]]; then
            echo "[INFO] Cray GNU Environment Detected"
            export compiler='gnu'
        fi
    elif [[ -n "$INTEL_COMPILER_VERSION" || -n "$I_MPI_ROOT" || -n "$MKLROOT" ]]; then
        echo "[INFO] Intel Compiler Detected"
        export compiler='intel'
    elif [[ -n "$PGI" || -n "$PGI_COMPILER" ]]; then
        echo "[INFO] PGI Compiler Detected"
        export compiler='pgi'
    elif [[ -n "$GCC_VERSION" || $(gcc --version 2>/dev/null) ]]; then
        echo "[INFO] GNU Compiler (GCC) Detected"
        export compiler='gnu'
    else
        echo "No recognized compiler detected"
    fi
}

#BOP
# !FUNCTION: detect_hpc_system
# !DESCRIPTION:
#    Identifies the HPC system and sets environment variables accordingly.
#
# !ARGUMENTS:
#    None.
#
# !RETURN VALUE:
#    None (sets environment variables hpc_name and WRAPPER).
#
#EOP
detect_hpc_system() {
    local sys_info=$(uname -a)

    if echo "$sys_info" | grep -q "cray_ari_s"; then
        export hpc_system="cray"
        export hpc_name="xc50"
        export WRAPPER="ftn"
        echo "[INFO] Detected: Cray XC50"
    
    elif echo "$sys_info" | grep -q "headnode.egeon.cptec.inpe.br"; then
        export hpc_system="linux"
        export hpc_name="egeon"
        export WRAPPER="mpif90"
        export LC_ALL="en_US.UTF-8"
        echo "[INFO] Detected: EGEON Cluster"
    
    elif echo "$sys_info" | grep -q "Linux ian"; then
        export hpc_system="linux"
        export hpc_name="jaci"
        export WRAPPER="ftn"
        export LC_ALL="en_US.UTF-8"
        echo "[INFO] Detected: JACI Cluster"
    
    else
        echo "[ERROR] Unknown machine: $(hostname)"
        return 1
    fi
}

#BOP
# !FUNCTION: load_env_system
# !DESCRIPTION:
#    Loads the environment setup required to compile the BAM model.
#
# !ARGUMENTS:
#    None.
#
# !RETURN VALUE:
#    None (loads necessary modules based on detected HPC system).
#
#EOP
load_env_system() {
    echo "[INFO] Loading system components for ${hpc_system}_${hpc_name}..."


    if [ "${hpc_name}" == "egeon" ]; then
        module -q purge
        module load intel/2021.4.0
        module load mpi/2021.4.0
        module load impi/2021.4.0
        module load netcdf/4.7.4
        module load pnetcdf/1.12.2 netcdf-fortran/4.5.3

    elif [ "${hpc_name}" == "jaci" ]; then
        #module reset
	#module swap PrgEnv-cray PrgEnv-intel
        module load PrgEnv-intel
        module load cray-mpich
        module load cray-netcdf
        #module load cray-pnetcdf

	# Mapeamento com base nas variáveis nativas do Cray PE
        export NETCDF_DIR="${CRAY_NETCDF_DIR}"
        export NETCDF_C_DIR="${CRAY_NETCDF_DIR}"
        export NETCDF_FORTRAN_DIR="${CRAY_NETCDF_DIR}"
        export NetCDF_Fortran_DIR="${CRAY_NETCDF_DIR}"
        
        export PNETCDF_DIR="${CRAY_PNETCDF_DIR}"
    	
    elif [ "${hpc_name}" == "xc50" ]; then
        . /opt/modules/default/etc/modules.sh
#        module -s purge
#        module load pbs
#        module load craype
        module load craype-x86-skylake
        module load craype-network-aries
        module load cray-netcdf
        module swap PrgEnv-cray PrgEnv-gnu
        export NETCDF_FORTRAN_DIR="${NETCDF_DIR}"
    fi
}

#BOP
# !FUNCTION: compile_component
# !DESCRIPTION:
#    Compiles a specified component of the BAM model.
#
# !CALLING SEQUENCE:
#    compile_component "component_name" "subdirectory"
#
# !ARGUMENTS:
#    - component: The name of the component to compile.
#    - subdir: The subdirectory containing the component source code.
#
# !RETURN VALUE:
#    None (executes compilation commands).
#
# !REMARKS:
#    - The function ensures the component directory exists before compiling.
#    - Supports "pre", "pos", and "model" components.
#
#EOP
compile_component() {
    local component=$1
    local subdir=$2
    
    local mkname=${compiler}_${hpc_name} 
    echo "[INFO] Compiling ${component}..."

    local home_component_bam_var="home_${component}_bam"
    local home_component_bam="${!home_component_bam_var}"

    if [ -z "$home_component_bam" ]; then
        echo "[ERROR] Variable ${home_component_bam_var} is not defined!" >&2
        exit 1
    fi

    local component_dir="${home_component_bam}/${subdir}"

    if [ ! -d "$component_dir" ]; then
        echo "[ERROR] Directory ${component_dir} does not exist!" >&2
        exit 1
    fi

    cd "${component_dir}" || exit 1
    echo "[INFO] Running make clean for ${mkname}..."
    make clean "${mkname}"
    echo "[INFO] Running make for ${mkname}..."
    make "${mkname}"

    if [[ "$component" == "pre" || "$component" == "model" ]]; then
        echo "[INFO] Running make install for ${component}..."
        make install
    fi

    echo "[INFO] Compilation of ${component} completed successfully."
}

#BOP
# !FUNCTION: compile
# !DESCRIPTION:
#    Compiles all BAM model components.
#
# !ARGUMENTS:
#    None.
#
# !RETURN VALUE:
#    None (executes compilation for each component).
#
#EOP
compile() {
    compile_component "pre" "build"
    compile_component "pos" "source"

    component="model"
    echo "[INFO] Compiling ${component}...${WRAPPER}"

    home_component_bam_var="home_${component}_bam"
    home_component_bam="${!home_component_bam_var}"

    if [ -z "$home_component_bam" ]; then
        echo "[ERROR] Variable ${home_component_bam_var} is not defined!" >&2
        exit 1
    fi

    if [ -e "${home_component_bam}/build" ]; then
        rm -rf "${home_component_bam}/build"
    fi

    mkdir -p "${home_component_bam}/build"
    cd "${home_component_bam}/build" || exit 1

    cmake -DCMAKE_Fortran_COMPILER=${WRAPPER} ..
    make -j$(nproc)
    make install
}

#BOP
# !FUNCTION: main
# !DESCRIPTION:
#    Main execution logic of the script, handling environment setup and compilation.
#
# !ARGUMENTS:
#    None.
#
# !RETURN VALUE:
#    None.
#
#EOP
main() {
    detect_hpc_system
    load_env_system
    detect_compiler
    compile

    cp -pfr "${home_pre_bam}/exec/"* "${subt_pre_bam}/exec"
    cp -pfr "${home_model_bam}/exec/"* "${subt_model_bam}/exec"
    cp -pfr "${home_pos_bam}/exec/"* "${subt_pos_bam}/exec"
}

main
#EOC
#-----------------------------------------------------------------------------#

