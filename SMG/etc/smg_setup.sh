#!/bin/bash -x
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
# !SCRIPT: smg_setup.sh
#
# !DESCRIPTION:
#   This script sets up the Global Modeling System (SMG) on the HPC, including
#   exporting environment variables, copying fixed files, configuring the
#   environment, compiling, and running test cases.
#
# !USAGE:
#   ./smg_setup.sh <option>
#
# !OPTIONS:
#   configure    - Configures SMG by creating directories, copying files, and modifying scripts.
#   compile      - Compiles the entire system, including BAM, GSI, and utilities.
#   testcase     - Prepares the environment and unpacks files for a test run.
#   help         - Displays this documentation.
#
# !REVISION HISTORY:
#   20 Dec 2017 - J. G. de Mattos - Initial version
#   07 Feb 2025 - Updates for robustness and efficiency
#EOP
#-----------------------------------------------------------------------------#

# Ensure hpc_name is defined before using it
if [[ -z "${hpc_name}" ]]; then
    echo "[WARNING]  Warning: hpc_name was not set before loading smg_setup.sh!"
    echo "Using default value: 'unknown'"
    export SUB="unknown"
fi

#BOP
# !FUNCTION: use_local_cmake
# !INTERFACE: use_local_cmake
# !DESCRIPTION:
#  Adds the project's local CMake binary directory to the PATH environment variable,
#  ensuring the consistent use of a specific CMake version throughout the build process.
#
#  This function dynamically detects the location of the local CMake installation based on
#  its relative position to the calling script and prepends its path to the existing PATH variable.
#
#  Usage:
#    use_local_cmake
#
# !REVISION HISTORY:
#  2024-03-05 - João Gerd - Initial implementation.
#EOP
#BOC
function use_local_cmake() {
    local script_path="$(realpath "${BASH_SOURCE[0]}")"
    local script_dir="$(dirname "$script_path")"
    local project_dir="$(dirname "$script_dir")"
    export PATH="$project_dir/utils/cmake/bin:$PATH"
    echo "[INFO] Using local CMake at: $project_dir/util/cmake/bin/cmake"
}
#EOC

#BOP
#  !FUNCTION: assign
#  !INTERFACE: assign <variable_name> <value>
#  !DESCRIPTION:
#   Defines and exports an environment variable from the provided arguments.
#EOP
#BOC
assign(){
  eval export $1=$2
}
#EOC

#BOP
#  !FUNCTION: vars_export
#  !INTERFACE: vars_export
#  !DESCRIPTION:
#   Exports necessary environment variables from a configuration file specific to the HPC system.
#   Ensures the configuration file exists before proceeding.
#EOP
#BOC
vars_export(){
  FilePaths=$(dirname ${BASH_SOURCE})/mach/${hpc_name}_paths.conf
  
  if [[ ! -f ${FilePaths} ]]; then
    echo "[FAIL] Error: Configuration file ${FilePaths} not found!" >&2
    exit 1
  fi

  while read line; do
    assign $line
  done < <(sed -r 's/^[ \t]*(.*[^ \t])[ \t]*$/\1/;/^$/d;/^#/d' ${FilePaths})
}
#EOC

#BOP
# !FUNCTION: disable_conda
# !DESCRIPTION:
#   Checks if Conda is active and deactivates it if necessary.
#EOP
#BOC
disable_conda() {
    if [[ -n "$CONDA_PREFIX" ]]; then
        echo "[WARNING]  Conda environment detected: $CONDA_PREFIX"
        echo "[ACTION] Deactivating Conda..."
        
        # Tenta desativar Conda
        if command -v conda &> /dev/null; then
            conda deactivate 2>/dev/null || source deactivate 2>/dev/null
        else
            echo "[WARNING]  Warning: Conda command not found, but CONDA_PREFIX is set!"
        fi
        
        unset CONDA_PREFIX
        unset CONDA_DEFAULT_ENV
        unset CONDA_PROMPT_MODIFIER
        echo "[ OK ] Conda has been disabled."
    else
        echo "[ OK ] No active Conda environment detected."
    fi
}
#EOC

#BOP
#  !FUNCTION: copy_fixed_files
#  !INTERFACE: copy_fixed_files
#  !DESCRIPTION:
#   Copies required fixed files necessary for running the model.
#   Uses a structured array to keep track of all required files, making it easier to manage.
#EOP
#BOC
copy_fixed_files(){
  vars_export

  if [ ${HOSTNAME:0:1} = 'e' ] || [ ${hpc_name} = "egeon" ] || [ ${hpc_name} = "jaci" ]; then
     echo "[INFO] Copying fixed files..."
     
     filesDataIn=(
       "AeroVar.Tab" 
       "ETAMPNEW_DATA" 
       "F_nwvl200_mu20_lam50_res64_t298_c080428.bin"
       "iceoptics_c080917.bin" 
       "ocnalbtab24bnd.bin"
     )
     for file in "${filesDataIn[@]}";do
        cp -pf ${public_bam}/MODEL/datain/$file ${subt_model_bam}/datain/
     done
   
     filesDataOut=(
       "HPRIME.dat"
     )
     for file in "${filesDataOut[@]}";do
        cp -pf ${public_bam}/PRE/dataout/$file ${subt_pre_bam}/dataout/
     done
   
     filesDataBC=(
       "sib2soilms.form" 
       "FluxCO2.bin" 
       "FluxCO2.ctl"
       "claymsk.form" 
       "clmt.form" 
       "deltat.form" 
       "ersst.bin"
       "ibismsk.form" 
       "ndviclm.form" 
       "sandmsk.form" 
       "sib2msk.form" 
       "soiltext.form"
     ) 
     for file in "${filesDataBC[@]}"; do
       cp -pf ${public_bam}/PRE/databcs/$file ${subt_pre_bam}/databcs/
     done
   
#     cp -pf ${home_model_bam}/datain/* ${subt_model_bam}/datain/
#     cp -pf ${home_pre_bam}/datain/* ${subt_pre_bam}/datain/
#     cp -pf ${home_pre_bam}/dataout/* ${subt_pre_bam}/dataout/
 
     cp -pf ${public_bam}/PRE/datain/2019/sst/* ${subt_pre_bam}/datain/
     cp -pf ${public_bam}/PRE/datain/2019/ncep_anl/* ${subt_pre_bam}/datain/
     cp -pf ${public_bam}/PRE/datain/2019/smc/* ${subt_pre_bam}/datain/

  elif [ ${HOSTNAME:0:1} = 'c' ]; then
     echo "[INFO] Linking fixed files..."
     ln -s /cray_home/joao_gerd/BAMFIX/model/datain/* ${subt_model_bam}/datain/
     ln -s /cray_home/joao_gerd/BAMFIX/pre/datain/* ${subt_pre_bam}/datain/
     ln -s /cray_home/joao_gerd/BAMFIX/pre/dataout/* ${subt_pre_bam}/dataout/
     ln -s /cray_home/joao_gerd/BAMFIX/pre/databcs/* ${subt_pre_bam}/databcs/
     ln -s /cray_home/joao_gerd/BAMFIX/pre/dataco2/* ${subt_pre_bam}/dataco2/
     ln -s /cray_home/joao_gerd/BAMFIX/pre/datasst/* ${subt_pre_bam}/datasst/
     ln -s /cray_home/joao_gerd/BAMFIX/pre/dataTop/* ${subt_pre_bam}/dataTop/ 
  fi

}
#EOC

#BOP
#  !FUNCTION: configure
#  !INTERFACE: configure
#  !DESCRIPTION:
#   Configures the SMG environment by creating necessary directories, copying fixed files,
#   adjusting required scripts, and ensuring the local CMake is correctly set up.
#   Ensures all essential paths are properly set up.
#
#  !CALLING SEQUENCE:
#   configure
#
#  !REMARKS:
#   - Prompts before proceeding unless AUTO_ACCEPT is set to "yes".
#   - Automatically extracts the local CMake binary if it isn't already present.
#   - Assumes CMake tarball is placed at ${SMG_ROOT}/utils/
#EOP
configure(){
  vars_export
  
  echo "[INFO] Configuring SMG..."
  
  # Confirm before proceeding unless AUTO_ACCEPT is enabled
  if [[ "${AUTO_ACCEPT}" != "yes" ]]; then
    read -p "Do you want to continue? (Y/N) " -n 1 -r
    echo
    [[ $REPLY =~ ^[Nn]$ ]] && exit 0
  fi


  # List of directories to be created
  dirs=(
    "${subt_dataout}"
    "${subt_bam}"
    "${subt_gsi}"
    "${subt_obs_dataout}"
    "${subt_run}" 
    "${subt_run_gsi}"
    "${subt_obs_run}"
    "${subt_pre_bam}/datain" 
    "${subt_pre_bam}/dataout" 
    "${subt_pre_bam}/databcs"
    "${subt_pre_bam}/datasst"
    "${subt_pre_bam}/dataco2"
    "${subt_pre_bam}/dataTop"
    "${subt_pre_bam}/exec"
    "${subt_model_bam}/datain" 
    "${subt_model_bam}/dataout"
    "${subt_model_bam}/exec"
    "${subt_pos_bam}/datain" 
    "${subt_pos_bam}/dataout"
    "${subt_pos_bam}/exec"
    "${subt_grh_bam}/datain"
    "${subt_grh_bam}/dataout"
    "${subt_gsi}/datain"
    "${subt_gsi}/dataout"
  )
  
  echo "Creating necessary directories..."
  for dir in "${dirs[@]}"; do
    if [[ ! -d $dir ]]; then
      mkdir -p $dir && echo "Created: $dir"
    else
      echo "Exists: $dir"
    fi
  done

  # Extract local CMake if not already extracted
  local project_dir="${RootDir}"
  local cmake_tarball="${project_dir}/utils/cmake-3.31.6-linux-x86_64.tar.gz"
  local cmake_dir="${project_dir}/utils/cmake"
  local cmake_exec="${cmake_dir}/bin/cmake"
  if [[ ! -x "${cmake_exec}" ]]; then
    echo "Extracting local CMake..."
    if [[ -f "${cmake_tarball}" ]]; then
      tar -xzf "${cmake_tarball}" \
          -C "${cmake_dir}" --strip-components=1
      echo "[INFO] Local CMake extracted successfully."
    else
      echo "[FAIL] Error: CMake tarball not found in utils directory."
      exit 1
    fi
  else
    echo "[INFO] Local CMake already extracted."
  fi
  
  # Ensure essential files are copied
  copy_fixed_files
  
  # Modify scripts and Makefiles
  modify_scripts

  echo "[ OK ] Configuration completed successfully."
}
#EOC

#BOP
#  !FUNCTION: modify_scripts
#  !INTERFACE: modify_scripts
#  !DESCRIPTION:
#   Modifies SMG scripts and Makefiles to ensure correct paths and environment variables.
#
#  !CALLING SEQUENCE:
#   modify_scripts
#
#  !REMARKS:
#   - Updates paths in multiple configuration files.
#   - Ensures environment variables are correctly referenced.
#EOP
#BOC
modify_scripts(){
  echo "[INFO] Modifying SMG scripts and Makefiles..."
  
  # List of scripts to update
  smg_scripts=(
    "${run_smg}/run_cycle.sh"
    "${scripts_smg}/runGSI"
    "${scripts_smg}/run_model.sh"
    "${scripts_smg}/run_obsmake.sh"
    "${scripts_smg}/run_blsdas.sh"
  )

  # Updating environment variable loading
  for script in "${smg_scripts[@]}"; do
    sed -i "/# Loading system variables/{n;d}" "$script"
    sed -i "/# Loading system variables/a\\source \"${home_smg}\"/config_smg.ksh vars_export" "$script"
  done

  # List of Makefiles to update
  makefiles=(
    "${home_model_bam}/source/Makefile"
    "${home_pos_bam}/source/Makefile"
  )

  # Updating executable paths
  for makefile in "${makefiles[@]}"; do
    sed -i "/# Path where the Model executable should be located/{n;d}" "$makefile"
    sed -i "/# Path where the Model executable should be located/a\\PATH2=\"${makefile%/*}/exec\"" "$makefile"
  done

  # Updating configurations in PostGridHistory.nml
  declare -A nml_replacements=(
    ["DirInPut"]="\"${subt_model_bam}/dataout\"/TQ0042L028"
    ["DirOutPut"]="\"${subt_grh_bam}/dataout\"/TQ0042L028"
    ["DirMain"]="\"${subt_bam}\""
  )

  for key in "${!nml_replacements[@]}"; do
    sed -i "/${key}='\//,1d" "${home_run_bam}/PostGridHistory.nml"
    sed -i "/!${key}=subt_bam/a\\${key}='${nml_replacements[$key]}'" "${home_run_bam}/PostGridHistory.nml"
  done

  # Updating environment variables in EnvironmentalVariables
  declare -A env_replacements=(
      ["HOMEBASE"]="${home_bam}"
      ["SUBTBASE"]="${subt_bam}"
      ["WORKBASE"]="${work_bam}"  # Adicionando WORKBASE se necessário
  )
  
  for key in "${!env_replacements[@]}"; do
      # Substitui a linha que começa com "export ${key}="
      sed -i "s|^export ${key}=.*|export ${key}=\"${env_replacements[$key]}\"|" "${home_run_bam}/EnvironmentalVariables"
      
      # Se a variável não existir, insere a linha abaixo do comentário "# BAM path in HOME"
      if ! grep -q "export ${key}=" "${home_run_bam}/EnvironmentalVariables"; then
          sed -i "/# BAM path in HOME/a\\export ${key}=\"${env_replacements[$key]}\"" "${home_run_bam}/EnvironmentalVariables"
      fi
  done


  echo "[ OK ] Script modifications completed."
}
#EOC



#BOP
#  !FUNCTION: compile
#  !INTERFACE: compile
#  !DESCRIPTION:
#   Compiles the entire SMG system, including BAM, GSI, and utilities.
#
#  !CALLING SEQUENCE:
#   compile
#
#  !REMARKS:
#   - Ensures all dependencies are available before compiling.
#   - Handles compilation of different components separately.
#EOP
#BOC
compile(){
  vars_export
  echo "[INFO] Starting compilation..."
  
  # Verify necessary directories
  [[ ! -d ${home_cptec}/bin ]] && mkdir -p ${home_cptec}/bin

  #if [[ ${HOSTNAME:0:1} == 'e' ]] && [[ ${HOSTNAME} != "eslogin01" && ${HOSTNAME} != "eslogin02" ]]; then
  #  echo "Please login to eslogin01 or eslogin02 before proceeding with the installation."
  #  exit 1
  #fi


  if [[ ${compgsi} -eq 1 ]]; then
    cd ${home_gsi}
    echo "[INFO] Compiling GSI..."
    echo "[INFO] PATH ${home_gsi}"

    . compile.sh -C ${compiler} 2>&1 | tee ${home_gsi}/compile.log
    if [[ ! -e ${home_gsi_src}/gsi.x ]]; then
      echo "[FAIL] Error: GSI compilation failed. Check compile.log."
      exit 1
    else  
      cp -pvfr ${home_gsi_src}/gsi.x ${home_gsi_bin}/      
      cp -pvfr ${home_gsi_src}/gsi.x ${home_cptec}/bin/      
    fi
  fi

  if [[ ${compang} -eq 1 ]]; then
    echo "[INFO] Compiling GSI bias correction utility..."

    source ./env.sh ${hpc_name} ${compiler}

    cd ${home_gsi}/util/global_angupdate
    ln -sf Makefile.conf.${hpc_name}-${compiler} Makefile.conf
    make -f Makefile clean
    make -f Makefile
    if [[ ! -e ${home_gsi}/util/global_angupdate/global_angupdate ]]; then
      echo "[FAIL] Error: GSI bias correction utility compilation failed."
      exit 1
    else  
      cp -pvfr ${home_gsi}/util/global_angupdate/global_angupdate ${home_cptec}/bin/      
    fi
    cp -pfvr ${home_gsi}/util/global_angupdate/global_angupdate ${home_cptec}/bin/global_angupdate
  fi

  if [[ ${compbam} -eq 1 ]]; then
    cd ${home_bam}
    echo "[INFO] Compiling BAM ..."
    echo "[INFO] PATH ${home_bam}"
    . compile.sh  2>&1 | tee ${home_bam}/compile.log
  fi


  # comp inctime como estava antes de instalacao na Jaci
  #if [[ ${compinctime} -eq 1 ]]; then
  #  echo "[INFO] Compiling inctime utility ..."
  #  echo "[INFO] PATH ${home_bam}"
  #  # This is just to ensure the intel env is loaded
  #  module swap gnu9/9.4.0 intel/2021.4.0
  #  cd ${util_inctime}/src
  #  export ARCH=Darwin_intel
  #  make
  #  if [[ ! -e ${util_inctime}/src/inctime ]]; then
  #    echo "[FAIL] Error: inctime utility compilation failed."
  #    exit 1
  #  else
  #    cp -pvfr ${util_inctime}/src/inctime ${home_cptec}/bin/
  #  fi
  #fi

  ### Adaptacao comp inctime considerando a maquina Jaci -> Carol 05AGT2026
  if [[ ${compinctime} -eq 1 ]]; then
    echo "[INFO] Compiling inctime utility ..."
    echo "[INFO] PATH ${home_bam}"
    
    cd ${util_inctime}/src
    make clean

    # Detecta se está na Jaci (por hostname ou variável de ambiente do Cray)
    if [[ -n "${CRAY_PRGENVINTEL}" ]] || [[ ${hpc_name} = "jaci" ]]; then
        echo "[INFO] Ambiente detectado: JACI (Cray PE)"
        # Na Jaci: Não troca módulos e força o wrapper 'ftn'
	export ARCH=jaci
	make
        #make FC=ftn F90=ftn
    else
        echo "[INFO] Ambiente nao eh JACI"
        # Na Egeon: Mantém os módulos legados e a chamada padrão
        module swap gnu9/9.4.0 intel/2021.4.0
        export ARCH=Darwin_intel
        make
    fi

    if [[ ! -e ${util_inctime}/src/inctime ]]; then
      echo "[FAIL] Error: inctime utility compilation failed."
      exit 1
    else
      mkdir -p ${home_cptec}/bin
      cp -pvfr ${util_inctime}/src/inctime ${home_cptec}/bin/
    fi
  fi

  echo "[ OK ] Compilation completed successfully."
}
#EOC

#BOP
#  !FUNCTION: testcase
#  !INTERFACE: testcase
#  !DESCRIPTION:
#   Prepares and unpacks files for a test case.
#EOP
#BOC
testcase(){
  echo -e "\e[34;1m Choose one of the available options for the testcase:\e[m"
  i=1
  for line in $(ls -1 ${public_bam}/PRE/datain); do
    year[$i]=${line}
    opts[$i]=${i}
    echo -e "\e[31;1m [$i]\e[m\e[37;1m - Testcase for \e[m\e[32;1m${year[$i]}\e[m"
    i=$((i+1))
  done
  read answer

  anlfile=$(ls -1 ${public_bam}/PRE/datain/${year[$answer]}/ncep_anl/gdas*|head -n 1)
  cp -pvfrL ${anlfile} ${subt_pre_bam}/datain/
  cp -pvfrL ${public_bam}/PRE/datain/${year[${answer}]}/sst/* ${subt_pre_bam}/datain/
#  cp -pvfrL ${public_bam}/PRE/datain/${year[${answer}]}/sno/* ${subt_pre_bam}/datain/
  cp -pvfrL ${public_bam}/PRE/datain/${year[${answer}]}/smc/*.vfm ${subt_pre_bam}/datain/
  cp -pvfrL ${public_bam}/PRE/datain/${year[${answer}]}/HybridLevels* ${subt_pre_bam}/datain/

  cp -pvfrL ${public_bam}/PRE/dataout/* ${subt_pre_bam}/dataout/
  cp -pvfrL ${public_bam}/PRE/databcs/* ${subt_pre_bam}/databcs/
  cp -pvfrL ${public_bam}/PRE/datasst/* ${subt_pre_bam}/datasst/

  cp -pvfrL ${public_bam}/MODEL/datain/* ${subt_model_bam}/datain/

  ln -sfv ${home_pos_bam}/datain/* ${subt_pos_bam}/datain/
}
#EOC

#BOP
#  !FUNCTION: help
#  !INTERFACE: help
#  !DESCRIPTION:
#   Provides help and usage instructions for the script.
#   Lists all available functions and their descriptions.
#EOP
#BOC
help(){
  echo "\nUsage: ${0##*/} <option>\n"
  echo "Available options:"
  grep -oP '^[a-z_]+(?=\(\)\{)' ${BASH_SOURCE} | while read function; do
    dsc=$(sed -n "/${function}(){$/ {n;p}" ${BASH_SOURCE} | sed "s/#*[DdEeSsCcRrIiPpTtIiOoNn]*://g")
    printf " * \e[1;31m%s\e[m --> \e[1;34m%s\e[m\n" "$function" "$dsc"
  done
}

#BOC
# Executes the selected option
#if [[ -n "$1" ]]; then
#  "$1" || echo "Invalid option! Use 'help' to see available options."
#else
#  help
#fi

#EOC

