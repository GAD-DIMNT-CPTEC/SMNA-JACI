#---------------------------------------------------------------------------------------------------#
#                           Library and GSI Compilation Script - Version 1.3                        #
#---------------------------------------------------------------------------------------------------#
#BOP                                                                                                #
# Description:                                                                                      #
#     Script to compile various libraries and the GSI model.                                        #
#                                                                                                   #
# Usage:                                                                                            #
#     ./compileLibraries.sh [options] -C <compiler>                                                 #
#                                                                                                   #
# Options:                                                                                          #
#        * -e  Compile ecbuild                                                                      #
#        * -b  Compile bacio                                                                        #
#        * -s  Compile sigbam                                                                       #
#        * -g  Compile sigio                                                                        #
#        * -f  Compile sfcio                                                                        #
#        * -n  Compile ncio                                                                         #
#        * -u  Compile bufr                                                                         #
#        * -w  Compile w3emc                                                                        #
#        * -m  Compile nemsio                                                                       #
#        * -p  Compile sp                                                                           #
#        * -i  Compile ip                                                                           #
#        * -r  Compile wrfio                                                                        #
#        * -c  Compile crtm                                                                         #
#        * -G  Compile GSI                                                                          #
#        * -a  Compile all (default)                                                                #
#        * --no-GSI  Compile all libraries except GSI                                               #
#        * -C  Specify compiler (intel or gnu) [required]                                           #
#        * -h  Display this help message                                                            #
#                                                                                                   #
# Notes:                                                                                            #
#     * The default behavior is to compile all libraries and the GSI model.                         #
#     * The compiler option (-C) is required.                                                       #
#     * Optional flags --nproc and --verbose provide more control over the compilation.             #
#                                                                                                   #
# Revisions:                                                                                        #
#     * 25-06-2024: de Mattos, J. G. - initial script version                                       #
#     * 26-06-2024: de Mattos, J. G. - added documentation and option parsing                       #
#     * 27-06-2024: de Mattos, J. G. - added compiler option                                        #
#     * 06-11-2024: de Mattos, J. G. - Refactored entire script                                     #
#     * 08-02-2025: de Mattos, J. G. - added Function to check basic system tools                   #
#     * 21-02-2025: de Mattos, J. G. - improved compile_library() to accept --nproc and --verbose   #
#     * 21-02-2025: de Mattos, J. G. - fixed argument parsing to correctly handle CMake options     #
#     * 07-03-2025: de Mattos, J. G. - added --no-GSI option to exclude GSI compilation             #
#EOP                                                                                                #
#---------------------------------------------------------------------------------------------------#


# Default: Compile everything
all=1

# Compilation flags
declare -A components=(
    [ecbuild]=0 [bacio]=0 [sigbam]=0 [sigio]=0 [sfcio]=0
    [ncio]=0 [bufr]=0 [w3emc]=0 [nemsio]=0 [sp]=0
    [ip]=0 [wrfio]=0 [crtm]=0 [GSI]=0
)

# Minimum required versions for each tool
declare -A requirements=(
    [cmake]="3.21.3"
    [make]="4.0"
    [git]="2.20.0"
)

# Local CMake version and directory name
LOCAL_CMAKE_VERSION="3.20.5-linux-x86_64"

# Compiler variable
compiler=""

# Function to display help message
usage() {
    echo "This script compiles the GSI system and its dependencies."
    echo "Usage: $0 [options] -C <compiler>"
    echo "Options:"
    echo "  -e  Compile ecbuild"
    echo "  -b  Compile bacio"
    echo "  -s  Compile sigbam"
    echo "  -g  Compile sigio"
    echo "  -f  Compile sfcio"
    echo "  -n  Compile ncio"
    echo "  -u  Compile bufr"
    echo "  -w  Compile w3emc"
    echo "  -m  Compile nemsio"
    echo "  -p  Compile sp"
    echo "  -i  Compile ip"
    echo "  -r  Compile wrfio"
    echo "  -c  Compile crtm"
    echo "  -G  Compile GSI"
    echo "  -a  Compile all (default)"
    echo "  -C  Specify compiler (intel or gnu) [required]"
    echo "  --no-GSI  Compile all libraries except GSI"
    echo "  -h  Display this help message"
}

# Parse command line options
while (( $# )); do
    opt=$1
    case "$opt" in
        -e) components[ecbuild]=1; all=0 ;;
        -b) components[bacio]=1; all=0 ;;
        -s) components[sigbam]=1; all=0 ;;
        -g) components[sigio]=1; all=0 ;;
        -f) components[sfcio]=1; all=0 ;;
        -n) components[ncio]=1; all=0 ;;
        -u) components[bufr]=1; all=0 ;;
        -w) components[w3emc]=1; all=0 ;;
        -m) components[nemsio]=1; all=0 ;;
        -p) components[sp]=1; all=0 ;;
        -i) components[ip]=1; all=0 ;;
        -r) components[wrfio]=1; all=0 ;;
        -c) components[crtm]=1; all=0 ;;
        -G) components[GSI]=1; all=0 ;;
        -a) all=1 ;;
        -C) compiler=$2; shift ;;
  --no-GSI) NO_GSI=true;;
        -h) usage; exit 0 ;;
        *) echo "[WARNING] Unknown argument: $opt" ;;
    esac
    shift
done

# Check if compiler option is set
if [[ -z "$compiler" ]]; then
    echo "Error: Compiler not specified. Use -C <compiler>."
    usage
    exit 1
fi

echo "[INFO] Compiler selected: $compiler"
echo "[INFO] Compilation options: $(for key in "${!components[@]}"; do [[ ${components[$key]} -eq 1 ]] && echo "$key"; done)"

# Function to detect the machine type and set the HPC name
detect_machine() {
    local sys_info=$(uname -a)

    if echo "$sys_info" | grep -q "cray_ari_s"; then
        export hpc_name="XC50"
        export SUB="cray"
        echo "[INFO] Detected: Cray XC50"
    
    elif echo "$sys_info" | grep -q "headnode.egeon.cptec.inpe.br"; then
        export hpc_name="egeon"
        export SUB="egeon"
        echo "[INFO] Detected: EGEON Cluster"

    elif echo "$sys_info" | grep -q "Linux ian"; then
        export hpc_name="jaci"
        export SUB="jaci"
        echo "[INFO] Detected: JACI Super"
    	
    else
        echo "[ERROR] Unknown machine: $(hostname)"
        echo "[ACTION] 1) Add the machine to the defined systems in etc/mach/"
        echo "[ACTION] 2) Define an option for it in the function copy_fixed_files inside etc/smg_setup.sh"
        return 1
    fi
}


# Function to load the environment setup script
load_env_script() {
    local machine=$1
    local compiler=$2

    if [ -f "env.sh" ]; then
        source env.sh "$machine" "$compiler"
    else
        echo "[ERROR] env.sh not found!"
        return 1
    fi
}

# Generic function to compile most libraries with optional parallelism and verbosity
compile_code() {
    local codename=$1
    local codedir=$2
    local install_dir=$3
    shift 3  # Remove the first three arguments, leaving the rest as CMake options

    # Default values 
    local nproc=$(nproc)
    local verbose=0
    local cmake_options=""
    
    # Extract optional arguments safely
    #  * --nproc=<N>  Set number of threads for compilation (default: all cores)                  #
    #  * --verbose     Enable verbose compilation output
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --nproc=*) nproc="${1#*=}"; shift ;;     # Example: --nproc=4
            --verbose) verbose=1; shift ;;           # Enable verbosity if flag is present
            *) cmake_options+=("$1"); shift ;;       # Capture all remaining arguments as CMake options
        esac
    done

    if [ ${components[$codename]} -eq 1 ] || [ $all -eq 1 ]; then
        echo "[INFO] Compiling $codename with $nproc threads..."

        if [ -z "$codedir" ] || [ -z "$install_dir" ]; then
            echo "[WARNING] Skipping $codename: Path not defined!"
            return 1
        fi

        #rm -rf "$install_dir"

        # Ensure we can enter the directory
        cd "$codedir" || { echo "[ERROR] Directory $codedir not found!"; exit 1; }
        
        mkdir -p build && cd build || { echo "[ERROR] Failed to create or enter build directory!"; exit 1; }

        # Enable verbose CMake if requested
        local cmake_cmd="cmake -DCMAKE_INSTALL_PREFIX=$install_dir ${cmake_options[*]} .."
        [ "$verbose" -eq 1 ] && cmake_cmd+=" -DCMAKE_VERBOSE_MAKEFILE=ON"

        echo "[INFO] Running: $cmake_cmd"
        eval "$cmake_cmd" || { echo "[ERROR] CMake failed for $codename!"; exit 1; }

        # Run make with verbosity if enabled
        local make_cmd="make -j$nproc"
        [ "$verbose" -eq 1 ] && make_cmd+=" VERBOSE=1"

        echo "[INFO] Running: $make_cmd"
        eval "$make_cmd" || { echo "[ERROR] Compilation failed for $codename!"; exit 1; }

        # Run make install and abort if it fails
        echo "[INFO] Installing $codename..."
        make install || { echo "[ERROR] Installation failed for $codename!"; exit 1; }

        echo "[INFO] $codename compiled successfully!"
    fi
}


# Specific compilation function for ecbuild
compile_ecbuild() {
    if [ ${components[ecbuild]} -eq 1 ] || [ $all -eq 1 ]; then
        echo "[INFO] Compiling ecbuild ..."

        rm -rf ${DIRGSI}/util/ecbuild
        cd ${DIRGSI}/util/ecbuild-develop || return 1
        mkdir -p bootstrap && cd bootstrap

        ../bin/ecbuild --prefix=${DIRGSI}/util/ecbuild ..
        make install

        echo "[INFO] ecbuild compiled successfully!"
    fi
}

# Specific compilation function for GSI
compile_gsi() {
    if [ ${components[GSI]} -eq 1 ] || [ $all -eq 1 ]; then
        echo "[INFO] Compiling GSI ..."

        cd ${DIRGSI} || return 1
        mkdir -p build && cd build

        cmake ..
        make

        echo "[INFO] GSI compiled successfully!"
    fi
}


# Fallback to local CMake if the system CMake is missing or too old
use_local_cmake() {
    # Resolve the directory of the current script
    local script_path
    script_path="$(realpath "${BASH_SOURCE[0]}")"
    local script_dir
    script_dir="$(dirname "$script_path")"

    # Build the path to the local cmake based on the version variable
    local local_cmake_dir="$script_dir/util/cmake-${LOCAL_CMAKE_VERSION}/bin"

    # Check if the local cmake binary exists and is executable
    if [[ -x "$local_cmake_dir/cmake" ]]; then
        # Prepend the local cmake directory to PATH
        export PATH="$local_cmake_dir:$PATH"
        echo "[INFO] System CMake version is below required. Using local CMake at: $local_cmake_dir/cmake"
    else
        echo "[ERROR] Local CMake not found at: $local_cmake_dir/cmake"
        return 1
    fi
}

# Check for existence and minimum version of basic tools
check_basic_tools() {
    local missing=()
    local wrong_version=()

    # Iterate over each required tool and its minimum version
    for tool in "${!requirements[@]}"; do
        local min_ver=${requirements[$tool]}

        # 1) Verify tool is installed
        if ! command -v "$tool" >/dev/null 2>&1; then
            missing+=("$tool (>= $min_ver)")
            continue
        fi

        # 2) Retrieve installed version by parsing --version or version output
        local ver_output
        ver_output=$("$tool" --version 2>/dev/null || "$tool" version 2>/dev/null)
        local inst_ver="0.0.0"
        if [[ $ver_output =~ ([0-9]+(\.[0-9]+)+) ]]; then
            inst_ver=${BASH_REMATCH[1]}
        fi

        # 3) Compare installed version against minimum using sort -V
        if [[ "$(printf '%s\n%s\n' "$min_ver" "$inst_ver" | sort -V | head -n1)" != "$min_ver" ]]; then
            wrong_version+=("$tool (installed: $inst_ver, required: >= $min_ver)")
        fi
    done

    # If any tool is missing or has too low a version
    if [[ ${#missing[@]} -gt 0 ]] || [[ ${#wrong_version[@]} -gt 0 ]]; then
        # Handle fallback for cmake only
        for msg in "${wrong_version[@]}" "${missing[@]}"; do
            if [[ $msg == cmake* ]]; then
                use_local_cmake || return 1
                # Re-validate cmake after setting local PATH
                command -v cmake >/dev/null 2>&1 || return 1
                break
            fi
        done

        # Re-check make and git (after potential CMake fallback)
        missing=(); wrong_version=()
        for tool in make git; do
            local min_ver=${requirements[$tool]}

            if ! command -v "$tool" >/dev/null 2>&1; then
                missing+=("$tool (>= $min_ver)")
                continue
            fi

            local ver_output
            ver_output=$("$tool" --version 2>/dev/null || "$tool" version 2>/dev/null)
            local inst_ver="0.0.0"
            if [[ $ver_output =~ ([0-9]+(\.[0-9]+)+) ]]; then
                inst_ver=${BASH_REMATCH[1]}
            fi

            if [[ "$(printf '%s\n%s\n' "$min_ver" "$inst_ver" | sort -V | head -n1)" != "$min_ver" ]]; then
                wrong_version+=("$tool (installed: $inst_ver, required: >= $min_ver)")
            fi
        done

        # Report any remaining issues
        if [[ ${#missing[@]} -gt 0 ]] || [[ ${#wrong_version[@]} -gt 0 ]]; then
            [[ ${#missing[@]} -gt 0 ]] && echo "[ERROR] Missing tools:" "${missing[*]}"
            [[ ${#wrong_version[@]} -gt 0 ]] && echo "[ERROR] Tools below minimum version:" "${wrong_version[*]}"
            return 1
        fi
    fi

    echo "[OK] All required tools are present with satisfactory versions."
    return 0
}



# Detect machine type
detect_machine || exit 1

# Load the environment setup script
load_env_script "$hpc_name" "$compiler" || exit 1

# Run basic tool check before proceeding
check_basic_tools || exit 1

# Compile all selected components
compile_ecbuild

echo "[INFO] Install directory: ${install_dir}"

# Libraries using the generic compilation function
compile_code "bufr" "${DIRLIB}/NCEPLIBS-bufr-develop" "${install_dir}" -DTEST_FILE_DIR="${DIRTEST}" -DKINDS="4:8:d:" --nproc=1  
compile_code "bacio" "${DIRLIB}/NCEPLIBS-bacio-develop" "${install_dir}"  --nproc=1
compile_code "sigbam" "${DIRLIB}/CPTECLIBS-sigioBAMMod-develop" "${install_dir}" --nproc=1
compile_code "sigio" "${DIRLIB}/NCEPLIBS-sigio-develop" "${install_dir}" --nproc=1
compile_code "sfcio" "${DIRLIB}/NCEPLIBS-sfcio-develop" "${install_dir}" --nproc=1
compile_code "ncio" "${DIRLIB}/NCEPLIBS-ncio-develop" "${install_dir}" --nproc=1
compile_code "w3emc" "${DIRLIB}/NCEPLIBS-w3emc-develop" "${install_dir}" -DBUILD=ON --nproc=1
compile_code "nemsio" "${DIRLIB}/NCEPLIBS-nemsio-develop" "${install_dir}" --nproc=1
compile_code "sp" "${DIRLIB}/NCEPLIBS-sp-develop" "${install_dir}" --nproc=1
compile_code "ip" "${DIRLIB}/NCEPLIBS-ip-develop" "${install_dir}" --nproc=1
compile_code "wrfio" "${DIRLIB}/NCEPLIBS-wrf_io-develop" "${install_dir}" --nproc=1
compile_code "crtm" "${DIRLIB}/crtm-master" "${install_dir}" --nproc=1

if [ ! ${NO_GSI} ];then
   compile_code "GSI" "${DIRGSI}" "${install_dir}" --nproc=1
fi

echo "[INFO] Compilation Process Finished!"

