#!/usr/bin/env bash

# Set number of processes to use
nprocs=1

# Exit on first error
set -o errexit

# All the directories containing source files
source_dirs="source"

# Get the build directory and optional source directories from the command line
if [[ $# -ge 1 ]]
then
    # Get absolute path of build_dir from first argument.
    build_dir="$(cd "$1" && pwd)"
    shift
    if [[ $# -ge 1 ]]
    then
        source_dirs="$(cd "$1" && pwd)"
        shift
    fi
else
    echo
    echo "Error: wrong number of arguments"
    echo
    echo "usage: $(basename $0) BUILD_DIRECTORY [SOURCE_DIRECTORIES]"
    echo

    exit 1
fi

echo
echo "BUILD_DIRECTORY: ${build_dir}"
echo "SOURCE_DIRECTORIES: ${source_dirs}"
echo

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source_dir="$( cd "${script_dir}/.." && pwd)"

# Run lint.sh on every source file in pion.
#
# Find all the .cpp files in the project and run clang-tidy on them.
# Need to cd into source_dir because source_dirs are defined relative to that.
#
cd "${source_dir}"
find ${source_dirs} \( \
    -iname "*.cpp" -o -iname "*.c" \
\) -print0 \
    | xargs -0 -n 1 -P "${nprocs}" -I TARGET_FILE "${script_dir}/lint.sh" "${build_dir}" TARGET_FILE "$@"
