#!/usr/bin/env bash

# Assume current script is in pion/tools 
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source_dir="${script_dir}/.."

pushd "${script_dir}"

cp "${script_dir}/Dockerfile_format" "${script_dir}/Dockerfile"

docker build -t pion:format .
docker run -it -v "${source_dir}":/home pion:format 

rm "${script_dir}/Dockerfile"

popd # script_dir
