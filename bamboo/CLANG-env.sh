# No module environment on the Mac
[[ $(uname) == "Darwin" ]] && return
# Initialise module environment if it is not
if [[ ! $(command -v module > /dev/null 2>&1) ]]; then
  . /usr/local/apps/module/init/bash
fi
module unload grib_api
module unload eccodes
#module switch gnu clang/3.6.2
#module switch gnu clang/3.9.1
module switch gnu clang
