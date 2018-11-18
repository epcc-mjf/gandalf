#!/bin/bash

# Source the file modules.bash that is in the same directory as this
# file.  Only works if this file is run, not sourced.  Does not follow
# symlinks.
script_dir=$(dirname $0)
source $script_dir/modules.bash
module list &> build_module.log

make clean
make executable 2>&1 | tee make.log
