#!/bin/bash

# This doesn't work quite well enough - the object files are cleaned
# between versions, so, for example, building with MEASURE, then
# plain, will cause Advisor/Amplifier/ITAC to notice that the files
# are out of date compared with the executable.

# Source the file modules.bash that is in the same directory as this
# file.  Only works if this file is run, not sourced.  Does not follow
# symlinks.
script_dir=$(dirname $0)
source $script_dir/modules.bash
module list &> make.log
echo >> make.log

#make clean
if [[ $MEASURE = true ]]; then
    if [[ $ITAC = true ]]; then
	rm -f make_measure_itac.log bin/gandalf_measure_itac
    else
	rm -f make_measure.log bin/gandalf_measure
    fi
else
    rm -f make_plain.log bin/gandalf_plain
fi

make executable 2>&1 | tee -a make.log
if [[ $MEASURE = true ]]; then
    if [[ $ITAC = true ]]; then
	mv bin/gandalf bin/gandalf_measure_itac
	mv make.log make_measure_itac.log
    else
	mv bin/gandalf bin/gandalf_measure
	mv make.log make_measure.log
    fi
else
    mv bin/gandalf bin/gandalf_plain
    mv make.log make_plain.log
fi
