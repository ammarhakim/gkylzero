#
# Get cppcheck from github and build it as instructed there. In brief
# do:
#
# git clone https://github.com/danmar/cppcheck.git
# cd cppcheck
# mkdir build; cd build
# cmake -DCMAKE_INSTALL_PREFIX=$HOME/gkylsoft/cppcheck -DCMAKE_BUILD_TYPE=Release ..
# make install -j
#
# You may need sudo permissions to install. Please use the "Release"
# build type as the non-optimized version of cppcheck is very slow.
#
# Errors are written to cppcheck-err.txt. Please ensure all code is
# "cppcheck clean"

$HOME/gkylsoft/cppcheck/bin/cppcheck -j8 --suppressions-list=cppcheck-suppress.txt -Izero -Iminus -Iapps  --enable=warning,performance,portability,information zero apps --output-file=cppcheck-err.txt
