# This assumes cppcheck is installed in the gkylsoft directory. Please
# download latest from http://cppcheck.sourceforge.net and use cmake
# to build it:
#
# cmake -DCMAKE_INSTALL_PREFIX=$HOME/gkylsoft/cppcheck -DCMAKE_BUILD_TYPE=Release .
# make -j4 install
#
# Errors are written to cppcheck-err.txt. Please ensure all code is
# "cppcheck clean"

$HOME/gkylsoft/cppcheck/bin/cppcheck --enable=warning,performance,portability,information zero apps regression 2>cppcheck-err.txt
