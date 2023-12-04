# This assumes cppcheck is in the PATH
#
# On Linux use apt install to get cppchecl. If needed, download latest
# from http://cppcheck.sourceforge.net and use cmake to build it:
#
# cmake -DCMAKE_INSTALL_PREFIX=$HOME/gkylsoft/cppcheck -DCMAKE_BUILD_TYPE=Release .
# make -j4 install
#
# Errors are written to cppcheck-err.txt. Please ensure all code is
# "cppcheck clean"

cppcheck --suppressions-list=cppcheck-suppress.txt -Izero -Iminus -Iapps  --enable=warning,performance,portability,information zero apps 2>cppcheck-err.txt
