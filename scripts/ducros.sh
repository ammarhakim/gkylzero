#!/bin/sh

# Computes modified Ducros sensor from Euler data

FILE_NAME=
LX=
NX=
USE_RHOU=0

# ----------------------------------------------------------------------------
# Function definitions
# ----------------------------------------------------------------------------

show_help() {
cat <<EOF

./ducros.sh -f=file_name -l=domain_size -n=num_cells -c

Run modified Ducros sensor. Use the -c switch to use rhou and not u.

EOF
}

# Helper functions

die() {
   echo "$*"
   echo
   echo "Dependency builds failed."
   echo
   exit 1
}

# ----------------------------------------------------------------------------
# MAIN PROGRAM
# ----------------------------------------------------------------------------

# Parse options

while [ -n "$1" ]
do
   value="`echo $1 | sed 's/[^=]*.\(.*\)/\1/'`"
   key="`echo $1 | sed 's/=.*//'`"
   if `echo "$value" | grep "~" >/dev/null 2>/dev/null`
   then
      echo
      echo '*WARNING*: the "~" sign is not expanded in flags.'
      echo 'If you mean the home directory, use $HOME instead.'
      echo
   fi
   case "$key" in
   -h)
      show_help
      exit 0
      ;;
   -c)
      USE_RHOU=1
      ;;   
   -f)
      [ -n "$value" ] || die "Missing value in flag $key."
      FILE_NAME="$value"
      ;;
   -l)
      [ -n "$value" ] || die "Missing value in flag $key."
      LX="$value"
      ;;
   -n)
      [ -n "$value" ] || die "Missing value in flag $key."
      NX="$value"
      ;;
   *)
      die "Error: Unknown flag: $1"
      ;;
   esac
   shift
done

if [ "$USE_RHOU" = "0" ]
then
    pgkyl ${FILE_NAME} -t inp sel -u inp -c1,2,3 -t rhou sel -u inp -c0 -t rho ev -t vel "rhou rho /" ev -t div2 "vel div sq" ev -t Omega "vel vel dot sqrt 0.1 * ${LX} ${NX} / /" ev -t Theta "div2 div2 Omega sq + 1e-10 + /" ev -t curl "vel curl" ev -t curl2 "curl curl dot" ev -t D0 "div2 div2 curl2 + 1e-10 + /" ev -t alpha "D0 Theta *"  pl -u alpha -a
else
    pgkyl ${FILE_NAME} -t inp sel -u inp -c1,2,3 -t rhou sel -u inp -c0 -t rho ev -t vel "rhou rho /" ev -t div2 "rhou div rho / sq" ev -t Omega "vel vel dot sqrt 0.1 * ${LX} ${NX} / /" ev -t Theta "div2 div2 Omega sq + 1e-10 + /" ev -t curl "vel curl" ev -t curl2 "curl curl dot" ev -t D0 "div2 div2 curl2 + 1e-10 + /" ev -t alpha "D0 Theta *"  pl -u alpha -a    
fi
