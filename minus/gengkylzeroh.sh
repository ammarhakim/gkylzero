#!/bin/sh

echo "Generating gkylzero.h amalgamated header file ..."

# list of header files, excluding private headers
header_list=`ls gkyl_*.h | grep -v "priv" | sort`

cat <<EOF
#pragma once

/* Amalgamated include. Generated automatically during Gkeyll "make install" */

#ifdef __cplusplus
extern "C" {
#endif
EOF

for head in $header_list
do
echo "#include <${head}>"
done

cat <<EOF1

#ifdef __cplusplus
}
#endif
EOF1


