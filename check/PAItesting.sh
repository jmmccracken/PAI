#!/bin/sh
#Required files for CCM test: PAItesting_A3.00_B2.60.dat & PAItestingAB3_out.dat
#Required files for PAI test: PAItesting_A3.00_B2.60.dat & PAItestingAB3_out_PAI3.dat
n=0
echo "Testing CCM ..."
./PAI -E 3 -t 1 -L 2000 -n 1 -Ey 3 -ty 1 -f check/PAItesting_A3.00_B2.60.dat -o check/PAItestingAB3_out.test -eY temp_eY.dat
DIFF=$(diff check/PAItestingAB3_out.dat check/PAItestingAB3_out.test)
if [ "$DIFF" = "" ] 
then
    echo "  --PASSED--"
    n=$((n+1))
else
    echo "  !! FAILED"
fi

echo "Testing PAI ..."
./PAI -E 3 -t 1 -L 2000 -n 1 -Ey 3 -ty 1 -f check/PAItesting_A3.00_B2.60.dat -o check/PAItestingAB3_out_PAI3.test -eY check/PAItestingAB3_eY_PAI3.dat -PAI
DIFF=$(diff check/PAItestingAB3_out_PAI3.dat check/PAItestingAB3_out_PAI3.test)
if [ "$DIFF" = "" ] 
then
    echo "  --PASSED--"
    n=$((n+1))
else
    echo "  !! FAILED"
fi

if [ $n = 2 ]
then
    exit 0
else
    exit 1
fi
