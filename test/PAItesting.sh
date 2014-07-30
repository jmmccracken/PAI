#!/bin/sh
#Required files for CCM test: PAItesting_A3.00_B2.60.dat & PAItestingAB3_out.dat
#Required files for PAI test: PAItesting_A3.00_B2.60.dat & PAItestingAB3_out_PAI3.dat

echo "Testing CCM ..."
./../PAI -E 3 -t 1 -L 2000 -n 1 -Ey 3 -ty 1 -f PAItesting_A3.00_B2.60.dat -o PAItestingAB3_out.test -eY temp_eY.dat
DIFF=$(diff PAItestingAB3_out.dat PAItestingAB3_out.test)
if [ "$DIFF" = "" ] 
then
    echo "  --PASSED--"
else
    echo "  !! FAILED"
fi

echo "Testing PAI ..."
./../PAI -E 3 -t 1 -L 2000 -n 1 -Ey 3 -ty 1 -f PAItesting_A3.00_B2.60.dat -o PAItestingAB3_out_PAI3.test -eY PAItestingAB3_eY_PAI3.dat -PAI
DIFF=$(diff PAItestingAB3_out_PAI3.dat PAItestingAB3_out_PAI3.test)
if [ "$DIFF" = "" ] 
then
    echo "  --PASSED--"
else
    echo "  !! FAILED"
fi

