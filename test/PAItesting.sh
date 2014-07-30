#!/bin/sh
#Required files for CCM test: PAItesting_A3.00_B2.60.dat & PAItestingAB3_out.dat
#Required files for PAI test: PAItesting_A3.00_B2.60.dat & PAItestingAB3_out_PAI3.dat

echo "Testing CCM ..."
./PAI -E 3 -t 1 -L 2000 -n 1 -Ey 3 -ty 1 -f test/PAItesting_A3.00_B2.60.dat -o test/PAItestingAB3_out.test -eY temp_eY.dat
DIFF=$(diff test/PAItestingAB3_out.dat test/PAItestingAB3_out.test)
if [ "$DIFF" = "" ] 
then
    echo "  --PASSED--"
else
    echo "  !! FAILED"
fi

echo "Testing PAI ..."
./PAI -E 3 -t 1 -L 2000 -n 1 -Ey 3 -ty 1 -f test/PAItesting_A3.00_B2.60.dat -o test/PAItestingAB3_out_PAI3.test -eY test/PAItestingAB3_eY_PAI3.dat -PAI
DIFF=$(diff test/PAItestingAB3_out_PAI3.dat test/PAItestingAB3_out_PAI3.test)
if [ "$DIFF" = "" ] 
then
    echo "  --PASSED--"
else
    echo "  !! FAILED"
fi

