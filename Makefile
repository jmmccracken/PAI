
all: PAI

ccmsrc: ccmsrc.cpp ccmsrc.h
	g++ -c ccmsrc.cpp -o ccmsrc -std=c++11 -lstdc++ -lm -pthread -Wl,--no-as-needed

ccmth: ccmth.cpp ccmsrc.h
	g++ -c ccmth.cpp -o ccmth -std=c++11 -lstdc++ -lm -pthread -Wl,--no-as-needed

PAI: ccmth ccmsrc
	g++ ccmsrc ccmth -o PAI -std=c++11 -lstdc++ -lm -pthread -Wl,--no-as-needed
	
test: PAI
	make testCCM
	make testPAI
	echo "--- All tests passed ---"

testCCM:
	./PAI -E 3 -t 1 -L 2000 -n 1 -Ey 3 -ty 1 -f check/PAItesting_A3.00_B2.60.dat -o check/PAItestingAB3_out.test -eY temp_eY.dat
	diff check/PAItestingAB3_out.dat check/PAItestingAB3_out.test

testPAI:
	./PAI -E 3 -t 1 -L 2000 -n 1 -Ey 3 -ty 1 -f check/PAItesting_A3.00_B2.60.dat -o check/PAItestingAB3_out_PAI3.test -eY check/PAItestingAB3_eY_PAI3.dat -PAI
	diff check/PAItestingAB3_out_PAI3.dat check/PAItestingAB3_out_PAI3.test

clean:
	rm -f PAI ccmth ccmsrc check/*.test check/PAItestingAB3_eY_PAI3.dat temp_eY.dat
	rm -f *~
