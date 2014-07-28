
all: PAI

ccmsrc: ccmsrc.cpp ccmsrc.h
	g++ -c ccmsrc.cpp -o ccmsrc -std=c++11 -lstdc++ -lm -pthread -Wl,--no-as-needed

ccmth: ccmth.cpp ccmsrc.h
	g++ -c ccmth.cpp -o ccmth -std=c++11 -lstdc++ -lm -pthread -Wl,--no-as-needed

PAI: ccmth ccmsrc
	g++ ccmsrc ccmth -o PAI -std=c++11 -lstdc++ -lm -pthread -Wl,--no-as-needed

clean:
	rm -f PAI ccmth ccmsrc
	rm -f *~