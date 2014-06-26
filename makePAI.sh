#!/bin/sh
echo "Make ccmsrc"
g++ -c ccmsrc.cpp -o ccmsrc -std=c++11 -lstdc++ -lm -pthread -Wl,--no-as-needed

echo "Make ccmth"
g++ -c ccmth.cpp -o ccmth -std=c++11 -lstdc++ -lm -pthread -Wl,--no-as-needed

echo "Make PAI"
g++ ccmsrc ccmth -o PAI -std=c++11 -lstdc++ -lm -pthread -Wl,--no-as-needed
