rm -r build
mkdir build
cmake -H. -Bbuild
cmake --build build

./build/frobenius cuww1.txt