#!/bin/sh
echo "Building liggghts_contactDem..."
# Or you can use g++ build
echo "Starting icpx build..."
icpx -std=c++17 main.cpp -lstdc++fs -o liggghts_contactDem
echo "Done"
