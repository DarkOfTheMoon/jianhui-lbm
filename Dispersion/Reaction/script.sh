#!/bin/sh

#a="hello world"

#echo "A is:"
#echo $a

#num=2
#echo "this is the ${num}nd"

rm *.bin
rm -rf output
rm -rf input
rm -rf working


mpic++ LB_RW.cpp -o lb
g++ Reaction.cpp -lpthread -o reaction

mkdir input
mkdir output
mkdir working

mpirun -np 6 lb INPUT_LB_RW.dat


for (( i=0; i<10; i++))
{

mv geo.bin ./input
mv vel.bin ./input

./reaction

mv ./input/geo.bin ./
mv ./input/vel.bin ./
cp ./working/GeometryUpdate.txt ./
mpirun -np 6 lb INPUT_LB_RW.dat 1

}
