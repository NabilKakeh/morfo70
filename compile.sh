################################################################################################
######
######  Compile Script for Morfo70: 
######  ~$sh compile.sh param1 param2 param3
######        param1: Compiler -> icc ; gcc
######        param2: Debug -> debug ; release
######        param3: OpenMP -> omp ; nomp
###### 
################################################################################################




if [ $1 = icc ]
then
	COMPILER=icc
	SUFIX=icc
else
	COMPILER=g++
	SUFIX=gcc
fi

if [ $2 = debug ]
then
	COMPILER=g++
	FLAG=-g
	SUFIX=debug
else	
	FLAG=-O3
fi

if [ $3 = omp ]
then
	OMP=1
	SUFIX=$SUFIX_omp
else	
	OMP=0
fi


$COMPILER -Wall -std=c++11 $FLAG  src/*.cpp -I include -o bin/morfo70_$SUFIX -DUSE_OMP=$OMP -lnetcdf -lm

