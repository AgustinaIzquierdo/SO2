#Makefile de Trabajo Practico N°1 de Sistemas Operativos
all: netcdf

netcdf : readNetcdf_2.c
	icc -xCORE_AVX512 -std=gnu99 -Wall -Werror -pedantic -L/home/Estudiante2/Development/SO2/hpc/libs/netcdf/lib/ -lnetcdf -O3 -o eNetcdf readNetcdf_2.c -fopenmp
	#gcc -Wall -Werror -pedantic -o eNetcdf readNetcdf_2.c `nc-config --cflags --libs`
clean:
	rm eNetcdf \

