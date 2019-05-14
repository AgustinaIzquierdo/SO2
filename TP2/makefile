#Makefile de Trabajo Practico N°1 de Sistemas Operativos

all: netcdf

netcdf : readNetcdf_2.c
        icc -xCORE_AVX512 -std=gnu99 -Wall -Werror -pedantic -o eNetcdf readNetcdf_2.c -L/home/Estudiante2/Develo$
        #gcc -Wall -Werror -pedantic -o eNetcdf readNetcdf_2.c `nc-config --cflags --libs`

clean:
        rm eNetcdf \

