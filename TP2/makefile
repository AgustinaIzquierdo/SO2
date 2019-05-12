#Makefile de Trabajo Practico N°1 de Sistemas Operativos

all: netcdf

netcdf : readNetcdf_2.c
	gcc -Wall -Werror -pedantic -o eNetcdf readNetcdf_2.c `nc-config --cflags --libs`

clean:
	rm eNetcdf \
