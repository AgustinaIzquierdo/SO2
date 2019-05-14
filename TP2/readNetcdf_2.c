#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include "/home/anij/Development/SO2/hpc/libs/netcdf/include/netcdf.h"
#include <omp.h>
#include <time.h>

//ejecutar con: gcc -o convolv readNetcdf.c -fopenmp `nc-config --cflags --libs`
//usuario y contraseñita: Estudiante2 -  3axzTWtX


#ifdef NAN
/* NAN is supported */
#endif

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {   printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

// /* nombre del archivo a leer */
#define FILE_NAME "OR_ABI-L2-CMIPF-M6C02_G16_s20191011800206_e20191011809514_c20191011809591.nc"
#define FILE_NAME2 "data.nc"

/* matriz de 21696 x 21696 */
#define NX 21696
#define NY 21696
/*Matriz para filtro*/
#define XY 3

/*Parametros globales*/
static int status,ncid2,id;
static size_t start[2]={0};
static size_t conteo[2]={0};


/*Declaracion de funciones*/
void convolve(float *imag, float filtro[XY][XY], float *imag_filt);
void save_data(int r, int s, float *imag_filt);

/*Funciones*/
void convolve(float *imag, float filtro[XY][XY], float *imag_filt)
{
    int row,col,x,y;
    row=NX;
    col=NY;
    
    omp_set_num_threads(4);
    #pragma omp parallel
    {
        for(x=0; x<row-1; x++)
        {   
            for(y=0; y<col-1; y++)
            {
                imag_filt[(x+1)*row+(y+1)]= (imag[(x)*row + (y)]       *filtro[0][0]
                                            +imag[(x)*row + (y+1)]     *filtro[0][1] 
                                            +imag[(x)*row + (y+2)]     *filtro[0][2]
                                            +imag[(x+1)*row + (y)]     *filtro[1][0]
                                            +imag[(x+1)*row + (y+1)]   *filtro[1][1] 
                                            +imag[(x+1)*row + (y+2)]   *filtro[1][2]
                                            +imag[(x+2)*row + (y)]     *filtro[2][0]
                                            +imag[(x+2)*row + (y+1)]   *filtro[2][1] 
                                            +imag[(x+2)*row + (y+2)]   *filtro[2][2]) *(float)0.00031746;
            }
        }
    }    
}

void save_data(int r,int s,float *imag_filt)
{
    /**/
    status = nc_open(FILE_NAME2, NC_WRITE, &ncid2);
    if (status != NC_NOERR) 
        ERR(status);

    // /*dimensiones de datos que voy a guardar en data.nc*/
    // status = nc_def_dim(ncid2, "filas", 21696, &filas);
    // if (status != NC_NOERR) 
    //     ERR(status);
    // status = nc_def_dim(ncid2, "columas", 21696, &columnas);
    // if (status != NC_NOERR) 
    //     ERR(status);

    // /*asignar esas dimensiones y id*/
    // dims[0] = filas;
    // dims[1] = columnas;
    // status = nc_def_var (ncid2, "CMI", NC_FLOAT, 2, dims, &id);
    // if (status != NC_NOERR) 
    //     ERR(status);

    /* Obtenemos elvarID de la variable CMI. */
    if ((status = nc_inq_varid(ncid2, "CMI", &id)))
        ERR(status);

    // status = nc_enddef(ncid2);
    // if (status != NC_NOERR) 
    //     ERR(status);

    if ((status = nc_put_vara_float(ncid2,id,start,conteo,imag_filt)))
        ERR(status);
    if ((status = nc_close(ncid2)))
        ERR(status);
}

int main()
{
    int ncid, varid;
    float *imag_in=malloc(2*NX*NY*sizeof(float));
    float *imag_filt=malloc(2*NX*NY*sizeof(float)); //matriz que almacena imagen filtrada
    float filtro[XY][XY] = {{-1.0, -1.0, -1.0},{-1.0, 8.0, -1.0},{-1.0, -1.0, -1.0}}; //filtro laplaciano
    int retval;

    clock_t inicio_convol , fin_convol;
	double time_convol;

    conteo[0]=NX;//cant filas
    conteo[1]=NY;//cant columnas

    start[0]=0;//posicion desde donde va a arrancar a guardar valores
    start[1]=0;//de imagen ya filtrada.

    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        ERR(retval);

    /* Obtenemos elvarID de la variable CMI. */
    if ((retval = nc_inq_varid(ncid, "CMI", &varid)))
        ERR(retval);
    
    /* Leemos la matriz. */
    // if ((retval = nc_get_var_float(ncid, varid, &imag_in[0][0])))
    //     ERR(retval);
    if ((retval = nc_get_vara_float(ncid, varid,start,conteo, imag_in)))
        ERR(retval);

    /* el desarrollo acá */
    /*-------------------------------------------------------------------------*/
    omp_set_num_threads(2);
    #pragma omp parallel
    {
        for(int n;n<(21696*21696);n++)
        {
            if(imag_in[n]==(-1))
            {
                imag_in[n]=(float)(0.0/0.0);
            }
        }
    }

	inicio_convol=clock();//tomo el tiempo en el que se inicia la convolucion

    convolve(imag_in, filtro, imag_filt);
    save_data(start[0],start[1], imag_filt);//arma matriz con datos de imagen
	
    fin_convol=clock();/*tomo el tiempo en el que finaliza la convolucion
                       y ya se tiene la imagen procesada.*/
    /*Calculo e imprimo tiempo total que le llevo filtrar imagen*/                       
    time_convol=(double)(fin_convol-inicio_convol)/CLOCKS_PER_SEC;
    printf("Tiempo de demora total (seg.): %.16g \n", time_convol*10);

    //--------------------------------------------------------------------------

    /* Se cierra el archivo y liberan los recursos*/
    if ((retval = nc_close(ncid)))
        ERR(retval);

    return 0;
}