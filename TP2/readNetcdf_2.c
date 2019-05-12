#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <omp.h>


#ifdef NAN
/* NAN is supported */
#endif

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {   printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

// /* nombre del archivo a leer */
#define FILE_NAME "pepe.nc"
#define FILE_NAME2 "pepe2.nc"

/* Lectura de una matriz de 21696 x 21696 */
#define NX 21696
#define NY 21696
#define XY 3


//ejecutar con: gcc -o convolv readNetcdf.c -fopenmp `nc-config --cflags --libs`
//usuario y contraseñita: Estudiante2 -  3axzTWtX

/*Parametros globales*/
static int status,ncid2,id;//filas,columnas;
// static int dims[2];

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
    printf("%f", imag[NX*NX]);
    printf("%f", imag_filt[NX*NX]);
}

void save_data(int r,int s,float *imag_filt)
{
    /**/
    status = nc_open(FILE_NAME2, NC_WRITE, &ncid2);//nc_create("/home/anij/facu/2019/1erSemestre/SO2/S02/TP2/data.nc", NC_CLOBBER, &ncid2);
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

    if ((status = nc_put_vara_float(ncid2,id,start,conteo,imag_filt)))
        ERR(status);
    if ((status = nc_close(ncid2)))
        ERR(status);

}

int main()
{
    int ncid, varid;
    float *imag_in=malloc(NX*NY*sizeof(float));
    float *imag_filt=malloc(NX*NY*sizeof(float)); //matriz que almacena imagen ya filtrada
    float filtro[XY][XY] = {{-1.0, -1.0, -1.0},{-1.0, 8.0, -1.0},{-1.0, -1.0, -1.0}}; //filtro laplaciano
    int retval;

    conteo[0]=NX;//cant filas
    conteo[1]=NY;//cant columnas
    start[0] = 0;
    start[1] = 0;

    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        ERR(retval);

    /* Obtenemos elvarID de la variable CMIsave_data. */
    if ((retval = nc_inq_varid(ncid, "CMI", &varid)))
        ERR(retval);
    
    /* Leemos la matriz. */
    // if ((retval = nc_get_var_float(ncid, varid, &imag_in[0][0])))
    //     ERR(retval);
    if ((retval = nc_get_vara_float(ncid, varid,start,conteo, imag_in)))
        ERR(retval);

    /* el desarrollo acá */
    //-------------------------------------------------------------------------
    
    convolve(imag_in, filtro, imag_filt);
    
    save_data(start[0],start[1], imag_filt);//arma matriz con datos de imagen
           
    //--------------------------------------------------------------------------

    /* Se cierra el archivo y liberan los recursos*/
    if ((retval = nc_close(ncid)))
        ERR(retval);

    return 0;
}