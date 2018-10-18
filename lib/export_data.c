/*
 * export_data.c
 *
 *  Created on: Sep 25, 2018
 *      Author: Derek W. Harrison
 */

#include <stdio.h>

#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
void export_temperature_field(char* file_name,
                              char* type_data,
                              grid_size_t grid_size,
                              double ***ptr)
{
    int i,j,k;
    int nx,ny,nz;

    nx = grid_size.nx;
    ny = grid_size.ny;
    nz = grid_size.nz;

    FILE *file;
    file = fopen(file_name,"w");
    if (file != NULL)
    {
        fprintf(file,"%s", type_data);

        for(j=1;j<=ny;j++)
        {
            fprintf(file,"\n\nj = %i\n\n", j);
            for(k=1;k<=nz;k++)
            {
                for(i=1;i<=nx;i++)
                    fprintf(file,"%f\t",ptr[i][j][k]);
                fprintf(file,"\n");
            }
        }
        fclose(file);
    }
    else
    {
        printf("Could not open file");
    }

}


/*-----------------------------------------------------------------------------------------------*/
void export_T_data_alongz_atxy(char* file_name,
                               grid_size_t grid_size,
                               grid_coordinates_t* grid_coordinates,
                               int node_x,
                               int node_y,
                               double ***T)
{
    int k;
    int nz;
    double ***X, ***Y, ***Z;

    X = grid_coordinates->X;
    Y = grid_coordinates->Y;
    Z = grid_coordinates->Z;

    nz = grid_size.nz;

    FILE *file;
    file = fopen(file_name,"w");
    if (file != NULL)
    {
        fprintf(file, "Temperature along Z at:\n");
        fprintf(file, "X = %f\tY = %f\n", X[node_x][node_y][1], Y[node_x][node_y][1]);
        fprintf(file, "Z\tT\n");

        for(k=1;k<=nz;k++)
        {
            fprintf(file,"%f\t%f\n", Z[node_x][node_y][k], T[node_x][node_y][k]);
        }
        fclose(file);
    }
    else
    {
        printf("Could not open file");
    }

}


/*-----------------------------------------------------------------------------------------------*/
void export_T_data_alongy_atxz(char* file_name,
                               grid_size_t grid_size,
                               grid_coordinates_t* grid_coordinates,
                               int node_x,
                               int node_z,
                               double ***T)
{
    int j;
    int ny;
    double ***X, ***Y, ***Z;

    X = grid_coordinates->X;
    Y = grid_coordinates->Y;
    Z = grid_coordinates->Z;

    ny = grid_size.ny;

    FILE *file;
    file = fopen(file_name,"w");
    if (file != NULL)
    {
        fprintf(file, "Temperature along Y at:\n");
        fprintf(file, "X = %f\tZ = %f\n", X[node_x][1][node_z], Z[node_x][1][node_z]);
        fprintf(file, "Y\tT\n");

        for(j = 1; j <= ny; j++)
        {
            fprintf(file,"%f\t%f\n", Y[node_x][j][node_z], T[node_x][j][node_z]);
        }
        fclose(file);
    }
    else
    {
        printf("Could not open file");
    }

}


/*-----------------------------------------------------------------------------------------------*/
void export_T_data_alongx_atyz(char* file_name,
                               grid_size_t grid_size,
                               grid_coordinates_t* grid_coordinates,
                               int node_y,
                               int node_z,
                               double ***T)
{
    int i;
    int nx;
    double ***X, ***Y, ***Z;

    X = grid_coordinates->X;
    Y = grid_coordinates->Y;
    Z = grid_coordinates->Z;

    nx = grid_size.nx;

    FILE *file;
    file = fopen(file_name,"w");
    if (file != NULL)
    {
        fprintf(file, "Temperature along X at:\n");
        fprintf(file, "Y = %f\tZ = %f\n", Y[1][node_y][node_z], Z[1][node_y][node_z]);
        fprintf(file, "X\tT\n");

        for(i = 1; i <= nx; i++)
        {
            fprintf(file,"%f\t%f\n", X[i][node_y][node_z], T[i][node_y][node_z]);
        }
        fclose(file);
    }
    else
    {
        printf("Could not open file");
    }

}


/*-----------------------------------------------------------------------------------------------*/
void write_data_to_vtk(char* file_name,
                       grid_size_t grid_size,
                       grid_coordinates_t* grid_coordinates,
                       double ***ptr)
{
    /*
     * Export data to vtk file for visualization in paraview.
     * The temperature field data must be printed in an order conforming
     * to the numbering scheme used:
     *
     * nn = i + (j-1)*nx + (k-1)*nx*ny
     *
     * where nn is the node number, nx the number of nodes used in the x direction,
     * ny the number of nodes used in the y direction and i, j, k are indices ranging
     * from 1 to nt where nt is the total number of nodes in the 3D grid, i.e.: nt = nx*ny*nz
     * with nz the number of nodes in the z direction.
     *
     * input    file_name
     * input    grid_size
     * input    grid_coordinates
     * input    ptr
     */

    int i,j,k;
    int nx,ny,nz,nt;
    double ***X, ***Y, ***Z;

    X = grid_coordinates->X;
    Y = grid_coordinates->Y;
    Z = grid_coordinates->Z;

    nx = grid_size.nx;
    ny = grid_size.ny;
    nz = grid_size.nz;

    nt = nx*ny*nz;

    FILE *file;
    file = fopen(file_name,"w");
    if (file != NULL)
    {
        fprintf(file,"# vtk DataFile Version 1.0\n");
        fprintf(file,"Temperature field\n");
        fprintf(file,"ASCII\n\n");
        fprintf(file,"DATASET RECTILINEAR_GRID\n");
        fprintf(file,"DIMENSIONS %i %i %i\n", nx, ny, nz);

        fprintf(file,"X_COORDINATES %i float\n", nx);
        for (i = 1; i <= grid_size.nx ; i++)
            fprintf(file,"%f ", X[i][1][1]);
        fprintf(file,"\nY_COORDINATES %i float\n", ny);
        for (j = 1; j <= grid_size.ny ; j++)
            fprintf(file,"%f ", Y[1][j][1]);
        fprintf(file,"\nZ_COORDINATES %i float\n", nz);
        for (k = 1; k <= grid_size.nz ; k++)
            fprintf(file,"%f ", Z[1][1][k]);

        fprintf(file,"\nPOINT_DATA %i\n", nt);
        fprintf(file,"FIELD FieldData 1\n");
        fprintf(file,"temperature 1 %i float\n", nt);
        for (k = 1; k <= grid_size.nz ; k++)
        for (j = 1; j <= grid_size.ny ; j++)
        for (i = 1; i <= grid_size.nx; i++)
            fprintf(file,"%f ",ptr[i][j][k]);

        fclose(file);
    }
    else
    {
        printf("Could not open file");
    }

}


/*-----------------------------------------------------------------------------------------------*/
void export_data(grid_size_t grid_size,
                 grid_coordinates_t* grid_coordinates,
                 double ***T)
{
    /*
     * Export data for visualization and analysis
     *
     * input    grid_size
     * input    grid_coordinates
     * input    T
     */

    export_temperature_field("heat3D.txt", "Temperature profile", grid_size, T);
    export_temperature_field("GridX.txt", "X grid coordinates", grid_size, grid_coordinates->X);
    export_temperature_field("GridY.txt", "Y grid coordinates", grid_size, grid_coordinates->Y);
    export_temperature_field("GridZ.txt", "Z grid coordinates", grid_size, grid_coordinates->Z);

    export_T_data_alongz_atxy("T_alongz_atxy.gnumeric", grid_size, grid_coordinates,
                             ((grid_size.nx - 1)/2)+1, ((grid_size.ny - 1)/2)+1, T);
    export_T_data_alongy_atxz("T_alongy_atxz.gnumeric", grid_size, grid_coordinates,
                             ((grid_size.nx - 1)/2)+1, ((grid_size.nz - 1)/2)+1, T);
    export_T_data_alongx_atyz("T_alongx_atyz.gnumeric", grid_size, grid_coordinates,
                             ((grid_size.ny - 1)/2)+1, ((grid_size.nz - 1)/2)+1, T);

    write_data_to_vtk("visualization.vtk", grid_size, grid_coordinates, T);
}
