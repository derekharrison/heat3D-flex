/*
 * export_data.c
 *
 *  Created on: Sep 25, 2018
 *      Author: Derek W. Harrison
 */

#include <stdio.h>

#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
void export_data(char* file_name,
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
