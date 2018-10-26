/*
 * main.c
 *
 *  Created on: September 25 2018
 *      Author: Derek W. Harrison
 *
 *      This code solves the 3D transient heat equation:
 *
 *      gammax*d2T/dx2 + gammay*d2T/dy2 + gammaz*d2T/dz2 + q(x,y,z,t) = rho*Cp*dT/dt
 *
 *      On a rectangular grid. The code handles mixed boundary conditions, i.e.
 *      both Neumann and Dirichlet boundry conditions
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../inc/boundary_functions.h"
#include "../inc/export_data.h"
#include "../inc/heat3D.h"
#include "../inc/main_utils.h"
#include "../inc/mappings.h"
#include "../inc/memory_functions.h"
#include "../inc/parameters_and_boundary_types.h"
#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
    bool                           exportData = FALSE;
    domain_size_t                 domain_size = {0};
    grid_size_t                     grid_size = {0};
    time_data_t                     time_data = {0};
    grid_coordinates_t*      grid_coordinates = NULL;
    double***                               T = NULL;
    boundary_conditions_t boundary_conditions = {0};
    boundary_type_faces_t boundary_type_faces = {0};
    physical_paramaters_t physical_parameters = {0};


    /* Set parameters and boundary conditions */
    set_parameters_and_boundary_types(&domain_size,
                                      &grid_size,
                                      &time_data,
                                      &physical_parameters,
                                      &boundary_type_faces,
                                      &exportData);


    /* Allocating memory for input and output of poisson solver */
    grid_coordinates = allocate_mem_grid_coordinates(grid_size.nx+1, grid_size.ny+1, grid_size.nz+1);
    T                = matrix3D(grid_size.nx+1, grid_size.ny+1, grid_size.nz+1);


    /* Setting boundary types and conditions */
    set_boundary_conditions(boundary_type_faces,
                            &boundary_conditions);


    /* Calling 3D heat conduction solver */
    heat3D(domain_size,
           grid_size,
           boundary_conditions,
           time_data,
           physical_parameters,
           &source_equation,
           grid_coordinates,
           T);


    /* Exporting data */
    if(exportData)
    {
        export_data(grid_size, grid_coordinates, T);
    }


    /* Freeing memory */
    free_grid_coordinates(grid_coordinates, grid_size.nx, grid_size.ny);
    free_memory_3D(T, grid_size.nx, grid_size.ny);

    return 0;

}

