/*
 * parameters_and_boundary_types.c
 *
 *  Created on: Oct 25, 2018
 *      Author: derek
 */

#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
void set_parameters_and_boundary_types(domain_size_t* domain_size,
                                       grid_size_t* grid_size,
                                       time_data_t* time_data,
                                       physical_paramaters_t* physical_parameters,
                                       boundary_type_faces_t* boundary_type_faces,
                                       bool* exportData)
{
    /*
     * Specify parameters and boundary types for 3D poisson solver
     *
     * output   domain_size
     * output   grid_size
     * output   time_data
     * output   physical_parameters
     * output   boundary_type_faces
     * output   exportData
     */

    /* Set parameters and boundary conditions */
    domain_size->Lx = 1.0;                               //length of domain along x coordinate
    domain_size->Ly = 1.0;                               //length of domain along y coordinate
    domain_size->Lz = 1.0;                               //length of domain along z coordinate

    grid_size->nx = 10;                                   //amount of nodes along x coordinate
    grid_size->ny = 10;                                   //amount of nodes along y coordinate
    grid_size->nz = 10;                                   //amount of nodes along z coordinate

    time_data->timesteps = 1;                          //number of timesteps
    time_data->ti        = 0.0;                          //initial time
    time_data->tf        = 0.25;                          //final time
    time_data->Tinitial  = 0.0;                         //inital temperature of system

    physical_parameters->conductivity.gammax = 1.0;     //conductivity along x coordinate
    physical_parameters->conductivity.gammay = 1.0;     //conductivity along y coordinate
    physical_parameters->conductivity.gammaz = 1.0;     //conductivity along z coordinate
    physical_parameters->Cp                  = 0.0;     //heat capacity
    physical_parameters->rho                 = 1.0;      //density

    boundary_type_faces->west_boundary   = DIRICHLET;    //west face boundary type
    boundary_type_faces->east_boundary   = DIRICHLET;    //east face boundary type
    boundary_type_faces->south_boundary  = DIRICHLET;    //bottom face boundary type
    boundary_type_faces->north_boundary  = DIRICHLET;    //north face boundary type
    boundary_type_faces->bottom_boundary = DIRICHLET;    //bottom face boundary type
    boundary_type_faces->top_boundary    = DIRICHLET;    //top face boundary type

    *exportData = TRUE;                                  //export data guard

}
