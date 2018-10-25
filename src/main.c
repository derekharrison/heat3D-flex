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

#include "../inc/export_data.h"
#include "../inc/heat3D.h"
#include "../inc/memory_functions.h"
#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
static double fixed_boundary_west(double y, double z, double t)
{
    /*
     * Specify the west face temperature distribution for the transient heat conduction equation
     * The distribution is of the form f(y,z,t).
     *
     * input    y
     * input    z
     * input    t
     *
     * return   temperature_west
     */

    double temperature_west = 1.0;

    return temperature_west;

}


/*-----------------------------------------------------------------------------------------------*/
static double fixed_boundary_east(double y, double z, double t)
{
    /*
     * Specify the east face temperature distribution for the transient heat conduction equation
     * The distribution is of the form f(y,z,t).
     *
     * input    y
     * input    z
     * input    t
     *
     * return   temperature_east
     */

    double temperature_east = 2.0;

    return temperature_east;

}


/*-----------------------------------------------------------------------------------------------*/
static double fixed_boundary_south(double x,double z, double t)
{
    /*
     * Specify the south face temperature distribution for the transient heat conduction equation
     * The distribution is of the form f(x,z,t).
     *
     * input    x
     * input    z
     * input    t
     *
     * return   temperature_south
     */

    double temperature_south = 3.0;

    return temperature_south;

}


/*-----------------------------------------------------------------------------------------------*/
static double fixed_boundary_north(double x,double z, double t)
{
    /*
     * Specify the north face temperature distribution for the transient heat conduction equation
     * The distribution is of the form f(x,z,t).
     *
     * input    x
     * input    z
     * input    t
     *
     * return   temperature_north
     */

    double temperature_north = 0.5;

    return temperature_north;

}


/*-----------------------------------------------------------------------------------------------*/
static double fixed_boundary_bottom(double x,double y, double t)
{
    /*
     * Specify the bottom face temperature distribution for the transient heat conduction equation
     * The distribution is of the form f(x,y,t).
     *
     * input    x
     * input    y
     * input    t
     *
     * return   temperature_bottom
     */

    double temperature_bottom = 1.3;

    return temperature_bottom;

}


/*-----------------------------------------------------------------------------------------------*/
static double fixed_boundary_top(double x,double y, double t)
{
    /*
     * Specify the top face temperature distribution for the transient heat conduction equation
     * The distribution is of the form f(x,y,t).
     *
     * input    x
     * input    y
     * input    t
     *
     * return   temperature_top
     */

    double temperature_top = 2.1;

    return temperature_top;

}


/*-----------------------------------------------------------------------------------------------*/
static double flux_boundary_west(double y, double z, double t)
{
    /*
     * Specify the west face flux equation for the transient heat conduction equation
     * The flux equation is of the form f(y,z,t).
     *
     * input    y
     * input    z
     * input    t
     *
     * return   flux_west
     */

    double flux_west = 800.0;

    return flux_west;

}


/*-----------------------------------------------------------------------------------------------*/
static double flux_boundary_east(double y, double z, double t)
{
    /*
     * Specify the east face flux equation for the transient heat conduction equation
     * The flux equation is of the form f(y,z,t).
     *
     * input    x
     * input    y
     * input    t
     *
     * return   flux_east
     */

    double flux_east = 800.0;

    return flux_east;

}


/*-----------------------------------------------------------------------------------------------*/
static double flux_boundary_south(double x,double z, double t)
{
    /*
     * Specify the south face flux equation for the transient heat conduction equation
     * The flux equation is of the form f(x,z,t).
     *
     * input    x
     * input    y
     * input    t
     *
     * return   flux_south
     */

    double flux_south = 800.0;

    return flux_south;

}


/*-----------------------------------------------------------------------------------------------*/
static double flux_boundary_north(double x,double z, double t)
{
    /*
     * Specify the north face flux equation for the transient heat conduction equation
     * The flux equation is of the form f(x,z,t).
     *
     * input    x
     * input    y
     * input    t
     *
     * return   flux_north
     */

    double flux_north = 800.0;

    return flux_north;

}


/*-----------------------------------------------------------------------------------------------*/
static double flux_boundary_bottom(double x,double y, double t)
{
    /*
     * Specify the bottom face flux equation for the transient heat conduction equation
     * The flux equation is of the form f(x,y,t).
     *
     * input    x
     * input    y
     * input    t
     *
     * return   flux_bottom
     */

    double flux_bottom = 800.0;

    return flux_bottom;

}


/*-----------------------------------------------------------------------------------------------*/
static double flux_boundary_top(double x,double y, double t)
{
    /*
     * Specify the top face flux equation for the transient heat conduction equation
     * The flux equation is of the form f(x,y,t).
     *
     * input    x
     * input    z
     * input    t
     *
     * return   flux_top
     */

    double flux_top = 800.0;

    return flux_top;

}


/*-----------------------------------------------------------------------------------------------*/
static double source_equation(double x, double y, double z, double t)
{
    /*
     * Specify the source equation for the transient heat conduction equation:
     *
     * gammax*d2T/dx2 + gammay*d2T/dy2 + gammaz*d2T/dz2 + q(x,y,z,t) = rho*Cp*dT/dt
     *
     * The source equation is of the form q(x,y,z,t).
     *
     * input    x
     * input    y
     * input    z
     *
     * return   q
     */

    double q = 100.0;

    return q;

}


/*-----------------------------------------------------------------------------------------------*/
static void set_boundary_conditions(boundary_type_faces_t boundary_type_faces,
                                    boundary_conditions_t* boundary_conditions)
{
    /*
     * Set input boundary conditions and types for heat3D solver
     *
     * input    boundary_type_faces
     * output   boundary_conditions
     */

    boundary_conditions->boundary_type_faces = boundary_type_faces;
    boundary_conditions->fixed_boundary_funcs.fixed_boundary_west = &fixed_boundary_west;
    boundary_conditions->fixed_boundary_funcs.fixed_boundary_east = &fixed_boundary_east;
    boundary_conditions->fixed_boundary_funcs.fixed_boundary_south = &fixed_boundary_south;
    boundary_conditions->fixed_boundary_funcs.fixed_boundary_north = &fixed_boundary_north;
    boundary_conditions->fixed_boundary_funcs.fixed_boundary_bottom = &fixed_boundary_bottom;
    boundary_conditions->fixed_boundary_funcs.fixed_boundary_top = &fixed_boundary_top;

    boundary_conditions->flux_boundary_funcs.flux_boundary_west = &flux_boundary_west;
    boundary_conditions->flux_boundary_funcs.flux_boundary_east = &flux_boundary_east;
    boundary_conditions->flux_boundary_funcs.flux_boundary_south = &flux_boundary_south;
    boundary_conditions->flux_boundary_funcs.flux_boundary_north = &flux_boundary_north;
    boundary_conditions->flux_boundary_funcs.flux_boundary_bottom = &flux_boundary_bottom;
    boundary_conditions->flux_boundary_funcs.flux_boundary_top = &flux_boundary_top;

}


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
    domain_size.Lx = 5.0;                               //length of domain along x coordinate
    domain_size.Ly = 1.0;                               //length of domain along y coordinate
    domain_size.Lz = 5.0;                               //length of domain along z coordinate

    grid_size.nx = 4;                                   //amount of nodes along x coordinate
    grid_size.ny = 5;                                   //amount of nodes along y coordinate
    grid_size.nz = 6;                                   //amount of nodes along z coordinate

    time_data.timesteps = 100;                          //number of timesteps
    time_data.ti        = 0.0;                          //initial time
    time_data.tf        = 0.1;                          //final time
    time_data.Tinitial  = 10.0;                         //inital temperature of system

    physical_parameters.conductivity.gammax = 15.1;     //conductivity along x coordinate
    physical_parameters.conductivity.gammay = 15.1;     //conductivity along y coordinate
    physical_parameters.conductivity.gammaz = 15.1;     //conductivity along z coordinate
    physical_parameters.Cp                  = 10.0;     //heat capacity
    physical_parameters.rho                 = 3.0;      //density

    boundary_type_faces.west_boundary   = DIRICHLET;    //west face boundary type
    boundary_type_faces.east_boundary   = DIRICHLET;    //east face boundary type
    boundary_type_faces.south_boundary  = DIRICHLET;    //bottom face boundary type
    boundary_type_faces.north_boundary  = DIRICHLET;    //north face boundary type
    boundary_type_faces.bottom_boundary = DIRICHLET;    //bottom face boundary type
    boundary_type_faces.top_boundary    = DIRICHLET;    //top face boundary type

    exportData = TRUE;                                  //export data guard


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

