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


/*-------------------------------------------------------------------------------*/
static double fixed_boundary_west(double y, double z, double t) {

	return 1.0;

}


/*-------------------------------------------------------------------------------*/
static double fixed_boundary_east(double y, double z, double t)
{

	return 2.0;

}


/*-------------------------------------------------------------------------------*/
static double fixed_boundary_south(double x,double z, double t)
{

	return 3.0;

}


/*-------------------------------------------------------------------------------*/
static double fixed_boundary_north(double x,double z, double t)
{

	return 0.5;

}


/*-------------------------------------------------------------------------------*/
static double fixed_boundary_bottom(double x,double y, double t)
{

	return 1.3;

}


/*-------------------------------------------------------------------------------*/
static double fixed_boundary_top(double x,double z, double t)
{

	return 2.1;

}


/*-------------------------------------------------------------------------------*/
static double flux_boundary_west(double y, double z, double t) {

	return 1.0;

}


/*-------------------------------------------------------------------------------*/
static double flux_boundary_east(double y, double z, double t)
{

	return 2.0;

}


/*-------------------------------------------------------------------------------*/
static double flux_boundary_south(double x,double z, double t)
{

	return 3.0;

}


/*-------------------------------------------------------------------------------*/
static double flux_boundary_north(double x,double z, double t)
{

	return 0.5;

}


/*-------------------------------------------------------------------------------*/
static double flux_boundary_bottom(double x,double y, double t)
{

	return 1.3;

}


/*-------------------------------------------------------------------------------*/
static double flux_boundary_top(double x,double z, double t)
{

	return 2.1;

}


/*-------------------------------------------------------------------------------*/
static double source_equation(double x, double y, double z, double t)
{
	/*
	 * Specify the source equation for the transient heat conduction equation:
	 *
	 * gammax*d2T/dx2 + gammay*d2T/dy2 + gammaz*d2T/dz2 + q(x,y,z,t) = rho*Cp*dT/dt
	 *
	 * The source equation is of the form q(x,y,z,t).
	 *
	 * input	x
	 * input	y
	 * input 	z
	 *
	 * return	q
	 */

	double q = 100.0;

	return q;

}


/*-------------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
	bool                      	   exportData = FALSE;
	domain_size_t                 domain_size = {0};
	grid_size_t                	    grid_size = {0};
	time_dep_input_t           time_dep_input = {0};
	gammas_t                           gammas = {0};
	grid_coordinates_t* 	 grid_coordinates = NULL;
	double***                               T = NULL;
	boundary_conditions_t boundary_conditions = {0};
	boundary_type_faces_t boundary_type_faces = {0};

	/* Set parameters and boundary conditions */
	domain_size.Lx = 5.0;								//length of domain along x coordinate
	domain_size.Ly = 1.0;								//length of domain along y coordinate
	domain_size.Lz = 5.0;								//length of domain along z coordinate

	grid_size.nx = 4;									//amount of nodes along x coordinate
	grid_size.ny = 4;									//amount of nodes along y coordinate
	grid_size.nz = 4;									//amount of nodes along z coordinate

	time_dep_input.timesteps = 100;						//number of timesteps
	time_dep_input.ti 		 = 0.0;						//initial time
	time_dep_input.tf 		 = 0.1;						//final time
	time_dep_input.rho 		 = 3.0;						//density
	time_dep_input.Cp 		 = 1.0;						//heat capacity
	time_dep_input.Tinitial  = 10.0;					//inital temperature of system

	gammas.gammax = 2.0;								//conductivity along x coordinate
	gammas.gammay = 2.0;								//conductivity along y coordinate
	gammas.gammaz = 2.0;								//conductivity along z coordinate

	boundary_type_faces.west_boundary   = DIRICHLET;	//west face boundary type
	boundary_type_faces.east_boundary   = DIRICHLET;	//east face boundary type
	boundary_type_faces.south_boundary  = DIRICHLET;	//bottom face boundary type
	boundary_type_faces.north_boundary  = DIRICHLET;	//north face boundary type
	boundary_type_faces.bottom_boundary = DIRICHLET;	//bottom face boundary type
	boundary_type_faces.top_boundary    = DIRICHLET;	//top face boundary type

	exportData = TRUE;						//export data guard


	/* Allocating memory for input and output of poisson solver */
	grid_coordinates = allocate_mem_grid_coordinates(grid_size.nx+1, grid_size.ny+1, grid_size.nz+1);
	T 				 = matrix3D(grid_size.nx+1, grid_size.ny+1, grid_size.nz+1);


	/* Setting boundary types and functions */
	boundary_conditions.fixed_boundary_funcs.fixed_boundary_west = &fixed_boundary_west;
	boundary_conditions.fixed_boundary_funcs.fixed_boundary_east = &fixed_boundary_east;
	boundary_conditions.fixed_boundary_funcs.fixed_boundary_south = &fixed_boundary_south;
	boundary_conditions.fixed_boundary_funcs.fixed_boundary_north = &fixed_boundary_north;
	boundary_conditions.fixed_boundary_funcs.fixed_boundary_bottom = &fixed_boundary_bottom;
	boundary_conditions.fixed_boundary_funcs.fixed_boundary_top = &fixed_boundary_top;

	boundary_conditions.flux_boundary_funcs.flux_boundary_west = &flux_boundary_west;
	boundary_conditions.flux_boundary_funcs.flux_boundary_east = &flux_boundary_east;
	boundary_conditions.flux_boundary_funcs.flux_boundary_south = &flux_boundary_south;
	boundary_conditions.flux_boundary_funcs.flux_boundary_north = &flux_boundary_north;
	boundary_conditions.flux_boundary_funcs.flux_boundary_bottom = &flux_boundary_bottom;
	boundary_conditions.flux_boundary_funcs.flux_boundary_top = &flux_boundary_top;

	boundary_conditions.boundary_type_faces = boundary_type_faces;

	/* Calling 3D heat conduction solver */
	heat3D(domain_size,
		   grid_size,
		   boundary_conditions,
		   time_dep_input,
		   gammas,
		   &source_equation,
		   grid_coordinates,
		   T);


	/* Exporting data */
	if(exportData)
	{
		export_data("heat3D.txt", "Temperature profile", grid_size, T);
		export_data("GridX.txt", "X grid coordinates", grid_size, grid_coordinates->X);
		export_data("GridY.txt", "Y grid coordinates", grid_size, grid_coordinates->Y);
		export_data("GridZ.txt", "Z grid coordinates", grid_size, grid_coordinates->Z);
	}


	/* Freeing memory */
	free_grid_coordinates(grid_coordinates, grid_size.nx, grid_size.ny);
	free_memory_3D(T, grid_size.nx, grid_size.ny);

	return 0;

}

