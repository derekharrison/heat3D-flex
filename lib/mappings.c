/*
 * mappings.c
 *
 *  Created on: Oct 25, 2018
 *      Author: derek
 */


#include <stdlib.h>

#include "../inc/boundary_functions.h"
#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
double* boundary_switch_mapper(double den)
{
    double* boundary_switch = malloc(sizeof(double) * 2);
    boundary_switch[NEUMANN] = -1.0/den;
    boundary_switch[DIRICHLET] = 1.0;

    return boundary_switch;
}


/*-----------------------------------------------------------------------------------------------*/
boundary_type_t* boundary_type_mapper(boundary_type_faces_t boundary_type_faces)
{
    /*
     * Create boundary type map
     *
     * input    boundary_type_faces
     *
     * return   boundary_type_map
     */

    boundary_type_t* boundary_type_map = malloc(sizeof(boundary_type_t) * N_BOUNDARIES);
    boundary_type_map[WEST] = boundary_type_faces.west_boundary;
    boundary_type_map[EAST] = boundary_type_faces.east_boundary;
    boundary_type_map[SOUTH] = boundary_type_faces.south_boundary;
    boundary_type_map[NORTH] = boundary_type_faces.north_boundary;
    boundary_type_map[BOTTOM] = boundary_type_faces.bottom_boundary;
    boundary_type_map[TOP] = boundary_type_faces.top_boundary;

    return boundary_type_map;

}


/*-----------------------------------------------------------------------------------------------*/
boundary_func_type_t* boundary_func_type_mapper()
{
    /*
     * Create boundary function type map
     *
     * input    boundary_type_faces
     *
     * return   boundary_func_type_map
     */

    boundary_func_type_t* boundary_func_type_map = malloc(sizeof(boundary_func_type_t) * N_BOUNDARIES);
    boundary_func_type_map[WEST].fixed_boundary = &fixed_boundary_west;
    boundary_func_type_map[WEST].flux_boundary = &flux_boundary_west;
    boundary_func_type_map[EAST].fixed_boundary = &fixed_boundary_east;
    boundary_func_type_map[EAST].flux_boundary = &flux_boundary_east;
    boundary_func_type_map[SOUTH].fixed_boundary = &fixed_boundary_south;
    boundary_func_type_map[SOUTH].flux_boundary = &flux_boundary_south;
    boundary_func_type_map[NORTH].fixed_boundary = &fixed_boundary_north;
    boundary_func_type_map[NORTH].flux_boundary = &flux_boundary_north;
    boundary_func_type_map[BOTTOM].fixed_boundary = &fixed_boundary_bottom;
    boundary_func_type_map[BOTTOM].flux_boundary = &flux_boundary_bottom;
    boundary_func_type_map[TOP].fixed_boundary = &fixed_boundary_top;
    boundary_func_type_map[TOP].flux_boundary = &flux_boundary_top;

    return boundary_func_type_map;
}


/*-----------------------------------------------------------------------------------------------*/
func_pointer* boundary_func_mapper(boundary_conditions_t* boundary_conditions)
{
    /*
     * Create boundary function map
     *
     * input    boundary_conditions
     *
     * return    boundary_func_map
     */

    func_pointer* boundary_func_map = malloc(sizeof(func_pointer) * N_BOUNDARIES);
    boundary_func_map[WEST] = boundary_conditions->boundary_funcs.boundary_west;
    boundary_func_map[EAST] = boundary_conditions->boundary_funcs.boundary_east;
    boundary_func_map[SOUTH] = boundary_conditions->boundary_funcs.boundary_south;
    boundary_func_map[NORTH] = boundary_conditions->boundary_funcs.boundary_north;
    boundary_func_map[BOTTOM] = boundary_conditions->boundary_funcs.boundary_bottom;
    boundary_func_map[TOP] = boundary_conditions->boundary_funcs.boundary_top;

    return boundary_func_map;
}


/*-----------------------------------------------------------------------------------------------*/
void boundary_mapper(func_pointer* boundary_map,
                     boundary_conditions_t* boundary_conditions)
{
    /*
     * Create boundary map
     *
     * input    boundary_map
     * output   boundary_conditions
     */

    boundary_conditions->boundary_funcs.boundary_west = boundary_map[WEST];
    boundary_conditions->boundary_funcs.boundary_east = boundary_map[EAST];
    boundary_conditions->boundary_funcs.boundary_south = boundary_map[SOUTH];
    boundary_conditions->boundary_funcs.boundary_north = boundary_map[NORTH];
    boundary_conditions->boundary_funcs.boundary_bottom = boundary_map[BOTTOM];
    boundary_conditions->boundary_funcs.boundary_top = boundary_map[TOP];

}
