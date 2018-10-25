/*
 * mappings.c
 *
 *  Created on: Oct 25, 2018
 *      Author: derek
 */


#include <stdlib.h>

#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
boundary_type_t* boundary_type_mapper(boundary_type_faces_t boundary_type_faces)
{
    /*
     * Create boundary type mapping
     *
     * input    boundary_type_faces
     *
     * return    m
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
func_pointer* boundary_func_mapper(boundary_conditions_t* boundary_conditions)
{
    /*
     * Create boundary func mapping
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
     * Create boundary mapping
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
