/*
 * mappings.c
 *
 *  Created on: Oct 25, 2018
 *      Author: derek
 */


#include <stdlib.h>

#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
boundary_type_t* boundary_type_map(boundary_type_faces_t boundary_type_faces)
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
