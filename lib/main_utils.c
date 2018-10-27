/*
 * main_utils.c
 *
 *  Created on: Oct 26, 2018
 *      Author: derek
 */

#include <stdlib.h>

#include "../inc/boundary_functions.h"
#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
void set_boundary_conditions(boundary_type_faces_t boundary_type_faces,
                             boundary_conditions_t* boundary_conditions)
{
    /*
     * Set input boundary conditions and types for heat3D solver
     *
     * input    boundary_type_faces
     * output   boundary_conditions
     */

    /* Set boundary types */
    boundary_conditions->boundary_type_faces = boundary_type_faces;


    /* Set boundary functions */
    boundary_conditions->boundary_funcs.boundary_west = &boundary_west;
    boundary_conditions->boundary_funcs.boundary_east = &boundary_east;
    boundary_conditions->boundary_funcs.boundary_south = &boundary_south;
    boundary_conditions->boundary_funcs.boundary_north = &boundary_north;
    boundary_conditions->boundary_funcs.boundary_bottom = &boundary_bottom;
    boundary_conditions->boundary_funcs.boundary_top = &boundary_top;

}
