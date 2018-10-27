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
    /*
     * Create boundary switch map
     *
     * input    den
     *
     * return   boundary_switch_map
     */

    double* boundary_switch_map = malloc(sizeof(double) * 2);
    boundary_switch_map[NEUMANN] = -1.0/den;
    boundary_switch_map[DIRICHLET] = 1.0;

    return boundary_switch_map;

}
