/*
 * main_utils.c
 *
 *  Created on: Oct 26, 2018
 *      Author: derek
 */

#include <stdlib.h>

#include "../inc/mappings.h"
#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
func_pointer set_boundary_funcs(boundary_type_t boundary_type,
                                func_pointer boundary_func_fixed,
                                func_pointer boundary_func_flux)
{
    /*
     * Set boundary function types for faces
     *
     * input    boundary_type
     * input    boundary_func_fixed
     * input    boundary_func_flux
     *
     * return   boundary_func
     */

    func_pointer boundary_func;

    if(boundary_type == DIRICHLET)
    {
        boundary_func = boundary_func_fixed;
    }
    else
    {
        boundary_func = boundary_func_flux;
    }

    return boundary_func;

}


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

    boundary_t boundary_type;


    /* Create array mappings */
    boundary_type_t* boundary_type_map = boundary_type_mapper(boundary_type_faces);
    boundary_func_type_t* boundary_func_type_map = boundary_func_type_mapper();
    func_pointer* boundary_func_map = boundary_func_mapper(boundary_conditions);


    /* Set boundary types */
    boundary_conditions->boundary_type_faces = boundary_type_faces;


    /* Set boundary functions */
    for(boundary_type = WEST; boundary_type <= TOP; boundary_type++)
    {
        boundary_func_map[boundary_type] = set_boundary_funcs(boundary_type_map[boundary_type],
                                                              boundary_func_type_map[boundary_type].fixed_boundary,
                                                              boundary_func_type_map[boundary_type].flux_boundary);
    }

    boundary_mapper(boundary_func_map,
                    boundary_conditions);


    /* Free array mappings */
    free(boundary_type_map);
    free(boundary_func_type_map);
    free(boundary_func_map);

}
