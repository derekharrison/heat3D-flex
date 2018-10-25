/*
 * main_utils.h
 *
 *  Created on: Oct 26, 2018
 *      Author: derek
 */

#ifndef MAIN_UTILS_H_
#define MAIN_UTILS_H_

#include "user_types.h"


func_pointer set_boundary_funcs(boundary_type_t boundary_type,
                                func_pointer boundary_func_fixed,
                                func_pointer boundary_func_flux);

void set_boundary_conditions(boundary_type_faces_t boundary_type_faces,
                             boundary_conditions_t* boundary_conditions);

#endif /* MAIN_UTILS_H_ */
