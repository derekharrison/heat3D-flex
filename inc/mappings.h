/*
 * mappings.h
 *
 *  Created on: Oct 25, 2018
 *      Author: derek
 */

#ifndef MAPPINGS_H_
#define MAPPINGS_H_


#include "user_types.h"

double* boundary_switch_mapper(double den);
boundary_type_t* boundary_type_mapper(boundary_type_faces_t boundary_type_faces);
boundary_func_type_t* boundary_func_type_mapper();
func_pointer* boundary_func_mapper(boundary_conditions_t* boundary_conditions);
void boundary_mapper(func_pointer* boundary_map,
                  boundary_conditions_t* boundary_conditions);

#endif /* MAPPINGS_H_ */
