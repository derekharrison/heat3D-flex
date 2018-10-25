/*
 * parameters_and_boundary_types.h
 *
 *  Created on: Oct 26, 2018
 *      Author: derek
 */

#ifndef PARAMETERS_AND_BOUNDARY_TYPES_H_
#define PARAMETERS_AND_BOUNDARY_TYPES_H_

void set_parameters_and_boundary_types(domain_size_t* domain_size,
                                       grid_size_t* grid_size,
                                       time_data_t* time_data,
                                       physical_paramaters_t* physical_parameters,
                                       boundary_type_faces_t* boundary_type_faces,
                                       bool* exportData);

#endif /* PARAMETERS_AND_BOUNDARY_TYPES_H_ */
