/*
 * heat3D.h
 *
 *  Created on: Sep 25, 2018
 *      Author: Derek W. Harrison
 */

#ifndef heat3D_H_
#define heat3D_H_

#include "user_types.h"


void heat3D(domain_size_t domain_size,
		    grid_size_t grid_size,
		    boundary_conditions_t boundary_conditions,
		    time_dep_input_t time_dep_input,
		    gammas_t gammas,
			double (*source)(double x,double y,double z,double t),
			grid_coordinates_t* grid_coordinates,
			double ***T);

#endif /* heat3D_H_ */
