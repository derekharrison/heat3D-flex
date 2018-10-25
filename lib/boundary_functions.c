/*
 * boundary_functions.c
 *
 *  Created on: Oct 25, 2018
 *      Author: derek
 */


/*-----------------------------------------------------------------------------------------------*/
double fixed_boundary_west(double y, double z, double t)
{
    /*
     * Specify the west face temperature distribution for the transient heat conduction equation
     * The distribution is of the form f(y,z,t).
     *
     * input    y
     * input    z
     * input    t
     *
     * return   temperature_west
     */

    double temperature_west = 1.0;

    return temperature_west;

}


/*-----------------------------------------------------------------------------------------------*/
double fixed_boundary_east(double y, double z, double t)
{
    /*
     * Specify the east face temperature distribution for the transient heat conduction equation
     * The distribution is of the form f(y,z,t).
     *
     * input    y
     * input    z
     * input    t
     *
     * return   temperature_east
     */

    double temperature_east = 2.0;

    return temperature_east;

}


/*-----------------------------------------------------------------------------------------------*/
double fixed_boundary_south(double x,double z, double t)
{
    /*
     * Specify the south face temperature distribution for the transient heat conduction equation
     * The distribution is of the form f(x,z,t).
     *
     * input    x
     * input    z
     * input    t
     *
     * return   temperature_south
     */

    double temperature_south = 3.0;

    return temperature_south;

}


/*-----------------------------------------------------------------------------------------------*/
double fixed_boundary_north(double x,double z, double t)
{
    /*
     * Specify the north face temperature distribution for the transient heat conduction equation
     * The distribution is of the form f(x,z,t).
     *
     * input    x
     * input    z
     * input    t
     *
     * return   temperature_north
     */

    double temperature_north = 0.5;

    return temperature_north;

}


/*-----------------------------------------------------------------------------------------------*/
double fixed_boundary_bottom(double x,double y, double t)
{
    /*
     * Specify the bottom face temperature distribution for the transient heat conduction equation
     * The distribution is of the form f(x,y,t).
     *
     * input    x
     * input    y
     * input    t
     *
     * return   temperature_bottom
     */

    double temperature_bottom = 1.3;

    return temperature_bottom;

}


/*-----------------------------------------------------------------------------------------------*/
double fixed_boundary_top(double x,double y, double t)
{
    /*
     * Specify the top face temperature distribution for the transient heat conduction equation
     * The distribution is of the form f(x,y,t).
     *
     * input    x
     * input    y
     * input    t
     *
     * return   temperature_top
     */

    double temperature_top = 2.1;

    return temperature_top;

}


/*-----------------------------------------------------------------------------------------------*/
double flux_boundary_west(double y, double z, double t)
{
    /*
     * Specify the west face flux equation for the transient heat conduction equation
     * The flux equation is of the form f(y,z,t).
     *
     * input    y
     * input    z
     * input    t
     *
     * return   flux_west
     */

    double flux_west = 800.0;

    return flux_west;

}


/*-----------------------------------------------------------------------------------------------*/
double flux_boundary_east(double y, double z, double t)
{
    /*
     * Specify the east face flux equation for the transient heat conduction equation
     * The flux equation is of the form f(y,z,t).
     *
     * input    x
     * input    y
     * input    t
     *
     * return   flux_east
     */

    double flux_east = 800.0;

    return flux_east;

}


/*-----------------------------------------------------------------------------------------------*/
double flux_boundary_south(double x,double z, double t)
{
    /*
     * Specify the south face flux equation for the transient heat conduction equation
     * The flux equation is of the form f(x,z,t).
     *
     * input    x
     * input    y
     * input    t
     *
     * return   flux_south
     */

    double flux_south = 800.0;

    return flux_south;

}


/*-----------------------------------------------------------------------------------------------*/
double flux_boundary_north(double x,double z, double t)
{
    /*
     * Specify the north face flux equation for the transient heat conduction equation
     * The flux equation is of the form f(x,z,t).
     *
     * input    x
     * input    y
     * input    t
     *
     * return   flux_north
     */

    double flux_north = 800.0;

    return flux_north;

}


/*-----------------------------------------------------------------------------------------------*/
double flux_boundary_bottom(double x,double y, double t)
{
    /*
     * Specify the bottom face flux equation for the transient heat conduction equation
     * The flux equation is of the form f(x,y,t).
     *
     * input    x
     * input    y
     * input    t
     *
     * return   flux_bottom
     */

    double flux_bottom = 800.0;

    return flux_bottom;

}


/*-----------------------------------------------------------------------------------------------*/
double flux_boundary_top(double x,double y, double t)
{
    /*
     * Specify the top face flux equation for the transient heat conduction equation
     * The flux equation is of the form f(x,y,t).
     *
     * input    x
     * input    z
     * input    t
     *
     * return   flux_top
     */

    double flux_top = 800.0;

    return flux_top;

}


/*-----------------------------------------------------------------------------------------------*/
double source_equation(double x, double y, double z, double t)
{
    /*
     * Specify the source equation for the transient heat conduction equation:
     *
     * gammax*d2T/dx2 + gammay*d2T/dy2 + gammaz*d2T/dz2 + q(x,y,z,t) = rho*Cp*dT/dt
     *
     * The source equation is of the form q(x,y,z,t).
     *
     * input    x
     * input    y
     * input    z
     *
     * return   q
     */

    double q = 100.0;

    return q;

}
