 // $Id$
 /**

 * @file ./EDIT

 * @brief Definitions of Grid class  
 
 * @author Mayank, Sachin, Shivang, Vinayak 
 * @author Mat. No. : 22040622
 */
 
 

///@brief Grid Class Declarations in header file Grid_int.h.
///@brief Grid object contains 3 vectors, local grid level and number of grid points information.*/
/// Class definitions in Grid_int.cpp


#include "Grid.h"

///@brief Grid Constructor Definition.

Grid::Grid(int l_level)
{
level=l_level;
double elem = pow(2,level); ///@param elem = number of elements (space between n grid points)///
ngp= elem+1; ///@param ngp = Number of Grid points.
h=1.0/elem; ///@param h = Grid spacing.

frc.resize((ngp)*(ngp));
u_app.resize((ngp)*(ngp),(0));
res.resize((ngp)*(ngp));
}

double Grid::get_hValue()
{
double temp =h;
return temp;
}


double Grid::get_ngpValue()
{
double temp =ngp;
return temp;
}

std::vector<double> Grid::get_u_app()
{
return u_app;
}


