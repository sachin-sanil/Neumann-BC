 // $Id$
 /**

 * @file ./EDIT

 * @brief Definitions of Grid_int class  
 
 * @author Mayank, Sachin, Shivang, Vinayak 
 * @author Mat. No. : 22040622
 */
 
 
///
	/*@brief Grid_int Class Declarations in header file Grid_int.h. Grid generation and applying 
	 boundary conditions.*/
/// 	Class definitions in Grid_int.cpp




#include "Grid_int.h"

///@brief Grid_int Constructor Definition.
Grid_int::Grid_int(int l_level)
{
level=l_level;
double elem = pow(2,level); 	///@param elem = number of elements (space between n grid points)///
ngp= elem+1;
h=1.0/elem;			///@param h = Grid spacing.
X.resize(ngp*ngp,0);		///Stores intial approximation with BC
}


/////////////////////////////// To access private members //////////////////////////////////////////

///@return h = Grid spacing.
double Grid_int::get_hValue()
{
double temp =h;
return temp;
}

///@return ngp = number of grid points in one row.
double Grid_int::get_ngpValue()
{
double temp =ngp;
return temp;
}

/////////////////////////////// To access X vector /////////////////////////////////////////////////

///@return X= Vector which stores initial U
std::vector<double> Grid_int::get_Xvalue(){
std::vector<double> u = X;
return u;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// To display X ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
 
void Grid_int::display_grid_int()
{ 
for(size_t i=0; i<ngp; ++i)
{
	{	
	for(size_t j=0; j<ngp;++j)
	std::cout<< X[i*ngp+j] << "\t";
	}
std::cout<<std::endl;
}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Functions to Apply Boundary conditions  ////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

///@brief function to apply boundary conditions 
///@param s= stores value of function sin(pi*x)*sinh(pi*y)
///@param h_= stores h value
///@param ngp_= stores ngp value



void Grid_int::boundary_con()
{

double s,pi=3.141592653589793238;
double h_=get_hValue();
double ngp_=get_ngpValue();
double temp = h_*pi;

for(size_t i=0; i<ngp_; ++i)
{
	for(size_t j=0; j<ngp_;++j)
	{		
		if(i==0 || i==ngp_-1) /// Traverses zero and n-1 row
			{
			s= sin(j*temp)*sinh(i*temp);

			//std::cout<<"loop i=3 s= "<< s << " i= " << i << " j " << j << std::endl;
			//std::cout << "i*ngp_+j " << i*ngp_+j<<std::endl;
			
			X[i*ngp_+j]=s; ///i*ngp_+j = maps 2D into 1D array X.
			}
		else if(j==0 || j==ngp_-1) ///Traverses zero and n-1 column
			{
			s= sin(j*temp)*sinh(i*temp);
			
			//std::cout<<"loop j=3 s= "<< s << " i= " << i << " j " << j << std::endl;
			//std::cout << "i*ngp_+j " << i*ngp_+j<<std::endl;

			X[i*ngp_+j]=s;
			}
	}
}
}	



//------------------------------------------------------U_exact----------------------------------------------------------------//

std::vector<double> Grid_int::U_exact(){

double pi=3.141592653589793238;
double h_=get_hValue();
double ngp_=get_ngpValue();
double temp = h_*pi;
std::vector<double> u_ex;

	for(size_t i=0;i<ngp_;++i){
	for(size_t j=0;j<ngp_;++j){

		u_ex.push_back(sin(j*temp)*sinh(i*temp));
				}
		}	

return u_ex;

}

//-------------------------------------------------------------------------------
void Grid_int::boundary_neu()
{

double s;
double h_=get_hValue();
double ngp_=get_ngpValue();



for(size_t i=0; i<ngp_; ++i)
{
	for(size_t j=0; j<ngp_;++j)
	{		
		if(i==0 || i==ngp_-1) /// Traverses zero and n-1 row
			{
			s= ((j*h_))*(1-(j*h_));

			//std::cout<<"loop i=3 s= "<< s << " i= " << i << " j " << j << std::endl;
			//std::cout << "i*ngp_+j " << i*ngp_+j<<std::endl;
			
			X[i*ngp_+j]=s; ///i*ngp_+j = maps 2D into 1D array X.
			}
		else if(j==0 || j==ngp_-1) ///Traverses zero and n-1 column
			{
			
			
			//std::cout<<"loop j=3 s= "<< s << " i= " << i << " j " << j << std::endl;
			//std::cout << "i*ngp_+j " << i*ngp_+j<<std::endl;
			std::cout<< "h value "<< h_ << std::endl;
			X[i*ngp_+j]= (-1.0)*h_;
			}
	}
}
}












