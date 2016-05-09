#include "Solver.h"
#include "Grid.h"
#include "Grid_int.h"
#include "Stencil.h"
#include <fstream>

Solver::Solver(int l_level, int n_Vcycle,char c)
{
    level= l_level;
    Vcycle = n_Vcycle;
    n_d= c;

    ngp_= pow(2,level)+1;
}

Solver::Solver(int l_level, int n_Vcycle, std::vector<double> &u)
{

level= l_level;
Vcycle = n_Vcycle;

ngp_= pow(2,level)+1; 
u_initial = u;

for(int i= level;i>0;--i)
	{
	Grid *g =  new Grid(i);
	lev_Vec.push_back(g);
	}
lev_Vec[0]->u_app = u_initial;
}

//------------------------------getu_app----------------------------------------------------------//

std::vector<double> Solver::get_u_app(int l_level)
{
 int x=level- l_level;
 std::vector<double> z= lev_Vec[x]->u_app;
 return z;
}

//------------------------------ get_res -------------------------------------//

std::vector<double> Solver::get_res(int l_level){
    int x=level- l_level;
   
    std::vector<double> z= lev_Vec[x]->res;
    return z;
}

///*********************************MAPPING FUNCTION**********************************************//

int map(int i,int j,int ngl)
{
return  (i*ngl)+j;
}

///*********************************PRINT FUNCTIONS**********************************************//

void Solver::display_u()
{
std::cout<< "level "<< level << std::endl; 
//ngp_= pow(2,level)+1; 
for(size_t i=0; i<ngp_; ++i)
	{
			{	
			for(size_t j=0; j<ngp_;++j)
			std::cout<< u_initial[i*ngp_+j] << "\t";
			}
	std::cout<<std::endl;
	}
}


void Solver::display_u_app(int l_level)
{
std::cout<< '\n' <<"U approximation at level " << l_level << std::endl;
int x = level-l_level;

double ngl_=lev_Vec[x]->get_ngpValue();

//double ngl_= pow(2,l_level)+1;

for(size_t i=0; i<ngl_; ++i)
	{
		{	
		for(size_t j=0; j<ngl_;++j)
		std::cout<< lev_Vec[x]->u_app[i*ngl_+j] << "\t";
		}
	std::cout<<std::endl;
	}
}




void Solver::display_frc(int l_level)
{
std::cout<< '\n'<< "Force vector at level " << l_level << std::endl;

int x = level-l_level;
double ngl_= lev_Vec[x]->Grid::get_ngpValue();


for(size_t i=0; i<ngl_; ++i)
	{
		{	
		for(size_t j=0; j<ngl_;++j)
		std::cout<< lev_Vec[x]->frc[i*ngl_+j] << "\t";
		}
	std::cout<<std::endl;
	}
}


void Solver::display_res(int l_level)
{
std::cout<< '\n' << "Residual vector at level " << l_level << std::endl;

int x = level-l_level;
double ngl_= lev_Vec[x]->Grid::get_ngpValue();


for(size_t i=0; i<ngl_; ++i)
	{
		{	
		for(size_t j=0; j<ngl_;++j)
		std::cout<< lev_Vec[x]->res[i*ngl_+j] << "\t";
		}
	std::cout<<std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///********************************* RED-BLACK G.S. *********************************************///
////////////////////////////////////////////////////////////////////////////////////////////////////

void Solver::RBGS(int l_level)
{
int x=level-l_level;
double ngl_= lev_Vec[x]->Grid::get_ngpValue();
double h2_= 1.0/((ngl_-1)*(ngl_-1));
double h_ = 1.0/(ngl_-1);

///---------------------------------- RED UPDATE -------------------------------------------------//

	 for(int i=1;i<ngl_-1;++i){
	
		if(i & 1)
		{
            for(int j=1;j<ngl_-1;j+=2){

			lev_Vec[x]->u_app[map(i,j,ngl_)] = 0.25*(lev_Vec[x]->u_app[map(i,j-1,ngl_)]	+
								lev_Vec[x]->u_app[map(i,j+1,ngl_)]	+
								lev_Vec[x]->u_app[map(i-1,j,ngl_)]	+
								lev_Vec[x]->u_app[map(i+1,j,ngl_)]	+
                                (h2_ * lev_Vec[x]->frc[map(i,j,ngl_)]));
		}
	    }
		else{

			for(int j=2;j<ngl_-1;j+=2){
		
			lev_Vec[x]->u_app[map(i,j,ngl_)] = 0.25*(lev_Vec[x]->u_app[map(i,j-1,ngl_)]	+
								lev_Vec[x]->u_app[map(i,j+1,ngl_)]	+
								lev_Vec[x]->u_app[map(i-1,j,ngl_)]	+
								lev_Vec[x]->u_app[map(i+1,j,ngl_)]	+
                                (h2_ * lev_Vec[x]->frc[map(i,j,ngl_)]));
		}
		}

	 }

///---------------------------------- BLACK UPDATE -----------------------------------------------//

	for(int i=1;i<ngl_-1;++i){
	
		if(i & 1)
        {
	  	for(int j=2;j<ngl_-1;j+=2)
			{
			lev_Vec[x]->u_app[map(i,j,ngl_)] = 0.25*(lev_Vec[x]->u_app[map(i,j-1,ngl_)]	+
								lev_Vec[x]->u_app[map(i,j+1,ngl_)]	+
								lev_Vec[x]->u_app[map(i-1,j,ngl_)]	+
								lev_Vec[x]->u_app[map(i+1,j,ngl_)]	+
								(h2_*lev_Vec[x]->frc[map(i,j,ngl_)]));
			}
	    	}
		else
		{

			for(int j=1;j<ngl_-1;j+=2)
			{
		
			lev_Vec[x]->u_app[map(i,j,ngl_)] = 0.25*(lev_Vec[x]->u_app[map(i,j-1,ngl_)]	+
								lev_Vec[x]->u_app[map(i,j+1,ngl_)]	+
								lev_Vec[x]->u_app[map(i-1,j,ngl_)]	+
								lev_Vec[x]->u_app[map(i+1,j,ngl_)]	+
								(h2_*lev_Vec[x]->frc[map(i,j,ngl_)]));
			}
		}

	 }

if (n_d=='n')
	{
		if (level==l_level)
		{
		for(int i=1;i<ngl_-1;++i)
		{
			for(int j=0;j<ngl_;++j)
			{
			if(j==0)	
			lev_Vec[x]->u_app[map(i,j,ngl_)] = lev_Vec[x]->u_app[map(i,1,ngl_)] -h_;
			else if(j==ngl_-1)
			lev_Vec[x]->u_app[map(i,j,ngl_)] = lev_Vec[x]->u_app[map(i,j-1,ngl_)]+h_;
			}
		}
		}
		else
		{
		for(int i=1;i<ngl_-1;++i){
			for(int j=0;j<ngl_;++j)
			{
			if(j==0)	
			lev_Vec[x]->u_app[map(i,j,ngl_)] = lev_Vec[x]->u_app[map(i,1,ngl_)];
			else if(j==ngl_-1)
			lev_Vec[x]->u_app[map(i,j,ngl_)] = -lev_Vec[x]->u_app[map(i,j-1,ngl_)];
			}
		}
		} 
	}

}


///***************************** SMOOTHING FUNCTIONS *********************************************//

void Solver::pre_smoothing(int l_level)
{
    for(int i=0;i<2;++i)
        this->RBGS(l_level);
}

void Solver::post_smoothing(int l_level)
{
	this->RBGS(l_level);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///*********************************** RESIDUAL **************************************************//
////////////////////////////////////////////////////////////////////////////////////////////////////

void Solver::residual(int l_level)
{

int x=level-l_level;
double ngl_= lev_Vec[x]->Grid::get_ngpValue();
double h2_ = ((ngl_-1)*(ngl_-1)); // h2_ =(h^2)
//std::cout<< "h sqr "<< h2_ <<std::endl ;
//double norm=0;
//double temp=0;
for(int i=1;i<ngl_-1;++i)
	{
	for(int j=1;j<ngl_-1;++j)
		{
		//std::cout<< "i "<< i <<" j " << j << " to map " << map(i,j,ngl_) << std::endl ;
		lev_Vec[x]->res[map(i,j,ngl_)]=	lev_Vec[x]->frc[map(i,j,ngl_)]+	(h2_*(lev_Vec[x]->u_app[map(i,j-1,ngl_)]		+
										lev_Vec[x]->u_app[map(i,j+1,ngl_)]			+
										lev_Vec[x]->u_app[map(i-1,j,ngl_)]			+
										lev_Vec[x]->u_app[map(i+1,j,ngl_)]			-
										(4*lev_Vec[x]->u_app[map(i,j,ngl_)])));
		
		}
	
	}

}


////////////////////////////////////////////////////////////////////////////////////////////////////
///*********************************** RESTRICTION ***********************************************//
////////////////////////////////////////////////////////////////////////////////////////////////////


void Solver::restriction(int l_level)
{
int x=level-l_level;
double ngl_= lev_Vec[x]->get_ngpValue();
double ngl_1dwn = lev_Vec[x+1]->get_ngpValue();
Stencil_rest r;


for(int i=0;i<ngl_;++i)
{
	if(i%2==0)
	{
	for(int j=0;j<ngl_;++j)
	{
		if(j%2==0)
		{	
			lev_Vec[x+1]->u_app[map(i/2,j/2,ngl_1dwn)]=0;	
		}
	}
}
}



for(int i=2;i<ngl_-2;i+=2)
{
	for(int j=2;j<ngl_-2;j+=2)
	{
	lev_Vec[x+1]->frc[map(i/2,j/2,ngl_1dwn)]=    	r.W*lev_Vec[x]->res[map(i,j-1,ngl_)] 		+ 
							r.C*lev_Vec[x]->res[map(i,j,ngl_)]   		+ 
							r.E*lev_Vec[x]->res[map(i,j+1,ngl_)]  		+ 
							r.SW*lev_Vec[x]->res[map(i-1,j-1,ngl_)] 	+ 
							r.S*lev_Vec[x]->res[map(i-1,j,ngl_)] 		+ 
							r.SE*lev_Vec[x]->res[map(i-1,j+1,ngl_)] 	+ 
							r.N*lev_Vec[x]->res[map(i+1,j,ngl_)] 		+ 
							r.NW*lev_Vec[x]->res[map(i+1,j-1,ngl_)] 	+ 
							r.NE*lev_Vec[x]->res[map(i+1,j+1,ngl_)];
	}
}

}


////////////////////////////////////////////////////////////////////////////////////////////////////
///*********************************** PROLONGATION **********************************************//
////////////////////////////////////////////////////////////////////////////////////////////////////

void Solver::prolongation(int l_level)
{
int x=level-l_level; 
double ngl_= lev_Vec[x]->get_ngpValue();   //coarse grid
double ngl_up = lev_Vec[x-1]->get_ngpValue();  //fine grid


Stencil_prol p;

for(size_t i=1;i<ngl_-1;i++)
	{
	for(size_t j=1;j<ngl_-1;j++)
		{

		lev_Vec[x-1]->u_app[map(2*i-1,2*j-1,ngl_up)]	+= 	p.SW*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i-1,2*j,ngl_up)]	+= 	p.S*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i-1,2*j+1,ngl_up)]  	+= 	p.SE*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i,2*j-1,ngl_up)] 	+= 	p.W*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i,2*j,ngl_up)]         	+= 	p.C*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i,2*j+1,ngl_up)] 	+= 	p.E*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i+1,2*j-1,ngl_up)] 	+=    	p.NW*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i+1,2*j,ngl_up)] 	+= 	p.N*lev_Vec[x]->u_app[map(i,j,ngl_)];
		lev_Vec[x-1]->u_app[map(2*i+1,2*j+1,ngl_up)] 	+= 	p.NE*lev_Vec[x]->u_app[map(i,j,ngl_)];
		}
	}

}

void store(double ngp_,std::vector<double>u)
{
double hgl_ = 1.0 / (ngp_-1);

std::ofstream mg;
mg.open("solution.txt");
for(size_t i=0; i<ngp_; ++i)
{
    {
    for(size_t j=0; j<ngp_;++j)
    mg<< i*hgl_<<"\t"<< j*hgl_<<"\t" <<u[i*ngp_+j] << "\n";
    }
}
mg.close();
}

//------------------------------------  Residual Norm --------------------//

double Solver::normResidual(std::vector<double> f_res,double previous_res){
static int count=1;
double res=0,div_res,norm_res;


for(size_t i=0;i<f_res.size();i++){

        res+=(f_res[i] * f_res[i]);

    }

    div_res = res / ((ngp_-2)*(ngp_-2));

    norm_res = sqrt(div_res);
   std::cout<<"Norm of residual after V-cycle "<< norm_res<<std::endl;
	
	if(count>1)
		{
			double q = norm_res / previous_res;
            std::cout<<"Convergence q rate after V-cycle "<< q <<std::endl;
		}


   
	count++;

return norm_res;
}

//--------------------------------------  Error task-5 ----------------------------------------------//

void error(std::vector<double>u_h,std::vector<double>u_exact){

std::ofstream mg;
mg.open("E-Norm.txt",std::fstream::in | std::fstream::out | std::fstream::app);

double norm=0.;


for(size_t k=0;k<u_h.size();k++){
	
	norm+=((u_h[k]-u_exact[k])*(u_h[k]-u_exact[k]));	
}
    std::cout<<"Error Norm is "<<sqrt(norm)<<std::endl;
  	mg<<u_h.size()<<"\t"<<sqrt(norm)<<"\n";

mg.close();	
}

//-----------------------------------------------------------------------------------------------------------------------------------------------//


//////////////////////////////Simulation///////////////////////////////////////////////////

void Solver::Simulation()
{
int l_Level = this-> level;
int n_Vcycle = this-> Vcycle;

double r[n_Vcycle];   // store residual norm value 
/// Initialise grid.
Grid_int v(l_Level);

/// Apply Boundary conditions.
if(n_d=='d')
v.boundary_con();
else 
v.boundary_neu();
v.display_grid_int();

std::vector<double> u_exact=v.U_exact();


std::vector<double> u = v.get_Xvalue();


for(int i=1; i<=n_Vcycle; ++i)
    {
  std::cout<< "Current V-cycle " << i <<std::endl;

    Solver S(l_Level,i,u);

          for (int j =l_Level; j>0; --j) // Pre Smoothning -> U Print -> Residual -> Residual print -> Restriction -> Force Print
        {
        S.pre_smoothing(j); //2*S.RBGS(j);

 	 S.residual(j);
  
        if(j!=1)
        {
        S.restriction(j);
       }
    }

        S.post_smoothing(1);
   
  
    for (int k =1; k<l_Level; ++k) // Prolongation -> U Print +1 Level -> Residual -> Residual print -> Restriction -> Force Print
    {

        S.prolongation(k);
 
        if(k!=1)
        {
            S.post_smoothing(k+1);
  
        }
    

}

 	u=S.get_u_app(l_Level);
	

	S.residual(l_Level);
	
	std::vector<double> f_res=S.get_res(l_Level);
	
	
	r[i]=S.normResidual(f_res,r[i-1]);



/*
std::cout<<"\n U exact solution ...."<<std::endl;
	
for(size_t i=0; i<ngp_; ++i)
{
    {
    for(size_t j=0; j<ngp_;++j)
    std::cout<< u_exact[i*ngp_+j] << "\t";
    }
std::cout<<std::endl;
}

std::cout<<"\n U Appro. solution ...."<<std::endl;


for(size_t i=0; i<ngp_; ++i)
{
    {
    for(size_t j=0; j<ngp_;++j)
    std::cout<< u[i*ngp_+j] << "\t";
    }
std::cout<<std::endl;
}*/


}
if(l_Level>=3 && l_Level<=8){ error(u,u_exact);  }
store(ngp_,u);

}



Solver::~Solver()
{
}

