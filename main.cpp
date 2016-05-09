#include<iostream>
#include<vector>
#include"Grid_int.h"
#include"Grid.h"
#include"Solver.h"
#include<sys/time.h>

int main(int argc, char *argv[])
{

int l_Level = std::stoi(argv[1]);
int n_Vcycle= std::stoi(argv[2]);
std::cout<< argv[3] <<std::endl;
//char dir_neu = argv[3];


std::cout<< "l = " << l_Level << " n = " << n_Vcycle << std::endl;


struct timeval x;
gettimeofday(&x,NULL);



Solver S(l_Level,n_Vcycle,*argv[3]);

double start=x.tv_sec + 1e-6*x.tv_usec;

S.Simulation();

gettimeofday(&x,NULL);
double end=x.tv_sec + 1e-6*x.tv_usec;
std::cout<<"Run time is  ";
std::cout<<(end-start)<<" seconds"<<std::endl;



return 0;

}
