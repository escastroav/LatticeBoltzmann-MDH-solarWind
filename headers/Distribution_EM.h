/*
Distribution.h:
Encargada de definir las funciones de distribucion en un arreglo tipo <vector> unidimensional.
El indice va de menor indice a mayor indice (r,j,p,i,Lx,Ly,Lz) y se accede al índice mediante:
ix + iy*Lx + iz*Lx*Ly + i*Lx*Ly*Lz + p*Lx*Ly*Lz*4 + j*Lx*Ly*Lz*4*3 + r*Lx*Ly*Lz*4*3*2
*/
#include <iostream>
#include <vector>


const int jS = 2;
const int empS = 3;
const int emiS = 4;

/*
lattice:
int Lx
int Ly
int Lz
*/

using namespace std;
//Class definition for distribution
class hDistribution
{
 public:
  int size=0;
  bool zero = false;
  double Lx=0,Ly=0,Lz=0;
  vector<double> array;
 
  //constructor
  hDistribution(int Lx0, int Ly0, int Lz0, bool zero0);
  ~hDistribution(void);
  //acces function element
  double function(int j, int p, int i, int ix, int iy, int iz);
  double function(int ix, int iy, int iz);
  //set function element
  double & operator()(int j, int p, int i, int ix, int iy, int iz);
  double & operator()(int ix, int iy, int iz);
  

  
 private:
};
//Implementating functions
hDistribution::hDistribution(int Lx0, int Ly0, int Lz0, bool zero0)
{
  zero = zero0;
  Lx=Lx0; Ly=Ly0; Lz=Lz0;
  size = zero? Lx*Ly*Lz : jS*empS*emiS*Lx*Ly*Lz;
  array.resize(size);
}

hDistribution::~hDistribution(void)
{
  array.clear();
}

double hDistribution::function(int j, int p, int i, int ix, int iy, int iz)
{  
  double function = array.at(ix + iy*Lx + iz*Lx*Ly + i*Lx*Ly*Lz + p*Lx*Ly*Lz*emiS + j*Lx*Ly*Lz*emiS*empS);

  return function;
}

double hDistribution::function(int ix, int iy, int iz)
{  
  double function = array.at(ix + iy*Lx + iz*Lx*Ly);

  return function;
}

double & hDistribution::operator()(int j, int p, int i, int ix, int iy, int iz)
{
  return array.at(ix + iy*Lx + iz*Lx*Ly + i*Lx*Ly*Lz + p*Lx*Ly*Lz*emiS + j*Lx*Ly*Lz*emiS*empS);
}

double & hDistribution::operator()(int ix, int iy, int iz)
{
  return array.at(ix + iy*Lx + iz*Lx*Ly);   
}
