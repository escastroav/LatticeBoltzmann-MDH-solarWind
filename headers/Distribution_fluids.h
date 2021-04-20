/*
Distribution.h:
Encargada de definir las funciones de distribucion en un arreglo tipo <vector> unidimensional.
El indice va de menor indice a mayor indice (r,j,p,i,Lx,Ly,Lz) y se accede al índice mediante:
ix + iy*Lx + iz*Lx*Ly + i*Lx*Ly*Lz + p*Lx*Ly*Lz*4 + j*Lx*Ly*Lz*4*3 + r*Lx*Ly*Lz*4*3*2
*/
#include <iostream>
#include <vector>

const int sS = 2;
const int pS = 3;
const int iS = 6;

/*
lattice:
int Lx
int Ly
int Lz
*/

using namespace std;
//Class definition for distribution
class fDistribution
{
 public:
  int size=0;
  bool zero = false;
  double Lx=0,Ly=0,Lz=0;
  vector<double> array;
 
  //constructor
  fDistribution(int Lx0, int Ly0, int Lz0, bool zero0);
  ~fDistribution(void);
  //acces function element
  double function(int s, int p, int i, int ix, int iy, int iz);
  double function(int s, int ix, int iy, int iz);
  //set function element
  double & operator()(int s, int p, int i, int ix, int iy, int iz);
  double & operator()(int s, int ix, int iy, int iz);
  

  
 private:
};
//Implementating functions
fDistribution::fDistribution(int Lx0, int Ly0, int Lz0, bool zero0)
{
  zero = zero0;
  Lx=Lx0; Ly=Ly0; Lz=Lz0;
  size = zero? sS*Lx*Ly*Lz : sS*pS*iS*Lx*Ly*Lz;
  array.resize(size);
}

fDistribution::~fDistribution(void)
{
  array.clear();
}

double fDistribution::function(int s, int p, int i, int ix, int iy, int iz)
{  
  double function = array.at(ix + iy*Lx + iz*Lx*Ly + i*Lx*Ly*Lz + p*Lx*Ly*Lz*iS + s*Lx*Ly*Lz*iS*pS);

  return function;
}

double fDistribution::function(int s, int ix, int iy, int iz)
{  
  double function = array.at(ix + iy*Lx + iz*Lx*Ly + s*Lx*Ly*Lz);

  return function;
}

double & fDistribution::operator()(int s, int p, int i, int ix, int iy, int iz)
{
  return array.at(ix + iy*Lx + iz*Lx*Ly + i*Lx*Ly*Lz + p*Lx*Ly*Lz*iS + s*Lx*Ly*Lz*iS*pS);
}

double & fDistribution::operator()(int s, int ix, int iy, int iz)
{
  return array.at(ix + iy*Lx + iz*Lx*Ly + s*Lx*Ly*Lz);   
}
