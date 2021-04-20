#include "headers/LB_MHD_D3Q19_F.h"

int Lx0=1;
int Ly0=64;
int Lz0=1;

class Hartmann : public LatticeBoltzmann
{
protected:
  double H=5,Bo=1,eta=1,rh0=1,rh1=1,g=1,ff=0;
public:
  Hartmann(void);
  void Inicie(void);
  void ImponerB(void);
  void Imprimase(const char* fileName,bool useNew);
};
Hartmann::Hartmann(void)
{
  g=0.0001;
  rh0=1.0;
  rh1=1.0;
  eta=1.0/3.0;
  H=13;
  Bo=H*sqrt(rh0*eta)/Ly0;
  Fext0[0]=g;Fext1[0]=g;
  //Fext0.cargue(g,0,0);
  //Fext1.cargue(g,0,0);
  //xi[0]=xi[1]=0;
}
void Hartmann::Inicie(void)
{
  double rho00=0,rho10=0,rhoc0=0;
  vector3D Q00,Q10,E0,B0,F00,F10,J0,Jp0,Ep0,Q0p0,Q1p0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	  {
	    rho00=rh0;	rho10=rh1;	rhoc0=0;
	    Q00.cargue(0,0,0);	Q10.cargue(0,0,0);	J0.cargue(0,0,0);
	    E0.cargue(0,0,0);		B0.cargue(0,Bo,0);	    
	    F00=F0(E0,B0,Q00,Q10,rho00,rho10);
	    F10=F1(E0,B0,Q00,Q10,rho00,rho10);	
	    Jp0=Jp(F00,F10,J0);		Ep0=Ep(E0,Jp0);
	    Q0p0=Q0p(Q00,F00);		Q1p0=Q1p(Q10,F10);

	    for(int p=0;p<3;p++)
	      for(int i=0;i<6;i++)
		{ 
		  f(0,p,i,ix,iy,iz)=f0eq(rho00,Q00,i,p);
		  f0(0,ix,iy,iz)=f0eq0(rho00,Q00,i,p);		  
		  f(1,p,i,ix,iy,iz)=f1eq(rho10,Q10,i,p);
		  f0(1,ix,iy,iz)=f1eq0(rho10,Q10,i,p);		  
		}
	  
	    for(int p=0;p<3;p++)
	      for(int i=0;i<4;i++)
		for(int j=0;j<2;j++)
		  { 
		    h(j,p,i,ix,iy,iz)=heq(Ep0,B0,j,i,p);
		    h0(ix,iy,iz)=0;
		  }
	  }
}
void Hartmann::ImponerB(void)
{
  double rho00=0,rho10=0,rhoc0=0;
  vector3D Q00,Q10,E0,B0,F00,F10,J0,Jp0,Ep0,Q0p0,Q1p0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	  {
	    B0[1]=Bo;	B0[2]=0;
	    if(iy==0||iy==Ly0-1)
	      {
		B0[0]=0;
		Q00[0]=0;
		Q10[0]=0;
	      }
	    else
	      {
		B0[0]=B(ix,iy,iz,false).x();
		Q00[0]=Q0(ix,iy,iz,false).x();
		Q10[0]=Q1(ix,iy,iz,false).x();
	      }
	    Q00[1]=0;	Q00[2]=0;
	    Q10[1]=0;	Q10[2]=0;	
	    rho00=rho0(ix,iy,iz,false);	rho10=rho1(ix,iy,iz,false);	rhoc0=rhoc();
	    J0=J(Q00,Q10);
	    E0.cargue(0,0,0);//=E(ix,iy,iz,false);
	    F00=F0(E0,B0,Q00,Q10,rho00,rho10);
	    F10=F0(E0,B0,Q00,Q10,rho00,rho10);	
	    Jp0=Jp(F00,F10,J0);		Ep0=Ep(E0,Jp0);
	    Q0p0=Q0p(Q00,F00);		Q1p0=Q1p(Q10,F10);

	    for(int p=0;p<3;p++)
	      for(int i=0;i<6;i++)
		{ 
		  fnew(0,p,i,ix,iy,iz)=f0eq(rho00,Q0p0,i,p);
		  f0new(0,ix,iy,iz)=f0eq0(rho00,Q0p0,i,p);		  
		  fnew(1,p,i,ix,iy,iz)=f1eq(rho10,Q1p0,i,p);
		  f0new(1,ix,iy,iz)=f1eq0(rho10,Q1p0,i,p);		  
		}
	  
	    for(int p=0;p<3;p++)
	      for(int i=0;i<4;i++)
		for(int j=0;j<2;j++)
		  { 
		    hnew(j,p,i,ix,iy,iz)=heq(Ep0,B0,j,i,p);
		    h0new(ix,iy,iz)=0;
		  }
	  }
}
void Hartmann::Imprimase(const char* fileName,bool useNew)
{
  ofstream outputFile(fileName);
  outputFile.precision(4);
  vector3D Q00, Q10, E0, B0, J0;					     
  int ix0=0,iz0 =0;
    for(int iy=0;iy<Ly;iy++)
    {
      Q00=Q0(ix0,iy,iz0,useNew);Q10=Q1(ix0,iy,iz0,useNew);	
      E0=E(ix0,iy,iz0,useNew);B0=B(ix0,iy,iz0,useNew);
      outputFile	
	<< iy << " "
	<< Q00.x() << " "
	<< Q00.y() << " "
	<< Q00.z() << " "
	<< Q10.x() << " "
	<< Q10.y() << " "
	<< Q10.z() << " "
	<< E0.x() << " "
	<< E0.y() << " "
	<< E0.z() << " "
	<< B0.x() << " "
	<< B0.y() << " "
	<< B0.z() << "\n";	      
    }
   outputFile.close();
}
int main()
{
  Hartmann hartmann = Hartmann();
  hartmann.ResizeDomain(Lx0,Ly0,Lz0);
  const char * fileName = "testHartmann_first13.dat";
  int t=0, tmax=1000;

  hartmann.Inicie();
  hartmann.Imprimase("time0013.dat",false);
  for(t=0;t<tmax;t++)
    {
      hartmann.Colisione();
      hartmann.ImponerB();
      hartmann.Adveccione();
    }
  hartmann.Imprimase(fileName,true);

  return 0;
}
