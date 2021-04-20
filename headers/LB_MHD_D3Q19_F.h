/*
Reproduccion LBM para MHD, primer intento
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include "Distribution_EM.h"
#include "Distribution_fluids.h"
#include "Vector.h"
//using namespace std;
/*
const int Lx=1;
const int Ly=1;
const int Lz=100;
*/
//Constates e0 y mu0 (INTOCABLE).
const double C=1.0;
const double C2 = 1.0;
const double eps0=1.0;
const double mu0=1.0;
const double nu=100;

const double tau_s=1.0;
const double tau_2=0.5;
const double tau_k=1.0;

const double gmma=1.0;

const double UtauS=1.0/tau_s;
const double UmUtauS=1-UtauS;
const double Utau2=1.0/tau_2;
const double UmUtau2=1-Utau2;
const double Utauk=1.0/tau_k;
const double UmUtauk=1-Utauk;

const double Ks = 1-0.5*UtauS;
const double K2 = 1-0.5*Utau2;
	       
class LatticeBoltzmann{
public:
  int Lx=1,Ly=1,Lz=1;
  double V0[3]={0.0,0.0,0.0};
  double e0[3]={0.0,0.0,0.0};
  double b0[3]={0.0,0.0,0.0};
  
  double V[3][6][3]; vector3D v[3][6],v0;// V[0][i][p]=V^p_ix,  V[1][i][p]=V^p_iy, V[2][i][p]=V^p_iz
  vector3D e[2][4][3]; //e[0][j][i][p]=e^p_ijx,e[1][j][i][p]=e^p_ijy,e[2][j][i][p]=e^p_ijz,
  vector3D b[2][4][3]; //b[0][j][i][p]=b^p_ijx,b[1][j][i][p]=b^p_ijy,b[2][j][i][p]=b^p_ijz,
  double w0=(1.0/(6*C2));
  double w[6] = {0,0,0,0,0,0};
  double xi[2] = {0,0};

  fDistribution f= fDistribution(Lx,Ly,Lz,false);		fDistribution fnew = fDistribution(Lx,Ly,Lz,false);
  fDistribution f0 = fDistribution(Lx,Ly,Lz,true);	fDistribution f0new = fDistribution(Lx,Ly,Lz,true);
  hDistribution h= hDistribution(Lx,Ly,Lz,false);		hDistribution hnew = hDistribution(Lx,Ly,Lz,false);
  hDistribution h0 = hDistribution(Lx,Ly,Lz,true);	hDistribution h0new = hDistribution(Lx,Ly,Lz,true);
  
  //double f[Lx][Ly][Lz][2][2][4][3], fnew[Lx][Ly][Lz][2][2][4][3]; // f[ix][iy][iz][r][j][i][p]
  //double f0[Lx][Ly][Lz][2], f0new[Lx][Ly][Lz][2]; // f[ix][iy][iz][r]
  //public:
  vector3D Fext0; vector3D Fext1;
  LatticeBoltzmann(void);
  void ResizeDomain(int Lx0, int Ly0, int Lz0);
  //fluidos
  double q[2]={1.0,1.0};	double m[2]={1.0,1.0};	double n[2]={1.0,1.0};
  double rho0(int ix, int iy, int iz, bool useNew);
  vector3D Q0(int ix, int iy, int iz, bool useNew);
  double rho1(int ix, int iy, int iz, bool useNew); 
  vector3D Q1(int ix, int iy, int iz, bool useNew);
  double rhoc(void); 
  //campos
  vector3D E(int ix, int iy, int iz, bool useNew);
  vector3D B(int ix, int iy, int iz, bool useNew);
  //Forzamientos
  vector3D F0(vector3D&E0,vector3D&B0,vector3D&Q0,vector3D&Q1,double&rho0,double&rho1);
  vector3D F1(vector3D&E0,vector3D&B0,vector3D&Q0,vector3D&Q1,double&rho0,double&rho1);    
  //campos auxiliares
  vector3D J(vector3D&Q0, vector3D&Q1);
  vector3D Jp(vector3D&F0, vector3D&F1, vector3D&J0);
  vector3D Ep(vector3D&E0, vector3D&Jp0);
  vector3D Q0p(vector3D&Q0,vector3D&F0);
  vector3D Q1p(vector3D&Q1,vector3D&F1);
  //constantes dielectricas relativas
  double epsr(int ix, int iy, int iz){return 1.0;}; 
  double mur(int ix, int iy, int iz){return 1.0;};
  double sigma(int ix, int iy, int iz){return 0.0;};
  //funciones de equilibrio
  double f0eq(double&rho0,vector3D&Q0p,int i,int p);
  double f0eq0(double&rho0,vector3D&Q0p,int i,int p);
  double f1eq(double&rho1,vector3D&Q1p,int i,int p);
  double f1eq0(double&rho1,vector3D&Q1p,int i,int p);
  double heq(vector3D&Ep0,vector3D&B0,int j, int i, int p);
  void Colisione(void);
  void Adveccione(void);
  //void Inicie(void);
  //void ImponerCampos(double&rho0,double*&D0,double*&B0,double*&H0,double*&E0,double*&J0,double*&Jp0,double*&Ep0,int t);
  //void Imprimase(const char* fileName,bool useNew);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //double C2=C*C;
  //w[0]=1.0/(6*C2);
  w[0]=1.0/(72*C2);
  w[1]=1.0/(72*C2);
  w[2]=1.0/(72*C2);
  w[3]=1.0/(72*C2);
  w[4]=1.0/(36*C2);
  w[5]=1.0/(36*C2);
  xi[0]=1/3;xi[1]=1/3;
  m[0]=1;m[1]=1;
  q[0]=1;q[1]=1;
  n[0]=1;n[1]=1;
  Fext0.cargue(0,0,0);	Fext1.cargue(0,0,0);
  double c12=C*sqrt(2.0), csq2=0.5;
  int imp=0,rec=0,rec1=0,rec2=0;
   
  V[0][0][0]=1;	V[1][0][0]=1;	V[2][0][0]=0;    //Plano XY p=0
  V[0][1][0]=-1;V[1][1][0]=1;	V[2][1][0]=0;	 //Plano XY p=0
  V[0][2][0]=-1;V[1][2][0]=-1;	V[2][2][0]=0;	 //Plano XY p=0
  V[0][3][0]=1;	V[1][3][0]=-1;	V[2][3][0]=0;	 //Plano XY p=0

  V[0][0][1]=1; V[1][0][1]=0; 	V[2][0][1]=1;	//Plano XZ p=1
  V[0][1][1]=-1;V[1][1][1]=0; 	V[2][1][1]=1;	//Plano XZ p=1
  V[0][2][1]=-1;V[1][2][1]=0; 	V[2][2][1]=-1;	//Plano XZ p=1
  V[0][3][1]=1; V[1][3][1]=0; 	V[2][3][1]=-1;	//Plano XZ p=1

  V[0][0][2]=0;	V[1][0][2]=1;	V[2][0][2]=1;	//Plano YZ p=2
  V[0][1][2]=0;	V[1][1][2]=-1;	V[2][1][2]=1;	//Plano YZ p=2
  V[0][2][2]=0;	V[1][2][2]=-1;	V[2][2][2]=-1;	//Plano YZ p=2
  V[0][3][2]=0;	V[1][3][2]=1;	V[2][3][2]=-1;	//Plano YZ p=2

  V[0][4][0]=-1;V[1][4][0]=0;	V[2][4][0]=0;	//Plano XY p=0
  V[0][5][0]=1;	V[1][5][0]=0;	V[2][5][0]=0;	//Plano XY p=0
  
  V[0][4][1]=0;	V[1][4][1]=-1;	V[2][4][1]=0;	//Plano XZ p=1
  V[0][5][1]=0;	V[1][5][1]=1;	V[2][5][1]=0;	//Plano XZ p=1
  
  V[0][4][2]=0;	V[1][4][2]=0;	V[2][4][2]=-1;	//Plano YZ p=2
  V[0][5][2]=0;	V[1][5][2]=0;	V[2][5][2]=1;	//Plano YZ p=2
 
  
  //Los vectores Velocidad V[p][i]=V^p_i como vectores
  v0.cargue(V0[0],V0[1],V0[2]);
  for(int p=0;p<3;p++)
    for(int i=0;i<6;i++){
      if(i<=3)
	v[p][i].cargue(c12*V[0][i][p],c12*V[1][i][p],c12*V[2][i][p]);
      else
	v[p][i].cargue(c12*V[0][i][p],c12*V[1][i][p],c12*V[2][i][p]);
  }
  
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++)
      {
	rec=i%4; rec1=(i+1)%4; rec2=(i+2)%4;
	e[0][rec1][p]=v[p][rec2]*0.5;
	e[1][rec1][p]=v[p][rec]*0.5;
      }
  for(int i=0;i<4;i++)
    {
      b[0][i][0]=v[0][i]^e[0][i][0]*csq2;		b[1][i][0]=v[0][i]^e[1][i][0]*csq2;

      b[0][i][1]=v[1][i]^e[0][i][1]*csq2;       	b[1][i][1]=v[1][i]^e[1][i][1]*csq2;
      
      b[0][i][2]=v[2][i]^e[0][i][2]*csq2;		b[1][i][2]=v[2][i]^e[1][i][2]*csq2;
      /* b[0][i][0].cargue(0,0,csq2);		b[1][i][0].cargue(0,0,-csq2);

      b[0][i][1].cargue(0,-csq2,0);       	b[1][i][1].cargue(0,csq2,0);
      
      b[0][i][2].cargue(csq2,0,0);		b[1][i][2].cargue(-csq2,0,0);*/
    }
 
}
void LatticeBoltzmann::ResizeDomain(int Lx0, int Ly0, int Lz0)
{
  Lx=Lx0;	Ly=Ly0;		Lz=Lz0;
  f = fDistribution(Lx,Ly,Lz,false);	fnew = fDistribution(Lx,Ly,Lz,false);
  f0 = fDistribution(Lx,Ly,Lz,true);	f0new = fDistribution(Lx,Ly,Lz,true);
  h = hDistribution(Lx,Ly,Lz,false);	hnew = hDistribution(Lx,Ly,Lz,false);
  h0 = hDistribution(Lx,Ly,Lz,true);	h0new = hDistribution(Lx,Ly,Lz,true);
}
double LatticeBoltzmann::rho0(int ix, int iy, int iz, bool useNew)
{
  double sum=0;
  
  for(int i=0;i<6;i++)
    for(int p=0;p<3;p++){	
      sum += useNew? fnew.function(0,p,i,ix,iy,iz) : f.function(0,p,i,ix,iy,iz);}
	sum += f0(0,ix,iy,iz);
  return sum;
}
vector3D LatticeBoltzmann::Q0(int ix, int iy, int iz, bool useNew)
{
  vector3D sum; sum.cargue(0,0,0);
  
    for(int i=0;i<6;i++)
      for(int p=0;p<3;p++)	
	  sum += useNew? fnew.function(0,p,i,ix,iy,iz)*v[p][i] : f.function(0,p,i,ix,iy,iz)*v[p][i];
  return sum;
}
double LatticeBoltzmann::rho1(int ix, int iy, int iz, bool useNew)
{
  double sum=0;
  
  for(int i=0;i<6;i++)
    for(int p=0;p<3;p++){	
      sum += useNew? fnew.function(1,p,i,ix,iy,iz) : f.function(1,p,i,ix,iy,iz);}
	sum += f0(1,ix,iy,iz);
  return sum;
}
vector3D LatticeBoltzmann::Q1(int ix, int iy, int iz, bool useNew)
{
  vector3D sum; sum.cargue(0,0,0);
  
    for(int i=0;i<6;i++)
      for(int p=0;p<3;p++)	
	  sum += useNew? fnew.function(1,p,i,ix,iy,iz)*v[p][i] : f.function(1,p,i,ix,iy,iz)*v[p][i];
  return sum;
}
double LatticeBoltzmann::rhoc()
{
  return (q[0]/m[0] + q[1]/m[1]);
}
vector3D LatticeBoltzmann::E(int ix, int iy, int iz, bool useNew)
{
  vector3D sum; sum.cargue(0,0,0);
  
    for(int i=0;i<4;i++)
      for(int p=0;p<3;p++)
      	for(int j=0;j<2;j++)	
	  sum += useNew? hnew.function(j,p,i,ix,iy,iz)*e[j][i][p] : h.function(j,p,i,ix,iy,iz)*e[j][i][p];
  return sum;
  
}
vector3D LatticeBoltzmann::B(int ix, int iy, int iz, bool useNew)
{
  vector3D sum; sum.cargue(0,0,0);
  
    for(int i=0;i<4;i++)
      for(int p=0;p<3;p++)      
	for(int j=0;j<2;j++)	
	  sum += useNew? hnew.function(j,p,i,ix,iy,iz)*b[j][i][p] : h.function(j,p,i,ix,iy,iz)*b[j][i][p];
  return sum;
  
}
vector3D LatticeBoltzmann::J(vector3D&Q0, vector3D&Q1)
{
  vector3D sum; sum.cargue(0,0,0);
  sum = Q0*(q[0]/m[0])+Q1*(q[1]/m[1]);
  return sum; 
}
vector3D LatticeBoltzmann::F0(vector3D&E0,vector3D&B0,vector3D&Q0,vector3D&Q1,double&rho0,double&rho1)
{
  vector3D sum;	sum.cargue(0,0,0);
  if(rho1==0){return sum;};
  double qm=q[0]/m[0], rho01=rho0/rho1;
  vector3D VxB=Q0^B0,V0V1=Q0-(rho01*Q1);
  sum = qm*(rho0*E0 + VxB) - nu*V0V1 + Fext0;
  return sum;
}
vector3D LatticeBoltzmann::F1(vector3D&E0,vector3D&B0,vector3D&Q0,vector3D&Q1,double&rho0,double&rho1)
{
  vector3D sum;	sum.cargue(0,0,0);
  if(rho0==0){return sum;};
  double qm=q[1]/m[1], rho10=rho1/rho0;
  vector3D VxB=Q1^B0,V1V0=Q1-(rho10*Q0);
  sum = qm*(rho1*E0 + VxB) - nu*V1V0 + Fext1;
  return sum;
}
vector3D LatticeBoltzmann::Jp(vector3D&F0, vector3D&F1, vector3D&J0)
{
  double qm0=q[0]/m[0], qm1=q[1]/m[1];
  vector3D sum; sum.cargue(0,0,0);

  sum = J0+0.5*(qm0*F0+qm1*F1);
  return sum; 
}
vector3D LatticeBoltzmann::Ep(vector3D&E0,vector3D&Jp0)
{
  vector3D sum; sum.cargue(0,0,0);

  sum = E0-Jp0*mu0*C2*0.5;
  return sum;
}
vector3D LatticeBoltzmann::Q0p(vector3D&Q0,vector3D&F0)
{
  vector3D sum; sum.cargue(0,0,0);
  sum = Q0 + 0.5*F0;
  return sum;
}
vector3D LatticeBoltzmann::Q1p(vector3D&Q1,vector3D&F1)
{
  vector3D sum; sum.cargue(0,0,0);
  sum = Q1 + 0.5*F1;
  return sum;
}

double LatticeBoltzmann::f0eq(double&rho0,vector3D&Q0p,int i,int p)
{
  double f=0, vQp=0, Q02=0,rho_g=pow(rho0,gmma-1);
  vQp = (v[p][i]*Q0p)/rho0;
  Q02 = norma2(Q0p)/(rho0*rho0);
  
  //f = 3*w[i]*(xi[0]*rho_g + vQp + (0.75/(C*C*rho0))*vQp*vQp - (0.5/rho0)*Q02);
  f=w[i]*rho0*(3*xi[0]*rho_g+3*vQp+(2.25/C2)*vQp*vQp-1.5*Q02);
  return f;  
}

double LatticeBoltzmann::f0eq0(double&rho0,vector3D&Q0p,int i,int p)
{
  double f=0, Q02=0,rho_g=pow(rho0,gmma-1),rho02=rho0;
  Q02 = norma2(Q0p)/(rho0*rho0);
  
  //f = 6*rho0*(C*C - xi[0]*rho_g + (0.25/rho02)*Q02);
  //f=6*rho0*C2*(1-(0.25/C2)*(4*xi[0]*rho_g+Q02));
  f=rho0*(1-(0.25/C2)*(4*xi[0]*rho_g+Q02));
  return f;  
}

double LatticeBoltzmann::f1eq(double&rho1,vector3D&Q1p,int i,int p)
{
  double f=0, vQp=0, Q12=0,rho_g=pow(rho1,gmma-1);
  vQp = (v[p][i]*Q1p)/rho1;
  Q12 = norma2(Q1p)/(rho1*rho1);
  
  //f = 3*w[i]*(xi[1]*rho_g + vQp + rho1*(0.75/(C*C*rho1))*vQp*vQp - (0.5/rho1)*Q12);
  f=w[i]*rho1*(3*xi[1]*rho_g+3*vQp+(2.25/C2)*vQp*vQp-1.5*Q12);
  return f;  
}


double LatticeBoltzmann::f1eq0(double&rho1,vector3D&Q1p,int i,int p)
{
  double f=0, Q12=0,rho_g=pow(rho1,gmma-1),rho12=rho1;
  Q12 = norma2(Q1p)/(rho1*rho1);
  
  //f = 6*rho1*(C*C - xi[1]*rho_g + (0.25/rho12)*Q12);
  //f=6*rho1*C2*(1-(0.25/C2)*(4*xi[1]*rho_g+Q12));
  f=rho1*(1-(0.25/C2)*(4*xi[1]*rho_g+Q12));
  return f;  
}


double LatticeBoltzmann::heq(vector3D&Ep0,vector3D&B0,int j,int i,int p)
{
  double f=0, eEp=0, bB=0;
  eEp = (e[j][i][p]*Ep0);	bB = (b[j][i][p]*B0);
  
  //f =  (0.125/C*C)*eEp + 0.125*bB;
  f =  (0.125/C2)*eEp + 0.125*bB;
  
  return f;  
}
void LatticeBoltzmann::Colisione(void)
{
  int ix=0,iy=0,iz=0,s=0,j=0,i=0,p=0;
  double rho00=0,rho10=0,rhoc0=0;
  vector3D Q00,Q10,E0,B0,F00,F10,J0,Jp0,Ep0,Q0p0,Q1p0;
  double vF0 = 0, vF1 = 0, eJ = 0, forz0 = 0, forz1 = 0, forz2 = 0;
  
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
	{
	  rho00=rho0(ix,iy,iz,false);	rho10=rho1(ix,iy,iz,false);	rhoc0=rhoc();
	  Q00=Q0(ix,iy,iz,false);	Q10=Q1(ix,iy,iz,false);		J0=J(Q00,Q10);
	  E0=E(ix,iy,iz,false);		B0=B(ix,iy,iz,false);
	  F00=F0(E0,B0,Q00,Q10,rho00,rho10);
	  F10=F0(E0,B0,Q00,Q10,rho00,rho10);	
	  Jp0=Jp(F00,F10,J0);		Ep0=Ep(E0,Jp0);
	  Q0p0=Q0p(Q00,F00);		Q1p0=Q1p(Q10,F10);

	  for(p=0;p<3;p++)
	    for(i=0;i<6;i++)
		{
		  vF0=v[p][i]*F00;	vF1=v[p][i]*F10;
		  forz0=Ks*vF0/(20*C2);	forz1=Ks*vF1/(20*C2);
		  fnew(0,p,i,ix,iy,iz)=UmUtauS*f.function(0,p,i,ix,iy,iz)+UtauS*f0eq(rho00,Q0p0,i,p)+forz0;
		  f0new(0,ix,iy,iz)=UmUtauk*f0.function(0,ix,iy,iz)+Utauk*f0eq0(rho00,Q0p0,i,p);		  
		  fnew(1,p,i,ix,iy,iz)=UmUtauS*f.function(1,p,i,ix,iy,iz)+UtauS*f1eq(rho10,Q1p0,i,p)+forz1;
		  f0new(1,ix,iy,iz)=UmUtauk*f0.function(1,ix,iy,iz)+Utauk*f1eq0(rho10,Q1p0,i,p);		  
		}
	  
	  for(p=0;p<3;p++)
	    for(i=0;i<4;i++)
	      for(j=0;j<2;j++)
		{
		  eJ=e[j][i][p]*Jp0;	forz2=K2*mu0*0.125*eJ;
		  hnew(j,p,i,ix,iy,iz)=UmUtau2*h.function(j,p,i,ix,iy,iz)+Utau2*heq(Ep0,B0,j,i,p)-forz2;
		  h0new(ix,iy,iz)=UmUtau2*h.function(ix,iy,iz);
		}
	}
}
void LatticeBoltzmann::Adveccione(void)
{
  int jx=0,jy=0,jz=0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	{
	  for(int p=0;p<3;p++)
	    for(int i=0;i<6;i++)
	      {
		jx=(ix+(int)v[p][i][0]+Lx)%Lx;jy=(iy+(int)v[p][i][1]+Ly)%Ly;jz=(iz+(int)v[p][i][2]+Lz)%Lx;	    
		f(0,p,i,jx,jy,jz) = fnew.function(0,p,i,ix,iy,iz);
		f0(0,ix,iy,iz) = f0new.function(0,ix,iy,iz);

		f(1,p,i,jx,jy,jz) = fnew.function(1,p,i,ix,iy,iz);
		f0(1,ix,iy,iz) = f0new.function(1,ix,iy,iz);
	      }
	  for(int p=0;p<3;p++)
	    for(int i=0;i<4;i++)
	      for(int j=0;j<2;j++)
		{
		  jx=(ix+(int)v[p][i][0]+Lx)%Lx;jy=(iy+(int)v[p][i][1]+Ly)%Ly;jz=(iz+(int)v[p][i][2]+Lz)%Lx;
		  h(j,p,i,jx,jy,jz) = hnew.function(j,p,i,ix,iy,iz);
		  h0(ix,iy,iz) = hnew.function(ix,iy,iz);
		}
	}
}/*
void LatticeBoltzmann::Inicie(void)
{
  double rho00=0,rho10=0,rhoc=0,H=5;
  double eta=
  vector3D Q00,Q10,E0,B0,F00,F10,J0,Jp00,Ep00,Q0p0,Q1p0;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
	  {
	    rho00=m[0]*n[0];	rho10=m[1]*n[1];	rhoc0=rhoc();
	    Q00.cargue(0,0,0);	Q10.cargue(0,0,0);	J0.cargue(0,0,0);
	    E0.cargue(0,0,0);		B0.cargue(0,Bo,0);
	    F00=F0(E0,B0,Q00,Q10,rho00,rho10);
	    F10=F0(E0,B0,Q00,Q10,rho00,rho10);	
	    Jp0=Jp(F00,F10,J0);		Ep0=Ep(E0,Jp0);
	    Q0p0=Q0p(Q00,F00);		Q1p0=Q1p(Q10,F10);

	    for(p=0;p<3;p++)
	      for(i=0;i<6;i++)
		{ 
		  f(0,p,i,ix,iy,iz)=f0eq(rho0,Q0p0,i,p);
		  f0(0,ix,iy,iz)=f0eq0(rho0,Q0p0,i,p);		  
		  f(1,p,i,ix,iy,iz)=f1eq(rho1,Q1p0,i,p);
		  f0(1,ix,iy,iz)=f1eq0(rho1,Q1p0,i,p);		  
		}
	  
	    for(p=0;p<3;p++)
	      for(i=0;i<4;i++)
		for(j=0;j<2;j++)
		  { 
		    h(j,p,i,ix,iy,iz)=heq(Ep0,B0,r,j,i,p);
		    h0(ix,iy,iz)=0;
		  }
	  }
}

void Hartmann::Imprimase(const char* fileName,bool useNew)
{
  ofstream outputFile(fileName);
  vector3D Q00, B00;					     
  int ix0=0,iz0 =0;
    for(int iy=0;iy<Ly;iy++)
    {
      Q00=Q0(ix,iy,iz,useNew);	
      B00=B0(ix,iy,iz,useNew);
      outputFile	
	<< iy << "\t"
	<< Q00.x() << "\t"
	<< B00.x() << "\n";	      
    }
   outputFile.close();
   }*/
