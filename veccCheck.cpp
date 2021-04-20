#include "headers/LB_MHD_D3Q19_F.h"

int main(void){
  LatticeBoltzmann LB = LatticeBoltzmann();
  vector3D one; one.cargue(1,1,1);
  cout << fixed;cout.precision(16);
  for(int p=0;p<3;p++)
    for(int i=0;i<4;i++){cout<<"V"<<p<<","<<i<<":\n";
    cout << LB.v[p][i][0] << "\t"
	 << LB.v[p][i][1] << "\t"
	 << LB.v[p][i][2] << "\n";cout<<"e:\n";
    cout << LB.e[0][i][p][0] << "\t"
	 << LB.e[0][i][p][1] << "\t"
	 << LB.e[0][i][p][2] << "|\t"
	 << LB.e[1][i][p][0] << "\t"
	 << LB.e[1][i][p][1] << "\t"
	 << LB.e[1][i][p][2] << "\n";cout<<"b:\n";
    cout << LB.b[0][i][p][0] << "\t"
	 << LB.b[0][i][p][1] << "\t"
	 << LB.b[0][i][p][2] << "|\t"
	 << LB.b[1][i][p][0] << "\t"
	 << LB.b[1][i][p][1] << "\t"
	 << LB.b[1][i][p][2] << "\n\n";      
	 }
  double sumV[3]={0,0,0}, sumE[3]={0,0,0}, sumB[3]={0,0,0}, sumVW[3]={0,0,0}, sumV6[3]={0,0,0};
  double sumW=0;
  double sumVVW[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  double sumVE[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  double sumEE[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  double sumBB[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  double sumVEE[3][3][3]=
    {
      {{0,0,0},{0,0,0},{0,0,0}},
      {{0,0,0},{0,0,0},{0,0,0}},
      {{0,0,0},{0,0,0},{0,0,0}}
    },
    sumVBB[3][3][3]=
    {
      {{0,0,0},{0,0,0},{0,0,0}},
      {{0,0,0},{0,0,0},{0,0,0}},
      {{0,0,0},{0,0,0},{0,0,0}}
    },
    sumVEB[3][3][3]=
    {
      {{0,0,0},{0,0,0},{0,0,0}},
      {{0,0,0},{0,0,0},{0,0,0}},
      {{0,0,0},{0,0,0},{0,0,0}}
    },
     sumVVVW[3][3][3]=
    {
      {{0,0,0},{0,0,0},{0,0,0}},
      {{0,0,0},{0,0,0},{0,0,0}},
      {{0,0,0},{0,0,0},{0,0,0}}
    };
  for(int i=0;i<6;i++)
    {
      sumW += LB.w[i];
    }
  for(int x=0;x<3;x++)
    for(int p=0;p<3;p++)
      for(int i=0;i<6;i++)
	{
	  sumVW[x] += LB.v[p][i][x]*LB.w[i];
	  sumV6[x] += LB.v[p][i][x];
	}
  for(int x=0;x<3;x++)
    for(int y=0;y<3;y++)
    for(int p=0;p<3;p++)
      for(int i=0;i<6;i++)
	{
	  sumVVW[x][y] += LB.v[p][i][x]*LB.v[p][i][y]*LB.w[i];
	}
  for(int x=0;x<3;x++)
    for(int p=0;p<3;p++)
      for(int i=0;i<4;i++)
	for(int j=0;j<2;j++)
	  {
	    sumV[x] += LB.v[p][i][x];
	    sumE[x] += LB.e[j][i][p][x];
	    sumB[x] += LB.b[j][i][p][x];
	  }
  for(int x=0;x<3;x++)
    for(int y=0;y<3;y++)
      for(int p=0;p<3;p++)
	for(int i=0;i<4;i++)
	  for(int j=0;j<2;j++)
	    {
	      sumVE[x][y] +=LB.v[p][i][x]*LB.e[j][i][p][y];
	      sumEE[x][y] +=LB.e[j][i][p][x]*LB.e[j][i][p][y];
	      sumBB[x][y] +=LB.b[j][i][p][x]*LB.b[j][i][p][y];
	    }
  for(int x=0;x<3;x++)
    for(int y=0;y<3;y++)
      for(int z=0;z<3;z++)	
	for(int p=0;p<3;p++)
	  for(int i=0;i<4;i++)
	    for(int j=0;j<2;j++)
	      {
		sumVEE[x][y][z] +=LB.v[p][i][x]*LB.e[j][i][p][y]*LB.e[j][i][p][z];
		sumVEB[x][y][z] +=LB.v[p][i][x]*LB.e[j][i][p][y]*LB.b[j][i][p][z];
		sumVBB[x][y][z] +=LB.v[p][i][x]*LB.b[j][i][p][y]*LB.b[j][i][p][z];
	      }
    for(int x=0;x<3;x++)
    for(int y=0;y<3;y++)
      for(int z=0;z<3;z++)	
	for(int p=0;p<3;p++)
	  for(int i=0;i<6;i++)
	    {
	      sumVVVW[x][y][z] += LB.v[p][i][x]*LB.v[p][i][y]*LB.v[p][i][z]*LB.w[i];
	    }
  cout << fixed;cout.precision(4);
  cout
    << "sum W = " << sumW << "\n";
   cout
    << "sum VW = " << sumVW[0]<<","<<sumVW[1]<<","<<sumVW[2]<<"\n"
    << "sum V6 = " << sumV6[0]<<","<<sumV6[1]<<","<<sumV6[2]<<"\n"
    << "sum V = " << sumV[0]<<","<<sumV[1]<<","<<sumV[2]<<"\n"
    << "sum B = " << sumB[0]<<","<<sumB[1]<<","<<sumB[2]<<"\n"
    << "sum E = " << sumE[0]<<","<<sumE[1]<<","<<sumE[2]<<"\n\n";
  cout
    << "sum VVW =\n"
    << sumVVW[0][0]<<","<<sumVVW[0][1]<<","<<sumVVW[0][2]<<"\n"
    << sumVVW[1][0]<<","<<sumVVW[1][1]<<","<<sumVVW[1][2]<<"\n"
    << sumVVW[2][0]<<","<<sumVVW[2][1]<<","<<sumVVW[2][2]<<"\n\n";
  cout
    << "sum VE =\n"
    << sumVE[0][0]<<","<<sumVE[0][1]<<","<<sumVE[0][2]<<"\n"
    << sumVE[1][0]<<","<<sumVE[1][1]<<","<<sumVE[1][2]<<"\n"
    << sumVE[2][0]<<","<<sumVE[2][1]<<","<<sumVE[2][2]<<"\n\n";
  cout
    << "sum EE =\n"
    << sumEE[0][0]<<","<<sumEE[0][1]<<","<<sumEE[0][2]<<"\n"
    << sumEE[1][0]<<","<<sumEE[1][1]<<","<<sumEE[1][2]<<"\n"
    << sumEE[2][0]<<","<<sumEE[2][1]<<","<<sumEE[2][2]<<"\n\n";
  cout
    << "sum BB =\n"
    << sumBB[0][0]<<","<<sumBB[0][1]<<","<<sumBB[0][2]<<"\n"
    << sumBB[1][0]<<","<<sumBB[1][1]<<","<<sumBB[1][2]<<"\n"
    << sumBB[2][0]<<","<<sumBB[2][1]<<","<<sumBB[2][2]<<"\n\n";
  
  cout
    << "sum VEE =\n"
    << sumVEE[0][0][0]<<","<<sumVEE[0][1][0]<<","<<sumVEE[0][2][0]<<"\n"
    << sumVEE[1][0][0]<<","<<sumVEE[1][1][0]<<","<<sumVEE[1][2][0]<<"\n"
    << sumVEE[2][0][0]<<","<<sumVEE[2][1][0]<<","<<sumVEE[2][2][0]<<"\n\n"
    
    << sumVEE[0][0][1]<<","<<sumVEE[0][1][1]<<","<<sumVEE[0][2][1]<<"\n"
    << sumVEE[1][0][1]<<","<<sumVEE[1][1][1]<<","<<sumVEE[1][2][1]<<"\n"
    << sumVEE[2][0][1]<<","<<sumVEE[2][1][1]<<","<<sumVEE[2][2][1]<<"\n\n"
    
    << sumVEE[0][0][2]<<","<<sumVEE[0][1][2]<<","<<sumVEE[0][2][2]<<"\n"
    << sumVEE[1][0][2]<<","<<sumVEE[1][1][2]<<","<<sumVEE[1][2][2]<<"\n"
    << sumVEE[2][0][2]<<","<<sumVEE[2][1][2]<<","<<sumVEE[2][2][2]<<"\n\n";

  cout
    << "sum VEB =\n"
    << sumVEB[0][0][0]<<","<<sumVEB[0][1][0]<<","<<sumVEB[0][2][0]<<"\n"
    << sumVEB[1][0][0]<<","<<sumVEB[1][1][0]<<","<<sumVEB[1][2][0]<<"\n"
    << sumVEB[2][0][0]<<","<<sumVEB[2][1][0]<<","<<sumVEB[2][2][0]<<"\n\n"
    
    << sumVEB[0][0][1]<<","<<sumVEB[0][1][1]<<","<<sumVEB[0][2][1]<<"\n"
    << sumVEB[1][0][1]<<","<<sumVEB[1][1][1]<<","<<sumVEB[1][2][1]<<"\n"
    << sumVEB[2][0][1]<<","<<sumVEB[2][1][1]<<","<<sumVEB[2][2][1]<<"\n\n"
    
    << sumVEB[0][0][2]<<","<<sumVEB[0][1][2]<<","<<sumVEB[0][2][2]<<"\n"
    << sumVEB[1][0][2]<<","<<sumVEB[1][1][2]<<","<<sumVEB[1][2][2]<<"\n"
    << sumVEB[2][0][2]<<","<<sumVEB[2][1][2]<<","<<sumVEB[2][2][2]<<"\n\n";
  cout
    << "sum VBB =\n"
    << sumVBB[0][0][0]<<","<<sumVBB[0][1][0]<<","<<sumVBB[0][2][0]<<"\n"
    << sumVBB[1][0][0]<<","<<sumVBB[1][1][0]<<","<<sumVBB[1][2][0]<<"\n"
    << sumVBB[2][0][0]<<","<<sumVBB[2][1][0]<<","<<sumVBB[2][2][0]<<"\n\n"
    
    << sumVBB[0][0][1]<<","<<sumVBB[0][1][1]<<","<<sumVBB[0][2][1]<<"\n"
    << sumVBB[1][0][1]<<","<<sumVBB[1][1][1]<<","<<sumVBB[1][2][1]<<"\n"
    << sumVBB[2][0][1]<<","<<sumVBB[2][1][1]<<","<<sumVBB[2][2][1]<<"\n\n"
    
    << sumVBB[0][0][2]<<","<<sumVBB[0][1][2]<<","<<sumVBB[0][2][2]<<"\n"
    << sumVBB[1][0][2]<<","<<sumVBB[1][1][2]<<","<<sumVBB[1][2][2]<<"\n"
    << sumVBB[2][0][2]<<","<<sumVBB[2][1][2]<<","<<sumVBB[2][2][2]<<"\n\n";
  cout
    << "sum VVVW =\n"
    << sumVVVW[0][0][0]<<","<<sumVVVW[0][1][0]<<","<<sumVVVW[0][2][0]<<"\n"
    << sumVVVW[1][0][0]<<","<<sumVVVW[1][1][0]<<","<<sumVVVW[1][2][0]<<"\n"
    << sumVVVW[2][0][0]<<","<<sumVVVW[2][1][0]<<","<<sumVVVW[2][2][0]<<"\n\n"
    
    << sumVVVW[0][0][1]<<","<<sumVVVW[0][1][1]<<","<<sumVVVW[0][2][1]<<"\n"
    << sumVVVW[1][0][1]<<","<<sumVVVW[1][1][1]<<","<<sumVVVW[1][2][1]<<"\n"
    << sumVVVW[2][0][1]<<","<<sumVVVW[2][1][1]<<","<<sumVVVW[2][2][1]<<"\n\n"
    
    << sumVVVW[0][0][2]<<","<<sumVVVW[0][1][2]<<","<<sumVVVW[0][2][2]<<"\n"
    << sumVVVW[1][0][2]<<","<<sumVVVW[1][1][2]<<","<<sumVVVW[1][2][2]<<"\n"
    << sumVVVW[2][0][2]<<","<<sumVVVW[2][1][2]<<","<<sumVVVW[2][2][2]<<"\n\n";
  

  return 0;
}
