#include "crystal.h"

#include <stdlib.h>
#include <string>
#include <cstring>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <stdio.h>

using namespace std;
void input(ifstream &in, DATA vt[], DATA2 bt[][MOL2], LOG& logdt, int typejudge[]){
		
	int N=0; 
	double u,w,w1,m1,n1,m2,n2;
	string m,n,e; 
	in>>m>>n;
	in>>u;
	in>>m>>n>>m>>n;
	in>>N;
	
	int ty[N+2];
	int id[N+2];
	int molid[N+5];	
	double x[N+2];
	double y[N+2];
	double z[N+2];
	int X[N+2]; 
	int Y[N+2]; 
	int Z[N+2]; 
	
	in>>m>>n>>m>>n>>m>>n>>m>>n>>m;
	in>>w>>w1>>e;
	in>>m1>>n1>>e;
	in>>m2>>n2>>e;
	in>>m>>n>>m>>n>>m>>n>>m>>n>>m>>n>>m;
	logdt.box_xlo=w;
	logdt.box_xhi=w1;
	logdt.box_ylo=m1;
	logdt.box_yhi=n1;
	logdt.box_zlo=m2;
	logdt.box_zhi=n2;
	logdt.box_x=w1-w;
	logdt.box_y=n1-m1;
	logdt.box_z=n2-m2;
	molid[-1]=-1;
	molid[-2]=-1;
	molid[N]=-1;
	ty[-1]=-1;
	logdt.Ncen=0;
	logdt.Mol=-1;
	
	for(int i=0;i<N;++i){
		in>>id[i]>>molid[i]>>ty[i]>>x[i]>>y[i]>>z[i]>>X[i]>>Y[i]>>Z[i];
		if(typejudge[ty[i]]==0||typejudge[ty[i-1]]==0){
			if(molid[i-1]==molid[i]){
				if(molid[i-2]==molid[i]){
					vt[logdt.Ncen].cenid=id[i-1];
					vt[logdt.Ncen].ty=ty[i-1];
					vt[logdt.Ncen].cenx=x[i-1];
					vt[logdt.Ncen].ceny=y[i-1];
					vt[logdt.Ncen].cenz=z[i-1];
					vt[logdt.Ncen].cenix=X[i-1];
					vt[logdt.Ncen].ceniy=Y[i-1];
					vt[logdt.Ncen].ceniz=Z[i-1];
					vt[logdt.Ncen].x=x[i]+X[i]*logdt.box_x-x[i-2]-X[i-2]*logdt.box_x;
					vt[logdt.Ncen].y=y[i]+Y[i]*logdt.box_y-y[i-2]-Y[i-2]*logdt.box_y;
					vt[logdt.Ncen].z=z[i]+Z[i]*logdt.box_z-z[i-2]-Z[i-2]*logdt.box_z;
					vt[logdt.Ncen].mol=molid[i-1];
					
					logdt.Ncen++;
				}if(molid[i-2]!=molid[i]){
					vt[logdt.Ncen].cenid=id[i-1];
					vt[logdt.Ncen].ty=ty[i-1];
					vt[logdt.Ncen].cenx=x[i-1];
					vt[logdt.Ncen].ceny=y[i-1];
					vt[logdt.Ncen].cenz=z[i-1];
					vt[logdt.Ncen].cenix=X[i-1];
					vt[logdt.Ncen].ceniy=Y[i-1];
					vt[logdt.Ncen].ceniz=Z[i-1];
					vt[logdt.Ncen].x=x[i-1]+X[i-1]*logdt.box_x-x[i]-X[i]*logdt.box_x;
					vt[logdt.Ncen].y=y[i-1]+Y[i-1]*logdt.box_y-y[i]-Y[i]*logdt.box_y;
					vt[logdt.Ncen].z=z[i-1]+Z[i-1]*logdt.box_z-z[i]-Z[i]*logdt.box_z;
					vt[logdt.Ncen].mol=molid[i-1];
					
					logdt.Ncen++;
				}
				
				if(i==N-1){
					vt[logdt.Ncen].cenid=id[i];
					vt[logdt.Ncen].ty=ty[i];
					vt[logdt.Ncen].cenx=x[i];
					vt[logdt.Ncen].ceny=y[i];
					vt[logdt.Ncen].cenz=z[i];
					vt[logdt.Ncen].cenix=X[i];
					vt[logdt.Ncen].ceniy=Y[i];
					vt[logdt.Ncen].ceniz=Z[i];
					vt[logdt.Ncen].x=x[i]+X[i]*logdt.box_x-x[i-1]-X[i-1]*logdt.box_x;
					vt[logdt.Ncen].y=y[i]+Y[i]*logdt.box_y-y[i-1]-Y[i-1]*logdt.box_y;
					vt[logdt.Ncen].z=z[i]+Z[i]*logdt.box_z-z[i-1]-Z[i-1]*logdt.box_z;
					vt[logdt.Ncen].mol=molid[i];
					
					logdt.Ncen++;
				}
			}
			if(molid[i-1]!=molid[i]){
				logdt.Mol++;
				logdt.Nmol[logdt.Mol]=-1;
				
				if(i>0){
					vt[logdt.Ncen].cenid=id[i-1];
					vt[logdt.Ncen].ty=ty[i-1];
					vt[logdt.Ncen].cenx=x[i-1];
					vt[logdt.Ncen].ceny=y[i-1];
					vt[logdt.Ncen].cenz=z[i-1];
					vt[logdt.Ncen].x=x[i-1]+X[i-1]*logdt.box_x-x[i-2]-X[i-2]*logdt.box_x;
					vt[logdt.Ncen].y=y[i-1]+Y[i-1]*logdt.box_y-y[i-2]-Y[i-2]*logdt.box_y;
					vt[logdt.Ncen].z=z[i-1]+Z[i-1]*logdt.box_z-z[i-2]-Z[i-2]*logdt.box_z;
					vt[logdt.Ncen].mol=molid[i-1];
					
					logdt.Ncen++;
				}
				
			}
			logdt.Nmol[logdt.Mol]++;
			bt[logdt.Mol][logdt.Nmol[logdt.Mol]].x=x[i];
			bt[logdt.Mol][logdt.Nmol[logdt.Mol]].y=y[i];
			bt[logdt.Mol][logdt.Nmol[logdt.Mol]].z=z[i];
			bt[logdt.Mol][logdt.Nmol[logdt.Mol]].ix=X[i];
			bt[logdt.Mol][logdt.Nmol[logdt.Mol]].iy=Y[i];
			bt[logdt.Mol][logdt.Nmol[logdt.Mol]].iz=Z[i];
			bt[logdt.Mol][logdt.Nmol[logdt.Mol]].realx=x[i]+X[i]*logdt.box_x;
			bt[logdt.Mol][logdt.Nmol[logdt.Mol]].realy=y[i]+Y[i]*logdt.box_y;
			bt[logdt.Mol][logdt.Nmol[logdt.Mol]].realz=z[i]+Z[i]*logdt.box_z;
			bt[logdt.Mol][logdt.Nmol[logdt.Mol]].ty=ty[i];
			bt[logdt.Mol][logdt.Nmol[logdt.Mol]].dttid=0;
		}
	}
  vt[logdt.Ncen].mol=0;

}