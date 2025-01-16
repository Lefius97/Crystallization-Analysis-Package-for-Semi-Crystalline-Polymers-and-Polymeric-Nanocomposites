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

#define AT 1000000

using namespace std;
void cluster_1(double sop_j, double cluster_p[], DATA vt[], LOG& logdt){
	int nmtotal=0,temp;
	int totcluster=1;
	int kcluster,clustertrue;
	for(int i=0;i<AT;i++){
		logdt.ncluster[i]=0;
	}
	logdt.ncluster[0]=1;
	int* clustnb=new int[logdt.Ncen];
	int* cluster2=new int[logdt.Ncen];
	for(int s=0;s<logdt.Ncen;++s){
		cluster2[s]=0;
	}
	cluster2[0]=1;
	double rx,ry,rz,vec_AB,vec_A,vec_B,CosD,p2;
	
	double Rdf=cluster_p[1];
	double cluster_judge=cluster_p[2];
	
	for(int k=1;k<logdt.Ncen;++k){
		nmtotal++;
		vec_A=vt[nmtotal].x*vt[nmtotal].x + vt[nmtotal].y*vt[nmtotal].y + vt[nmtotal].z*vt[nmtotal].z;
		vec_A=sqrt(vec_A);
		int nbtotal=-1;
		int clustmin=10000000;
		for(int s=0;s<logdt.Ncen;++s){
			clustnb[s]=-1;
		}
		for(int i=0;i<nmtotal;++i){
			rx=vt[nmtotal].cenx-vt[i].cenx;
			ry=vt[nmtotal].ceny-vt[i].ceny;
			rz=vt[nmtotal].cenz-vt[i].cenz;
			if(rx>logdt.box_x/2){rx=rx-logdt.box_x;}
			if(rx<-logdt.box_x/2){rx=rx+logdt.box_x;}
			if(ry>logdt.box_y/2){ry=ry-logdt.box_y;}
			if(ry<-logdt.box_y/2){ry=ry+logdt.box_y;}
			if(rz>logdt.box_z/2){rz=rz-logdt.box_z;}
			if(rz<-logdt.box_z/2){rz=rz+logdt.box_z;}
			rx=rx*rx+ry*ry+rz*rz;
			rx=sqrt(rx);
			if(rx<=Rdf && vt[nmtotal].sop>=sop_j && vt[i].sop>=sop_j){
				vec_AB=vt[nmtotal].x*vt[i].x + vt[nmtotal].y*vt[i].y + vt[nmtotal].z*vt[i].z;
				vec_B=vt[i].x*vt[i].x + vt[i].y*vt[i].y + vt[i].z*vt[i].z;
				vec_B=sqrt(vec_B);
				CosD=vec_AB/(vec_A*vec_B);
				p2=(3*CosD*CosD-1)/2;
				if(p2>=cluster_judge){
					kcluster=cluster2[i];
					clustertrue=findproclust(kcluster,logdt);
					for(int s=0;s<=nbtotal;++s){
						if(clustertrue==clustnb[s]){
							goto part;
						}
					}
					nbtotal++;
					clustnb[nbtotal]=clustertrue;
					clustmin=min(clustertrue, clustmin);
				}
			}
			part:continue;
		} 
		
		if(nbtotal==-1){
			totcluster++;
			cluster2[nmtotal]=totcluster;
			logdt.ncluster[totcluster]=1;
			goto part2;
		}else{
			cluster2[nmtotal]=clustmin;
			temp=0;
			for(int m=0;m<=nbtotal;++m){
				clustertrue=clustnb[m];
				temp+=logdt.ncluster[clustertrue];
			}
			logdt.ncluster[clustmin]=temp;
			logdt.ncluster[clustmin]+=1;
			for(int m=0;m<=nbtotal;++m){
				if(clustnb[m]!=clustmin){
					int kk=clustnb[m];
					logdt.ncluster[kk]=-clustmin;
				}
			}
		}
		part2:continue;
	}
	
	int nl[nmtotal];
	int site[nmtotal];

	double id_centroidx0[nmtotal];
	double id_centroidy0[nmtotal];
	double id_centroidz0[nmtotal];
	double id_centroidx1[nmtotal];
	double id_centroidy1[nmtotal];
	double id_centroidz1[nmtotal];
	double id_vx[nmtotal];
	double id_vy[nmtotal];
	double id_vz[nmtotal];	
	double id_gyration[nmtotal];
	double id_gyration_gyx2[nmtotal];
	double id_gyration_gyy2[nmtotal];
	double id_gyration_gyz2[nmtotal];
	for(int i=0;i<nmtotal;++i){
		nl[i]=0;
		site[i]=0;

		id_centroidx0[i]=0.0;
		id_centroidy0[i]=0.0;
		id_centroidz0[i]=0.0;
		id_centroidx1[i]=0.0;
		id_centroidy1[i]=0.0;
		id_centroidz1[i]=0.0;
		id_gyration[i]=0.0;
		id_gyration_gyx2[i]=0.0;
		id_gyration_gyy2[i]=0.0;
		id_gyration_gyz2[i]=0.0;
	}
	for(int i=0;i<nmtotal;++i){
		double xx=vt[i].cenx;
		double yy=vt[i].ceny;
		double zz=vt[i].cenz;
		
		kcluster=cluster2[i];
		int idc=findproclust(kcluster,logdt);
		site[i]=idc;
		nl[idc]++;
		if(id_centroidx0[idc]==0.0) id_centroidx0[idc]=xx;
		if(id_centroidy0[idc]==0.0) id_centroidy0[idc]=yy;
		if(id_centroidz0[idc]==0.0) id_centroidz0[idc]=zz;
		
		if((xx-id_centroidx0[idc])>0.5*logdt.box_x){
			xx=xx-logdt.box_x;
		}
		if((xx-id_centroidx0[idc])<-0.5*logdt.box_x){
			xx=xx+logdt.box_x;
		}
		if((yy-id_centroidy0[idc])>0.5*logdt.box_y){
			yy=yy-logdt.box_y;
		}
		if((yy-id_centroidy0[idc])<-0.5*logdt.box_y){
			yy=yy+logdt.box_y;
		}
		if((zz-id_centroidz0[idc])>0.5*logdt.box_z){
			zz=zz-logdt.box_z;
		}
		if((zz-id_centroidz0[idc])<-0.5*logdt.box_z){
			zz=zz+logdt.box_z;
		}
		id_centroidx1[idc]+=xx;
		id_centroidy1[idc]+=yy;
		id_centroidz1[idc]+=zz;
		
		double rr=sqrt(vt[i].x*vt[i].x+vt[i].y*vt[i].y+vt[i].z*vt[i].z);
		id_vx[idc]=vt[i].x/rr;
		id_vy[idc]=vt[i].y/rr;
		id_vz[idc]=vt[i].z/rr;
	}
	for(int i=0;i<nmtotal;++i){
		if(nl[i]>0){
			
			id_centroidx1[i]=id_centroidx1[i]/nl[i];
			id_centroidy1[i]=id_centroidy1[i]/nl[i];
			id_centroidz1[i]=id_centroidz1[i]/nl[i];

			if(id_centroidx1[i]<logdt.box_xlo) id_centroidx1[i]+=logdt.box_x;
			if(id_centroidx1[i]>logdt.box_xhi) id_centroidx1[i]-=logdt.box_x;
			
			if(id_centroidy1[i]<logdt.box_ylo) id_centroidy1[i]+=logdt.box_y;
			if(id_centroidy1[i]>logdt.box_yhi) id_centroidy1[i]-=logdt.box_y;
			
			if(id_centroidz1[i]<logdt.box_zlo) id_centroidz1[i]+=logdt.box_z;
			if(id_centroidz1[i]>logdt.box_zhi) id_centroidz1[i]-=logdt.box_z;
			
		}
	}
	double gy_x, gy_y, gy_z;
	for(int i=0;i<nmtotal;++i){	
		
		gy_x=vt[i].cenx-id_centroidx1[site[i]];
		gy_y=vt[i].ceny-id_centroidy1[site[i]];
		gy_z=vt[i].cenz-id_centroidz1[site[i]];
		if(gy_x>logdt.box_x/2){gy_x=gy_x-logdt.box_x;}
		if(gy_x<-logdt.box_x/2){gy_x=gy_x+logdt.box_x;}
		if(gy_y>logdt.box_y/2){gy_y=gy_y-logdt.box_y;}
		if(gy_y<-logdt.box_y/2){gy_y=gy_y+logdt.box_y;}
		if(gy_z>logdt.box_z/2){gy_z=gy_z-logdt.box_z;}
		if(gy_z<-logdt.box_z/2){gy_z=gy_z+logdt.box_z;}
		id_gyration_gyx2[site[i]]=id_gyration_gyx2[site[i]]+(gy_x*gy_x)/nl[site[i]];
		id_gyration_gyy2[site[i]]=id_gyration_gyy2[site[i]]+(gy_y*gy_y)/nl[site[i]];
		id_gyration_gyz2[site[i]]=id_gyration_gyz2[site[i]]+(gy_z*gy_z)/nl[site[i]];
		id_gyration[site[i]]=id_gyration[site[i]]+(gy_x*gy_x+gy_y*gy_y+gy_z*gy_z)/nl[site[i]];
		if(nl[site[i]]>1){
			vt[i].site=site[i]+1;
		}else{
			vt[i].site=0;
		}
		
		vt[i].site_size=nl[site[i]];
	}
	int cl_last=0;
	int cl_last2=0;
	int cl_all=0;
	for(int i=0;i<nmtotal;++i){
		if(nl[i]>1){
			logdt.sop_cluster_gy[cl_last]=sqrt(id_gyration[i]);
			logdt.sop_cluster_gyx[cl_last]=sqrt(id_gyration_gyx2[i]);
			logdt.sop_cluster_gyy[cl_last]=sqrt(id_gyration_gyy2[i]);
			logdt.sop_cluster_gyz[cl_last]=sqrt(id_gyration_gyz2[i]);
			logdt.sop_cluster_size[cl_last]=nl[i];
			logdt.sop_cluster_cx[cl_last]=id_centroidx1[i];
			logdt.sop_cluster_cy[cl_last]=id_centroidy1[i];
			logdt.sop_cluster_cz[cl_last]=id_centroidz1[i];
			logdt.sop_cluster_vx[cl_last]=id_vx[i];
			logdt.sop_cluster_vy[cl_last]=id_vy[i];
			logdt.sop_cluster_vz[cl_last]=id_vz[i];
			cl_last++;
			cl_all+=nl[i];
		}
		if(nl[i]>10){
			cl_last2++;
		}
	}
	logdt.sop_cluster_num=cl_last;
	logdt.sop_cluster_num2=cl_last2;
	logdt.Ncl_all_sop=cl_all;
}

void cluster_2(double dtt_j, double dtt_sj, double cluster_p[], DATA vt[], DATA2 bt[][MOL2], LOG& logdt){
	int nmtotal=0,temp;
	int totcluster=0;
	int kcluster,clustertrue;
	for(int i=0;i<AT;i++){
		logdt.ncluster[i]=0;
	}
	logdt.ncluster[0]=1;
	int* clustnb=new int[logdt.Ncen+1];
	int* cluster2=new int[logdt.Ncen+1];
	for(int s=0;s<logdt.Ncen;++s){
		cluster2[s]=0;
	}
	cluster2[0]=1;
	double rx,ry,rz,vec_AB,vec_A,vec_B,CosD,p2;
	
	double Rdf=cluster_p[1];
	double cluster_judge=cluster_p[2];
	
	int k1,k2,s1,s2,id1,id2;
	
	for(int k=1;k<logdt.Ncen;++k){
		nmtotal++;
		k1=vt[nmtotal].btid_k;
		s1=vt[nmtotal].btid_s;
		id1=bt[k1][s1].dttid;
		vec_A=vt[nmtotal].x*vt[nmtotal].x + vt[nmtotal].y*vt[nmtotal].y + vt[nmtotal].z*vt[nmtotal].z;
		vec_A=sqrt(vec_A);
		int nbtotal=-1;
		int clustmin=10000000;
		for(int s=0;s<logdt.Ncen;++s){
			clustnb[s]=-1;
		}
		for(int i=0;i<nmtotal;++i){
			k2=vt[i].btid_k;
			s2=vt[i].btid_s;
			id2=bt[k2][s2].dttid;
			rx=vt[nmtotal].cenx-vt[i].cenx;
			ry=vt[nmtotal].ceny-vt[i].ceny;
			rz=vt[nmtotal].cenz-vt[i].cenz;
			if(rx>logdt.box_x/2){rx=rx-logdt.box_x;}
			if(rx<-logdt.box_x/2){rx=rx+logdt.box_x;}
			if(ry>logdt.box_y/2){ry=ry-logdt.box_y;}
			if(ry<-logdt.box_y/2){ry=ry+logdt.box_y;}
			if(rz>logdt.box_z/2){rz=rz-logdt.box_z;}
			if(rz<-logdt.box_z/2){rz=rz+logdt.box_z;}
			rx=rx*rx+ry*ry+rz*rz;
			rx=sqrt(rx);
			if(rx<=Rdf && double(logdt.dtt_num[id1])>dtt_sj && double(logdt.dtt_num[id2])>dtt_sj){
				vec_AB=vt[nmtotal].x*vt[i].x + vt[nmtotal].y*vt[i].y + vt[nmtotal].z*vt[i].z;
				vec_B=vt[i].x*vt[i].x + vt[i].y*vt[i].y + vt[i].z*vt[i].z;
				vec_B=sqrt(vec_B);
				CosD=vec_AB/(vec_A*vec_B);
				p2=(3*CosD*CosD-1)/2;
				if(p2>=cluster_judge){
					kcluster=cluster2[i];
					clustertrue=findproclust(kcluster,logdt);
					for(int s=0;s<=nbtotal;++s){
						if(clustertrue==clustnb[s]){
							goto part;
						}
					}
					nbtotal++;
					clustnb[nbtotal]=clustertrue;
					clustmin=min(clustertrue, clustmin);
				}
			}
			part:continue;
		}
	
		if(nbtotal==-1){
			totcluster++;
			cluster2[nmtotal]=totcluster;
			logdt.ncluster[totcluster]=1;
			goto part2;
		}else{
			cluster2[nmtotal]=clustmin;
			temp=0;
			for(int m=0;m<=nbtotal;++m){
				clustertrue=clustnb[m];
				temp+=logdt.ncluster[clustertrue];
			}
			logdt.ncluster[clustmin]=temp;
			logdt.ncluster[clustmin]+=1;
			for(int m=0;m<=nbtotal;++m){
				if(clustnb[m]!=clustmin){
					int kk=clustnb[m];
					logdt.ncluster[kk]=-clustmin;
				}
			}
		}
		part2:continue;
	}
	
	int nl[nmtotal];
	int site[nmtotal];

	double id_centroidx0[nmtotal];
	double id_centroidy0[nmtotal];
	double id_centroidz0[nmtotal];
	double id_centroidx1[nmtotal];
	double id_centroidy1[nmtotal];
	double id_centroidz1[nmtotal];
	double id_vx[nmtotal];
	double id_vy[nmtotal];
	double id_vz[nmtotal];	
	double id_gyration[nmtotal];
	double id_gyration_gyx2[nmtotal];
	double id_gyration_gyy2[nmtotal];
	double id_gyration_gyz2[nmtotal];
	for(int i=0;i<nmtotal;++i){
		nl[i]=0;
		site[i]=0;

		id_centroidx0[i]=0.0;
		id_centroidy0[i]=0.0;
		id_centroidz0[i]=0.0;
		id_centroidx1[i]=0.0;
		id_centroidy1[i]=0.0;
		id_centroidz1[i]=0.0;
		id_gyration[i]=0.0;
		id_gyration_gyx2[i]=0.0;
		id_gyration_gyy2[i]=0.0;
		id_gyration_gyz2[i]=0.0;
	}
	for(int i=0;i<nmtotal;++i){
		double xx=vt[i].cenx;
		double yy=vt[i].ceny;
		double zz=vt[i].cenz;
		
		kcluster=cluster2[i];
		int idc=findproclust(kcluster,logdt);
		site[i]=idc;
		nl[idc]++;
		if(id_centroidx0[idc]==0.0) id_centroidx0[idc]=xx;
		if(id_centroidy0[idc]==0.0) id_centroidy0[idc]=yy;
		if(id_centroidz0[idc]==0.0) id_centroidz0[idc]=zz;
		
		if((xx-id_centroidx0[idc])>0.5*logdt.box_x){
			xx=xx-logdt.box_x;
		}
		if((xx-id_centroidx0[idc])<-0.5*logdt.box_x){
			xx=xx+logdt.box_x;
		}
		if((yy-id_centroidy0[idc])>0.5*logdt.box_y){
			yy=yy-logdt.box_y;
		}
		if((yy-id_centroidy0[idc])<-0.5*logdt.box_y){
			yy=yy+logdt.box_y;
		}
		if((zz-id_centroidz0[idc])>0.5*logdt.box_z){
			zz=zz-logdt.box_z;
		}
		if((zz-id_centroidz0[idc])<-0.5*logdt.box_z){
			zz=zz+logdt.box_z;
		}
		id_centroidx1[idc]+=xx;
		id_centroidy1[idc]+=yy;
		id_centroidz1[idc]+=zz;
		
		double rr=sqrt(vt[i].x*vt[i].x+vt[i].y*vt[i].y+vt[i].z*vt[i].z);
		id_vx[idc]=vt[i].x/rr;
		id_vy[idc]=vt[i].y/rr;
		id_vz[idc]=vt[i].z/rr;
	}
	for(int i=0;i<nmtotal;++i){
		if(nl[i]>0){
			
			id_centroidx1[i]=id_centroidx1[i]/nl[i];
			id_centroidy1[i]=id_centroidy1[i]/nl[i];
			id_centroidz1[i]=id_centroidz1[i]/nl[i];

			if(id_centroidx1[i]<logdt.box_xlo) id_centroidx1[i]+=logdt.box_x;
			if(id_centroidx1[i]>logdt.box_xhi) id_centroidx1[i]-=logdt.box_x;
			
			if(id_centroidy1[i]<logdt.box_ylo) id_centroidy1[i]+=logdt.box_y;
			if(id_centroidy1[i]>logdt.box_yhi) id_centroidy1[i]-=logdt.box_y;
			
			if(id_centroidz1[i]<logdt.box_zlo) id_centroidz1[i]+=logdt.box_z;
			if(id_centroidz1[i]>logdt.box_zhi) id_centroidz1[i]-=logdt.box_z;
			
		}
	}
	double gy_x, gy_y, gy_z;
	for(int i=0;i<nmtotal;++i){	
		
		gy_x=vt[i].cenx-id_centroidx1[site[i]];
		gy_y=vt[i].ceny-id_centroidy1[site[i]];
		gy_z=vt[i].cenz-id_centroidz1[site[i]];
		if(gy_x>logdt.box_x/2){gy_x=gy_x-logdt.box_x;}
		if(gy_x<-logdt.box_x/2){gy_x=gy_x+logdt.box_x;}
		if(gy_y>logdt.box_y/2){gy_y=gy_y-logdt.box_y;}
		if(gy_y<-logdt.box_y/2){gy_y=gy_y+logdt.box_y;}
		if(gy_z>logdt.box_z/2){gy_z=gy_z-logdt.box_z;}
		if(gy_z<-logdt.box_z/2){gy_z=gy_z+logdt.box_z;}
		id_gyration_gyx2[site[i]]=id_gyration_gyx2[site[i]]+(gy_x*gy_x)/nl[site[i]];
		id_gyration_gyy2[site[i]]=id_gyration_gyy2[site[i]]+(gy_y*gy_y)/nl[site[i]];
		id_gyration_gyz2[site[i]]=id_gyration_gyz2[site[i]]+(gy_z*gy_z)/nl[site[i]];
		id_gyration[site[i]]=id_gyration[site[i]]+(gy_x*gy_x+gy_y*gy_y+gy_z*gy_z)/nl[site[i]];
		if(nl[site[i]]>1){
			bt[vt[i].btid_k][vt[i].btid_s].site=site[i]+1;
		}else{
			bt[vt[i].btid_k][vt[i].btid_s].site=0;
		}
		
		bt[vt[i].btid_k][vt[i].btid_s].site_size=nl[site[i]];
		
	}
	int cl_last=0;
	int cl_last2=0;
	int cl_all=0;
	for(int i=0;i<nmtotal;++i){
		if(nl[i]>1){
			logdt.cluster_gy[cl_last]=sqrt(id_gyration[i]);
			logdt.cluster_gyx[cl_last]=sqrt(id_gyration_gyx2[i]);
			logdt.cluster_gyy[cl_last]=sqrt(id_gyration_gyy2[i]);
			logdt.cluster_gyz[cl_last]=sqrt(id_gyration_gyz2[i]);
			logdt.cluster_size[cl_last]=nl[i];
			logdt.cluster_cx[cl_last]=id_centroidx1[i];
			logdt.cluster_cy[cl_last]=id_centroidy1[i];
			logdt.cluster_cz[cl_last]=id_centroidz1[i];
			logdt.cluster_vx[cl_last]=id_vx[i];
			logdt.cluster_vy[cl_last]=id_vy[i];
			logdt.cluster_vz[cl_last]=id_vz[i];
			cl_last++;
			cl_all+=nl[i];
		}
		if(nl[i]>10){
			cl_last2++;
		}
	}
	logdt.cluster_num=cl_last;
	logdt.cluster_num2=cl_last2;
	logdt.Ncl_all=cl_all;

}

int findproclust(int sn, LOG& logdt){
	int r,t,k;
	
	r=sn;
	t=-logdt.ncluster[r];
	if(t<0) goto find1;
	
	find2:
	r=t;
	k=t;
	t=-logdt.ncluster[k];
	if(t<0){
		logdt.ncluster[sn]=-r;
	}else{
		goto find2;
	}
	
	find1:
		return r;

}
