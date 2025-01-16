#include <stdio.h>
#include <fstream>

#define AT 1000000
#define MOL1 10000
#define MOL2 1000
using namespace std;

struct DATA
{
	double cenx,ceny,cenz,x,y,z,sop;
	int cenid,mol,ty,site,site_size,cenix,ceniy,ceniz,moltype_sop,moltype_dtt;
	int btid_k,btid_s;
};
extern struct DATA vector[AT];

struct DATA2
{
	double x,y,z,realx,realy,realz,vx,vy,vz;
	int ty,dttid,site,site_size,ix,iy,iz;
};
extern struct DATA2 bondatom[MOL1][MOL2];

struct LOG
{
	double box_x,box_y,box_z,sop_avg,dtt_avg,dttc_avg;
	double box_xlo,box_xhi,box_ylo,box_yhi,box_zlo,box_zhi;
	int Ncen;
	
	int Ncount;
	
	int Ntail_s,Ntie_s,Nloop_s;
	double Ltail_s,Ltie_s,Lloop_s;
	int Ntail_d,Ntie_d,Nloop_d;
	double Ltail_d,Ltie_d,Lloop_d;
	
	double sop_crystal_pro,sop_crystal_pro_v;	
	double sop_cluster_gy_max,sop_cluster_gy_avg,sop_cluster_gy_x_max,sop_cluster_gy_x_avg,sop_cluster_gy_y_max,sop_cluster_gy_y_avg,sop_cluster_gy_z_max,sop_cluster_gy_z_avg,sop_cluster_size_avg;	
	double sop_cluster_size_avgW;
	int sop_crystal_num,sop_cluster_num,sop_cluster_size_max;
	int sop_cluster_num2,Ncl_all_sop;
	int sop_cluster_size[AT];
	double sop_cluster_gy[AT];
	double sop_cluster_gyx[AT];
	double sop_cluster_gyy[AT];
	double sop_cluster_gyz[AT];
	double sop_cluster_cx[AT];
	double sop_cluster_cy[AT];
	double sop_cluster_cz[AT];
	double sop_cluster_vx[AT];
	double sop_cluster_vy[AT];
	double sop_cluster_vz[AT];
	
	int ncluster[AT];
	
	double crystal_pro,crystal_pro_v;	
	double cluster_gy_max,cluster_gy_avg,cluster_gy_x_max,cluster_gy_x_avg,cluster_gy_y_max,cluster_gy_y_avg,cluster_gy_z_max,cluster_gy_z_avg,cluster_size_avg;	
	double cluster_size_avgW;
	int crystal_num,cluster_num,cluster_size_max;
	int cluster_num2,Ncl_all;
	int cluster_size[AT];
	double cluster_gy[AT];
	double cluster_gyx[AT];
	double cluster_gyy[AT];
	double cluster_gyz[AT];
	double cluster_cx[AT];
	double cluster_cy[AT];
	double cluster_cz[AT];
	double cluster_vx[AT];
	double cluster_vy[AT];
	double cluster_vz[AT];
	int dtt_num[AT];
	int Mol,Ndtt;
	int Nmol[AT];
	double probe;
};
extern struct LOG logdata;

void input(ifstream &in, DATA vt[], DATA2 bt[][MOL2], LOG& logdt, int typejudge[]);
void calculate_1(double sop_r, double sop_j, int mode_c, double cluster_p[], ofstream& out4, DATA vt[], LOG& logdt, int Nt, int mode_v, int v_dump, ofstream& out7, int mode_chain_sop);
void calculate_2(double dtt_j, double dtt_sj, int mode_c, double cluster_p[], ofstream& out6, DATA vt[], DATA2 bt[][MOL2], LOG& logdt, int Nt, int mode_v, int v_dump, ofstream& out8, int mode_chain_dtt);
void volume_cry1(DATA vt[], LOG& logdt, double sop_j, int v_dump, ofstream& out7, int Nt);
void volume_cry2(DATA vt[], DATA2 bt[][MOL2], LOG& logdt, double dtt_sj, int v_dump, ofstream& out8, int Nt);
void cluster_1(double sop_j, double cluster_p[], DATA vt[], LOG& logdt);
void cluster_2(double dtt_j, double dtt_sj, double cluster_p[], DATA vt[], DATA2 bt[][MOL2], LOG& logdt);
int findproclust(int sn, LOG& logdt);
void output(int Ntemp, int mode_sop, int mode_dtt, int mode_c_sop, int mode_c_dtt, ofstream& out1, ofstream& out2, ofstream& out3, ofstream& out5, DATA vt[], DATA2 bt[][MOL2], LOG& logdt, int mode_v, int mode_chain_sop, int mode_chain_dtt);
void conformation_sop(DATA vt[], LOG& logdt, double sop_j);
void conformation_dtt(DATA2 bt[][MOL2], DATA vt[], LOG& logdt, double dtt_sj);