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
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cstdio>

using namespace std;
int main(int argc, char *argv[]) {

	struct DATA vector[AT];
	struct DATA2 bondatom[MOL1][MOL2];
	struct LOG logdata;
	
	int com=1;
	string input_file_name;
	string output_folder_name;
	int mode_clsuter_sop=0;
	int mode_clsuter_dtt=0;
	int mode_chain_sop=0;
	int mode_chain_dtt=0;
	int mode_dtt=0;
	int mode_sop=0;
	int mode_volume=0;
	double sop_r=0.0,sop_j=0.0,dtt_j=0.0,dtt_size_j=0.0;
	double cluster_p_sop[3];
	double cluster_p_dtt[3];
	int typejudge[20]={0};
	int volume_dump=0;
	int dt=1;
	
	if (argc<2) {
		cout << endl;
		cout << "Crystallization Analysis Package for Semi-crystalline Polymers & Polymeric Nanocomposites" << endl;
		cout << "Usage:" << endl 
			<< "Program -in Input_file -out Folder_name Parameter" << endl;
		cout << "Parameter:" << endl
			<< "      -sop R crystal_judge -c sop_R cluster_judge -cf; e.g. -sop 1.44 0.8 -c 1.05 0.95 -cf" << endl
			<< "      -dtt crystal_judge crystal_size_judge -c dtt_R cluster_judge -cf; e.g. -dtt 0.95 14 -c 2.0 0.625 -cf" << endl
			<< "      -ig typenumber" << endl
			<< "      -v Lx_probe (d)" << endl
			<< "      -ti time interval" << endl;
		cout << "Note: atom input file must follow this format sequence: id mol type x y z ix iy iz; and the box must be triclinic"<<endl;
		cout << "      -sop : choose crystallization analysis by sop; need to define range and judgement"<<endl;
		cout << "      -dtt : choose crystallization analysis by dtt; need to define judgement and stem length"<<endl;
		cout << "      -c : choose cluster analysis or not; need to define range and judgement"<<endl;
		cout << "      -cf : choose conformation analysis or not; only applicable to linear polymers"<<endl;
		cout << "      -v : choose calculate crystallinity by volume; need to define the probe size and choose to dump the snapshot file or not"<<endl;
		cout << "      -ig : choose ignore the atoms' type; should be the last command"<<endl;
		cout << "      -ti : time interval; default 1"<<endl;
		cout << endl;
		cout << "v6.0 Haoyu Wu Feb.23 2024"<<endl;
		cout << endl;
		return 1;
	}else{
		while(com < argc){
			if( strcmp(argv[com],"-in") == 0 ){
				com++;
				input_file_name = argv[com];
			}else if( strcmp(argv[com],"-out") == 0 ){
				com++;
				output_folder_name = argv[com];
			}else if( strcmp(argv[com],"-sop") == 0 ){
				mode_sop=1;
				com++;
				sop_r = atof(argv[com]);
				com++;
				sop_j = atof(argv[com]);
				com++;
				if( strcmp(argv[com],"-c") == 0 ){
					mode_clsuter_sop = 1;
					com++;
					cluster_p_sop[1] = atof(argv[com]);
					com++;
					cluster_p_sop[2] = atof(argv[com]);
				}else{
					com--;
				}
        com++;
        if( strcmp(argv[com],"-cf") == 0 ){
					mode_chain_sop = 1;
				}else{
					com--;
				}
			}else if( strcmp(argv[com],"-dtt") == 0 ){
				mode_dtt=1;
				com++;
				dtt_j = atof(argv[com]);
				com++;
				dtt_size_j = atof(argv[com]);
				com++;
				if( strcmp(argv[com],"-c") == 0 ){
					mode_clsuter_dtt = 1;
					com++;
					cluster_p_dtt[1] = atof(argv[com]);
					com++;
					cluster_p_dtt[2] = atof(argv[com]);
				}else{
					com--;
				}
        com++;
        if( strcmp(argv[com],"-cf") == 0 ){
					mode_chain_dtt = 1;
				}else{
					com--;
				}
			}else if( strcmp(argv[com],"-ig") == 0 ){
				com++;
				while(argv[com]!='\0'){
					typejudge[atoi(argv[com])]=1;
					com++;
				}
			}
			else if( strcmp(argv[com],"-ti") == 0 ){
				com++;
				dt=atoi(argv[com]);
			}else if( strcmp(argv[com],"-v") == 0 ){
				mode_volume = 1;
				com++;
				logdata.probe = atof(argv[com]);
				com++;
				if( strcmp(argv[com],"d") == 0 ){
					volume_dump=1;
				}else{
					com--;
				}
			}else {
				cout << "ERROR: Unrecognized option -- "<<argv[com]<<endl;
				exit(0);
			}
			com++;
		}
	}
	if (input_file_name.empty()) {
		cout << "ERROR: Wrong input file name"<<endl;
		return 1;
	}
	
	string dir="./";
	string output_folder=dir+output_folder_name;
	if(access(output_folder.c_str(), F_OK) == -1){
		int status = mkdir(output_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (status != 0) {
            cout << "ERROR: Wrong creating directory" << endl;
            return 1;
        }
	}
	string data1=output_folder+"/log.dat";
	ofstream out1;
	out1.open((data1).c_str());
	
	string data2=output_folder+"/crystal.xyz";
	ofstream out2;
	out2.open((data2).c_str());
	
	string data3=output_folder+"/cluster(sop)_log.dat";
	ofstream out3;
	out3.open((data3).c_str());
	
	string data4=output_folder+"/cluster(sop)_t.dat";
	ofstream out4;
	out4.open((data4).c_str());
	
	string data5=output_folder+"/cluster(dtt)_log.dat";
	ofstream out5;
	out5.open((data5).c_str());
	
	string data6=output_folder+"/cluster(dtt)_t.dat";
	ofstream out6;
	out6.open((data6).c_str());
	
	string data7=output_folder+"/volume(sop).xyz";
	ofstream out7;
	out7.open((data7).c_str());
	
	string data8=output_folder+"/volume(dtt).xyz";
	ofstream out8;
	out8.open((data8).c_str());

	if(mode_sop!=1&&mode_dtt!=1){
		cout << "ERROR: Wrong mode set" << endl;
		return 1;
	}
	
	ifstream in((input_file_name).c_str());
	if(!in.is_open()){cout << "Error: Input file is empty"<<endl;exit(0);}
	int Ntemp=0;
	logdata.Ncount=0;
	while(!in.eof()){	
		input(in,vector,bondatom,logdata,typejudge);
		bool tag=true;
		if(in.eof()){
			tag=false;
			break;
		}
		if(tag==false) break;
	
		Ntemp++;
		if(Ntemp>0&&Ntemp%dt==0){
			logdata.Ncount++;
			if(mode_sop==1){

				calculate_1(sop_r,sop_j,mode_clsuter_sop,cluster_p_sop,out4,vector,logdata,Ntemp,mode_volume,volume_dump,out7,mode_chain_sop);

			}if(mode_dtt==1){

				calculate_2(dtt_j,dtt_size_j,mode_clsuter_dtt,cluster_p_dtt,out6,vector,bondatom,logdata,Ntemp,mode_volume,volume_dump,out8,mode_chain_dtt);

			}

			output(Ntemp,mode_sop,mode_dtt,mode_clsuter_sop,mode_clsuter_dtt,out1,out2,out3,out5,vector,bondatom,logdata,mode_volume,mode_chain_sop,mode_chain_dtt);
			
		}
	}
	in.close();
	out1.close();
	out2.close();
	out3.close();
	out4.close();
	out5.close();
	out6.close();
	out7.close();
	out8.close();
	
	int rm01,rm02;
	if(volume_dump==0){
		rm01=remove((data7).c_str());
		rm02=remove((data8).c_str());
	}
	
	int rm1,rm2,rm3,rm4;
	if(mode_clsuter_sop==0){
		rm1=remove((data3).c_str());
		rm2=remove((data4).c_str());
	}if(mode_clsuter_dtt==0){
		rm3=remove((data5).c_str());
		rm4=remove((data6).c_str());
	}
	
	if(rm1==-1 || rm2==-1 || rm3==-1 || rm4==-1 || rm01==-1 || rm02==-1){
		cout << "ERROR: Can not remove the cluster dat" << endl;
		return 1;
	}
	
	cout<<"Crystal analysis is finished"<<endl;	
}