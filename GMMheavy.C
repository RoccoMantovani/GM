#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <bitset>
#include <Eigen/SVD>
#include <algorithm>
#include <string>
#include </home/rocco/Desktop/finalGM/GMM/GMengine.h>
#define REGULARIZER 10e-200
#define OPTNSTEPS 400
#define NREL 1
#define NSEQ 1000
using namespace std;


int GMM (string filename, char synth, int j, int nvertices, 
	int nsample, int nconf, int nQ, double**&SumLogPostRatio, 
	double a, int iter, double*&teams)
{
int pgen=2184;							// !
ofstream outs;
//***create phi matrix***
int**phi = new int*[nconf]; 
for (int g=0; g<nconf; g++){
	phi[g] = new int [nconf];
}
fillphimatrix(nconf, phi);
writematrix(phi,nconf,nconf,"phi.txt");

//***create & sort s matrix***
bool **s = new bool*[nsample]; 
for (int g=0; g<nsample; g++){
	s[g] = new bool [nvertices];
}
readboolmat(s, nvertices, nsample, filename);
quickspinsort(s, nvertices, nsample); 
writematrix(s,nvertices,nsample, "sortedstates.txt");

//***create frequency matrix***
int ndiff=0;
int*freq=new int[nsample];
countocc(s, freq, nsample, nvertices);
for(int i=0;i<nsample;i++){
	if(freq[i]!=0) {
		ndiff++;
	}
}
int**fqm=new int*[ndiff]; 
for(int g=0;g<ndiff;g++){
	fqm[g]=new int[nvertices+1];
}
freqmat(s,fqm, freq, nsample, nvertices, ndiff); 
quickfreqsort(fqm, nvertices, ndiff);
writematrix(fqm,nvertices+1,ndiff,"fqm.txt");

//***KPart barriers***
/* I take sorted fqm and assign a binary value "barrier"
to each pair of adiacent states in it:
states:		    s0    s1    s2    s3    s4    s5       ... 	      s(nC-1)
barriers:	(b0)   b1    b2    b3    b4    b5      ...	b(nC-1)
b0 will always be =0;
When bj=1 the corresponding barrier is considered.
Changing bj will identify a ("freq_ordering_preserving")
partition of the states.
*/
outs.open("results.txt");
int KpIndex=0;
int *Kbarrier=new int[nconf];
int nbarriers=0;
for(int q=0;q<nconf;q++){
	Kbarrier[q]=0;
}
for (int a=1;a<ndiff;a++){
	if(fqm[a][nvertices]!=fqm[a-1][nvertices]){
		Kbarrier[a]=1;
		nbarriers++;
	}
}
if(ndiff<nconf){
	Kbarrier[ndiff]=1;
	nbarriers++;
}
outs<<'\t';
for(int q=1;q<nconf;q++){
	KpIndex+=pow(2,(nconf-1-q))*Kbarrier[q];
	outs<<Kbarrier[q];
}
outs<<endl;
//***pPart barriers***
int *Pbarrier=new int[nconf];
for(int q=0;q<nconf;q++){
	Pbarrier[q]=0;
}

int* mjperm=new int[ndiff];
int* Mj=new int[ndiff];
int** mj=new int*[ndiff];
	for(int y=0;y<ndiff;y++){
	mj[y]=new int;	
	*mj[y]=0;
	mjperm[y]=y;
	Mj[y]=0;
}
int mjparts;

int** m0=new int*[nconf-ndiff];
int* m0perm=new int[nconf-ndiff];
int* M0=new int[nconf-ndiff];
for(int y=0;y<nconf-ndiff;y++){
	m0[y]=new int;	
	*m0[y]=0;
	m0perm[y]=y;
	M0[y]=0;
}
int m0parts;

int ncg=(int)(pow(2,nbarriers)); 
int nfine=(int)(pow(2,ndiff))-(int)(pow(2,nbarriers)); 
int nbroke0=(int)(pow(2,nconf-1))-(int)(pow(2,ndiff));

double * Ksplit=new double[nQ];
double * pmulti=new double[nQ];

int y; 
int last; 
int lost;
int legit;
int o;
ofstream pmul;
pmul.open("multi_p.txt");
for (int p=0;p<nQ;p++){
	cout<<"partition n° " << p << endl;	
	for (int l=0;l<nconf-1;l++){
		Pbarrier[nconf-1-l]=((p%int(pow(2,l+1)))>=int(pow(2,l)));
	}	
	y=0;
//	
	last=ndiff;
	Ksplit[p]=1;
	pmulti[p]=0;
	for(int y=0;y<nconf-ndiff;y++){
		M0[y]=1;
		*m0[y]=0;
	}
	for (int l=ndiff;l<nconf;l++){
		if(Kbarrier[l]==0&&Pbarrier[l]==1){
			Ksplit[p]=3;
			*m0[y]=l-last;
			y++;
			last=l;
		}
	}	
	*m0[y]=nconf-last;
	m0parts=y+1;
	last=1;	
	lost=0;
	sortint(m0, m0perm, nconf-ndiff, 0);
	outs << "populations of unobserved parts:\n";	
	for(int w=0;w<m0parts;w++){
		outs<<*m0[w]<<" ";
	}
	outs<<endl;
	while(last<m0parts){
		if(*m0[last]!=0){
			if (*m0[last]==*m0[last-1]){
				M0[lost]++;
			}
			else{
				lost++;
			}
			last++;
		}
		else{
			break;
		}
	}
	outs << "twin parts of unobserved:\n";	
	for(int y=0;y<lost+1;y++){
		outs<<M0[y]<<" ";
	}
	outs<<endl<<endl;
	pmulti[p]=tgamma(nconf-ndiff+1);
	for(int y=0;y<nconf-ndiff;y++){
		pmulti[p]/=tgamma(*m0[y]+1);
		pmulti[p]/=tgamma(M0[y]+1);
	}
//
	for(int y=0;y<ndiff;y++){
		Mj[y]=1;
		*mj[y]=0;
	}
	legit=0; 
	o=0; 
	last=0;
	lost=0;
	for (int l=1;l<ndiff+1;l++){
		if(Kbarrier[l]==1){
			if(o>0){
				*mj[o]=l-max(last,legit);
				sortint(mj,mjperm, ndiff, 0);
				pmulti[p]*=tgamma(l-legit+1);
				outs<<"*"<<tgamma(l-legit+1);
				for(int u=0;u<o+1;u++){		
					pmulti[p]/=tgamma(*mj[u]+1);
					outs<<"/"<<tgamma(*mj[u]+1);
					outs<<"("<<*mj[u]<<")";
					if(u==0){continue;}
					if (*mj[u]==*mj[u-1]){
						Mj[lost]++;
					}
					else{
						lost++;
					}
					*mj[u-1]=0;
				}
				*mj[o]=0;
				outs <<endl;
				for(int y=0;y<lost+1;y++){
					//pmulti[p]/=tgamma(Mj[y]+1);
					outs<<Mj[y]<<"M";
				}
				outs<<endl<<"--------"<<endl;
				for(int y=0;y<o+1;y++){
					Mj[y]=1;
				}
				o=0;
				lost=0;
			}
			legit=l;
			if(Pbarrier[l]==1){
				last=l;
			}
		}
		if(Kbarrier[l]==0&&Pbarrier[l]==1){
			if (Ksplit[p]!=3&&Ksplit[p]!=2){			
				Ksplit[p]=2;
				outs<<"found a fine!!\n";
			}	
			*mj[o]=l-max(last,legit);
			last=l;	
			o++;
		}
	}
	pmul<<p<<'\t'<<pmulti[p]<<'\t'<<Ksplit[p]<<endl;
	outs<<p<<'\t'<<pmulti[p]<<'\t'<<Ksplit[p]<<endl;
}
pmul.close();

//compute posteriors for all partitions
double * PostRatio=new double[nQ];
for (int p=0;p<nQ;p++){
	PostRatio[p]=0;
}
double * LogPostRatio=new double[nQ];
for (int p=0;p<nQ;p++){
	LogPostRatio[p]=0;
}

ofstream post;
string pref="latest/";
string mid="a"+to_string(a)+"_j"+to_string(j);
allpartitions(PostRatio, LogPostRatio, nvertices, nconf, nsample, ndiff, fqm, a);
post.open(pref+"postratios"+mid+".txt");
for (int p=0;p<nQ;p++){
	post<<p<<'\t'<<LogPostRatio[p]<<'\t'<<Ksplit[p]<<endl;
}
post.close();

int*Permutation=new int[nQ];
double**SortedLogPostRatio=new double* [nQ];
for (int p=0;p<nQ;p++){
	Permutation[p]=p;
	SortedLogPostRatio[p]=new double;
	*SortedLogPostRatio[p]=LogPostRatio[p];
	*SumLogPostRatio[p]+=LogPostRatio[p];
}
sortcouplings(SortedLogPostRatio,Permutation,nQ,0);
double**ReducedSortedLogPostRatio=new double* [nQ];
for (int p=0;p<nQ;p++){
	ReducedSortedLogPostRatio[p]=new double;
	*ReducedSortedLogPostRatio[p]=*SortedLogPostRatio[p]-*SortedLogPostRatio[nQ-1];
	//cout<<*ReducedSortedLogPostRatio[p]<<endl;
}
double broke0=0, fine=0, cg=0;
post.open(pref+"bestpartitions_"+mid+".txt");
for(int p=nQ-1;p>0;p--){
	post<<*SortedLogPostRatio[p]<<'\t';
	post<<Permutation[p]<<'\t'<<toBinary(Permutation[p],nconf-1)<<'\t'<<Ksplit[Permutation[p]];
	if(Ksplit[Permutation[p]]==1){
		cg+=exp(*ReducedSortedLogPostRatio[p]); //cout << cg << " ";
	}
	else if(Ksplit[Permutation[p]]==2){
		fine+=pmulti[Permutation[p]]*exp(*ReducedSortedLogPostRatio[p]);
		nfine+=(pmulti[Permutation[p]]-1);
		outs<<pmulti[Permutation[p]]<<'\t';
	}
	else{
		broke0+=pmulti[Permutation[p]]*exp(*ReducedSortedLogPostRatio[p]); 
		nbroke0+=pmulti[Permutation[p]]-1;
	}
	if(Permutation[p]==pgen){
		post<<" ***TRUE ";
	}
	if(Permutation[p]==KpIndex){
		post<<" ***K(sample) ";
	}
	post<<endl;
}
post.close();
//cout <<"fine/(broke0+cg) = "<<(fine/(broke0+cg))<<endl;
//cout <<"(fine+broke0)/cg = "<<((fine+broke0)/cg)<<endl;
teams[0]+=(fine/(fine+broke0+cg));
teams[1]+=(broke0/(fine+broke0+cg));
teams[2]+=(cg/(fine+broke0+cg));
if(j==0){post.open("latest/teams.txt");}
else{post.open("latest/teams.txt",ios::app);}
post<<(fine/(fine+broke0+cg))<<'\t'<<(broke0/(fine+broke0+cg))<<'\t'<<(cg/(fine+broke0+cg))<<endl;
post.close();
cout 	<< "ncg "<<'\t'<<ncg<<"\tnfine "<<'\t'<<nfine<<"\tnbroke0 "<<'\t'<<nbroke0<<endl;
cout	<<'\t'<<(cg/((fine+broke0+cg)))
	<<'\t'<<(fine/((fine+broke0+cg)))
	<<'\t'<<(broke0/((fine+broke0+cg)))<<endl;
cout	<<'\t'<<(cg/((fine+broke0+cg)*ncg))
	<<'\t'<<(fine/((fine+broke0+cg)*nfine))
	<<'\t'<<(broke0/((fine+broke0+cg)*nbroke0))<<endl;
outs 	<< "ncg "<<'\t'<<ncg<<"\tnfine "<<'\t'<<nfine<<"\tnbroke0 "<<'\t'<<nbroke0<<endl;
outs	<<'\t'<<(cg/((fine+broke0+cg)))
	<<'\t'<<(fine/((fine+broke0+cg)))
	<<'\t'<<(broke0/((fine+broke0+cg)))<<endl;
outs	<<'\t'<<(cg/((fine+broke0+cg)*ncg))
	<<'\t'<<(fine/((fine+broke0+cg)*nfine))
	<<'\t'<<(broke0/((fine+broke0+cg)*nbroke0))<<endl;
outs.close();
if(j==iter-1){
	cout<<"LAST RUN"<<endl;
	for (int p=0;p<nQ;p++){
		Permutation[p]=p;
	}
	sortcouplings(SumLogPostRatio,Permutation,nQ,0);
	post.open("latest/ALLBestpartitions_a"+to_string(a)+".txt");
	for(int p=nQ-1;p>-1;p--){
		post<<*SumLogPostRatio[p]<<'\t';
		post<<Permutation[p]<<'\t'<<toBinary(Permutation[p],nconf-1)<<'\t'<<Ksplit[Permutation[p]];
		if(Permutation[p]==pgen){
			post<<" *** ";
		}
		post<<endl;
	}
	post.close();	

}


for (int g=0; g<nconf; g++){delete[]phi[g];}delete[]phi;
for (int g=0; g<nsample; g++){delete[]s[g];}delete[]s;
delete[]freq;
for(int g=0;g<ndiff;g++){delete[]fqm[g];}delete[]fqm;
delete[]Kbarrier;
delete[]Pbarrier;
delete[]Ksplit;
delete[]LogPostRatio;
delete[]PostRatio;
delete[]Permutation;
for (int p=0;p<nQ;p++){delete[]SortedLogPostRatio[p];}delete[]SortedLogPostRatio;
for (int p=0;p<nQ;p++){delete[]ReducedSortedLogPostRatio[p];}delete[]ReducedSortedLogPostRatio;

return 0;

}

//////////////////////
int main()
{
int iterations=1;
int pgen=2184;
double beta=2;
double a = 1;
double *teams = new double[3];
for (int t=0;t<3;t++){
	teams[t]=0;
}
int nsample=1000;
string filename;
char synth='n';
for (int j=0;j<iterations;j++){
	cout <<"\nDo you need to generate a synthetic dataset now?\n"
	     <<"(you need to have a proper \"decidedinteractions.txt\" file ready)\n"
	     <<"[y/n]\n";
	//cin>>synth;
	if(synth=='y'){
		cout<<"ok.\n\n";
		sampler(beta,nsample);
		cout<<"dataset generated correctly;\nI will now begin the analysis.\n\n";
		filename="synthetics/1ex0sgc.txt";
	}
	else{
		cout<< "Which file should I open?\n"; 
		cin >> filename; 
		cout << "I'll open: "<< filename << "!" << endl;
	}
	//***count things
	int *count = new int[2];
	countdata(count, filename);
	int nvertices=count[0];
	int nsample=count[1]; 
	clog<<"\nnvertices is "<<nvertices<<endl; 
	const int nconf=(int)pow(2, nvertices);
	int nQ=pow(2,(nconf-1));
	//***do
	double ** SumLogPostRatio=new double*[nQ];
	for (int p=0;p<nQ;p++){
		SumLogPostRatio[p]=new double;
		*SumLogPostRatio[p]=0;
	}	
	GMM(filename,synth,j,nvertices,nsample,nconf, nQ, SumLogPostRatio,a, iterations, teams);
}
for (int t=0;t<3;t++){teams[t]/=iterations;}

return 0;
}





