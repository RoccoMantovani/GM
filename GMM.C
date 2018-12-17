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


int GMM (string filename, char synth)
{
//***count things
int *count = new int[2];
countdata(count, filename);
int nvertices=count[0];
int nsample=count[1]; 
clog<<"\nnvertices is "<<nvertices<<endl; 
const int nconf=(int)pow(2, nvertices);

//***create phi matrix
int**phi = new int*[nconf]; 
for (int g=0; g<nconf; g++){
	phi[g] = new int [nconf];
}
fillphimatrix(nconf, phi);
writematrix(phi,nconf,nconf,"phi.txt");

//***create & sort s matrix
bool **s = new bool*[nsample]; 
for (int g=0; g<nsample; g++){
	s[g] = new bool [nvertices];
}
readboolmat(s, nvertices, nsample, filename);
quickspinsort(s, nvertices, nsample); 
writematrix(s,nvertices,nsample, "sortedstates.txt");
//***create frequency matrix
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
//***Kpart barriers
int *b=new int[nconf];
for(int q=0;q<nconf;q++){
	b[q]=0;
}
int *bb=new int[nconf];
for(int q=0;q<nconf;q++){
	bb[q]=0;
}
for (int a=1;a<ndiff;a++){
	if(fqm[a][nvertices]!=fqm[a-1][nvertices]){
		b[a]=1;
	}
}
cout<<'\t';
for(int q=1;q<nconf;q++){
	cout<<b[q];
}
cout<<endl;
//compute posteriors for all partitions
int nQ=pow(2,(nconf-1));
double * PostRatio=new double[nQ];
for (int p=0;p<nQ;p++){
	PostRatio[p]=0;
}
double * LogPostRatio=new double[nQ];
for (int p=0;p<nQ;p++){
	LogPostRatio[p]=0;
}
double * Ksplit=new double[nQ];
for (int p=0;p<nQ;p++){
	Ksplit[p]=1;
	for (int l=0;l<nconf-1;l++){
		bb[nconf-1-l]=((p%int(pow(2,l+1)))>=int(pow(2,l)));
		if(b[l]==0&&bb[l]==1){
			Ksplit[p]=2;
		}	
	}
}
ofstream post;
double a;

a=0.001;
allpartitions(PostRatio, LogPostRatio, nvertices, nconf, nsample, ndiff, fqm, a);
post.open("postratios_0001.txt");
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
}
sortcouplings(SortedLogPostRatio,Permutation,nQ,0);
post.open("bestpartitions.txt");
for(int p=nQ-1;p>nQ-1-30;p--){
	post<<*SortedLogPostRatio[p]<<'\t';
	post <<toBinary(Permutation[p],nconf-1)<<endl;
}
post.close();




/*
//***identify frequency partition:
int*TempIndex=new int[ndiff]; 				
int Ksize=countKsize(ndiff, nvertices, fqm, TempIndex);
int oldsize=Ksize;
int*ind=new int[oldsize+1];
int*oldind=new int[oldsize+1];
fpartition(ndiff, nvertices, fqm, TempIndex, ind, Ksize);
for (int p=0;p<oldsize+1;p++){
	oldind[p]=ind[p];
}
*/
return 0;
}

//////////////////////
int main()
{
string filename;
char synth;
cout <<"\nDo you need to generate a synthetic dataset now?\n"
     <<"(you need to have a proper \"decidedinteractions.txt\" file ready)\n"
     <<"[y/n]\n";
cin>>synth;
if(synth=='y'){
	cout<<"ok.\n\n";
	sampler();
	cout<<"dataset generated correctly;\nI will now begin the analysis.\n\n";
	filename="synthetics/1ex0sgc.txt";
}
else{
	cout<< "Which file should I open?\n"; 
	cin >> filename; 
	cout << "I'll open: "<< filename << "!" << endl;
}
GMM(filename,synth);
return 0;
}





