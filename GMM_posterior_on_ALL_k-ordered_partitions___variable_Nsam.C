/*19/12/2018
This program computes the posterior P(Q|sample) for ALL partitions obtainable by
considering the list of all possible states ordered by their frequency in 
the sample, and enforcing (nstates-1) binary "barriers" between adjacent states in this list.
I then have 2(^nstates - 1) possible "ordered" partitionings of the state space.

Our considerations should only be affected by the relative positions, in the list,
of states that appear in the sample a DIFFERENT number of times; thus, we correct our posterior 
probabilities by multiplying each of them by a multiplicity factor counting the number of different
possible partitions sharing the specified "barrier" structure, but with states in the same 
frequency classes scrambled.

We compare the cumulative probabilities of
1) P is a coarse graining of K
2) P is such that the "unobserved states" set is not broken, but at least one "k!=0" set of the K partition is;
3) P is such that no "k!=0" set of the K partition is broken, but the unobserved set is.
4) P is broken both in observet & in unobserved states.
computing all of them for different sample sizes, fixed a generative distribution.
*/


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <bitset>
#include </scratch/rmantova/Eigen/SVD>
#include <algorithm>
#include <string>
#include </scratch/rmantova/GM-exact-ordered-partitions-loglik-ratios/GMengine.h>
//#include </home/rocco/Desktop/finalGM/GMM/GMengine.h>
#define REGULARIZER 10e-200
#define OPTNSTEPS 400
#define NREL 1
#define NSEQ 1000
using namespace std;


int GMM (string filename, char synth, int j, int nvertices, 
	int nsample, int nconf, long long int nQ, double**&SumLogPostRatio, 
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
for(int q=1;q<nconf;q++){
	KpIndex+=pow(2,(nconf-1-q))*Kbarrier[q];
	outs<<Kbarrier[q];
}
outs<<endl<<endl;
//***pPart barriers***
int *Pbarrier=new int[nconf];
for(int q=0;q<nconf;q++){
	Pbarrier[q]=0;
}

long long int ncg=      (int)(pow(2,nbarriers)); 
long long int nf=       (int)(pow(2,ndiff))-(int)(pow(2,nbarriers)); 
long long int ncg0b=    (int)(pow(2,nbarriers+nconf-1-ndiff))-(int)(pow(2,nbarriers));
long long int nf0b=     (int)(pow(2,nconf-1))-ncg0b-nf-ncg;

double * Ksplit=new double[nQ];
double * pmulti=new double[nQ];

int y; 
int last; 
int first;
int lost;
int legit;
int o;
int Klast;
int R=0;
int*mj=new int[nconf];
ofstream pmul;
pmul.open("multi_p.txt");
for (int p=0;p<nQ;p++){
    outs<<endl;
    //prepare partition - draw barriers
	cout<<"partition n° " << p << endl;
	for (int l=0;l<nconf-1;l++){
		Pbarrier[nconf-1-l]=((p%int(pow(2,l+1)))>=int(pow(2,l)));
	}	
	y=0;
	Ksplit[p]=1;
	pmulti[p]=1;
    for (int t=0;t<ndiff;t++){
        mj[t]=0;
    }
	legit=0; 
	o=0; 
	last=0;
    Klast=1;
	lost=0;
    first=0;
    R=0;
	for (int l=1;l<nconf+2;l++){
		if(l==nconf||Kbarrier[l]==1){
//            outs<<"'\t'K1   ";
			if (l!=nconf&&Pbarrier[l]==0){                     //K1 P0
//                outs<<"'\t'P0   ";
                if (Klast==1){
                    if (last!=legit){
                        pmulti[p]*=tgamma(l-legit+1);
                        pmulti[p]/=tgamma(l-last+1);
                        pmulti[p]/=tgamma(last-legit+1);
                    }
                }
                else{
                    pmulti[p]*=tgamma(l-legit+1); 
                    //outs<<pmulti[p]<<"   ";
                    pmulti[p]/=tgamma(l-last+1);
                    //outs<<pmulti[p]<<"   ";
                    pmulti[p]/=tgamma(last-legit+1);
                    //outs<<pmulti[p]<<"   ";
                    pmulti[p]*=tgamma(first-legit+1);
                    //outs<<pmulti[p]<<"   ";
                    pmulti[p]*=tgamma(last-first+1);
                    //outs<<pmulti[p]<<"   ";
                    pmulti[p]/=tgamma(last-legit+1);
                    //outs<<pmulti[p]<<"   ";
                }
                Klast=0;
            }
            else{                                   //K1 P1
 //               outs<<"'\t'P1   ";
                if(Klast==0){
                    pmulti[p]*=tgamma(l-legit+1);
                    pmulti[p]/=tgamma(l-first+1);
                    pmulti[p]/=tgamma(first-legit+1);
                }
                Klast=1;
                mj[R]=1;
            }
            /*if (l==nconf){
                last=l;
            }*/
            pmulti[p]*=tgamma(last-first+1);        //K1
//            outs<<last-first<<"    \n";
            first=l; 
            last=l;
            legit=l;
            for(int t=0;t<R+1;t++){
                pmulti[p]/=tgamma(mj[t]+1);
                mj[t]=0;
            }
            pmulti[p]/=tgamma(R+1);
            //outs<<pmulti[p]<<endl;
            R=0;
            if(l==nconf){
                break;
            }
        }
        else{                                       //K0 P1
            //outs<<"'\t'K0   ";
            if(Pbarrier[l]==1){
                if(l<ndiff){
                    Ksplit[p]=2;
                }
                else{
                    Ksplit[p]+=2;
                }
 //               outs<<"'\t'P1   ";
                if(first==legit){
                    first=l;
                }
                else{
                    R++;
                }
                mj[R]=1;
                last=l;
                //outs<<"last   "<<last<<"    first    "<<first<<"     ";
                //outs<<pmulti[p]<<endl;
            }
            else{                                   //K0 P0
//                outs<<"'\t'P0   ";
                if(first!=legit){
                    mj[R]++;
                }
//                outs <<pmulti[p]<< "\tR "<<R<<"\tmj[R] "<<mj[R]<<'\t'<<endl;
            }
        }
    }
	pmul<<p;
    for(int w=1;w<nconf;w++){
        outs<<Pbarrier[w];
    }
    outs<<'\t'<<pmulti[p]<<'\t'<<Ksplit[p]<<endl;
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
double broke0=0, fine=0, cg=0, finebroke0=0;
post.open(pref+"bestpartitions_"+mid+".txt");
for(int p=nQ-1;p>0;p--){
	post<<*SortedLogPostRatio[p]<<'\t';
	post<<Permutation[p]<<'\t'<<toBinary(Permutation[p],nconf-1)<<'\t'<<Ksplit[Permutation[p]];
	if(Ksplit[Permutation[p]]==1){
		cg+=exp(*ReducedSortedLogPostRatio[p]); 
        //cout << cg << " ";
	}
	else if(Ksplit[Permutation[p]]==2){
		fine+=pmulti[Permutation[p]]*exp(*ReducedSortedLogPostRatio[p]);
		nf+=(pmulti[Permutation[p]]-1);
// 		outs<<pmulti[Permutation[p]]<<'\t';
	}
	else if(Ksplit[Permutation[p]]==3){
		broke0+=pmulti[Permutation[p]]*exp(*ReducedSortedLogPostRatio[p]); 
		ncg0b+=pmulti[Permutation[p]]-1;
//         cout<<ncg0b<<'\t';
	}
    else if(Ksplit[Permutation[p]]==4){
//         outs<<"444444   ";
		finebroke0+=pmulti[Permutation[p]]*exp(*ReducedSortedLogPostRatio[p]); 
		nf0b+=pmulti[Permutation[p]]-1;
//         cout<<nf0b<<'\t';
	}	
	if(Permutation[p]==pgen){
		post<<" ***TRUE ";
	}
	if(Permutation[p]==KpIndex){
		post<<" ***K(sample) ";
	}
	post<<endl;
}
cout<<endl;
post.close();
//cout <<"fine/(broke0+cg) = "<<(fine/(broke0+cg))<<endl;
//cout <<"(fine+broke0)/cg = "<<((fine+broke0)/cg)<<endl;
teams[0]=(cg/(fine+broke0+cg+finebroke0));
teams[1]=(fine/(fine+broke0+cg+finebroke0));
teams[2]=(broke0/(fine+broke0+cg+finebroke0));
teams[3]=(finebroke0/(fine+broke0+cg+finebroke0));
if(j==0){post.open("latest/teams.txt");}
else{post.open("latest/teams.txt",ios::app);}
post<<teams[0]<<'\t'<<teams[1]<<'\t'<<teams[2]<<'\t'<<teams[3]<<endl;
post.close();
teams[0]=(cg/(ncg)*(fine+broke0+cg+finebroke0));
teams[1]=(fine/(nf)*(fine+broke0+cg+finebroke0));
teams[2]=(broke0/(ncg0b)*(fine+broke0+cg+finebroke0));
teams[3]=(finebroke0/(nf0b)*(fine+broke0+cg+finebroke0));
if(j==0){post.open("latest/representatives.txt");}
else{post.open("latest/representatives.txt",ios::app);}
post<<teams[0]<<'\t'<<teams[1]<<'\t'<<teams[2]<<'\t'<<teams[3]<<endl;
post.close();

cout 	<< "ncg "<<'\t'<<ncg<<"\tnf "<<'\t'<<nf<<"\tncg0b "<<'\t'<<ncg0b<<"\tnf0b "<<'\t'<<nf0b<<endl;
cout	<<'\t'<<teams[0]<<'\t'<<teams[1]<<'\t'<<teams[2]<<'\t'<<teams[3]<<endl;
cout	<<'\t'<<teams[0]/ncg<<'\t'<<teams[1]/nf<<'\t'<<teams[2]/ncg0b<<'\t'<<teams[3]/nf0b<<endl;
 outs    <<"ncg\t\tnf\t\tncg0b\t\tnf0b\n";
 outs 	<< ncg<<'\t'<<nf<<'\t'<<ncg0b<<'\t'<<nf0b<<endl;
 outs	<<teams[0]<<'\t'<<teams[1]<<'\t'<<teams[2]<<'\t'<<teams[3]<<endl;
 outs	<<teams[0]/ncg<<'\t'<<teams[1]/nf<<'\t'<<teams[2]/ncg0b<<'\t'<<teams[3]/nf0b<<endl;
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
delete[]mj;
delete[]pmulti;
for (int p=0;p<nQ;p++){delete[]SortedLogPostRatio[p];}delete[]SortedLogPostRatio;
for (int p=0;p<nQ;p++){delete[]ReducedSortedLogPostRatio[p];}delete[]ReducedSortedLogPostRatio;

return 0;

}

//////////////////////
int main()
{
int iterations=500;
int pgen=2184;
double beta=1;
double a = 1;
double *teams = new double[3];
for (int t=0;t<3;t++){
	teams[t]=0;
}
/*double *runningteams = new double[3];
for (int t=0;t<3;t++){
	runningteams[t]=0;
}*/
int nsample=1000;
string filename;
char synth='y';
for (int j=0;j<iterations;j++){
	cout <<"\nDo you need to generate a synthetic dataset now?\n"
	     <<"(you need to have a proper \"decidedinteractions.txt\" file ready)\n"
	     <<"[y/n]\n";
	//cin>>synth;
	if(synth=='y'){
		cout<<"ok.\n\n";
		sampler(beta,5*(j+1));
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
	long long int nQ=pow(2,(nconf-1));
	//***do
	double ** SumLogPostRatio=new double*[nQ];
	for (int p=0;p<nQ;p++){
		SumLogPostRatio[p]=new double;
		*SumLogPostRatio[p]=0;
	}	
	GMM(filename,synth,j,nvertices,nsample,nconf, nQ, SumLogPostRatio,a, iterations, teams);
	for (int p=0;p<nQ;p++){delete[]SumLogPostRatio[p];}
	delete[]SumLogPostRatio;
}
//for (int t=0;t<3;t++){runningteams[t]/=iterations;}

return 0;
}





