//GMhead.h
using namespace std;
#define REGULARIZER 10e-200
#define STP 0.3
//----------------------------------------------------------
void gsampler(double*&g, int nv, double beta, int nsample, int &nfile, int&ex)
{
int mgtd=(int)pow(10,ex);
long double * p=new long double[(int)pow(2,nv)];
double expo=0.; 
int prod=1; 
int ix=0, iy=0; 
int x=0, y=0;
for (int s=0;s<pow(2,nv);s++){			//states s
	expo=g[0];
	for (int i=1;i<pow(2,nv);i++){		//couplings
		prod=1;	
		ix=i; 
		iy=s;
		for (int j=0; j<nv;j++){
			x=(ix>=(pow(2,nv-1-j)));
			y=(iy>=(pow(2,nv-1-j)));
			if (x==1){
				prod*=(y*2-1); 
				ix-=(pow(2,nv-1-j));
			}
			if (y==1){
				iy-=(pow(2,nv-1-j));
			}  	
		}		
		expo+=g[i]*prod;
	}
	p[s]=exp(beta*expo);
}
double Z=0;
for (int s=0;s<pow(2,nv);s++){
	Z+=p[s];
} 
for (int s=0;s<pow(2,nv);s++){
	p[s]=p[s]/Z; 
	//cout << s << '\t' << p[s]<<endl;
}

long double * dots=new long double[(int)pow(2,nv)];
dots[0]=0.;
for (int s=1;s<pow(2,nv);s++){
	dots[s] = dots[s-1] + p[s-1];
}

srand (time (0));

ofstream spinconfigs;
spinconfigs.open("synthetics/"+to_string(mgtd)+"ex"+to_string(nfile)+"sgc.txt");

for (int t=0;t<nsample;t++)
{
	long double x= ((double)rand())/RAND_MAX; 
	int i=0;
	while(x>dots[i]){
		i++;
		if(i==(int)pow(2,nv)){
			break;
		}
	} 
	i=i-1;
	if(t!=0){
		spinconfigs<<'\n';
	}	
	for (int bit=0;bit<nv;bit++)
	{
		if(i>=pow(2,nv-1-bit)){
			spinconfigs<<1;
			i-=pow(2,nv-1-bit);
		}
		else {
			spinconfigs<<0;
		}
	}	
}
spinconfigs<<'\n';
spinconfigs.close();
spinconfigs.open("synthetics/"+to_string(mgtd)+"ex"+to_string(nfile)+"g.txt");
for (int y=1;y<(int)pow(2,nv);y++){
spinconfigs<<g[y]<<endl;
}
spinconfigs.close();
nfile++;
delete [] p;
delete [] dots;
return;
}
//----------------------------------------------------------------
void gennext(double*&g, int*&indx, int &nv, int cur, int &npts, 
		double &gmax, double &gmin, int &nsample, int&nfile, int &ex)
{
if (indx[cur]!=0){
	for (int i=0;i<npts+1;i++){
		g[indx[cur]]=gmin+(gmax-gmin)*i/npts;
		if(indx[cur+1]==0){
			gsampler(g,nv,1,nsample,nfile, ex);
		}
		else{
			gennext(g,indx,nv,cur+1,npts,gmax,gmin,nsample,nfile, ex);
		}
	}
}
return;
}
//----------------------------------------------------------------

int explorer(int &nve, int npts, double gmin, double gmax, int nsample, int ex)
{
int ntotrun=0;
ifstream mod;
mod.open("decmod.txt");
int nv=0;
mod>>nv;
nve=nv;
double * g=new double[(int)pow(2,nv)]; 
for (int i=0;i<(int)pow(2,nv);i++){
	g[i]=0;
} 
int a=1; 
int*indx=new int[(int)pow(2,nv)];
for (int i=0;i<(int)pow(2,nv);i++){
	indx[i]=0;
} 
int i=0;
mod>>a; 
while (a!=111){
	while (a!=0){
		if (a>nv || a<0) {
			cout << "the model is not well "
			     << "written. Check it."; 
			return 0;
		}
		indx[i]+=(int)pow(2,nv-a);
		mod >> a; 
	}
	i++;
	mod>>a;
}
int nfile=0;
cout<<"\n GOT HERE\n";
gennext(g,indx,nv,0,npts,gmin,gmax,nsample,nfile,ex);
ntotrun=(int)pow(npts+1,i);
clog<<"\n[explorer] NTOTRUN " << ntotrun << endl;
cout<<"\n[explorer] NTOTRUN " << ntotrun << endl;
return ntotrun;
}
//----------------------------------------------------------------
void sampler(double beta=0, int nsample=0)
{
int nfile=0, ex=0;
ifstream cip; 
cip.open("decidedinteractions.txt");
int nv=0; 
cip>>nv; 
cout << "number of vertices is " 
     << nv 
     << endl;
if(beta==0){
	cout<<"What's the value of beta? ";
	cin>>beta;
}
cout<<"ok.\n";
if(nsample==0){
	cout <<"How many data points do you want to generate?\n";
	cin>>nsample;
}
//int nsample=895;
double * g=new double[(int)pow(2,nv)]; 
for (int i=0;i<(int)pow(2,nv);i++){
	g[i]=0;
}
int a=1, indx=0;
cip>>a; 
while (a!=111){
	indx=0;
	while (a!=0){
		if (a>nv || a<0) {
			cout << "the couplings list is not well "
			     << "written. Check it."; 
			return;
		}
		indx+=(int)pow(2,nv-a);
		cip >> a; 
	}
	cip>>g[indx]; 
	cip>>a;
}

//for (int i=0;i<pow(2,nv);i++){cip>>g[i];};
gsampler(g,nv,beta,nsample,nfile,ex);
/*
long double * p=new long double[(int)pow(2,nv)];
double expo=0.; 
int prod=1; 
int ix=0, iy=0; 
int x=0, y=0;
for (int s=0;s<pow(2,nv);s++){			//states s
	expo=g[0];
	for (int i=1;i<pow(2,nv);i++){		//couplings
		prod=1;	
		ix=i; 
		iy=s;
		for (int j=0; j<nv;j++){
			x=(ix>=(pow(2,nv-1-j)));
			y=(iy>=(pow(2,nv-1-j)));
			if (x==1){
				prod*=(y*2-1); 
				ix-=(pow(2,nv-1-j));
			}
			if (y==1){
				iy-=(pow(2,nv-1-j));
			}  	
		}		
		expo+=g[i]*prod;
	}
	p[s]=exp(beta*expo);
}
double Z=0;
for (int s=0;s<pow(2,nv);s++){
	Z+=p[s];
} 
for (int s=0;s<pow(2,nv);s++){
	p[s]=p[s]/Z; 
	//cout << s << '\t' << p[s]<<endl;
}

long double * dots=new long double[(int)pow(2,nv)];
dots[0]=0.;
for (int s=1;s<pow(2,nv);s++){
	dots[s] = dots[s-1] + p[s-1];
}

srand (time (0));

ofstream spinconfigs;
spinconfigs.open("sgc.txt");

for (int t=0;t<nsample;t++)
{
	long double x= ((double)rand())/RAND_MAX; 
	int i=0;
	while(x>dots[i]){
		i++;
	} 
	i=i-1;
	if(t!=0){
		spinconfigs<<'\n';
	}	
	for (int bit=0;bit<nv;bit++)
	{
		if(i>=pow(2,nv-1-bit)){
			spinconfigs<<1;
			i-=pow(2,nv-1-bit);
		}
		else {
			spinconfigs<<0;
		}
	}	
}
spinconfigs<<'\n';
spinconfigs.close();
*/
delete [] g;
return;
}
//----------------------------------
void countdata(int *&count,string filename)	
/*
>counts rows (#observations) and columns (#spins) of the dataset;
>INPUT:	> file "filename" containing the list of observed
 	  configurations; 
	> dummy "count";
>OUTPUT: stores #columns in count[0] and #rows in count[1].
*/
{
int i=0, j=0, nrows, ncolumns; 
string state;
bool col=false;
ifstream spinfil;
spinfil.open(filename); 
if(!spinfil){
	clog 	<< "[countdata] error. Abort here";
	cout 	<< "document not valid." 
		<< '\t' 
		<< endl;
	exit(-1);
}
spinfil>>state; 
ncolumns=state.length();
while(1){
	spinfil>>state; 
	i++;
	if (state.length()!=ncolumns||!spinfil){
		nrows=i;
		break;
	}
}
cout 	<< "[countdata] We have "
	<< (count[0] = ncolumns) 
	<< " simultaneously observed binary variables for " 
	<< (count[1] = nrows) 
	<< " observations.\n";
spinfil.close();
}


//-------------------------------------------------------
void readboolmat(bool **&s, int nvertices, 
		int nsample, string filename) 	 
/*
>reads the dataset;
>INPUT: > file "filename" containing the list of observed
 	  configurations;
	> #columns of "filename" ("nvertices")
	> #rows of "filename" ("nsample")
	> dummy "s" of dimensions [nvertices]X[nsample];
>OUTPUT: stores read dataset in "s[spin][observation]"
*/
{											
string line; 
ifstream spinfi;
spinfi.open(filename); 
if(!spinfi) {
	clog 	<< "[countdata] error. Abort here";
	cout << "document is not valid" 
	     << '\t' 
	     << endl;
}
for (int a=0;a<nsample;a++){
	spinfi>>line; 
	for (int b=0;b<nvertices;b++){
		s[a][b]=(bool)((int)(line.at(b) - '0'));
	}
}
spinfi.close();
clog<< "[readboolmat] all ok here\n";
return;
}

//---------------------------------------------------------
bool eqspin(bool*&bj, bool*&bk, int nvertices) 	
/*
> "==" operator for bool arrays;
>INPUT:	> integer "nvertices";
	> bool arrays "bj", "bk", both of length "nvertices"; 
>OUTPUT: returns "1" if strings "bj", "bk" are identical;
	 returns "0" otherwise.
*/
{
int m=0; 
bool ret=1;
while(m<nvertices){
	if (bj[m]!=bk[m]){
		ret=0;
		break;
	}
	m++;
}
if (m==nvertices) {
	ret=1;
	}
return ret;
}
//---------------------------------------------------------
bool mineqspin(bool*&bj, bool*&bk, int nvertices) 
/*
> "<=" operator for bool arrays;
>INPUT:	> integer "nvertices";
	> bool arrays "bj", "bk", both of length "nvertices"; 
>OUTPUT: returns "1" if the binary numbers encoded in "bj", "bk"
	 are such that bj<=bk;
	 returns "0" otherwise.
*/
{
int m=0; 
bool ret=0;
while(m<nvertices){
	if (bj[m]==0&&bk[m]==1){
		ret=1;
		break;
	} 
	else if(bj[m]==1&&bk[m]==0){
		ret=0;
		break;
	}
	m++;
}
if (m==nvertices) {
	ret=1;
}
return ret;
}
//---------------------------------------------------------
void quickss(bool**&b, int low, int high, int 
		nvertices,bool*&temp) 
/*
>quicksort recursive engine called by a boolean rows list 
 sorting function (quickspinsort);
>INPUT: > sample matrix "b" : [nvertices]X[nsample];
	> "nvertices";
	> "low", "high": boundaries of the present
	  window of quicksorting algorithm;
>OUTPUT: calls itself as many times as necessary; at the final
 	return the selected "room" of binary rows has been 
	rearranged, in "increasing" order, according to the 
	"mineqspin" comparator.  
>during the process of rearrangement the actual
 binary strings are never moved in memory,
 the function only swaps the order in which 
 they are pointed;
*/
{
if(low>=high){
	return;
}
int pivot = high, i=low-1;
for (int j=low;j<=high;j++){
	if (mineqspin(b[j],b[pivot], nvertices)){
		i++;
		temp=b[i];
		b[i]=b[j];
		b[j]=temp;
	}
}
quickss (b, low, i-1, nvertices,temp); 
quickss (b, i+1, high, nvertices,temp);
//cout<< "[quickss] all ok here\n";
return;
}
//---------------------------------------------------------
void quickspinsort(bool**&s, int nvertices, int nsample) 	
/*
>takes a bool matrix (spins) and uses a given row vs row 
comparator to arrange its rows in increasing order.
>INPUT: > sample matrix "s" : [nvertices]X[nsample];
	> "nvertices", "nsample";
>OUTPUT: at the final return the whole matrix of binary rows has
	 been rearranged, in "increasing" order, according to the 
	"mineqspin" comparator.  
*/
{
int low=0, high=nsample-1; bool*temp = s[0];
quickss(s, low, high, nvertices, temp);
clog<< "[quickspinsort] all ok here\n";
return;
}
//---------------------------------------------------------
void writematrix(int **&s, int nrows, int ncols, string name) 
/*
>writes any matrix of integers on file.
>INPUT: > integer matrix to be written "s" : [nrows]x[ncols];
	> "nrows", "ncolumns";
	> name of desired file "name";
>OUTPUT: a file [name].txt is created and filled.
*/		
{											
int i=0, j=0, copied=0;
ofstream spinfi;
spinfi.open(name);
for (j=0;j<ncols;j++){
	for (i=0;i<nrows;i++) {
		spinfi << s[j][i] << " ";
	} 
	spinfi<<endl;
}
spinfi.close();
clog<< "[writematrix(int)] all ok here\n";
return;
}

//---------------
void writematrix(bool **&s, int ncols, int nrows,
				string name) 	
/*
>writes any matrix of booleans on file.  
>INPUT: > bool matrix to be written "s" : [nrows]x[ncols];
	> "nrows", "ncolumns";
	> name of desired file "name";
>OUTPUT: a file [name].txt is created and filled.
*/	
{											
int i=0, j=0, copied=0;
ofstream spinfi;
spinfi.open(name);
for (j=0;j<nrows;j++){
	for (i=0;i<ncols;i++) {
		spinfi << s[j][i];
	} 
	spinfi<<endl;
}
spinfi.close();
clog<< "[writematrix(bool)] all ok here\n";
return;
}
//---------------------------------------------------------
string toBinary(int n, int nbits)
{
    string r;
    for(int i=0;i<nbits;i++) {r=(n%2==0 ?"0":"1")+r; n/=2;}
    return r;
}

//---------------------------------------------------------
int countocc(bool**&s,int*&freq, int nsample, int nvertices)
/*
>counts occurrencies of specific rows throughout a bool matrix; 
>INPUT: > bool matrix, with rows already rearranged in increasing
	  order, "s" : [nsample]X[nvertices];
	> "nsample", "nvertices";
>OUTPUT: returns an integer array "f" of length [nsample] of which		DESCRIPTION TO BE UPDATED
 	 only the first [ndiff] elements are different from 0, 
	 where ndiff is the number of different rows observed at 
	 least once; 
	 the first [ndiff] elements of "f" are the # of 
	 occurrecies of different rows in order of first
	 appearance (and thus, in increasing order). 
	 The other elements of "f" are garbage.
*/
{
int k=0;
for (int j=0;j<nsample;j++){
	freq[j]=0;
}
freq[0]=1;
for (int j=1; j<nsample;j++){
	if (eqspin(s[j], s[j-1], nvertices)) {
		freq[k]++;
	} 
	else {
		k=j;
		freq[k]++;
	}
}
clog<< "[countocc] all ok here\n";
return 0;
}

//---------------------------------------------------------
void freqmat(bool**&s, int**&fqmat, int*&freq, int nsample, 
			int nvertices, int ndiff)	
/*
>fills the "frequency matrix" of a list of binary rows.
>INPUT: > bool matrix (dataset) "s" : [nsample]X[nvertices] with 
	  rows already rearranged in increasing order;
	> array "freq" of length [nsample] whose first [ndiff]
	  elements are the # of occurrecies of the different rows
	  seen in the "s" matrix, in order of first appearance
	  (and thus, in increasing order); 
	> "nsample", "nvertices", "ndiff".
OUTPUT:  returns a matrix "freqmat" with dimensions
	 [ndiff] X [nvertices+1]; each row contains as first 
	 [nvertices] elements a specific row of bools (now 
	 integers 0, 1) observed in matrix "s", and as last 
	 element an integer counting the number of occurrencies, 
	 in matrix "s", of that particular row. The rows of 
	 "freqmat" are arranged by the increasing order of
 	 the (1,0) strings, NOT by the number of occurrences.
*/
{
for(int i=0, k=0;i<nsample;i++){
	if(freq[i]!=0){
		for(int j=0;j<nvertices;j++){
			fqmat[k][j]=(int)(s[i][j]);
		}
		fqmat[k][nvertices]=freq[i];
		k++;
	}
}
clog<< "[freqmat] all ok here\n";
return;
}
//---------------------------------------------------------
void quickfs(int**&fqm, int low, int high, 
			int nvertices, int*&temp) 	
/*
>quicksort recursive engine called by a "frequency matrix"
 sorting function;
>INPUT: > a matrix "fqm" with dimensions
	 [ndiff] X [nvertices+1]; each row contains as first 
	 [nvertices] elements a specific row of 0,1 digits, and 
	 as last element an integer representing the number of 
	 occurrencies, in the dataset, of that particular row.
	 The rows of "freqmat" are arranged by the increasing
	 order of the (1,0) strings, NOT by the number of
	 occurrencies.
	> "nvertices";
	> "low", "high": boundaries of the present
	  "window" looked at by the quicksorting algorithm;
	> useless dummy thing "temp" which will soon disappear.
>OUTPUT:> calls itself as many times as necessary; at the final
	  return the selected (low:high) rows of the matrix fqm 
	  have been rearranged, by the increasing order OF
	  OCCURRENCIES.  
	> useless dummy thing "temp" changes irrelevantly.
>during the process of rearrangement the actual integer strings
 are never moved in memory: the function only changes the order
 in which they are pointed at.
*/
{
if(low>=high){
	return;
}
int pivot = high, i=low-1;
for (int j=low;j<=high;j++){
	if (fqm[j][nvertices]>=fqm[pivot][nvertices]){
		i++;
		temp=fqm[i];
		fqm[i]=fqm[j];
		fqm[j]=temp;
	}
}
quickfs (fqm, low, i-1, nvertices, temp); 
quickfs (fqm, i+1, high, nvertices, temp);
return;
}
//---------------------------------------------------------
void quickfreqsort(int**&fqm, int nvertices, int ndiff)   
/*
>takes a "frequency matrix" and sorts its rows by #occurrencies.
>INPUT: > a matrix "fqm" with dimensions
	 [ndiff] X [nvertices+1]; each row contains as first 
	 [nvertices] elements a specific row of 0,1 digits, and 
	 as last element an integer representing the number of 
	 occurrencies, in the dataset, of that particular row.
	 The rows of "freqmat" are arranged by the increasing
	 order of the (1,0) strings, NOT by the number of
	 occurrencies.
	> "nvertices", "ndiff".
>OUTPUT: at return, the rows of the whole matrix fqm 
	 have been rearranged, by the increasing order of their
	 OCCURRENCIES.  
*/
{
int low=0, high=ndiff-1; 
int*temp = fqm[0];
quickfs(fqm, low, high, nvertices, temp);
clog<< "[quickfreqsort] all ok here\n";
return;
}

//----------------------------------------------------------
void quickcoupsort(double**&c, int*&o, int low, int high, 
		   double*&tempd, int&tempi, const bool abss)	
/*
>quicksort recursive engine called by a "coupling constants" 
 sorting function;
>INPUT: > "c" pointer to pointers to single doubles (coupling
 	  constants); it's the thing to be sorted;
	> "o" encodes a specific permutation of [2^nspins]
	  elements which already occurred BEFORE this call;
	> "low, high" boundaries of the present quicksort window;
	> useless dummies "tempd", "tempi";
	> if "abss" ==1 my comparator will work "by absolute 
	  value", if it's ==0 it will work with usual "<=".
>OUTPUT: calls itself as many times as necessary; at the final 
	 return:
	> the selected (low:high) pointers "c[j]" have 
	 been arranged by increasing order of the values they 
	 point at;
	> "o" contains the updated global permutation;
	> useless dummies change irrelevantly. 
*/
{
if(low>=high){
	return;
}
int pivot = high; 
int i = low-1;
if (abss){
	for (int j=low;j<=high;j++){
		if (abs(*c[j])<=abs(*c[pivot])){
			i++;
			tempd=c[i];
			c[i]=c[j];
			c[j]=tempd;
			tempi=o[i];
			o[i]=o[j];
			o[j]=tempi;
		}
	}
}
else{
	for (int j=low;j<=high;j++){
		if ((*c[j])<=(*c[pivot])){
			i++;
			tempd=c[i];
			c[i]=c[j];
			c[j]=tempd;
			tempi=o[i];
			o[i]=o[j];
			o[j]=tempi;
		}
	}
}
quickcoupsort (c, o, low, i-1, tempd, tempi, abss); 
quickcoupsort (c, o, i+1, high, tempd, tempi, abss);
return;
}
//----------------------------------------------------------
void quickintsort(int**&c, int*&o, int low, int high, 
		   int*&tempd, int&tempi, const bool abss)	
/*
>quicksort recursive engine called by a "coupling constants" 
 sorting function;
>INPUT: > "c" pointer to pointers to single doubles (coupling
 	  constants); it's the thing to be sorted;
	> "o" encodes a specific permutation of [2^nspins]
	  elements which already occurred BEFORE this call;
	> "low, high" boundaries of the present quicksort window;
	> useless dummies "tempd", "tempi";
	> if "abss" ==1 my comparator will work "by absolute 
	  value", if it's ==0 it will work with usual "<=".
>OUTPUT: calls itself as many times as necessary; at the final 
	 return:
	> the selected (low:high) pointers "c[j]" have 
	 been arranged by increasing order of the values they 
	 point at;
	> "o" contains the updated global permutation;
	> useless dummies change irrelevantly. 
*/
{
if(low>=high){
	return;
}
int pivot = high; 
int i = low-1;
if (abss){
	for (int j=low;j<=high;j++){
		if (abs(*c[j])>=abs(*c[pivot])){
			i++;
			tempd=c[i];
			c[i]=c[j];
			c[j]=tempd;
			tempi=o[i];
			o[i]=o[j];
			o[j]=tempi;
		}
	}
}
else{
	for (int j=low;j<=high;j++){
		if ((*c[j])>=(*c[pivot])){
			i++;
			tempd=c[i];
			c[i]=c[j];
			c[j]=tempd;
			tempi=o[i];
			o[i]=o[j];
			o[j]=tempi;
		}
	}
}
quickintsort (c, o, low, i-1, tempd, tempi, abss); 
quickintsort (c, o, i+1, high, tempd, tempi, abss);
return;
}
//----------------------------------------------------------
void sortcouplings(double**&c, int*&o, int nc, const bool abss)	
/*
>takes "c" pointer to pointers to coupling constants and
 sorts the pointers c[j] by the increasing order of 
 the values they point at.
>INPUT: > "c" pointer to nv^2 pointers to single double (coupling
 	  constants); it's the thing to be sorted;
	> dummy "o" of length [nv]; 
	> "nv";
	> if "abss" ==1 my comparator will work "by absolute 
	  value", if it's ==0 it will work with usual "<=".
>OUTPUT: calls itself as many times as necessary; at the final 
	 return:
	> the pointers "c[j]" have been arranged by increasing 
	  order of the values they point at;
	> "o" contains the global permutation used to sort;
*/
{
double *tempd = c[1]; 
int tempi=o[1];
int low =1;
int high=nc-1; 
quickcoupsort(c,o, low, high, tempd, tempi, abss);
clog<< "[sortcouplings] all ok here\n";
}
//--------------------------------------------------------------
void sortint(int**&c, int*&o, int nc, const bool abss)	
/*
>takes "c" pointer to pointers to coupling constants and
 sorts the pointers c[j] by the increasing order of 
 the values they point at.
>INPUT: > "c" pointer to nv^2 pointers to single double (coupling
 	  constants); it's the thing to be sorted;
	> dummy "o" of length [nv]; 
	> "nv";
	> if "abss" ==1 my comparator will work "by absolute 
	  value", if it's ==0 it will work with usual "<=".
>OUTPUT: calls itself as many times as necessary; at the final 
	 return:
	> the pointers "c[j]" have been arranged by increasing 
	  order of the values they point at;
	> "o" contains the global permutation used to sort;
*/
{
int *tempd = c[1]; 
int tempi=o[1];
int low =0;
int high=nc-1; 
quickintsort(c,o, low, high, tempd, tempi, abss);
clog<< "[sortcouplings] all ok here\n";
}
//----------------------------------------------------------

void writechi(int**&chi, int chisize, int nvertices, string prefix) 		
/*
>writes the chi matrix in file "chi.txt"
>INPUT:	> integer matrix "chi" of dimensions [2^nvertices]X[Ksize]
	 to be written on file;
	>"Ksize", "nvertices".
>OUTPUT: at return, "chi" is unchanged; a file "chi.txt" has been 
	 written. 
*/
{											
int i=0, j=0, copied=0;
ofstream ochi;
ochi.open(prefix+"chi.txt");
for (j=0;j<pow(2,nvertices);j++){
	for (i=0;i<chisize;i++) {
		ochi << chi[j][i]<<" ";
	} 
	ochi<<endl;
}
ochi.close();
clog<< "[writechi] all ok here\n";
return;
}

//----------------------------------------------------------
void readchi(double **&chi2, int Ksize, int nvertices) 	
/*
>called only if requested (part of "regularization" procedure).
>reads a modified chi matrix to be computed with external 
 tools during a specific idle time of main software's execution.
>INPUT:	>integer matrix "chi2" of dimensions [Ksize]X[nvertices]
	 to be filled;
	>"Ksize", "nvertices".
>OUTPUT: at return, "chi2" is filled with read data.  
*/
{											
int i=0, j=0, into; 
char copied; 
ifstream spinfi;
spinfi.open("chi2.txt");
cout << "opened!" 
     << '\n';
if(!spinfi) {
	clog 	<< "[countdata] error. Abort here";
	cout << "document is not valid" 
	     << '\t' 
	     << endl;
}
for (j=0;j<pow(2,nvertices);j++){
	for (i=0;i<Ksize;i++) {
		spinfi >> chi2[j][i];
	}
}
spinfi.close();
clog<< "[readchi] all ok here\n";
return;
}
//----------------------------------------------------------
void fillphimatrix(int nconf, int**&phi)
/*
>fills the phi[m][n] matrix for the desired graph size,
 the first index [m] denoting a specific monomial spin operator
 (=product of spins corresponding to "1"s in the binary
 representation of m), the second index [n] denoting a specific
 spin configuration (again, determined the trivial way by binary 
 rep of n) for the selected operator to be acted upon.
>INPUT:	> dummy "phi" matrix of dimensions [nconf]X[nconf].
	> "nconf" (-->==2^graphsize)
>OUTPUT: > at return, phi is filled with the outcomes of the 
	   operator-on-state actions described above.
*/
{
int ba, bb;
for (int ii=0, i=0, re=0; i<nconf; i++){
	ii=i;
	for (int jj=0, j=0; j<nconf; j++){
		ii=i;jj=j; re=1;
		for (int c=0;c<9;c++){
			ba=ii%2;
			bb=(jj%2)*2-1;
			ii/=2;
			jj/=2;
			if (ba==1)re*=bb;
		}					
		phi[i][j]=re;
	}
}
clog<< "[fillphimatrix] all ok here\n";
return;
}

//----------------------------------------------------------
int countKsize (int ndiff, int nvertices, int**fqm,
				 int*&TempIndex)
/*
>counts the number of sets constituting the "frequency partition" 
 for the dataset under study;
>INPUT:	> frequency matrix "fqm" (see output of "freqmat(...)"),
	  of dimensions [ndiff]X[nvertices+1];
	> "nvertices", "ndiff";
	> dummy int TempIndex[] of dimension [ndiff].
>OUTPUT:> returns the number of sets in the freq partition;
	> at return, dummy input TempIndex of dimension [ndiff] 
	  has been filled this way:
	  **first [Ksize] values --> a list of the indexes of 
	  the starting (fqm) rows for the different partitions.
 	  **[Ksize+1] to [ndiff] values --> garbage.
*/
{
int Ksize=1;
TempIndex[0]=0; 
for(int j=1;j<ndiff;j++){
	if (fqm[j][nvertices]!=fqm[j-1][nvertices]){
		TempIndex[Ksize]=j;
		Ksize++;
	}
}
cout <<"\n" 
     <<"[countKsize] frequency partition is composed of " 
     << Ksize 
     << "+"
     << (int)(ndiff!=(int)pow(2,nvertices))
     <<" sets.\n\n";
clog<< "[countKsize] all ok here\n";
return Ksize;
}

//----------------------------------------------------------
void fpartition(int ndiff, int nvertices, int**fqm,
		 int*&TempIndex, int*&ind, int Ksize, string prefix, char xplorer)
/*
>creates the correct "index list" to identify rows belonging to 
 different K-parition sets in the frequency matrix; computes and 
 writes on file the "merged" frequencies relative to K-sets.
>INPUT:	> frequency matrix "fqm" (see output of "freqmat(...)"),
	  of dimensions [ndiff]X[nvertices+1];
	> "TempIndex" as given by the function "countKsize";
	> dummy "ind" of length [Ksize+1] to be filled;
	> "nvertices", "ndiff", "Ksize".
>OUTPUT: at return: 
	> the not-garbage contents of "TempIndex" have 
	  been copied to constitute the whole contents of "ind"; 
	  last entry of "ind" is = [ndiff];
	> a file "freqss.txt" has been created and filled with an 
	  array of frequencies representing #occurrencies(rows 
	  from a specific set of the K partition).
*/
{
for (int l=0;l<Ksize;l++){
	ind[l]=TempIndex[l];
} 
ind[Ksize] = ndiff;
int*Kk=new int[Ksize]; 			
for (int u=0;u<Ksize;u++){
	Kk[u]=(fqm[ind[u]][nvertices])*(ind[u+1]-ind[u]);
}
ofstream kk; 
if (xplorer!='y') {
	kk.open(prefix+"freqss.txt");
	for (int y=0;y<Ksize;y++){
		kk<<Kk[y]<<endl;
	}
	kk.close();
	kk.open("k1_N.txt",ios_base::app);
	kk<<Kk[0]<<"\t";
	kk.close();
}

delete [] Kk;
clog<< "[fpartition] all ok here\n";
return;
}
//-----------------------------------------------------------------------
void allpartitions(double*&PostRatio, double*&LogPostRatio, int nv, int nconf, int N, int ndiff, int**&fqm, double a)
{
	ofstream testq;
	testq.open("testq.txt");
	int go,l,j;
	int*b=new int[nconf];
	for (int i=0;i<nconf;i++){
		b[i]=0;
	}
	int Qsize=1;
	int nQ=pow(2,(nconf-1));
	//***generate barrier-type partitions
	for (int p=1; p<nQ;p++){
		//***generate a partition	
		Qsize=1;
		for (int l=0;l<nconf-1;l++){
			b[nconf-1-l]=((p%int(pow(2,l+1)))>=int(pow(2,l)));
			if (b[nconf-1-l]!=0){
				Qsize++;
			}
		}
		for (int u=0;u<nconf;u++){testq<<b[u]<<" ";}
		testq<<endl;
		//***compute K[],m[] for selected partition
		int*K=new int[Qsize];
		int*m=new int[Qsize];
		for (int i=0;i<Qsize;i++){
			K[i]=0;
			m[i]=0;
		}

		l=0;
		j=0;
		while(1){
			if(b[l]==1){
				j++;
			}
			m[j]++;
			if(l<ndiff){
				K[j]+=fqm[l][nv];
			}
			l++;
			if(l==nconf){
				break;
			}
		}
		for (int i=0;i<Qsize;i++){testq<<K[i]<<'\t'<<m[i]<<endl;}
		LogPostRatio[p]=N*log(nconf)+lgamma(a*Qsize)-lgamma(a*Qsize+N);
		for (int j=0;j<Qsize;j++){
			LogPostRatio[p]+=lgamma(K[j]+a)-lgamma(a)-K[j]*log(m[j]);
		}
		LogPostRatio[p]-=Qsize;	
		/*PostRatio[p]=pow(nconf,N)*(tgamma(Qsize)/tgamma(Qsize+N));
		for (int j=0;j<Qsize;j++){
			PostRatio[p]*=(tgamma(K[j]+1)/pow(m[j],K[j]));
		}*/
		delete[]K;
		delete[]m;
	}
	return;
}
//-----------------------------------------------------------------------

int qstar (int nconf, int nsample, int*&ind, int &Ksize, int*&K, double a, int demo)
/*
*/
{
//find best merge
double*gain=new double[Ksize];
int stop=0;
int win=0;
double newk;
int nobs=nconf-ind[Ksize];
int allk=Ksize+(nobs!=0);
for(int i=0;i<Ksize-1;i++){
	gain[i]=(lgamma(K[i]+a)+lgamma(K[i+1]+a)-lgamma(K[i]+K[i+1]+a)-lgamma(a));
	gain[i]+=K[i]*log(1+(double)(ind[i+2]-ind[i+1])/(ind[i+1]-ind[i]));
	gain[i]+=K[i+1]*log(1+(double)(ind[i+1]-ind[i])/(ind[i+2]-ind[i+1]));
	gain[i]-=(lgamma(a*(allk-1))+lgamma(nsample+a*allk));
	gain[i]+=(lgamma(a*allk)+lgamma(nsample+a*(allk-1)));
	//gain[i] is computed
	if(gain[i]<=gain[win]){
		win=i;
	}
}
if(nobs>0){
	gain[Ksize-1]=0;
	gain[Ksize-1]+=K[0]*log(1+(double)nobs/(ind[Ksize]-ind[Ksize-1]));
	cout<<"kweak "<<K[0]<<"  nobs "<<nobs<<" frac "<<(ind[Ksize]-ind[Ksize-1]) <<endl;	
	cout<<"gain for unobs is "<<gain[Ksize-1]<<endl;
	gain[Ksize-1]-=(lgamma(a*(allk-1))+lgamma(nsample+a*allk));
	cout<<"gain for unobs is "<<gain[Ksize-1]<<endl;	
	gain[Ksize-1]+=(lgamma(a*allk)+lgamma(nsample+a*(allk-1)));
	cout<<"gain for unobs is* "<<gain[Ksize-1]<<endl;
	if(gain[Ksize-1]<=gain[win]){
		win=Ksize-1;
	}
	if (demo==1){
		if(gain[Ksize-1]<0){
			win=Ksize-1;
		}
	}
}
//merge
if (win!=Ksize-1){
	if(gain[win]<0){
		K[win]+=K[win+1];
		for (int i=win+1;i<Ksize;i++){
			ind[i]=ind[i+1];
			
		}
		for (int i=win+1;i<Ksize-1;i++){
			K[i]=K[i+1];
		}
		K[Ksize-1]=0;
		ind[Ksize]=0;
		Ksize--;
	}
	else{
		stop=1;
	}
}
else{
	if(gain[Ksize-1]<0){
		ind[Ksize]=nconf;
		cout <<"it's sample point n°"<<nsample<<" and I'm merging unobserved states\n";
	}
	else{
		stop=1;
	}
}
delete[]gain;	
return stop;
//
}

//----------------------------------------------------------
void fillchimatrix(int nconf, int Ksize, int**&chi, 
			int*ind, int nv, int**fqm, int ndiff)
/*
>INPUT:	> dummy "chi" matrix of dimensions [nconf]X[Ksize];
	> index list "ind" for specific rows of "fqm" (-->which 
	  correspond to starting points of K-partition sets),
	  of length [Ksize+1]; 
	> frequency matrix "fqm" (see output of "freqmat(...)"),
	  of dimensions [ind[Ksize]]X[nvertices+1];
	> "nv", "Ksize", "nconf".
>OUTPUT:> at return, "chi" has been filled; its first index 
	  denoting (via its nv binary rep digits) a specific spin 
	  operator, the second index denoting a specific set of 
	  the frequency partition; the spin operator acts on each 
	  of the configs contained in this set and the outcomes 
	  are summed to become the corresponding entry of "chi".
*/
{
int yy=0, prod, x=0;
for (int y=0;y<nconf;y++) {
	chi[y][Ksize]=0;
	for (int z=0;z<Ksize-1;z++){
		chi[y][z]=0;
		for(int k=ind[z]; k<ind[z+1];k++){
			prod=1;
			yy=y;			
			for (int n=0;n<nv;n++){
				x=(yy>=(pow(2,nv-1-n))); 
				if (x==1){
					prod*=(fqm[k][n]*2-1);
					yy-=(pow(2,nv-1-n));
				}  
			}
			if(yy!=y){
				chi[y][z] +=prod;
				chi[y][Ksize]-=prod;
			}
		}
		for(int k=ind[Ksize-1]; k<ndiff;k++){
			prod=1;
			yy=y;			
			for (int n=0;n<nv;n++){
				x=(yy>=(pow(2,nv-1-n))); 
				if (x==1){
					prod*=(fqm[k][n]*2-1);
					yy-=(pow(2,nv-1-n));
				}  
			}
			if(yy!=y){
				chi[y][Ksize-1]+=prod;
				chi[y][Ksize]-=prod;
			}
		}
	}
	if(ind[Ksize]>ndiff){				
		chi[y][Ksize-1]+=chi[y][Ksize];
		chi[y][Ksize]=0;
	}
}
for (int j=1;j<Ksize+1;j++){
	chi[0][j-1]=ind[j]-ind[j-1];
}
chi[0][Ksize]=nconf-ind[Ksize];
clog<< "[fillchimatrix] all ok here\n";
return;
}
//----------------------------------------------------------
void computecouplings (int nconf, int ndiff, double**&couplings, 
			int Ksize,int**chi, int*K, 
			int nvertices, int nsample, int*ind, 
			double**&couplingsa, double**&couplingsg,
			double**&couplingsga, string prefix)
/*
>uses the GM formula to compute the coupling constants g[mu];
>INPUT:	> dummies: "couplings", "couplingsg", "couplingsg", all 
	  pointers to pointers to single doubles ([nconf]X[1]), 
	  all to be filled;
	> frequency matrix "fqm" (see output of "freqmat(...)"),
	  of dimensions [ndiff]X[nvertices+1];
	> chi matrix "chi" (see output of "fillchimatrix(...)"),
	  of dimensions [nconf]X[Ksize];
	> fqm-index list for K-partition: "ind" (see output of
	  fpartition(...), of length [Ksize+1];  
	> "nconf", "ndiff", "Ksize", "nvertices", "nsample".
>OUTPUT:> at return:
	> computed g[mu] are stored identically in 
	  "*couplings[mu]" and in "*couplingsa[mu]";
	> special couplings gg[mu] (-->obtained faking 1 
	  observation for each unobserved state are stored in 
	  "*couplingsg[mu]");
>NB: as we're not allowed to take the log of 0, a global double 
 constant of very high value is used to simulate the divergence 
 to occur. 
*/
{

double xplode=log(REGULARIZER/nsample);
double q;
double corr;
double freqo;

ofstream coups;
//1.when unobserved explode:
coups.open(prefix+"mucouplings.txt"); 
for (int y=1;y<nconf;y++){
	*couplings[y]=0;
	for(int o=0; o<Ksize;o++){
		freqo=(double)K[o]/nsample;
		clog<<"[compcoups]freqo "<<freqo<<endl;			
		*couplings[y]+=(chi[y][o] * log(freqo))/nconf;
	} 
	corr=(chi[y][Ksize] * xplode/nconf);
	*couplings[y]+=corr;
	coups <<*couplings[y]<<endl;	
}
coups.close();
for (int mu=1;mu<nconf;mu++){
	*couplingsa[mu]=*couplings[mu];
}
//2.when unobserved are assigned 1 fake observation each:	
coups.open(prefix+"mucouplings_faked.txt"); 
double qg, corrg;
int nchigh=nsample+nconf-ndiff;
for (int y=1;y<nconf;y++){
	*couplingsg[y]=0;
	for(int o=0; o<Ksize;o++){
		freqo=(double)K[o]/(nchigh);
		*couplingsg[y]+=(chi[y][o] * log(freqo))/nconf;
	} 
	corrg=(chi[y][Ksize] * log(1./(nchigh)))/nconf;
	*couplingsg[y]+=corrg;
	coups <<*couplingsg[y]<<endl;	
}
coups.close();
for (int mu=1;mu<nconf;mu++){
	*couplingsga[mu]=*couplingsg[mu];
}
clog<< "[computecouplings] all ok here\n";
return;
}
//----------------------------------------------------------
template <typename T> 
double gmpm(double*g,int s, T**&phi, int nop)
{
//clog<<"[gmpm]called successfully: s="<<s<<";\n";
double gmuphimu=0;
for (int mu=1;mu<nop;mu++){
	gmuphimu+=((g[mu])*phi[mu][s]);
//	clog<<"[gmpm]g["<<mu<<"]*phi["<<mu<<"]["<<s<<"];\n";	
}
//clog<<"[gmpm] all ok here\n";
return gmuphimu;
}
//----------------------------------------------------------
template <typename T> 
double ensemble(int nconf, int nop, double*&muc,
		 T**&phi, double*&pos, double*&avgphi)
/*
> computes the canonical probabilities p(s|g), and the ensemble  
  averages of spin operators under that measure;
>INPUT:	> "muc" of length [nop] (-->coupling constants of the 
	  model);
	> dummy "pos" of length [nconf];
	> phi matrix "phi" (see output of "fillphimatrix(...)"),
	  of dimensions [nop]X[nconf];
	> "nconf", "nop";
	> dummy "avgphi", of length [nop], to be filled.
>OUTPUT:> returns the value f the partition function Z(g);
	> "pos" is filled with probabilities p(s|g) 
	  (elements are arranged, the trivial way, by '1's in the 
	  binary representation of the array index representing 
	  the "+" spins of the corresponding state "s");
*/
{
clog<<"[ensemble]*called successfully;\n";
double Zg=0;
double*gmuphimu = new double[nconf];
double tesst=0; 
clog<<"[ensemble]*declarations successful;\n";
for (int s=0;s<nconf;s++){
	gmuphimu[s]=gmpm(muc,s,phi,nop);
}
clog<<"[ensemble]*gmpms;\n";
for(int s=0;s<nconf;s++){
	Zg+=exp(gmuphimu[s]);
} 
clog<<"[ensemble]*Zg built;\n";
for (int s=0;s<nconf;s++){
	pos[s]=exp(gmuphimu[s])/Zg; 
	tesst+=pos[s];
}
clog<<"[ensemble]*pos built;\n";
//cout << "tesst: " 
//     << tesst 	//[test] this must be =1!!!
//     << endl; 
for (int mu=0;mu<nop;mu++){
	avgphi[mu]=0;
	for (int s=0;s<nconf;s++){
		avgphi[mu]+=phi[mu][s]*pos[s];
	}
}
clog<<"[ensemble]*avgs built;\n";
delete [] gmuphimu;
clog<< "[ensemble] all ok here\n";
return Zg;
}
//----------------------------------------------------------
template <typename T> 
void empiricalaverages(int nop, int nv, int ndiff, T**&phi,
			int**&fqm, int nsample, double*&empavg)
/*
>computes the empirical averages of all spin operators "phi[mu]".
>INPUT:	> frequency matrix "fqm" (see output of "freqmat(...)"),
	  of dimensions [ndiff]X[nv+1];
	> dummy "empavg", of length [nop], to be filled;
	> "nsample", whose value must coincide with the sum of 
	  the last column values of "fqm";
	> "ndiff", "nv", "nconf";
>OUTPUT:> at return, "empavg" has been filled with the average 
	  value of all the monomial spin operators across the 
	  dataset encoded in "fqm". 
*/
{
//THIS FUNCTION IS INCOMPLETE
for(int lam=0;lam<nop;lam++){
	empavg[lam]=0;
	for(int s=0;s<ndiff;s++){	
		int label=0;
		for(int c=0;c<nv;c++){
			label+=fqm[s][nv-1-c]*(int)pow(2,c);
		}
		empavg[lam]+=phi[lam][label]*(double)fqm[s][nv]/nsample;
	}
}
clog<< "[empiricalaverages] all ok here\n";
return;
}
//----------------------------------------------------------
void sortandcoutcoup(int nv, int nconf, int*&order, 
			double**&couplings, 
			bool abs, string filename)
/*
>sorts couplings and keeps track of the resulting global 
 permutation.
>INPUT:	> "couplings", pointer to pointers to single doubles 
	  ([nconf]X[1]), these doubles being the MLE coupling 
	  constants as given by GM formula;
	> "filename" of the file I want to write in the sorted 
	  couplings;
	> dummy "order", of length [nconf].
	> "abs" whose logical value, when positive, enforces 
	  sorting by absolute value;
	> "nconf".
>OUTPUT: at return:
	> the pointers "couplings[j]" have been rearranged by the 
	  ascending order of their pointed value;
	> "order" now contains, for each coupling, the place it 
	  occupied before the rearrangement took place;
	> sorted "couplings" & "order" are written on file; "order
	  [j] are printed in their binary rep, so that they can 
	  be read directly as identifying the corresponding 
 	  monomial spin operators.
*/
{
ofstream coups;
for (int i=0; i<nconf; i++){
	order[i]=i;
}
sortcouplings(couplings, order, nconf, abs);
coups.open(filename);
for (int y=1;y<nconf;y++){
	coups << *couplings[y] 
	      << "\t\t\t"; 
	for (int u=nv;u>0;u--){
		coups << (order[y]%(int)pow(2,u))/(int)pow(2,u-1);
	}
	coups << endl;
}
coups.close();	
clog<< "[sortandcoutcoup] all ok here\n";
return;
}
//----------------------------------------------------------
void findglam(	int nconf, int nrel, int nsample, double**&relg, double like,
			double**&psi, double*&remav,double*deltaloglik, double step,
			double*loglik, string prefix){
//---
int ppp=0;
int maxsteps=100000000;
nrel++;
int prevwin=0;
int prevsign=0;
double*avgpsi0=new double[nrel];
for (int lam=0;lam<nrel;lam++){
	avgpsi0[lam]=0;
}
double*delta=new double[nrel];
for (int l=0;l<nrel;l++){
	delta[l]=0;
}
int sign=0;
int*fix=new int[nrel];
for (int y=0;y<nrel;y++){
	fix[y]=0;
}
double*pos=new double[nconf];
for(int s=0;s<nconf;s++){
	pos[s]=0;
}
double*dot=new double[nrel];
double lik;
int win=0;
ofstream likk;
//---
likk.open("lik.txt");
int search;
int found =0;
int perc=10;
search=1; 
for (int lam=0;lam<nrel;lam++){
	dot[lam]=0;
} 
while(search){
	ppp++; 
	if(ppp>maxsteps){
		cout<< "TOOMANY";
		break;
	}
	if(ppp%10000==0)cout<<".";
	if(ppp%(maxsteps/100)==0){
		cout<<endl<<perc<<"%\n";
		perc+=10;
	}
	ensemble(nconf,nrel,dot,psi,pos,avgpsi0);
	delta[0]=0;
	prevwin=win;
	found=0;	
	for (int lam=1;lam<nrel;lam++){
		delta[lam]=remav[lam]-avgpsi0[lam]; 
		if(abs(delta[lam])>=abs(delta[win])){
			if (fix[lam]==0){
				if(ppp>maxsteps-1000){likk<<"\nupdatewin\n";}			
				win=lam;
				found=1;
			}
		}
	}
	if(ppp>maxsteps-1000){likk<<"\nfound is "<<found << endl;}
	prevsign=sign;		
	sign=(int)(delta[win]/abs(delta[win]));
	if(ppp>maxsteps-1000){likk<<prevwin<<'\t'<<win<<'\t'<<prevsign<<'\t'<<sign<<'\t'<<delta[win]<<endl;}
	if(abs(delta[win])<0.0000001){
		cout<<"WE'RE CLOSE";
		break;
	}
	dot[win]=dot[win]+sign*step;
	if(ppp>2&&win==prevwin&&sign!=prevsign){
		fix[win]=1;
	}
	if(found==0){
		cout<<"WE'RE THERE";
		break;
	}
}	
for (int lam=1;lam<nrel;lam++){
	*relg[lam]=dot[lam];
}
lik=0;
for (int lam=1;lam<nrel;lam++){
	lik+=dot[lam]*remav[lam];
}
lik-=log(ensemble(nconf,nrel,dot,psi,pos,avgpsi0));
lik*=nsample;
*loglik=lik;
*deltaloglik=lik-like;
likk<<"loglik under optimal glam is "<<*loglik<<endl;
likk<<"Lloss under optimal glam is "<<*deltaloglik<<endl;
likk.close();

//		free memory
delete[]delta;
delete[] avgpsi0;
delete []dot; 
delete[]pos;
delete[]fix;
clog<< "[findglam] all ok here\n";
return;
}



//----------------------------------------------------------
void rankrelcouplings(	int nconf, int nrel, double**&relg, double like,
			double**&psi, double*&remav,double**&deltalik, double step,
			double**&lik0, string prefix){
nrel++;
double**avgpsi0=new double*[nrel];
for(int k=0;k<nrel;k++){
	avgpsi0[k]=new double[nrel];
	for (int l=0;l<nrel;l++){
		avgpsi0[k][l]=0;
	}
}
double**delta=new double*[nrel];
for(int k=0;k<nrel;k++){
	delta[k]=new double[nrel];
	for (int l=0;l<nrel;l++){
		delta[k][l]=0;
	}
}
double*pos=new double[nconf];
for(int s=0;s<nconf;s++){
	pos[s]=0;
}
double*dot=new double[nrel];

double lik;
int win;
ofstream likk;
likk.open("likk.txt");
int search;
for (int k=1;k<nrel;k++){
	search=1;
	cout << "analyzing coupling n° "
	     << k
	     << "\n"; 
	for (int lam=0;lam<nrel;lam++){
		dot[lam]=0;
	} 
	dot[k]=0;
	while(search){
		ensemble(nconf,nrel,dot,psi,pos,avgpsi0[k]);
		win=0;
		delta[k][0]=0;
		delta[k][k]=0;
		for (int lam=1;lam<nrel;lam++){
			if(lam!=k){
				delta[k][lam]=remav[lam]-avgpsi0[k][lam]; 
			}
			if(abs(delta[k][lam])>abs(delta[k][win])){
				win=lam;
			}
		}
		if(abs(delta[k][win])<0.00001){
			search=0;
			break;
		}		
		dot[win]=dot[win]+step; 
	}	
	if (k==0) {
		for (int lam=1;lam<nrel;lam++){
			*relg[lam]=dot[lam];
		}
	}
	lik=0;
	for (int lam=1;lam<nrel;lam++){
		lik+=dot[lam]*remav[lam];
	}
	lik-=log(ensemble(nconf,nrel,dot,psi,pos,avgpsi0[k]));
	likk<<lik<<endl;
	*lik0[k]=lik;
	*deltalik[k]=lik-like;
}
likk.close();

//		free memory
for(int k=0;k<nrel;k++){
	delete[]avgpsi0[k];
}
delete[] avgpsi0;
delete []dot; 
delete[]pos;
clog<< "[rankrelcouplings] all ok here\n";
return;
}


//----------------------------------------------------------
void testdeltaphis(int nconf, double*&empavg, double*&avgphi, string prefix)
/*
>TEST
*/
{
ofstream coups;
coups.open(prefix+"deltaphis.txt");
double *dphi=new double [nconf]; 
for (int mu=0;mu<nconf;mu++){
	dphi[mu]=empavg[mu]-avgphi[mu];
	coups<<dphi[mu]<<endl;
} 
coups.close();
delete[]dphi;
}
//----------------------------------------------------------
void svdchi	(	int nconf, int**&phi, double**&psi, 
			double **&Umat, double **&Wmat,
			double*&LL, string prefix, int nspec)
{
char anyany;
cout << "I'm not good at performing SVD decompositions :("
     << endl;
cout << "take \"chi.txt\" and "
     << "go get your SVD decomposition elsewhere\n" 
     << "(I heard of a fine MATLAB guy named reg_ok.m...);\n" 	
     << "when ready, put the obtained matrices in this folder:\n"
     << "I'll look for files named \"Umat.txt\",\n"
     << "\"Wmat.txt\", and \"Lam.txt\".\n\n"
     << "When you're done and ready, hit a random key:\n";
cout<<"HEYYYY";
cin >> anyany;
cout<<"HEYYYY";
ifstream icou;
icou.open(prefix+"Umat.txt");
for (int a=0;a<nconf-1;a++){
	for(int b=0;b<nconf-1;b++){
		icou>>Umat[a][b];
	}
}
cout<<"HEYYYY";
icou.close();
icou.open(prefix+"Wmat.txt");
for (int a=0;a<nspec;a++){
	for(int b=0;b<nspec;b++){
		icou>>Wmat[a][b];
	}
}
icou.close();
icou.open(prefix+"LL.txt");
for (int a=0;a<nspec;a++){
	icou>>LL[a];
}
cout<<"HEYYYY";
icou.close();
for(int s=0;s<nconf;s++){
	psi[0][s]=1;
}
cout<<"HEYYYY";
for(int lam=1;lam<nspec+1;lam++){
	for(int s=0;s<nconf;s++){
		psi[lam][s]=0;
		for(int mu=1;mu<nconf;mu++){
			psi[lam][s]+=Umat[mu-1][lam-1]*phi[mu][s];
		}
	}
}
return;	
}

//----------------------------------------------------------
int readrealcouplings(int nvertices, int nconf, double*&greal, 
				int*&tru,char xplorer,string prefix)
/*
>to call if the dataset was sinthetic; this function reads from 
 file the real couplings that originated the sample, so that it 
 can judge the degree of inference success.
*/
{
ifstream inn; 
if (xplorer=='y'){
	inn.open(prefix+"g.txt");
	for(int u=0;u<nconf;u++){
		inn>>greal[u];
		tru[u]=(greal[u]*1000)/1;
	}
	return 0;
}
inn.open("decidedinteractions.txt");
int nv=0; 
inn>>nv; 
if(nv!=nvertices){
	cout << "sample I read does not match model given. Error."
	     << endl; 
	return 1;
}
for (int i=0;i<(int)pow(2,nv);i++){
	greal[i]=0;
}
for (int i=0;i<nconf;i++){
	tru[i]=0;
}
int a=1, indx=0;
inn>>a;
while (a!=111)
{
	indx=0;
	while (a!=0){
		if (a>nv || a<0) {
			cout << "the g list is not well written. "
			     << "Check it."; 
			return 2;
		}
		indx+=(int)pow(2,nv-a);
		inn >> a;
	}
	inn>>greal[indx]; 
	tru[indx]=((greal[indx]))/1;
	inn>>a;
}
inn.close();
clog<< "[readrealcouplings] all ok here\n";
return 0;
}
//----------------------------------------------------------
void coutfinalcoup(int nv, int nconf, double**&couplings, 
			int*&order, 
			int*&tru, double*&greal, bool dir, 
			string filename) 
/*
*/
{
ofstream coups;	
coups.open(filename);
for(int y=1, l=0;y<nconf;y++){
	l=((int)dir*(nconf-1))+((-1)*(((int)dir*2)-1)*(y-(int)dir));
	coups << fixed<<setprecision(7)<<*couplings[l];
	coups << "\t\t\t";
	for (int u=nv;u>0;u--){
		coups << (order[l]%(int)pow(2,u))/(int)pow(2,u-1);
	}
	coups << "\t\t\t"
	      << tru[order[l]];
	coups << "\t"
	      << fixed<<setprecision(2)<<greal[order[l]];
	coups << endl;	
}
coups.close();
clog<< "[coutfinalcoup**] all ok here\n";
return;
}
//----------------------------------------------------------
void coutfinalcoup(int nv,int nconf, double*&couplings, 
			int*&order,
			int*&tru, double*&greal, bool dir, 
			string filename) 
/*
*/
{
ofstream coups;	
coups.open(filename);
for(int y=1, l=0;y<nconf;y++){
	l=((int)dir*(nconf-1))+((-1)*(((int)dir*2)-1)*(y-(int)dir));
	coups << fixed<<setprecision(7)<<couplings[l] 
	      << "\t\t\t";
	for (int u=nv;u>0;u--){
		coups << (order[l]%(int)pow(2,u))/(int)pow(2,u-1);
	}
	coups << "\t\t\t"
	      << tru[order[l]]
	      << "\t"
	      << fixed<<setprecision(2)<<greal[order[l]]
	      << endl;	
}
coups.close();
clog<< "[coutfinalcoup*] all ok here\n";
}
//----------------------------------------------------------
//----------------------------------------------------------
int looper(int nv, string filename, bool vsr=0)
/*
>SUPER MESSY
*/
{
int l,i,co; 		//declarations 
int nc=pow(2,nv); clog<<"nc = "<<nc<<endl;
double *d=new double[nc-1]; 
double*gre=new double[nc-1];
int *tr=new int[nc-1];
string* state=new string[nc-1];
string* loopsum=new string[nc-1];
string zer = "0";
for(int y=0;y<nv-1;y++){
	zer+="0";
}
for(int i=0;i<nc-1;i++){
	loopsum[i]="";
}
loopsum[0]="G0";
int**gens=new int*[nv]; 
for(int i=0;i<nv;i++){
	gens[i]=new int[nv];
}
int**s=new int*[nc-1]; 
for(int i=0;i<nc-1;i++){
	s[i]=new int[nv+1];
}
int***bb=new int**[nv];
for(int i=0;i<nv;i++){
	bb[i]=new int*[(int)pow(2,i)];
	for(int j=0;j<(int)pow(2,i);j++){
		bb[i][j]=new int[nv];	
	}
  }	
ifstream fin;

fin.open(filename+".txt");	//read
for (int i=0;i<nc-1;i++){
	fin>>d[i];
	fin>>state[i]; 
	clog <<state[i]<<"\t\t"
	     <<state[i].length()<<endl;		//
	if (state[i]==zer){
		//cout << "discarding " <<zer<<endl;
		if (vsr){
			fin>>tr[i]>>gre[i];
		}
		i--;
		continue;
	}
	for (int b=0;b<nv;b++){
		s[nc-2-i][b]=state[i][b]-'0'; clog<<s[nc-2-i][b]; 
	}
	clog<<endl;
	if(vsr){
		fin>>tr[i]>>gre[i];
	}
}
fin.close();	

i=0;
int next=0; 
int found;
for (int l=0;l<nv;l++){			//find nv generators
	i=0;
	found=2;
	while (found!=0) {
		found=2;
		if(l==0){
			//cout<<"CREATE FIRST BUBBLE\n";
			s[next][nv]=0;			
			break;
		}
		for (int ll=0;ll<l;ll++){
			for (int a=0;a<pow(2,ll);a++){
				for (int b=0;b<nv;b++){
					if (s[next][b]!=bb[ll][a][b]){	//"if s[next]!=bb[ll][a]
						break;
					}
					if (b==nv-1){
						found=1;
						s[next][nv]=ll;
						loopsum[next]+=("G"+to_string(ll));
						co=0;
						for (int y=0;y<nv;y++){
							if (((a%(int)pow(2,y+1))/(int)pow(2,y))==1){
								loopsum[next]+=("+G"+to_string(y));
							}
						}
					}
				}
				if(found==1) break;	
			}
			if(found==1) {
				next++;
				break;
			}
			if(ll==l-1){
					clog << "NOT FOUND!\n\n";
				s[next][nv]=l;
				loopsum[next]+=("G"+to_string(l));
				found=0;
			}
		}
	}
	
	for (int b=0;b<nv;b++){
		gens[l][b]=s[next][b];
		bb[l][0][b]=gens[l][b];
	}
	i=1; 
	for (int ll=0;ll<l;ll++){
		for (int a=1;a<pow(2,ll)+1;a++){
			for (int b=0;b<nv;b++){
				bb[l][i][b]=(bb[l][0][b]!=bb[ll][a-1][b]); 
			}
			i++;	
		}
	}
	next++;
}
//cout << "classifichiamo"<<endl;				//	all the rest
found=2; int end=0;

while (found!=0 && end!=1) {
	found=2;
	for (int ll=0;ll<nv;ll++){
		for (int a=0;a<pow(2,ll);a++){
			for (int b=0;b<nv;b++){
				if (s[next][b]!=bb[ll][a][b]) {break;}
				if (b==nv-1){
					found=1;
					s[next][nv]=ll;
					loopsum[next]+=("G"+to_string(ll));
					for (int y=0;y<nv;y++){
						if (((a%(int)pow(2,y+1))/(int)pow(2,y))==1){
							loopsum[next]+=("+G"+to_string(y));co++;
						}						
					}
					}
			}
			if(found==1) break;	
		}
		
		if(found==1) {
			if(next==pow(2,nv)-2){
				//cout<<"abbiamo finito!"<<endl;
				end=1;
				break;
			}
			next++;
			break;
		}
		if(ll==nv-1){clog << "NOT FOUND! strange...\n\n";found=0;}
	}
}
ofstream add;
add.open(filename+"_looped.txt");
add <<"inferred g"<<"\t"<<"operator"<<"\t\t";
if(vsr){
	add<<"r1"<<"\t"<<"real coupling"<<"\t";
}
add <<"generators involved"<<"\t\t"<<"in family:";
add << endl;
for (int i=0;i<nc-1;i++){
	add<<fixed<<setprecision(7)<<d[i]<<"\t";
	add<<state[i]<<"\t";
	if(vsr){
		add<<tr[i]<<"\t"<<gre[i]<<"\t";
	}
	add <<loopsum[nc-2-i]<<"\t";
	if (loopsum[nc-2-i].length()<8){
		add<<"\t\t\t";
	}
	else if (loopsum[nc-2-i].length()<16){
		add<<"\t\t";
	}
	else if (loopsum[nc-2-i].length()<24){
		add<<"\t";
	}
	for (int k=0;k<s[nc-2-i][nv];k++){add<<"---";}
	add<<"---"<<s[nc-2-i][nv]<<endl;
}
add.close();
//cout << "execution successful";
delete[]d;
delete[]tr;
delete[]state;
delete []gre;
delete []loopsum;
for(int i=0;i<nv;i++){
	delete[]gens[i];
}
delete []gens;
for(int i=0;i<nc-1;i++){
	delete[]s[i];
}
delete []s;
for(int i=0;i<nv;i++){
	for(int j=0;j<(int)pow(2,i);j++){
		delete[]bb[i][j];	
	}
	delete[]bb[i];
  }
delete []bb;
clog<< "[looper] all ok here\n";
return 0;
}

//---------------------------------------------------------

void plotorders(int nv, string*&plots, int nplots, string prefix){
int nc=0, bb=0; 
char a;
ifstream gimu;
ofstream gemu;
clog <<"\n[plotorders] "<<(nc=(int)pow(2,nv));
double*gmu=new double[nc];
double*mrk=new double[nc]; 
int *o=new int[nc];
for (int p=0;p<nplots;p++){
	gimu.open(prefix+"gmu"+plots[p]);
	for (int i=0;i<nc;i++){
		gimu>>gmu[i];
		o[i]=0;
		for(int j=0;j<nv;j++){
			gimu>>a;
			if (a=='1'){
				o[i]++;
			}
		}
		if (o[i]==0) {
			gmu[i]=0;
		}
	}
	gimu.close();

	gemu.open(prefix+"orders"+plots[p]);
	for (int mu=0;mu<nc;mu++){
		gemu<<o[mu]<<"\t"<<gmu[mu]<<endl;
	}
	gemu.close();
}
/*
gimu.open("scaledrelevance.txt");
for (int i=0;i<nc;i++){
	o[i]=0; 
	gimu>>mrk[i];
	for(int j=0;j<nv;j++){
		gimu>>a;
		cout<<a << " ";
		if (a=='1'){
			o[i]++;
		}
	}
cout << o[i]<<endl;
}
//bb=i;for (int b=0;b<nv;b++){o[i]+=bb%2;bb/=2;}}
gimu.close();

gemu.open("markes.txt");
for (int mu=0;mu<nc;mu++){
	gemu<<o[mu]<<"\t"<<mrk[mu]<<endl;
}
gemu.close();
*/

/*
for (int k=0;k<nv;k++)
{
gemu.open("order"+std::to_string(k+1)+".txt");
for (int mu=0;mu<nc;mu++){if(o[mu]==k+1){gemu<<gmu[mu]<<endl;}}
gemu.close();
}
*/
delete [] o;
delete []gmu;
delete []mrk;
return;
}

///////////////////////////////////////////////

//----------------------------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------



