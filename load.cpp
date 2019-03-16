#define DB 0
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <math.h>
//#include "fft.h"
//
#include "waves.h"

using namespace std;

//------------------------------------------------------------
int Grid::load(char *fname){
	FILE *fp;
	char cbff[128];

	fp=fopen(fname,"r");
	int ndat=0;
	while(fgets(cbff,128,fp) != NULL) ndat++;
	fclose(fp);

	printf("ndat=%d\n",ndat);

	Xcod=(double *)malloc(sizeof(double)*ndat);
	Ycod=(double *)malloc(sizeof(double)*ndat);
	vals=(int *)malloc(sizeof(int)*ndat);
	fp=fopen(fname,"r");
	int i;
	for(i=0;i<ndat;i++){
		fscanf(fp,"%lf,%lf,%d\n",Xcod+i,Ycod+i,vals+i);
	}
	
	double y0=Ycod[0];
	Nx=1;
	while(Ycod[Nx]==y0){
		Nx++;
	}
	Ny=ndat/Nx;
	dx=Xcod[1]-Xcod[0];
	dy=Ycod[Ny]-Ycod[0];
	printf("Nx=%d, Ny=%d\n",Nx,Ny);
	printf("dx=%lf, dy=%lf\n",dx,dy);
	
	fclose(fp);
	return(ndat);
};
//------------------------------------------------------------
Array2D::Array2D(int nx,int ny){
	Nx=nx;
	Ny=ny;
	ndat=Nx*Ny;

	int i,j;
	A2=(double *)malloc(sizeof(double)*ndat);
	for(i=0;i<ndat;i++) A2[i]=0.0;

	A=(double **)malloc(sizeof(double*)*Nx);
	for( i=0;i<Nx;i++) A[i]=A2+Ny*i;

	for(i=0;i<2;i++){
		dx[i]=1.0;
		Xa[i]=0.0;
	}

};
Array2D::Array2D(){
	Nx=1;
	Ny=1;
	A2=(double *)malloc(sizeof(double)*ndat);
	A=(double **)malloc(sizeof(double*)*Nx);
	for(int i=0;i<2;i++){
		dx[i]=1.0;
		Xa[i]=0.0;
	}
};
/*
Array2D::~Array2D(){
	free(A);
	free(A2);
};
*/
void Array2D::set_Xa(double x, double y){
	Xa[0]=x; Xa[1]=y;
};
void Array2D::set_dx(double x, double y){
	dx[0]=x; dx[1]=y;
};
void Array2D::set_Wd(){

	Xb[0]=Xa[0]+(Nx-1)*dx[0];
	Xb[1]=Xa[1]+(Ny-1)*dx[1];
	Wd[0]=Xb[0]-Xa[0];
	Wd[1]=Xb[1]-Xa[1];
};

void Array2D::out(char *fn){
	FILE *fp=fopen(fn,"w");

	int j;
	fprintf(fp,"%d,%d\n",Nx,Ny);
	for( int i=0;i<Nx;i++){
	for( j=0;j<Ny;j++){
		fprintf(fp,"%le\n",A[i][j]);
	}
	}
	fclose(fp);
};
//------------------------------------------------------------
void Array3D::load(char *dir_name){
	FILE *fp;
	char fname[128];

	int i,j,k;

	sprintf(fname,"%s/scope_%d.csv",dir_name,0);
	awv.load(fname);

	k=0;
	for(j=0;j<Ny;j++){
	printf("j=%d\n",j);
	for(i=0;i<Nx;i++){
		sprintf(fname,"%s/scope_%d.csv",dir_name,k);
		awv.amp=A[i][j];
		awv.load(fname);
		k++;
	}
	}

	//for(j=0;j<awv.Nt;j++) printf("%le\n",A[0][0][j]);
};
Array3D::Array3D(int nx, int ny, int nz){
	Nx=nx;
	Ny=ny;
	Nz=nz;

	Nd[0]=Nx;
	Nd[1]=Ny;
	Nd[2]=Nz;
	
	ndat=Nx*Ny*Nz;
	int i,j,k;  

	A3=(double *)malloc(sizeof(double)*ndat);
	for(k=0;k<ndat;k++) A3[k]=0.0;

	A2=(double **)malloc(sizeof(double *)*nx*ny);	
	for(j=0;j<nx*ny; j++) A2[j]=A3+j*nz;

	A=(double ***)malloc(sizeof(double **)*nx);
	for(i=0;i<nx;i++) A[i]=A2+i*ny;

	for(i=0;i<3;i++){
		dx[i]=1.0;
		Xa[i]=0.0;
	}
}
Array3D::~Array3D(){
	free(A3);
	free(A2);
	free(A);
};
void Array3D::set_Xa(double x, double y, double z){
	Xa[0]=x; Xa[1]=y; Xa[2]=z;
};
void Array3D::set_dx(double x, double y, double z){
	dx[0]=x; dx[1]=y; dx[2]=z;
};
void Array3D::set_Wd(){

	//int N[3]={Nx,Ny,Nz};
	//printf("Xa=%lf %lf %lf\n",Xa[0],Xa[1],Xa[2]);
	//printf("dx=%lf %lf %lf\n",dx[0],dx[1],dx[2]);
	for(int i=0;i<3;i++){
		Xb[i]=Xa[i]+(Nd[i]-1)*dx[i];
		Wd[i]=Xb[i]-Xa[i];
	}
};
void Array3D::print_dim(){
	printf("Array size=(%d, %d, %d)\n",Nx,Ny,Nz);
};
void Array3D::set_val(double v){
	for(int i=0;i<ndat;i++) A3[i]=v;
};
void Array3D::Lp(int p){
	int i,j;
	Array2D Lp(Nx,Ny);

	Wv1D wv;
	wv=awv;
	wv.print_info();

	double t1=wv.t1;
	double t2=wv.t2;
	char fn[128];
	if(p==2){
		sprintf(fn,"L2.out");
		for(i=0;i<Nx;i++){
			printf("i=%d/%d\n",i+1,Nx);
		for(j=0;j<Ny;j++){
			wv.amp=A[i][j];
			Lp.A[i][j]=wv.L2(t1,t2);
		}
		}
	}else if(p==0){
		sprintf(fn,"Linf.out");
		for(i=0;i<Nx;i++){
			printf("i=%d/%d\n",i+1,Nx);
		for(j=0;j<Ny;j++){
			wv.amp=A[i][j];
			Lp.A[i][j]=wv.max(t1,t2);
		}
		}
	};
	Lp.out(fn);
};
Array2D Array3D::proj(){
	Array2D Bdat(Ny,Nz);
	int i,j,k;
	//printf("Nx,Ny,Nz=%d %d %d\n",Nx,Ny,Nz);
	for(i=0; i<Ny; i++){
		//printf("i=%d\n",i);
	for(j=0; j<Nz; j++){
		for(k=0; k<Nx; k++) Bdat.A[i][j]+=A[k][i][j]; 
		Bdat.A[i][j]/=Nx;
	}
	}
	return(Bdat);
};
void Array3D::CorrY(){

	int i,j;
	double tmax,Amax;
	Array2D Tm(Nx,Ny-1);
	Array2D Amp(Nx,Ny-1);
	char fn[128]="debug.dat";

	Wv1D wv1,wv2,cor;
	wv1=awv;
	wv2=awv;
	wv1.print_info();
	wv2.print_info();
	int tmp;
	for(i=0;i<Nx;i++){
		printf("i=%d/%d\n",i+1,Nx);
	for(j=0;j<Ny-1;j++){
		wv1.amp=A[i][j];
		wv2.amp=A[i][j+1];
		cor=corr(wv1,wv2,&tmax,&Amax);
		Tm.A[i][j]=tmax;
		Amp.A[i][j]=Amax;
	}
	}

	char fname[128];
	strcpy(fname,"tmax.out");
	Tm.out(fname);
	strcpy(fname,"amax.out");
	Amp.out(fname);

};
void Array3D::CorrX(){

	int i,j;
	double tmax,Amax;
	Array2D Tm(Nx-2,Ny);
	Array2D Amp(Nx-2,Ny);

	Wv1D wv1,wv2,cor;
	wv1=awv;
	wv2=awv;
	wv1.print_info();
	wv2.print_info();
	int tmp;
	for(i=0;i<Nx-2;i++){
		printf("i=%d/%d\n",i+1,Nx-1);
	for(j=0;j<Ny;j++){
		wv1.amp=A[i][j];
		wv2.amp=A[i+2][j];
		cor=corr(wv1,wv2,&tmax,&Amax);
		Tm.A[i][j]=tmax;
		Amp.A[i][j]=Amax;
	}
	}

	char fname[128];
	strcpy(fname,"tmax2.out");
	Tm.out(fname);
	strcpy(fname,"amax2.out");
	Amp.out(fname);
};
//------------------------------------------------------------
Wv1D::Wv1D(){
	amp=0;
	time=0;
	mllc=false;
	fft_stat=0;
};
Wv1D::Wv1D(int ndat){
	Nt=ndat;
	amp=(double *)malloc(sizeof(double)*Nt);
	time=(double *)malloc(sizeof(double)*Nt);
	mllc=true;
	fft_stat=0;
	strcpy(data_file,"");
	t1=0.0;
	t2=0.0;
	dt=0.0;
};
Wv1D::Wv1D(char *fname){
	amp=0;
	time=0;
	mllc=false;
	fft_stat=0;
	Wv1D::load(fname);
};
void Wv1D::print_info(){
	printf("File: %s\n",data_file);
	printf("(t1,t2)=%lf,%lf\n",t1,t2);
	printf("dt=%lf\n",dt);
	printf("Nt=%d\n",Nt);
}
int Wv1D::load(char *fname){
	FILE *fp=fopen(fname,"r");

	strcpy(data_file,fname);
	fscanf(fp,"%lf, %lf\n",&t1,&dt);
	t1*=1.e06;
	dt*=1.e06;
	fscanf(fp,"%d\n",&Nt);
	if(!mllc){
		amp=(double *)malloc(sizeof(double)*Nt);
		time=(double *)malloc(sizeof(double)*Nt);
		mllc=true;
	}
	double sum=0.0;
	for(int i=0;i<Nt;i++){
		fscanf(fp,"%lf\n",amp+i);
		time[i]=t1+dt*i;
		sum+=amp[i];
	};
	t2=time[Nt-1];
	sum/=Nt;
	for(int i=0;i<Nt;i++) amp[i]-=sum;
	fclose(fp);
}
int Wv1D::FFT(int isgn){
	int p=ceil(log2(Nt));
	Np=pow(2,p);

	if(fft_stat==0) Amp=(complex<double> *)malloc(sizeof(complex<double>)*Np);
	if(isgn==1){
		for(int i=0;i<Nt;i++) Amp[i]=complex<double>(amp[i],0.0);
		for(int i=Nt;i<Np;i++) Amp[i]=complex<double>(0.0,0.0);
	}

	fft(Amp,Np,isgn);
	fft_stat=isgn;

	return(Np);
};
double Wv1D::L2(double t1, double t2){
	int i1,i2;
	i1=int(t1/dt);
	i2=int(t2/dt);
	if(i1<0) i1=0;
	if(i2<0) i2=0;
	if(i1>=Nt-1) i1=Nt-1;
	if(i2>=Nt-1) i2=Nt-1;

	int i=0;
	double sum=0.0;
	for(i=i1;i<=i2;i++) sum+=(amp[i]*amp[i]);
	return(sqrt(sum)/(i2-i1+1));
};
double Wv1D::max(double t1, double t2){
	int i1,i2;
	i1=int(t1/dt);
	i2=int(t2/dt);
	if(i1<0) i1=0;
	if(i2<0) i2=0;
	if(i1>=Nt-1) i1=Nt-1;
	if(i2>=Nt-1) i2=Nt-1;

	int i=0;
	double amax=fabs(amp[i2]);
	for(i=i1;i<i2;i++){
	       	if(amax < fabs(amp[i])) amax=fabs(amp[i]);
	}
	return(amax);
};
void Wv1D::out_Amp(char *fn,int ofst){
	FILE *fp=fopen(fn,"w");
	int j;
	double xx,x1,dx;
	double df;
	if(fft_stat==1){
		df=1./dt/Np;
		dx=df;
		x1=0.0;
	}else{
		x1=t1;
		dx=dt;
	}
	for(int i=0;i<Np;i++){
		j=(ofst+i)%Np;
		xx=(i-ofst)*dx;
		fprintf(fp,"%le %le %le %le\n",xx,Amp[j].real(),Amp[j].imag(),abs(Amp[j]));
	}
	fclose(fp);
};
void Wv1D::out_amp(char *fn){
	FILE *fp=fopen(fn,"w");
	for(int i=0;i<Nt;i++) fprintf(fp,"%le\n",amp[i]);

	fclose(fp);
};
//------------------------------------------------------------
Wv1D corr(Wv1D wv1, Wv1D wv2, double *tmax, double *Amax){
	wv1.FFT(1);
	wv2.FFT(1);

	Wv1D wv3(wv1.Np);
	wv3.FFT(1);
	double A1=0.0,A2=0.0;
	complex <double > Z1,Z2;
	Z1=complex<double>(0.0,0.0);
	Z2=complex<double>(0.0,0.0);
	wv3=wv1;
	for(int i=0;i<wv1.Np;i++){
		Z1=wv1.Amp[i];
		Z2=conj(wv2.Amp[i]);
		wv3.Amp[i]=Z1*Z2;
		A1+=abs(Z1*Z1);
		A2+=abs(Z2*Z2);
	};
	wv3.FFT(-1);
	A1=sqrt(A1*A2);
	int imax=0;
	double A0=abs(wv3.Amp[0].real()/A1);
	for(int i=0;i<wv1.Np;i++){
		wv3.Amp[i]/=A1;
		if(A0 < abs(wv3.Amp[i].real())){
			A0=abs(wv3.Amp[i].real());
			imax=i;
		}
	}
	double t0=imax*wv1.dt;
	if(imax>wv1.Np/2) t0=-(wv1.Np-imax)*wv1.dt;
	(*tmax)=t0;
	(*Amax)=A0;
	//printf("%d %le %le\n",imax,t0,A0);
	return(wv3);
};
//------------------------------------------------------------
#if DB ==1
int main(){

	Grid Gd;	// Measurement Grid 
	char fname[128]="../W20H30_fine/xyl.csv";
	Gd.load(fname); // import grid info.


	int i;
/*
	for(i=5;i+41<=2500;i+=41){
	sprintf(fname,"../W20H30_fine/scope_%d.csv",i);
	Wv1D awv1(fname); // import wave data 1

	sprintf(fname,"../W20H30_fine/scope_%d.csv",i+41);
	Wv1D awv2(fname); // import wave data 2
	Wv1D cor;	// correlation function

	char fn[128];
	cor=corr(awv1,awv2);
	strcpy(fn,"cor.out");
	cor.out_Amp(fn,cor.Np/2); 
	//cor.print_info();
	}
*/


	sprintf(fname,"../W20H30_fine/scope_%d.csv",0);
	Wv1D awv1(fname); // import wave data 1
	Array3D WV(Gd.Nx,Gd.Ny,awv1.Nt);
	WV.print_dim();
	char dir_name[128]="../W20H30_fine";
	WV.load(dir_name);
	WV.awv.print_info();

	//WV.CorrY();
	//WV.CorrX();
	WV.Lp(2);

	return(0);
};
#endif
int main(){

	Grid Gd;	// Measurement Grid 
	char fname[128]="../W20H30_fine/xyl.csv";
	Gd.load(fname); // import grid info.

	//---------------------------------------
	sprintf(fname,"../W20H30_fine/scope_%d.csv",0);
	Wv1D awv1(fname); // import wave data 1
	Array3D WV(Gd.Nx,Gd.Ny,awv1.Nt);

	WV.set_Xa(Gd.Xcod[0],Gd.Ycod[0],awv1.t1);
	WV.set_dx(Gd.dx,Gd.dy,awv1.dt);
	WV.set_Wd();

	WV.print_dim();
	char dir_name[128]="../W20H30_fine";
	WV.load(dir_name);

	//---------------------------------------
	Array2D Bwv;
	Bwv=WV.proj();
	puts("done");
	char fn[128]="bwv.out";
	Bwv.out(fn);
	
	return(0);
};
