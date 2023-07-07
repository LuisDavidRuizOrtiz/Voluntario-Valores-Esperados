#include <math.h>
#include <stdio.h>
#include "complex.h" 
#define N 200
#define pi 3.141592


void derivada(fcomplex *input, fcomplex *output)
{
int k,n;
output[0]=RCmul((1.0/0.4),(Csub(input[1],input[0])));
for(k=1;k<N;k++)
output[k]=RCmul((1.0/0.4),(Csub(input[k],input[k-1])));
}

void derivadaR(double *input, double *output)
{
int k,n;
output[0]=((1.0/0.4)*(input[1]-input[0]));
for(k=1;k<N;k++)
output[k]=((1.0/0.4)*(input[k]-input[k-1]));
}


int main()
{
fcomplex phi[N],A0[N],b[N],alpha[N],beta[N],gamma[N],chi[N],phiconj[N],phip[N],phipp[N],expp,expp2,exppaux;
int i,j,k,l,q,nc=(N/10.0);
double h=0.4,lamb=0.0002,k0,s,V[N],Vp[N],x0,Am[N],Ap[N],norma2,norma,normap,expx,expx2,delx,delp,expVp,derivP,derivX,expxaux;

FILE *f1,*f2,*f3,*f4,*f5,*f6,*f7,*f8,*f9,*f10,*f11,*f12,*f13,*f14,*f15,*f16,*f17,*f18;
f1=fopen("modulo.txt","w");
f2=fopen("pot.txt","w");
f3=fopen("inicial.txt","w");
f4=fopen("Dt(expp).txt","w");
f5=fopen("chi.txt","w");
f6=fopen("norma.txt","w");
f7=fopen("VpEsperado.txt","w");
f8=fopen("alpha.txt","w");
f9=fopen("mxdot.txt","w");
f10=fopen("MomentoEsperado.txt","w");
f11=fopen("PosicionEsperada.txt","w");
f12=fopen("Esperadox2.txt","w");
f13=fopen("Esperadop2.txt","w");
f14=fopen("Indetermin.txt","w");
f15=fopen("IndetX.txt","w");
f16=fopen("IndetP.txt","w");
f17=fopen("Ehrenfest1.txt","w");
f18=fopen("Ehrenfest2.txt","w");
x0=50.0;
k0=(2*pi*nc/N);
s=(1.0/(4*k0*k0));
                         //Definir Potencial
for(j=0;j<N;j++)
	{
	if(((4.8*N/5.0)<j)&&((5.0*N/5.0)>j))
	V[j]=1000;
	else 
	V[j]=lamb*pow(((j)-(N/2.0)),2);
	
	}
derivadaR(V,Vp);	
	                //Funcion de onda inicial
for(j=1;j<N-1;j++)
phi[j]=Cgauss(k0*j,exp(-pow(((j*h)-x0),2)/(2*N*N*h*h/(16*16))));
phi[0]=Complex(0,0);
phi[N-1]=Complex(0,0);

for(j=0;j<N;j++)
{
fprintf(f3,"%lf\t%lf\n",j*h,Cabs(phi[j]));
}
                     //Norma de phi inicial
norma2=0.0;
for(i=0;i<N;i++)
{
norma2=norma2+(h*Cabs(phi[i])*Cabs(phi[i]));
} 
norma=sqrt(norma2);
printf("%lf",norma);
                     

                        //Calcular A+,A-,A0
for(j=1;j<N;j++)
	{
	Am[j]=1.0;
	Ap[j]=1.0;
	A0[j]=Complex(-2-V[j],2.0/s);
	b[j]=Cmul(Complex(0,4.0/s),phi[j]);
	}     
	
	               //Calcular alfa
alpha[N-1]=Complex(0,0);
for(j=N-2;j>=0;j--)
	{
	gamma[j+1]=Cdiv(Complex(1.0,0),Cadd(A0[j+1],RCmul(Ap[j+1],alpha[j+1])));
	alpha[j]=RCmul(-Am[j+1],gamma[j+1]);
	}
for(j=0;j<N;j++)
fprintf(f8,"%lf\t%lf\t%lf\n",j*h,alpha[j].r,alpha[j].i);
	
	                      //Set chi inicial =0
for(i=0;i<N;i++)
{
chi[i]=Complex(0,0);
}	
	
for(q=0;q<10000;q++)         //Inicio iteraciones
{                            
                            //Calculo b 
  for(j=0;j<N;j++)
  b[j]=Cmul(Complex(0,4.0/s),phi[j]);
  
  

                             //Calculo de beta en t
beta[N-1]=Complex(0,0);
for(j=N-2;j>=0;j--)
	{
	beta[j]=Cmul(gamma[j+1],Csub(b[j+1],RCmul(Ap[j+1],beta[j+1])));
	}	

chi[0]=Complex(0,0);      //Calculo de chi en t a partir de alpha y beta         
for(j=1;j<N;j++)
{
chi[j]=Cadd(Cmul(alpha[j-1],chi[j-1]),beta[j-1]);
}
chi[N-1]=Complex(0,0);

                        //Calculo de phi en t a partir de phi en t-h y chi
for(j=0;j<N;j++)
{
phi[j]=Csub(chi[j],phi[j]);
}
phi[N-1]=Complex(0,0);
              

                                                           //Calculo de valor esperado de p 
exppaux=expp;
derivada(phi, phip);
expp=Complex(0.0,0.0);
for(j=0;j<N;j++)
{
expp=Cadd(expp,(Cmul(Conjg(phi[j]),(Cmul(Complex(0.0,-h),phip[j])))));
}                                                        

                                                          //Calculo de d<p>/dt
derivP=(4.0*k0*k0/h*h)*(Cabs(expp)-Cabs(exppaux));


                                                           //Calculo del valor esperado de x
expxaux=expx;
expx=0.0;
for(j=0;j<N;j++)
{
expx=expx+(1.0/norma)*(j*h*h*Cabs(phi[j])*Cabs(phi[j]));
}  
                                                             //Calculo de d<x>/dt
                                                             
derivX=(4.0*k0*k0/h*h)*(expx-expxaux);


                                                            //Valor esp de x²
expx2=0.0;
for(j=0;j<N;j++)
{
expx2=expx2+(1.0/norma)*(j*h*h*j*h*Cabs(phi[j])*Cabs(phi[j]));
}                     
                                                      
                                                            //Valor esp de p² 
derivada(phip, phipp);
expp2=Complex(0.0,0.0);
for(j=0;j<N;j++)
{
expp2=Cadd(expp2,(Cmul(Conjg(phi[j]),(RCmul(-h,phipp[j])))));
}                                                      

                                                           //Calculo del valor esperado de V'
expVp=0.0;
for(j=0;j<N;j++)
{
expVp=expVp+(Vp[j]*h*Cabs(phi[j])*Cabs(phi[j]));
}    

 
                                 
                                                         //Indeterminaciones en p y x                 
delx=sqrt(fabs((expx2)-(pow(expx,2))));
delp=sqrt(fabs((Cabs(expp2))-(pow(Cabs(expp),2))));



                         

if(q%10==0) 
{ 
 fprintf(f11,"%i\t%lf\n",q,expx);
 fprintf(f10,"%i\t%lf\n",q,Cabs(expp)); 
 fprintf(f12,"%i\t%lf\n",q,expx2); 
 fprintf(f13,"%i\t%lf\n",q,Cabs(expp2)); 
 fprintf(f14,"%i\t%lf\n",q,delx*delp); 
 fprintf(f15,"%i\t%lf\n",q,delx);
 fprintf(f16,"%i\t%lf\n",q,delp); 
 fprintf(f7,"%i\t%lf\n",q,-expVp);   
 fprintf(f4,"%i\t%lf\n",q,derivP);  
 fprintf(f9,"%i\t%lf\n",q,0.5*derivX); 
 fprintf(f17,"%i\t%lf\n",q,Cabs(expp)-0.5*derivX);
 fprintf(f18,"%i\t%lf\n",q,derivP+expVp);                  
}
                       
                         //plot |phi|² y potencial
if(q%2==0)
{
for(j=0;j<(N-1);j++)
	{
	fprintf(f1,"%i\t%lf\n",j,(Cabs(phi[j])*Cabs(phi[j])));
	fprintf(f2,"%i\t%lf\n",j,V[j]);
	fprintf(f5,"%i\t%lf\t%lf\n",j,chi[j].r,chi[j].i);
	}
fprintf(f1,"%i\t%lf\n\n\n",(N-1),Cabs(phi[N-1])*Cabs(phi[N-1]));
fprintf(f2,"%i\t%lf\n\n\n",(N-1),V[N-1]);
normap=0.0;
for(j=0;j<N;j++)
	normap=normap+(h*Cabs(phi[j])*Cabs(phi[j]));
fprintf(f6,"%i\t%lf\n",q,(normap/norma));
}


}


fclose(f1);
fclose(f2);
fclose(f3);
fclose(f4);
fclose(f5);
fclose(f6);
fclose(f7);
fclose(f8);
fclose(f9);
fclose(f10);
fclose(f11);
fclose(f12);
fclose(f13);
fclose(f14);
fclose(f15);
fclose(f16);
fclose(f17);
fclose(f18);
return 0;
}
