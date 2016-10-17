#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define RL 10 
#define CL 10 
#define N  RL*CL
#define MC 100000 
#define MCDIS 50000 

extern double drand64(void);

int mod(int a,int L)
{
int t;
if (a==-1) 
t=L-1;
else
t=a%L;
return t;
}

void Metropolis(double s[][CL],double T)
{
int i,j,mm;
double dE;
for (mm=0;mm<N;mm++)
{
i=RL*drand64();
j=CL*drand64();

if((i%2!=0 && j%2!=0)||(i%2==0 && j%2==0))
dE=s[i][j]*2*(s[(i+1)%RL][j]+s[mod((i-1),RL)][j]+s[i][(j+1)%CL]);
else
dE=s[i][j]*2*(s[(i+1)%RL][j]+s[mod((i-1),RL)][j]+s[i][mod((j-1),CL)]);

if(dE <= 0 || drand64()<exp(-dE/T) ) s[i][j] = -s[i][j];
}
}

double Energy(double s[][CL])
{
double E=0; 
int i,j;
for (i=0;i<RL;i++){
for (j=0;j<CL;j++){
if((i%2!=0 && j%2!=0)||(i%2==0 && j%2==0))
E+=-s[i][j]*(s[(i+1)%RL][j]+s[mod((i-1),RL)][j]+s[i][(j+1)%CL]);
else
E+=-s[i][j]*(s[(i+1)%RL][j]+s[mod((i-1),RL)][j]+s[i][mod((j-1),CL)]);
}
}
E=E/(double)2;
return E;
}

int main(int argc, char *argv[]) 
{
double e_avr,es_sum,c,e,e_var;
double s[RL][CL]; 
int i,j,k,mc;
double  T;
FILE *file;
char *file_name;
file_name = "lab3_mc.txt";
file = fopen(file_name,"w");


for (T=0.5;T<=5;T+=0.01)
{
e_avr=0;
es_sum=0;
for (i=0;i<RL;i++)  /* initialize lattice */
{
for (j=0;j<CL;j++)
s[i][j]=1;
}

for (mc=0;mc<MC;mc++)  /* Monte Carlo */
{
    Metropolis(s,T);
if(mc>=MCDIS){
e=Energy(s)/(double)N;
e_avr+=e/(double)(MC-MCDIS);
es_sum+=pow(e,2)/(double)(MC-MCDIS);
}
}
c=(double)N*(es_sum-pow(e_avr,2))/pow(T,2); 
fprintf(file,"%f %f %f %f\n",T,e_avr,sqrt(es_sum-pow(e_avr,2)),c);
    }
fclose(file);
return 0;
}








