#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <time.h>

FILE *fp1;
FILE *fp2;
FILE *fp3;

#define MAX 130

int main(void)
{
clock_t begin=clock();
int i,j,k,n, ncorr, capn, ni=capn+1, nj=capn+1, maxsteps=500, maxcorrits=100000, work_total;
double x[MAX][MAX], y[MAX][MAX], u[MAX][MAX], ustar[MAX][MAX], unew[MAX][MAX], v[MAX][MAX], vstar[MAX][MAX], vnew[MAX][MAX], p[MAX][MAX], pprime[MAX][MAX], T[MAX][MAX], Tnew[MAX][MAX], Hnx[MAX][MAX], Holdx[MAX][MAX], Hny[MAX][MAX], Holdy[MAX][MAX], massbal[MAX][MAX], dx, dy ;
double length=1, height=1, mu, kcond=.0254, rho=1,cfl=1e6, dt,time=0, cp=1004, uinf=1, Tinf=100, Twall=0, dudxe, dudxw, dudyn, dudys, dvdxe, dvdxw, dvdyn, dvdys, up, vp, ue, uw, un, us, ve, vw, vn, vs, pp, pe, pw, pn, ps, dpdx, dpdy,dpdyn, dpdys, dpdxe, dpdxw, deltap, deltapold, deltau, deltav;
double convcheck=1e-8,corrcheck=1e-10, urf=1.9,sumcontinuity, resid, l2residnormu, l2residnormv, l2residnormp, l2pscale, l2uscale, l2vscale, maxu, maxv, maxp, poffset;
double Re;
double press, velx, vely, hx, hy, mb;
fp3 = fopen("./convergence.dat","w");



//Calculate mesh parameters
Re=1000;
capn=128;
mu=rho*uinf*length/Re;
dx=length/capn;
dy=height/capn;
double FoMax = 0.4;
double CoMax = 0.4;
double dtMaxFo = FoMax*dx*dx*rho/mu;
double dtMaxCo = CoMax*dx/uinf;
if (dtMaxFo < dtMaxCo) dt = dtMaxFo;
else dt = dtMaxCo;
double Fo = mu*dt/(dx*dx*rho);
double Co = uinf*dt/dx;



printf("dx = %lf\tdy %lf\tdt = %lf\tCo = %lf\tFo = %lf\n",dx, dy, dt,Co,Fo);

//Initialize solution
for(i=1;i<=capn;i++){
		for(j=1;j<=capn;j++){
	u[i][j]=0;
	v[i][j]=0;
	p[i][j]=0;
	T[i][j]=0;
}
//printf("\n");
}

for(n=0;n<=maxsteps;n++)
{
//Start Main Loop

//Predict X-Momentum
for(i=1;i<=capn-1;i++){
	for(j=1;j<=capn;j++){
up=u[i][j];
ue=.5*(up+u[i+1][j]);
	if (i==capn-1) ue=0.5*up;
uw=.5*(up+u[i-1][j]);
	if(i==1) uw=0.5*up;
un=.5*(up+u[i][j+1]);
	if(j==capn) un=uinf;
us=.5*(up+u[i][j-1]);
	if(j==1)us=0;
vn=.5*(v[i][j]+v[i+1][j]);
	if(j==capn) vn=0;
vs=.5*(v[i][j-1]+v[i+1][j-1]);
	if(j==1) vs=0;
dudxe=(u[i+1][j]-up)/dx;
	if(i==capn-1) dudxe=(-1.0*up)/dx;
dudxw=(up-u[i-1][j])/dx;
	if(i==1) dudxw = up/dx;
dudyn=(u[i][j+1]-up)/dy;
	if(j==capn) dudyn=(uinf-up)/(0.5*dy);
	//if(j==capn) dudyn=(8*uinf-9*up-u[i][j-1])/(3*dy);
dudys=(up-u[i][j-1])/dy;
	if(j==1) dudys=up/(0.5*dy);
	//if(j==1) dudys=(9*up-u[i][j+1])/(3*dy);
dpdx=(p[i+1][j]-p[i][j])/dx;
Holdx[i][j]=Hnx[i][j];
Hnx[i][j]=mu/(rho*dx)*(dudxe-dudxw)+mu/(rho*dy)*(dudyn-dudys)-(ue*ue-uw*uw)/dx-(vn*un-vs*us)/dy;
ustar[i][j]=u[i][j]+dt*(1.5*Hnx[i][j]-0.5*Holdx[i][j])-dt*dpdx/rho;
//printf("%lf\t",Hnx[i][j]);
} /*printf("\n");*/}

//Predict Y-Momenturm
for(i=1;i<=capn;i++){
	for(j=1;j<=capn-1;j++){
vp=v[i][j];
vn=.5*(vp+v[i][j+1]);
	if(j==capn-1) vn = 0.5*vp;
vs=.5*(vp+v[i][j-1]);
	if(j==1) vs = 0.5*vp;
ve=.5*(vp+v[i+1][j]);
	if(i==capn) ve = 0;
vw=.5*(vp+v[i-1][j]);
	if(i==1) vw = 0;
ue=.5*(u[i][j]+u[i][j+1]);
	if(i==capn) ue = 0;
uw=.5*(u[i-1][j]+u[i-1][j+1]);
	if(i==1) uw=0;
dvdyn=(v[i][j+1]-vp)/dy;
	if(j==capn-1) dvdyn=-0.5*vp/dy;
dvdys=(vp-v[i][j-1])/dy;
	if(j==1) dvdys=0.5*vp/dy;
dvdxe=(v[i+1][j]-vp)/dx;
	if(i==capn) dvdxe=-vp/(0.5*dx);
	//if(i==capn) dvdxe==(-9.0*vp-v[i-1][j])/(3*dx);
dvdxw=(vp-v[i-1][j])/dx;
	if(i==1) dvdxw = vp/(0.5*dx);
	//if(i==1) dudys=(9*vp-v[i+1][j])/(3*dy);
dpdy=(p[i][j+1]-p[i][j])/dy;
Holdy[i][j]=Hny[i][j];
Hny[i][j]=mu/(rho*dy)*(dvdyn-dvdys)+mu/(rho*dx)*(dvdxe-dvdxw)-(ue*ve-uw*vw)/dx-(vn*vn-vs*vs)/dy;
vstar[i][j]=v[i][j]+dt*(1.5*Hny[i][j]-0.5*Holdy[i][j])-dt*dpdy/rho;
}}


//Pressure Correction
for(i=1;i<=capn;i++){
	for(j=1;j<=capn;j++){
work_total++;
//pprime[i][j]=0;
}}

for(ncorr=0;ncorr<=maxcorrits;ncorr++){
l2residnormp=0;
	for(i=1;i<=capn;i++){
		for(j=1;j<=capn;j++){
pp=pprime[i][j];
pn=pprime[i][j+1];
if(j==capn)pn=pp;
ps=pprime[i][j-1];
if(j==1) ps=pp;
pe=pprime[i+1][j];
if(i==capn) pe=pp;
pw=pprime[i-1][j];
if(i==1) pw=pp;
ue=ustar[i][j];
if(i==capn) ue=0;
uw=ustar[i-1][j];
if(i==1) uw=0;
vn=vstar[i][j];
if(j==capn) vn=0;
vs=vstar[i][j-1];
if(j==1) vs=0;
deltap=.25*((pn+ps+pw+pe))-((.25/dt)*(ue*dy-uw*dy+vn*dx-vs*dx));
l2residnormp+=dx*dy*pow(deltap-pprime[i][j],2);
pprime[i][j]=pprime[i][j]*(1.0-urf)+urf*deltap;
}}

//Check Continuity
sumcontinuity=0;
for(i=1;i<=capn;i++){
	for(j=1;j<=capn;j++){
if(i!=capn){
dpdx=(pprime[i+1][j]-pprime[i][j])/dx;
unew[i][j]=ustar[i][j]-dt*dpdx;
}
if(j!=capn){
dpdy=(pprime[i][j+1]-pprime[i][j])/dy;
vnew[i][j]=vstar[i][j]-dt*dpdy;
}
ue=unew[i][j];
if(i==capn) ue=0;
uw=unew[i-1][j];
if(i==1) uw=0;
vn=vnew[i][j];
if(j==capn) vn=0;
vs=vnew[i][j-1];
if(j==1) vs=0;
sumcontinuity+=fabs(ue-uw+vn-vs);
}}

l2residnormp=sqrt(l2residnormp);
if(l2residnormp==0) printf("zero pressure norm\n");
if(ncorr==0 && n==0)l2pscale=l2residnormp;
l2residnormp=l2residnormp/l2pscale;


//printf("\t\t%g\t%g\n",l2residnormp, sumcontinuity);
if((l2residnormp)<corrcheck || ncorr>=maxcorrits) {
l2residnormu=0;
for(i=1;i<=capn;i++){
	for(j=1;j<=capn;j++){
if(i!=capn){
l2residnormu+=fabs(pow(unew[i][j]-u[i][j],2));
u[i][j]=unew[i][j];
}
if(j!=capn){
l2residnormv+=fabs(pow(vnew[i][j]-v[i][j],2));
v[i][j]=vnew[i][j];
}
}}
	break;
}
}
l2residnormu=sqrt(l2residnormu);
if(n==0) l2uscale=l2residnormu;
l2residnormu=l2residnormu/l2uscale;

l2residnormv=sqrt(l2residnormv);
if(n<=1) l2vscale=l2residnormv;
l2residnormv=l2residnormv/l2vscale;


//Update Pressure
for(i=1;i<=capn;i++){
	for(j=1;j<=capn;j++){
p[i][j]+=pprime[i][j];
massbal[i][j]=rho*(u[i][j]-u[i-1][j]+v[i][j]-v[i][j-1])/(1e-6+0.25*(u[i][j]+u[i-1][j]+v[i][j]+v[i][j-1]));
}}

maxp=0;
poffset=p[1][1];
//Relevel Pressure
for(i=1;i<=capn;i++){
	for(j=1;j<=capn;j++){
p[i][j]=p[i][j]-poffset;
if(fabs(p[i][j])>maxp)maxp=p[i][j];
}}

//Verify Continuity
sumcontinuity=0;
for(i=1;i<=capn;i++)
	for(j=1;j<=capn;j++){{
sumcontinuity+=rho*fabs(((u[i][j]-u[i-1][j])*dx+(v[i][j]-v[i][j-1])*dy));
}}

if(n%100==0 || n==0)printf("n\ttime\tu-res \tv-res \tp-resid \tp-steps \tucl \tvcl \tmaxp \tcontinuity-check\n");
if(n%10==0 || n<10)printf("%d\t%0.3lf\t%0.3g\t%0.3g\t%0.5g\t%d\t%0.5g\t%0.5g\t%0.5g\t%0.5g\n",n,time, l2residnormu,l2residnormv,l2residnormp,ncorr, u[capn/2][capn/2], v[capn/2][capn/2], maxp, sumcontinuity);
if(n%100==0) fprintf(fp3,"%g\t%g\n",time,l2residnormu);
time+=dt;
if(l2residnormu<convcheck) break;
//End Main Loop
}



//Print Solution
fp2 = fopen("./solution.dat","w");
fprintf(fp2,"TITLE = \"Solution Data\"\n");
fprintf(fp2,"variables=\"x(m)\"\"y(m)\"\"p(N/m^2)\"\"u(m/s)\"\"v(m/s)\"\"T(K)\"\"Hx(-)\"\"Hy(-)\"\"massbal(kg/s)\"\n");
/*
"\"T(K)\
*/
  fprintf(fp2, "zone T=\"n=%d\"\n",n);
  fprintf(fp2, "I= %d J= %d\n",capn+1, capn+1);
  fprintf(fp2, "DATAPACKING=POINT\n");
for(j=0;j<=capn;j++){
	for(i=0;i<=capn;i++){

		velx=0.5*(u[i][j]+u[i][j+1]);
		hx=0.5*(Hnx[i][j]+Hnx[i][j+1]);
		vely=0.5*(v[i][j]+v[i+1][j]);
		hy=0.5*(Hny[i][j]+Hny[i][j+1]);
		press=0.25*(p[i][j]+p[i+1][j]+p[i][j+1]+p[i+1][j+1]);
		if(j==0) {
		velx=0;
		vely=0;
			}
		if(i==0 || i==capn){
		velx=0;
		vely=0;
		}
		if(j==capn){
		velx=uinf;
		vely=0;
		}
		/*
		if(j==capn) press=0.5*(p[i][capn]+0.125*(p[i][capn]-p[i][capn-1])+p[i+1][capn]+0.125*(p[i+1][capn]-p[i+1][capn-1]));
		if(j==0) press=0.5*(p[i][1]+0.125*(p[i][1]-p[i][2])+p[i+1][1]+0.125*(p[i+1][1]-p[i+1][2]));
		if(i==capn) press=0.5*(p[capn][j]+0.125*(p[capn][j]-p[capn-1][j])+p[capn][j+1]+0.125*(p[capn][j+1]-p[capn-1][j+1]));
		if(i==0) press=0.5*(p[1][j]+0.125*(p[1][j]-p[2][j])+p[1][j+1]+0.125*(p[1][j+1]-p[2][j+1]));
		if(i==0 && j==0) press=p[1][1]+0.125*(p[1][1]-p[1][2])+0.125*(p[1][1]-p[2][1]);
		if(i==capn && j==capn) press=p[capn][capn]+0.125*(p[capn][capn]-p[capn][capn-1])+0.125*(p[capn][capn]-p[capn-1][capn]);
		if(i==capn && j==0) press=p[capn][1]+0.125*(p[capn][1]-p[capn][2])+0.125*(p[capn][1]-p[capn-1][1]);
		if(i==0 && j==capn) press=p[1][capn]+0.125*(p[1][capn]-p[1][capn])+0.125*(p[1][capn]-p[2][capn]);
		*/
		if(j==capn) press=0.5*(p[i][capn]+p[i+1][capn]);
		if(j==0) press=0.5*(p[i][1]+p[i+1][1]);
		if(i==capn) press=0.5*(p[capn][j]+p[capn][j+1]);
		if(i==0) press=0.5*(p[1][j]+p[1][j+1]);
		if(i==0 && j==0) press=p[1][1];
		if(i==capn && j==capn) press=p[capn-1][capn];
		if(i==capn && j==0) press=p[capn][1];
		if(i==0 && j==capn) press=p[1][capn];




		fprintf(fp2,"%e %e %e %e %e %e %e %e %e\n", (i)*dx, (j)*dy , press, velx, vely,T[i][j], Hnx[i][j],Hny[i][j], massbal[i][j]);
}}

printf("Check\n");

fp1 = fopen("./midlines.dat","w");
fprintf(fp1,"0.0\t0.0\t0.0\n");
for(j=1;j<=capn;j++)
{
fprintf(fp1,"%g\t%g\t%g\n",dy*(j-0.5),u[capn/2][j],v[j][capn/2]);
}
fprintf(fp1,"1.0\t1.0\t0.0\n");

//fclose(fp1);
fclose(fp2);
fclose(fp1);
fclose(fp3);
clock_t end=clock();
printf("Program Complete Elapsed Time = %lf seconds total sweeps =%d\n", (end-begin)/(1.0*CLOCKS_PER_SEC), work_total);

return 0;

}
