// odeint.cc
//
//

#include "nr.h"
#include "nrutil.h"
#include <stdio.h>
#include "math.h"

class Ode_Int {
public:
  int ignore, kount;
  double dxsav, minstep;
  void init(int n);
  void tidy(void);
  void go(double x1, double x2, double xstep,
	  double eps, void (*derivs)(double, double[],double[]));
  void go_simple(double x1, double x2, int nstep,
		 void (*derivs)(double, double[],double[]));
  void go_scale(double x1, double x2, double step,
		 void (*derivs)(double, double[],double[]));
  void set_bc(int n, double num);
  double get_x(int i);
  double get_y(int n, int i);
  double get_d(int n, int i);
  double *xp,**yp;

private:
  double **dydxp,*hstr,*ystart;
  int kmax,nok,nbad,nvar;
  void rkck(double y[], double dydx[], int n, double x, double h,
	    double yout[],
	    double yerr[], void (*derivs)(double, double[], double[]));
  void rkqs(double y[], double dydx[], int n, double *x, double htry, 
	    double eps,	double yscal[], double *hdid, double *hnext,
	    void (*derivs)(double, double[], double[]));
  void odeint(double ystart[], int nvar, double x1, double x2, double eps, 
	      double h1,double hmin, int *nok, int *nbad,
	      void (*derivs)(double, double [], double []));
#define float double
  void rk4(float y[], float dydx[], int n, float x, float h, float yout[],
	     void (*derivs)(float, float [], float []));
  void rkdumb(float vstart[], int nvar, float x1, float x2, int nstep,
	void (*derivs)(float, float [], float []));
  void rkscale(float vstart[], int nvar, float x1, float x2, float h1,
	void (*derivs)(float, float [], float []));
#undef float 

};

void Ode_Int::tidy(void)
{
  free_vector(this->ystart,1,this->nvar);
  free_vector(this->hstr,1,this->kmax);
  free_vector(this->xp,1,this->kmax);
  free_matrix(this->yp,1,this->nvar,1,this->kmax);
  free_matrix(this->dydxp,1,this->nvar,1,this->kmax);
}

void Ode_Int::init(int n)
{
  this->kmax=10000;
  this->nvar=n;
  this->ignore=0;
  this->dxsav=0.0;
  this->minstep=0.0;

  this->xp=vector(1,this->kmax);
  this->yp=matrix(1,this->nvar,1,this->kmax);
  //printf("HERE\n");
  this->dydxp=matrix(1,this->nvar,1,this->kmax);
  // printf("HERE\n");
  this->hstr=vector(1,this->kmax);
  this->ystart=vector(1,this->nvar);
}


void Ode_Int::set_bc(int n, double num)
{
  this->ystart[n]=num;
}


double Ode_Int::get_d(int n, int i)
{
  return this->dydxp[n][i];
}


double Ode_Int::get_x(int i)
{
  return this->xp[i];
}

double Ode_Int::get_y(int n, int i)
{
  return this->yp[n][i];
}


void Ode_Int::go(double x1, double x2, double xstep,
		 double eps, void (*derivs)(double, double[],double[]))
{
  if (this->dxsav == 0.0) this->dxsav=xstep;

  odeint(this->ystart,this->nvar,x1,x2,eps,xstep,this->minstep,&this->nok,
	 &this->nbad,derivs);
}


void Ode_Int::go_simple(double x1, double x2, int nstep,
		 void (*derivs)(double, double[],double[]))
{
  rkdumb(this->ystart,this->nvar,x1,x2,nstep,derivs);
  this->kount=nstep+1;
}

void Ode_Int::go_scale(double x1, double x2, double step, 
		       void (*derivs)(double, double[],double[]))
{
  rkscale(this->ystart,this->nvar,x1,x2,step,derivs);
}



#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define MAXSTP 100000
#define TINY 1.0e-30

void Ode_Int::odeint(double ystart[], int nvar, double x1, double x2, 
		     double eps, double h1, double hmin, int *nok, int *nbad,
		     void (*derivs)(double, double [], double []))
  // see NR pp 721.
{
  int nstp,i;
  double xsav,x,hnext,hdid,h;
  double *yscal,*y,*dydx;

  yscal=vector(1,nvar);
  y=vector(1,nvar);
  dydx=vector(1,nvar);
  x=x1;
  h=SIGN(h1,x2-x1);
  *nok = (*nbad) = kount =0;
  for (i=1;i<=nvar;i++) y[i]=ystart[i];

  if (kmax>0) xsav=x-this->dxsav*2.0;
  for (nstp=1;nstp<=MAXSTP;nstp++) {
    //    printf("t=%lg delt=%lg\n", x, h);
    (*derivs)(x,y,dydx);
    for (i=1;i<=nvar;i++)
      yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
      //      yscal[i]=1.0;
      //yscal[i]=fabs(y[i]);
      //yscal[i]=fabs(y[i]);

    if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(this->dxsav)) {
      this->xp[++kount]=x; this->hstr[kount]=h;
      for(i=1;i<=nvar;i++) {
	this->yp[i][kount]=y[i];
	this->dydxp[i][kount]=dydx[i];
      }
      xsav=x;
    }
    if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
    //printf("h=%lg\n", h);
    rkqs(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
    //printf("hdid=%lg, hnext=%lg\n", hdid, hnext);
    if (hdid == h) ++(*nok); else ++(*nbad);
    if ((x-x2)*(x2-x1) >= 0.0) {
      for (i=1;i<=nvar;i++) ystart[i]=y[i];
      if (kmax) {
	this->xp[++kount]=x;
	for (i=1;i<=nvar;i++) {
	  this->yp[i][kount]=y[i];
	  this->dydxp[i][kount]=dydx[i];
	}      
      }
      free_vector(dydx,1,nvar);
      free_vector(y,1,nvar);
      free_vector(yscal,1,nvar);
      return;
    }
    if (fabs(hnext) <= hmin) printf("Step size too small in odeint..\n");
    h=hnext;
  }
  printf("Too many steps in odeint...\n");
}


void Ode_Int::rkqs(double y[], double dydx[], int n, double *x, 
		   double htry, double eps, double yscal[], double *hdid, 
		   double *hnext,
		   void (*derivs)(double, double[], double[]))

  // Numerical recipes pp 719 

{
   int i;
   double errmax,h,htemp,xnew,*yerr,*ytemp;

   yerr=vector(1,n);
   ytemp=vector(1,n);
   h=htry;
   for (;;) {
      rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
      errmax =0.0;
      for (i=1;i<=(n-this->ignore);i++) 
	errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
      /*    for (i=1;i<=n;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));*/
      errmax /= eps;
      if (errmax <= 1.0) break;
      htemp=SAFETY*h*pow(errmax,PSHRNK);
      h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
      xnew=(*x)+h;
      if (xnew == *x) printf("stepsize underflow in rkqs!!!");
   }
   if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
   else *hnext=5.0*h;
   *x += (*hdid=h);
   for (i=1;i<=n;i++) y[i]=ytemp[i];
   free_vector(ytemp,1,n);
   free_vector(yerr,1,n);
}


void Ode_Int::rkck(double y[], double dydx[], int n, double x, 
		   double h, double yout[], double yerr[], 
		   void (*derivs)(double, double[], double[]))
{
   int i;
   static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
      b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42=-0.9,b43=1.2,
      b51=-11.0/54.0,b52=2.5,b53=-70.0/27.0,b54=35.0/27.0,
      b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
      b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
      c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
      dc5=-277.0/14336.0;
   double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
      dc4=c4-13525.0/55296.0,dc6=c6-0.25;
   double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

   ak2=vector(1,n);
   ak3=vector(1,n);
   ak4=vector(1,n);
   ak5=vector(1,n);
   ak6=vector(1,n);
   ytemp=vector(1,n);
   for (i=1;i<=n;i++)
      ytemp[i]=y[i]+b21*h*dydx[i];
   (*derivs)(x+a2*h,ytemp,ak2);
   for (i=1;i<=n;i++)
      ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
   (*derivs)(x+a3*h,ytemp,ak3);
   for (i=1;i<=n;i++)
      ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
   (*derivs)(x+a4*h,ytemp,ak4);
   for (i=1;i<=n;i++)
     ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
   (*derivs)(x+a5*h,ytemp,ak5);
   for (i=1;i<=n;i++)
     ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
   (*derivs)(x+a6*h,ytemp,ak6);
   for (i=1;i<=n;i++)
      yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
   for (i=1;i<=n;i++)
      yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
   free_vector(ytemp,1,n);
   free_vector(ak6,1,n);
   free_vector(ak5,1,n);
   free_vector(ak4,1,n);
   free_vector(ak3,1,n);
   free_vector(ak2,1,n);
}



#define float double

#define NRANSI

void Ode_Int::rkscale(float vstart[], int nvar, float x1, float x2, float h1,
	void (*derivs)(float, float [], float []))
{
	int i,k;
	float x,h;
	float *v,*vout,*dv;
	float hsum;

	v=vector(1,nvar);
	vout=vector(1,nvar);
	dv=vector(1,nvar);
	for (i=1;i<=nvar;i++) {
		v[i]=vstart[i];
		this->yp[i][1]=v[i];
	}
	this->xp[1]=x1;
	x=x1;
	h=h1; // user supplies initial step size
	k=1; // count the number of steps
	while ( x2>(x+h) ) {
		(*derivs)(x,v,dv);

		// choose a step size
		hsum=0.0;
		for (i=1;i<=nvar;i++) {
		  hsum+=fabs(dv[i]/v[i]);
		}
		h=0.1/hsum;
		//printf("x=%lg, h=%lg\n  ", x, h);

		rk4(v,dv,nvar,x,h,vout,derivs);
		if ((float)(x+h) == x) nrerror("Step size too small in routine rkdumb");
		x += h;
		this->xp[k+1]=x;
		for (i=1;i<=nvar;i++) {
			v[i]=vout[i];
			this->yp[i][k+1]=v[i];
			this->dydxp[i][k+1]=dv[i];
		}

		// choose a new step size
		//hsum=0.0;
		//for (i=1;i<=nvar;i++) {
		//  hsum+=fabs(dv[i]/v[i]);
		//	}
		//h=0.1/hsum;
		//printf("x=%lg, h=%lg\n  ", x, h);
		k++;
	}
	this->kount=k;
	free_vector(dv,1,nvar);
	free_vector(vout,1,nvar);
	free_vector(v,1,nvar);
}

void Ode_Int::rkdumb(float vstart[], int nvar, float x1, float x2, int nstep,
	void (*derivs)(float, float [], float []))
{
	int i,k;
	float x,h;
	float *v,*vout,*dv;

	v=vector(1,nvar);
	vout=vector(1,nvar);
	dv=vector(1,nvar);
	for (i=1;i<=nvar;i++) {
		v[i]=vstart[i];
		this->yp[i][1]=v[i];
	}
	this->xp[1]=x1;
	x=x1;
	h=(x2-x1)/nstep;
	for (k=1;k<=nstep;k++) {
		(*derivs)(x,v,dv);
		rk4(v,dv,nvar,x,h,vout,derivs);
		if ((float)(x+h) == x) nrerror("Step size too small in routine rkdumb");
		x += h;
		this->xp[k+1]=x;
		for (i=1;i<=nvar;i++) {
			v[i]=vout[i];
			this->yp[i][k+1]=v[i];
		}
	}
	free_vector(dv,1,nvar);
	free_vector(vout,1,nvar);
	free_vector(v,1,nvar);
}


void Ode_Int::rk4(float y[], float dydx[], int n, float x, float h, float yout[],
	void (*derivs)(float, float [], float []))
{
	int i;
	float xh,hh,h6,*dym,*dyt,*yt;

	dym=vector(1,n);
	dyt=vector(1,n);
	yt=vector(1,n);
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
	(*derivs)(xh,yt,dyt);
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym);
	for (i=1;i<=n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	free_vector(yt,1,n);
	free_vector(dyt,1,n);
	free_vector(dym,1,n);
}
#undef NRANSI

#undef float
