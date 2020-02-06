#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "odeint.h"

//#define LOWER 7.2589e+14 //SDSS u min
//#define UPPER 9.97645e+14 // SDSS u max
//#define LOWER 5.16438e+14 //SDSS g min
//#define UPPER 8.20226e+14 // SDSS g max
//#define LOWER 4.14651e+14 //SDSS r min
//#define UPPER 5.54658e+14 // SDSS r max
//#define LOWER 3.47384e+14 //SDSS i min
//#define UPPER 4.64434e+14// SDSS i max
//#define LOWER 2.68151e+14 //SDSS z min
//#define UPPER 3.8658e+14 // SDSS z max
//#define LOWER 9.95988e+14 //GALEX NUV min
//#define UPPER 1.77287e+15 // GALEX NUV  max
//#define LOWER 1.65631e+15 //GALEX FUV min
//#define UPPER 2.23559e+15 // GALEX FUV max
//#define LOWER 8.74669e+14 //Swift M2 min
//#define UPPER 1.87079e+15 // Swift M2 max
//#define LOWER 5.18448e+14 //Swift W1 min
//#define UPPER 1.87079e+15 // Swift W1 max
//#define LOWER 5.99286e+14 //Swift W2 min
//#define UPPER 1.87079e+15 // Swift W2 max
//#define LOWER 7.22391e+14 //Bessel U min
//#define UPPER 9.82926e+14 // Bessel U max
#define LOWER 5.35344e+14 //Bessel B min
#define UPPER 8.1025e+14 // Bessel B max
//#define LOWER 4.34482e+14 //Bessel V min
//#define UPPER 6.24568e+14 // Bessel V max
//#define LOWER 3.52697e+14 //Bessel R min
//#define UPPER 5.35344e+14 // Bessel R max
//#define LOWER 3.29442e+14 //Bessel I min
//#define UPPER 4.22243e+14 // Bessel I max

//#define VEGA 0 //AB mags
//#define VEGA 0.79 // Bessel U
#define VEGA -0.09 //Bessel B
//#define VEGA 0.02 // Bessel V
//#define VEGA 0.21 // Bessel R
//#define VEGA 0.45 //Bessel I

#define AREA "/Users/annahaynie/Desktop/Carnegie/CSM/Bessell_B_nu.txt"
#define TEMP "/Users/annahaynie/SNEC/SNEC-1.01/Data/T_eff.dat"
#define RADIUS "/Users/annahaynie/SNEC/SNEC-1.01/Data/rad_photo.dat"

#define h 6.62606885E-027 // units = erg*s
#define c 2.99792458E+010 // units = cm/s
#define d 3.08567758E+019 //10 pc in cm for absolute mag calculation
#define k_b 1.380658E-016 //boltzmann constant units erg/K

#define OUTPUT "/Users/annahaynie/SNEC/SNEC-1.01/Bmag_full.txt"

//  swift_mag2.c

//  Same goals as mag_array.c but using the integration method found in double_integration.c

//  gives same results as mag_array.c but is faster, but still too slow to loop through all time steps on my personal machine- my estimate is that it would take between 5-8 hours to do this so I have been manually changing the time at each step

//  finding magnitudes in Swift bands using correctly calculated F_nu values

//  DOES TAKE TIME DELAY EFFECTS INTO ACCOUNT

//  Created by Anna Haynie on 6/19/19.


Ode_Int ODE;
Ode_Int ODE2;

void readdata();

void derivs(double x, double y[], double dydx[]);
void derivs2(double x, double y[], double dydx[]);

//int get_filelines_area();
//int get_filelines_temp();

int FILELINES_area=19;
int FILELINES_temp=15905;

int get_index_temp(double tvalue);
int get_index_area(double nu);

double get_tr(double tvalue, double theta);
double get_area(double nu);
double get_Bnu(double tr, double nu);

double t_file[50000], temp_file[50000],t2_file[50000], rad_file[50000],lambda_file[50000], nu_file[50000], area_file[50000];

//using a single time until satisfied that file works/is fast enough
double tvalue;
double nu;
double error=1.0E+02;

int main()
{
    readdata();
    //for a given time either in the loop or specified above, do the "outer" integral over nu
    for (int i=5349; i<=FILELINES_temp; i++) {
        tvalue=t_file[i];
        //int j = get_index_temp(tvalue);
        //printf("%lf %i %i \n",tvalue, i, j-1);
        
        ODE.init(2);
        ODE.set_bc(1,0.0);
        ODE.set_bc(2,0.0);
        ODE.go(LOWER,UPPER,1.0E+04,5.0,derivs);
        
        FILE *results;
        results = fopen(OUTPUT,"a");
        
        //use output of 'outer' integrals to calculate mags
        double N,D,r,eta,m, temp;
        r=rad_file[i];
        temp=temp_file[i];
        eta=(r/d)*(r/d);
        N=ODE.get_y(1,ODE.kount);
        D=ODE.get_y(2,ODE.kount);
        m=-2.5*log10(eta*N/D)-48.60-VEGA;
        fprintf(results,"%lg %lg \n",tvalue,m);
        //printf("r=%lf T=%lf N=%lf D=%lf m=%lf \n", r, temp, N, D, m);
        
        ODE.tidy();
        fclose(results);
    }
}

void readdata()
{
    FILE *data, *data2, *data3;
    data = fopen(AREA,"r");
    data2 = fopen(TEMP,"r");
    data3 = fopen(RADIUS,"r");
    for(int i=0;i<FILELINES_area;i++) {
        fscanf(data, "%lf %lf %lf", &lambda_file[i], &nu_file[i] , &area_file[i]);
    }
    for(int i=0;i<FILELINES_temp;i++) {
        fscanf(data2, "%lf %lf", &t_file[i], &temp_file[i]);
        fscanf(data3, "%lf %lf", &t2_file[i], &rad_file[i] );
        //printf("%lf \n", t_file[i]);
    }
    fclose(data);
    fclose(data2);
    fclose(data3);
}

void derivs(double x, double y[], double dydx[])
{
  //set up "outer" integral over nu of F_nu*area and area
    double area;
    nu=x;
    area=get_area(nu);
    
    //do "inner" integral over theta for every given value of nu that the "outer" integral calls
    ODE2.init(1);
    ODE2.set_bc(1,0.0);
    ODE2.go(0.0,PI/2,PI*1.0e-3,5.0e-7,derivs2);
    
    //use output of "inner" integral to calculate F_nu to be used in "outer" integral
    double fnu;
    fnu=ODE2.get_y(1,ODE2.kount);
    
    ODE2.tidy();
    
    dydx[1]=fnu*area;
    dydx[2]=area;
    
    //printf("%lf %lf \n", area, fnu);
}

void derivs2(double x, double y[], double dydx[])
{
    //set up integral of L_nu' over theta
    double theta,tr,Bnu;
    theta=x;
    tr=get_tr(tvalue, theta);
    Bnu=get_Bnu(tr,nu);
    dydx[1]=Bnu*sin(theta)*cos(theta);
    //printf("%lf %lf %lf \n", theta, tr, Bnu);
}

//int get_filelines_area()
//{
    //count number of lines in area file
  //  FILE *fp;
  //  char ch;
  //  int linesCount=0;
  //  fp=fopen(AREA,"r");
  //  while((ch=fgetc(fp))!=EOF)
  //  {
    //    if(ch=='\n')
      //      linesCount++;
   // }
   // fclose(fp);
   // return linesCount;
//}

//int get_filelines_temp()
//{
    //count number of lines in temp file
  //  FILE *fp;
  //  char ch;
  //  int linesCount=0;
  //  fp=fopen(TEMP,"r");
  //  while((ch=fgetc(fp))!=EOF)
  //  {
    //    if(ch=='\n')
      //      linesCount++;
   // }
   // fclose(fp);
   // return linesCount;
//}

int get_index_temp(double tvalue)
{
    //find index of given nu in area file
    int j = 0;
    for(int i=0;i<FILELINES_temp;i++) {
        if(t_file[i]<tvalue){
            j++;
        }
        //printf("%i %i %lf \n", i,j, t_file[i]);
    }
    return j;
}

int get_index_area(double nu)
{
    //find index of given nu in area file
    int j = 1;
    for(int i=0;i<FILELINES_area;i++) {
        if(nu_file[i]>nu)
            j = i;
    }
    return j;
}

double get_tr(double tvalue, double theta)
{
    //taking time delay effects into account
    //calculate retardation time give a theta on the surface of the star to find the true r(tr) and T(tr)
    double tr, r, td, t;
    int i;
    i=get_index_temp(tvalue);
    r=rad_file[i-1];
    t=t_file[i-1];
    td=(r/c)*(1-cos(theta));
    tr=t-td;
    do {
        int j=get_index_temp(tr);
        double r2=rad_file[j-1];
        double td2=(r2/c)*(1-cos(theta));
        tr=t-td2;
        error=(fabs(r-r2))/r;
        r=r2;
    } while (error>1.0E-03);
    
    //printf("%lf %lf %i %lf %lf %lf \n", tvalue,theta,i,t,r,tr);
    
    return tr;
}

double get_area(double nu)
{
    //linear interpolation to find area for a given nu called by the "outer" integral
    double area, a, b;
    int i,j;
    j=get_index_area(nu);
    i=j+1;
    a=(nu_file[j]-nu)/(nu_file[j]-nu_file[i]);
    b=(nu-nu_file[i])/(nu_file[j]-nu_file[i]);
    area=a*area_file[i]+b*area_file[j];
    return area;
}

double get_Bnu(double tr, double nu)
{
    //calculate L_nu' to be integrated*sin(theta) over theta to get the correct F_nu
    int i;
    double temp, c2, nu3, exp_, Bnu;
    i=get_index_temp(tr);
    temp=temp_file[i-1];
    c2=c*c;
    nu3=nu*nu*nu;
    exp_=(h*nu)/(k_b*temp);
    Bnu=(4*PI*h/c2)*(nu3)/(exp(exp_)-1);
    
    //printf("%i %lf %lf %lf %lf \n", i, tr,temp,nu, Bnu);
    return Bnu;
}

#include "swift_mag2.h"
