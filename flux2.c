#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "odeint.h"
#define TEMP "T_eff.txt"
#define RADIUS "rad_photo.txt"
#define h 6.62606885E-027 // units = erg*s
#define c 2.99792458E+010 // units = cm/s
#define k_b 1.380658E-016 //boltzmann constant units erg/K
#define sigma_sb 5.6704E-05 //stefan boltzmann constant in erg⋅cm−2⋅s−1⋅K−4
#define OUTPUT1 "flux_fnu_test.txt"
#define OUTPUT2  "flux_integratedflux.txt"


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

int get_filelines_temp();

const int FILELINES_temp=get_filelines_temp();

int get_index_temp(double tvalue);

double get_tr(double tvalue, double theta);
double get_Bnu(double tr, double nu);

double t_file[50000], temp_file[50000],t2_file[50000], rad_file[50000];

double tvalue=2.3E+05;
double nu;
double error=1.0E+02;

int main()
{
    readdata();
    //for a given time either in the loop or specified above, do the "outer" integral over nu
    ODE.init(1);
    ODE.set_bc(1,0.0);
    ODE.go(1.0E+12,6.0E+16,1.0E+03,5.0,derivs);
        
    FILE *results;
    results = fopen(OUTPUT2,"a");
        
    //use output of 'outer' integrals to calculate mags
    double F, temp, Tf;
    int i=get_index_temp(tvalue);
    temp=temp_file[i];
    F=ODE.get_y(1,ODE.kount);
    Tf=pow((F/sigma_sb),0.25);
    //fprintf(results,"%lg %lg %lg %lg \n",tvalue,F, temp, Tf);
    printf("%lg %lg %lg %lg \n",tvalue,F, temp, Tf);
    ODE.tidy();
    fclose(results);

}

void readdata()
{
    FILE *data, *data2;
    data = fopen(TEMP,"r");
    data2 = fopen(RADIUS,"r");
    for(int i=0;i<FILELINES_temp;i++) {
        fscanf(data, "%lf %lf", &t_file[i], &temp_file[i]);
        fscanf(data2, "%lf %lf", &t2_file[i], &rad_file[i] );
    }
    fclose(data);
    fclose(data2);
}

void derivs(double x, double y[], double dydx[])
{
  //set up "outer" integral over nu of F_nu*area and area
    double f_nu;
    nu=x;
    
    //do "inner" integral over theta for every given value of nu that the "outer" integral calls
    ODE2.init(1);
    ODE2.set_bc(1,0.0);
    ODE2.go(0.0,PI/2,PI*1.0e-3,5.0e-7,derivs2);
    
    //use output of "inner" integral to calculate F_nu to be used in "outer" integral
    f_nu=ODE2.get_y(1,ODE2.kount);
    
    FILE *results;
    results = fopen(OUTPUT1,"a");
    fprintf(results,"%lf %lf \n", nu, f_nu);
    fclose(results);
    
    ODE2.tidy();
    
    dydx[1]=f_nu;
}

void derivs2(double x, double y[], double dydx[])
{
    //set up integral of B_nu over theta
    double theta,tr,Bnu;
    theta=x;
    tr=get_tr(tvalue, theta);
    Bnu=get_Bnu(tr,nu);
    dydx[1]=Bnu*sin(theta)*cos(theta);
}

int get_filelines_temp()
{
    //count number of lines in area file
    FILE *fp;
    char ch;
    int linesCount=0;
    fp=fopen(TEMP,"r");
    while((ch=fgetc(fp))!=EOF)
    {
        if(ch=='\n')
            linesCount++;
    }
    fclose(fp);
    return linesCount;
}

int get_index_temp(double tvalue)
{
    //find index of given time in temp file
    int j = 1;
    for(int i=0;i<FILELINES_temp;i++) {
        if(t_file[i]<tvalue)
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
    r=rad_file[i];
    t=t_file[i];
    td=(r/c)*(1-cos(theta));
    tr=t-td;
    do {
        int i=get_index_temp(tr);
        double r2=rad_file[i];
        double td2=(r2/c)*(1-cos(theta));
        tr=t-td2;
        error=(fabs(r-r2))/r;
        r=r2;
    } while (error>1.0E-03);
    
    return tr;
}


double get_Bnu(double tr, double nu)
{
    //calculate B_nu to be integrated*sin(theta) over theta to get the correct F_nu
    int i;
    double temp, r, r2, c2, nu3, exp_, Bnu;
    i=get_index_temp(tr);
    temp=temp_file[i];
    r=rad_file[i];
    r2=r*r;
    c2=c*c;
    nu3=nu*nu*nu;
    exp_=(h*nu)/(k_b*temp);
    Bnu=(4*PI*PI*h/c2)*(nu3)/(exp(exp_)-1); // already did the integral over phi giving an additional factor of 2pi to Planck's formula
    return Bnu;
} // extra pi rn

#include "swift_mag2.h"
