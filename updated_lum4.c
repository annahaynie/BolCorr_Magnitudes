#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "odeint.h"
#define RADFILE "rad_photo.txt"
#define LUMFILE "lum_observed.txt"
#define OUTPUT "test.txt"

//  updated_lum4.c
//  
//  Finding the integrated luminosity by taking time dilation effects into account 
//  the tr calculation relies on r(t), so a second iteraton would find tr' using r(tr) and then a third gets tr'' using r(tr') etc...
//  until the error between r(t) and r(tr) is satisfactorily small
//
//  Created by Anna Haynie on 5/31/19.

Ode_Int ODE;
void derivs(double x, double y[], double dydx[]);

void readdata();

int get_filelines();
const int FILELINES=get_filelines();

int get_index(double tvalue);
double t_r(double tvalue, double theta);
double get_lum(double t_r_value);

double tvalue_file[50000], lumvalue_file[50000],t2value_file[50000], radvalue_file[50000];

double tvalue=20;
double error=1.0E+005;

int main()
{
    readdata();
   
    for(int i=1;i<FILELINES;i++) {
        //find updated luminosity for each time step
        tvalue=tvalue_file[i];
        //do integral over theta
        ODE.init(1);
        ODE.set_bc(1,0.0);
        ODE.go(0.0,PI/2,PI*1.0e-2,5.0e-7,derivs);
        
        FILE *results;
        results = fopen(OUTPUT,"a");
        double updated_lum_result;
        updated_lum_result= ODE.get_y(1,ODE.kount);
        fprintf(results,"%lg %lg\n",tvalue_file[i],updated_lum_result);
        ODE.tidy();
        fclose(results);
    //printf("time= %lf lum= %lf \n", tvalue ,updated_lum_result);
    }
}

void readdata()
{
    FILE *data_file, *data_file2;
    data_file = fopen(LUMFILE,"r");
    data_file2 = fopen(RADFILE,"r");
    for(int i=0;i<FILELINES;i++) {
        fscanf(data_file, "%lf %lf", &tvalue_file[i], &lumvalue_file[i] );
        fscanf(data_file2, "%lf %lf", &t2value_file[i], &radvalue_file[i] );
    }
    fclose(data_file);
    fclose(data_file2);
}

void derivs(double x, double y[], double dydx[])
{
    //set up integral for L(tr)*sin(theta) over theta
    double theta, t_r_value, lum;
    theta = x;
    t_r_value= t_r(tvalue, theta);
    lum=get_lum(t_r_value);
    dydx[1] = sin(theta)*lum;
    printf("%lf %lf \n", t_r_value,lum);
}

int get_filelines()
{
    //count number of lines in luminosity file
    FILE *fp;
    char ch;
    int linesCount=0;
    fp=fopen(LUMFILE,"r");
    while((ch=fgetc(fp))!=EOF)
    {
        if(ch=='\n')
            linesCount++;
    }
    fclose(fp);
    return linesCount;
}

int get_index(double tvalue)
{
    //find index of a given time in the luminosity file
    int j = 1;
    for(int i=0;i<FILELINES;i++) {
        if(tvalue_file[i]<tvalue)
            j = i;
    }
    return j;
}

double t_r(double tvalue, double theta)
{
    //taking time delay effects into account
    // calculates retardation time given some theta on the star
    double t_r_value, r, r2, r3, d,t, d2;
    int i,j;
    i=get_index(tvalue);
    t=t2value_file[i];
    r=radvalue_file[i];
    d=(r/3.0E+10)*(1-cos(theta));
    t_r_value=t-d;
   // printf("%lg %lg %lg %lg \n", theta, r, d, t_r_value);
    //iterate through this same process until the error between r(t) and r(tr) is sufficiently small
    do {
        int i=get_index(t_r_value);
        double r2=radvalue_file[i];
        d2=(r2/3.0E+10)*(1-cos(theta));
        t_r_value=t-d2;
        error = (fabs(r-r2))/r;
        r=r2;
    } while (error>1.0E-03);
    
    //printf("%lf %lf %lf %lf %lf \n", t, r, theta, d2, t_r_value);
    return t_r_value;
}

double get_lum(double t_r_value)
{
    //linear interpolation to find luminosity for a given tr value calculated from a given theta called by the integral
    double lum, a, b;
    int i, j;
    i=get_index(t_r_value);
    j=i+1;
    a = (tvalue_file[j]-t_r_value)/(tvalue_file[j]-tvalue_file[i]);
    b = (t_r_value-tvalue_file[i])/(tvalue_file[j]-tvalue_file[i]);
    lum = (a*lumvalue_file[i]+b*lumvalue_file[j]);
    return lum;
}


#include "updated_lum4.h"
