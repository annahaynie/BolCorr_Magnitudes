#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "odeint.h"
#define AREA "/Users/annahaynie/Desktop/Carnegie/CSM/GALEX_NUV_nu.txt"
//#define LOWER 8.74669E+014 //min nu for M2
//#define UPPER 1.87079E+015 //max nu for M2
//#define LOWER 5.18448E+14//min nu for W1
//#define UPPER 1.87079E+15 //max nu for W1
//#define LOWER 5.99286E+014 //min nu for W2
//#define UPPER 1.87079E+015 //max nu for W2
#define LOWER 9.95988e+14//min nu for GALEX
#define UPPER 1.77287e+15 //max nu for GALEX
#define TEMP "/Users/annahaynie/SNEC/Short-runs/short4/Data/T_eff.dat"
#define RADIUS "/Users/annahaynie/SNEC/Short-runs/short4/Data/rad_photo.dat"
#define h 6.62606885E-027 // units = erg*s
#define c 2.99792458E+010 // units = cm/s
#define d 3.08567758E+019 //10 pc in cm for absolute mag calculation
#define k_b 1.380658E-016 //boltzmann constant units erg/K
#define OUTPUT "/Users/annahaynie/Desktop/no_time_delay.txt"
//
//  Swift_Magnitudes.c
//  Calculates magnitudes in Swift bands NOT TAKING TIME DELAY EFFECTS INTO ACCOUNT
//  Calculates F_nu using Planck's Law with r(t) and T(t) so only needs to do a single set of integrals over nu for a given time step
//
//  Created by Anna Haynie on 6/5/19.
//

Ode_Int ODE;
void derivs(double x, double y[], double dydx[]);

void readdata();

int get_filelines_nu();
const int FILELINES_nu=get_filelines_nu();

int get_filelines_temp();
const int FILELINES_temp=get_filelines_temp();

double lambda_file[50000], nu_file[50000], area_file[50000],t_file[50000], temp_file[50000], t2_file[50000], rad_file[50000];

int get_index_temp(double tvalue);
int get_index_nu(double nu);

double get_flux(double tvalue, double nu);
double get_area(double nu);

double tvalue;

int main()
{
    readdata();
    for (int i=1; i<=FILELINES_temp; i++) {
        //for each time step, do the steps to calculate the magnitude at that time
        tvalue=t_file[i];
        
        //do integral over nu
        ODE.init(2);
        ODE.set_bc(1,0.0);
        ODE.set_bc(2,0.0);
        ODE.go(LOWER,UPPER,1.0E+08,5.0,derivs);
        
        FILE *results;
        results = fopen(OUTPUT,"a");
        
        //use output of integrals to calculate magnitudes
        double N,D,r,eta,m, temp;
        r=rad_file[i];
        temp=temp_file[i];
        eta=(r/d)*(r/d);
        N=ODE.get_y(1,ODE.kount);
        D=ODE.get_y(2,ODE.kount);
        m=-2.5*log10(eta*N/D)-48.60;
        fprintf(results,"%lg %lg\n",tvalue,m);
        
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
    for(int i=0;i<FILELINES_nu;i++) {
        fscanf(data, "%lf %lf %lf", &lambda_file[i], &nu_file[i] , &area_file[i]);
    }
    for(int i=0;i<FILELINES_temp;i++) {
        fscanf(data2, "%lf %lf ", &t_file[i], &temp_file[i]);
        fscanf(data3, "%lf %lf", &t2_file[i], &rad_file[i] );
    }
    fclose(data);
    fclose(data2);
    fclose(data3);
}

void derivs(double x, double y[], double dydx[])
{
    //set up to integrate F_nu*area and area over nu
    double nu, flux, area;
    nu=x;
    flux=get_flux(tvalue, nu);
    area=get_area(nu);
    
    dydx[1]=flux*area;
    dydx[2]=area;
}

int get_filelines_nu()
{
    //count the number of lines in the area file
    FILE *fp;
    char ch;
    int linesCount=0;
    fp=fopen(AREA,"r");
    while((ch=fgetc(fp))!=EOF)
    {
        if(ch=='\n')
            linesCount++;
    }
    
    fclose(fp);
    return linesCount;
}

int get_filelines_temp()
{
    //count the number of lines in the temp file
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
    //find the index of a given time in the temp file
    int j = 1;
    for(int i=0;i<FILELINES_temp;i++) {
        if(t_file[i]<tvalue)
            j = i;
    }
    return j;
}

int get_index_nu(double nu)
{
    //find the index of a given nu in the area file
    int j = 1;
    for(int i=0;i<FILELINES_nu;i++) {
        if(nu_file[i]>nu)
            j = i;
    }
    return j;
}

double get_area(double nu)
{
    //linear interpolation to find the area for a given nu called by the integral
    double area, a, b;
    int i,j;
    j=get_index_nu(nu);
    i=j+1;
    a=(nu_file[j]-nu)/(nu_file[j]-nu_file[i]);
    b=(nu-nu_file[i])/(nu_file[j]-nu_file[i]);
    area=a*area_file[i]+b*area_file[j];
    return area;
}

double get_flux(double tvalue, double nu)
{
    //calculate F_nu using Planck's Law 
    int i;
    double temp,c2, flux, nu3, exp_;
    i= get_index_temp(tvalue);
    temp=temp_file[i];
    c2=pow(c,2);
    nu3=pow(nu,3);
    exp_=exp((h*nu)/(k_b*temp));
    flux=(2*h/c2)*(nu3/(exp_-1));
    return flux;
}
