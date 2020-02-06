#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "odeint.h"

#define h 6.62606885E-027 // units = erg*s
#define c 2.99792458E+010 // units = cm/s
#define d 3.08567758E+019 //10 pc in cm for absolute mag calculation
#define k_b 1.380658E-016 //boltzmann constant units erg/K
#define sigma_sb 5.6704E-05 //stefan boltzmann constant
#define M_bol_sun 4.74 // absolute bolometric mag of the sun
#define L_sun 3.839E+033 // bolometric luminosity of the sun in erg/s

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
//#define LOWER 5.35344e+14 //Bessel B min
//#define UPPER 8.1025e+14 // Bessel B max
//#define LOWER 4.34482e+14 //Bessel V min
//#define UPPER 6.24568e+14 // Bessel V max
//#define LOWER 3.52697e+14 //Bessel R min
//#define UPPER 5.35344e+14 // Bessel R max
#define LOWER 3.29442e+14 //Bessel I min
#define UPPER 4.22243e+14 // Bessel I max


//#define VEGA 0 //AB mags
//#define VEGA 0.79 // Bessel U
//#define VEGA -0.09 //Bessel B
//#define VEGA 0.02 // Bessel V
//#define VEGA 0.21 // Bessel R
#define VEGA 0.45 //Bessel I

#define AREA "/Users/annahaynie/Desktop/Carnegie/CSM/Bessell_I_nu.txt"
#define TEMP "/Users/annahaynie/Desktop/Carnegie/CSM/temp_list.txt"
#define OUTPUT "/Users/annahaynie/Desktop/Carnegie/CSM/BC_Bessel.I.txt"
//
//  BC_mag2.c
//
// This is a "remake" of the make_BC.c program to make sure that I am calculating the corrections in the right way so that it can reproduce the corrections already tabulated in SNEC
//
//  Created by Anna Haynie on 1/14/20.
//
// for a given band at a given temperature we have: BolCorr = M_bol- M_band
// M_bol= -2.5Log10(L/L_sun)+M_bol_sun where L is the luminosity of the star
// L = 4PI*R^2*sigma*T^4 as given by the Stefan-Boltzmann law - this assumption of a black body is a good approximation of the luminosity of the star before and after the peak of the explosion
// M_band= -2.5Log10(eta*N/D) -48.6
// eta=(R/d)^2 where we set d=10pc for absolute magnitude and R is radius
// N=Integral(A*F) over nu where A is the effective area for that band and F is the flux
// D= Integral(A) over nu
// F= (2*h*nu^3/c^2)*(exp(h*nu/kT)-1)^-1 as given by Planck's Law
// For a series of temperatures and frequencies that we expect our explosion to range over, we can calculate F and L and given a specific band we can calculate N and D letting us find the BolCorr at a specific temperature for each band

Ode_Int ODE;
void readdata();
void derivs(double x, double y[], double dydx[]);

int get_filelines_temp();
int temp_lines=get_filelines_temp();
//count number of lines in temp file

int get_filelines_area();
int area_lines=get_filelines_area();
//count number of lines in area file

int get_index_area(double nu); //given a frequency, identify the corresponding line in area file

double get_area(double nu); // given a frequency, find the corresponding area using linear interpolation
double get_fnu(double temp,double nu); // given a temp and frequency, solve for Flux using Planck's Law

double temp_file[40000],lambda_file[40000],nu_file[40000],area_file[40000]; //initialize arrays

double temp;
double nu;
double eta=(1/d)*(1/d); //because Log10R^2 shows up in both M_bol and M_band, it will cancel in their subtraction and won't matter, so I have removed it here and will leave it out when calculating L

int main()
{
    readdata();
    for (int i=0; i<temp_lines; i++) {
        temp=temp_file[i];
        double T4=temp*temp*temp*temp;
        
        ODE.init(2);
        ODE.set_bc(1,0.0);
        ODE.set_bc(2,0.0);
        ODE.go(LOWER,UPPER,1.0E+10,5.0,derivs);
        
        double N, D, R2, L, M_band, M_bol, BC;
        N=ODE.get_y(1,ODE.kount);
        D=ODE.get_y(2,ODE.kount);
        M_band=-2.5*log10(eta*N/D)-48.6-VEGA;
        
        L=4*PI*sigma_sb*T4;
        M_bol=-2.5*log10(L/L_sun)+M_bol_sun;
        
        BC=M_bol-M_band;
        
        ODE.tidy();
        
        FILE *results;
        results = fopen(OUTPUT,"a");
        fprintf(results,"%lg %lg \n",temp, BC);
        fclose(results);
        
    }
    
}
void readdata()
{
    FILE *data, *data2;
    data = fopen(AREA,"r");
    data2 = fopen(TEMP,"r");
    for(int i=0;i<area_lines;i++) {
        fscanf(data, "%lf %lf %lf", &lambda_file[i], &nu_file[i] , &area_file[i]);
    }
    for(int i=0;i<temp_lines;i++) {
        fscanf(data2, "%lf", &temp_file[i]);
        //printf("%lg %lg \n", temp_file[2], temp_file[3]);
    }
    fclose(data);
    fclose(data2);
}
void derivs(double x, double y[], double dydx[]) //set up integral calculation
{
    double area, fnu;
    nu=x;
    area=get_area(nu);
    fnu=get_fnu(temp,nu);
    
    dydx[1]=fnu*area;
    dydx[2]=area;
}

int get_filelines_area() // count number of lines in area file
{
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

int get_filelines_temp() // count number of lines in temp file
{
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

int get_index_area(double nu) // find line in area file given a frequency
{
    int j = 1;
    for(int i=0;i<area_lines;i++) {
        if(nu_file[i]>nu)
            j = i;
    }
    return j;
}
double get_area(double nu) //find area via linear interpolation given nu - used at every step of integrals for N and D
{
    double area, a, b;
    int i,j;
    j=get_index_area(nu);
    i=j+1;
    a=(nu_file[j]-nu)/(nu_file[j]-nu_file[i]);
    b=(nu-nu_file[i])/(nu_file[j]-nu_file[i]);
    area=a*area_file[i]+b*area_file[j];
    return area;
}

double get_fnu(double temp, double nu) //find flux given a frequency and temperature - used at every step of integral for N
{
    double c2, nu3, ex,fnu;
    c2=c*c;
    nu3=nu*nu*nu;
    ex=(h*nu)/(k_b*temp);
    
    fnu=(2*h*nu3)*PI*(1/c2)*(1/(exp(ex)-1));
    return fnu;
}

#include "BC_mag2.h"
