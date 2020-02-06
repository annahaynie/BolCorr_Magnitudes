#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "odeint.h"
#define RAD "/Users/annahaynie/SNEC/SNEC-1.01/Data/rad_photo.dat"
#define LUM "/Users/annahaynie/SNEC/SNEC-1.01/Data/lum_observed.dat"
#define TEMP "/Users/annahaynie/SNEC/SNEC-1.01/Data/T_eff.dat"
#define BCFILE "/Users/annahaynie/SNEC/BolCorr3.txt"
#define OUTPUT "/Users/annahaynie/SNEC/SNEC-1.01/Data/BC_mag_td.txt"


#define m_bol_sun 4.74 // absolute bolometric mag of the sun
#define L_sun 3.826E+033 // bolometric luminosity of the sun in ergs/s
#define sigma_sb 5.6704E-05 //stefan boltzmann constant units erg*cm^-2*K^-4
#define z 0.1665
#define c 2.99792458E+010 // units = cm/s

#define VEGA_U 0.79 // Bessel U
#define VEGA_B -0.09 //Bessel B
#define VEGA_V 0.02 // Bessel V
#define VEGA_R 0.21 // Bessel R
#define VEGA_I 0.45 //Bessel I

//
//  BC_mag_td.c
//  
//
//  Created by Anna Haynie on 2/3/20.
//

Ode_Int ODE;
void derivs(double x, double y[], double dydx[]);

void readdata();

int get_filelines();
const int FILELINES=get_filelines();

int get_filelines_BC();
const int FILELINES_BC=get_filelines_BC();

int get_index(double tvalue);
int get_index_BC(double temp);
double t_r(double tvalue, double theta);
double get_lum(double t_r_value);
double get_BC(double temp);

double t_file[50000], lum_file[50000],t2_file[50000], rad_file[50000],t3_file[50000], temp_file[50000], BC_array[5000], m[5000],temp2_file[50000], ptf_file[5000], u_file[50000],g_file[50000], r_file[50000], i_file[50000], z_file[50000], U_file[50000], B_file[50000], V_file[50000], R_file[50000], I_file[50000],  NUV_file[50000], FUV_file[50000], M2_file[50000], W1_file[50000], W2_file[50000];

double tvalue;
double temp;
double error=1.0E+005;

int main()
{
    readdata();
   
    for(int i=1;i<FILELINES;i++) {
        //find updated luminosity for each time step
        tvalue=t_file[i];
        //do integral over theta
        ODE.init(1);
        ODE.set_bc(1,0.0);
        ODE.go(0.0,PI/2,PI*1.0e-2,5.0e-7,derivs);
        
        FILE *results;
        results = fopen(OUTPUT,"a");
        double L, BC;
        L= ODE.get_y(1,ODE.kount);
        //int i = get_index(tvalue);
        temp = temp_file[i];
        //double lum = lum_file[i];
        
        get_BC(temp);
        
        for (int j=0; j<11; j++) {
            m[j]=m_bol_sun-2.5*log10(L/L_sun)-BC_array[j];
            //printf("%lf %lf %lf %lf %lf \n", tvalue, lum, L, L/L_sun, BC_array[j]);
        }
        
        ODE.tidy();
        //printf("%lf \n", temp);
        fprintf(results,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", tvalue, m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14]);
        //fprintf(results,"%lf %lf %lf \n", tvalue, lum, L);
        fclose(results);
    }
}

void readdata()
{
    FILE *data, *data2, *data3, *data4;
    data = fopen(LUM,"r");
    data2 = fopen(RAD,"r");
    data3 = fopen(TEMP,"r");
    data4 = fopen(BCFILE,"r");
    for(int i=0;i<FILELINES;i++) {
        fscanf(data, "%lf %lf", &t_file[i], &lum_file[i] );
        fscanf(data2, "%lf %lf", &t2_file[i], &rad_file[i] );
        fscanf(data3, "%lf %lf", &t3_file[i], &temp_file[i] );
        //printf("%i %lf \n",i,temp_file[i]);
    }
    for(int i=0;i<FILELINES_BC;i++) {
        fscanf(data4, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &temp2_file[i], &u_file[i], &g_file[i], &r_file[i], &i_file[i], &z_file[i], &U_file[i], &B_file[i], &V_file[i], &R_file[i], &I_file[i],  &NUV_file[i], &FUV_file[i], &M2_file[i], &W1_file[i], &W2_file[i]);
        //printf("%i %lf %lf %lf %lf %lf %lf \n",i,temp2_file[i], u_file[i],g_file[i],r_file[i],i_file[i],z_file[i]);
        //fscanf(data4, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", &temp2_file[i], &ptf_file[i], &u_file[i], &g_file[i], &r_file[i], &i_file[i], &z_file[i], &U_file[i], &B_file[i], &V_file[i], &R_file[i], &I_file[i]);
        //printf("%i %lf \n", i, temp2_file[i]);
    }
    fclose(data);
    fclose(data2);
    fclose(data3);
    fclose(data4);
}

void derivs(double x, double y[], double dydx[])
{
    //set up integral for L(tr)*sin(theta) over theta
    double theta, t_r_value, lum;
    theta = x;
    t_r_value= t_r(tvalue, theta);
    if (t_r_value<0) {
        lum=0;
    } else {
        lum=get_lum(t_r_value);
    }
    dydx[1] = sin(theta)*lum;
    //printf("%lf %lf %lf \n", tvalue,t_r_value,lum);
}

int get_filelines()
{
    //count number of lines in luminosity file
    FILE *fp;
    char ch;
    int linesCount=0;
    fp=fopen(LUM,"r");
    while((ch=fgetc(fp))!=EOF)
    {
        if(ch=='\n')
            linesCount++;
    }
    fclose(fp);
    return linesCount;
}

int get_filelines_BC()
{
    //count number of lines in luminosity file
    FILE *fp;
    char ch;
    int linesCount=0;
    fp=fopen(BCFILE,"r");
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
        if(t_file[i]<tvalue)
            j = i;
    }
    return j;
}

int get_index_BC(double temp)
{
    //find index of a given time in the luminosity file
    int j = 1;
    for(int i=0;i<FILELINES_BC;i++) {
        if(temp2_file[i]<temp)
            j=i;
    }
    //printf("%i %lf %lf \n", j, temp, temp2_file[j]);
    return j;
}

double t_r(double tvalue, double theta)
{
    //taking time delay effects into account
    // calculates retardation time given some theta on the star
    double t_r_value, r, r2, r3, d,t, d2;
    int i,j;
    i=get_index(tvalue);
    t=t2_file[i];
    r=rad_file[i];
    d=(r/c)*(1-cos(theta));
    t_r_value=t-d;
   // printf("%lg %lg %lg %lg \n", theta, r, d, t_r_value);
    //iterate through this same process until the error between r(t) and r(tr) is sufficiently small
    do {
        int i=get_index(t_r_value);
        double r2=rad_file[i];
        d2=(r2/c)*(1-cos(theta));
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
    a = (t_file[j]-t_r_value)/(t_file[j]-t_file[i]);
    b = (t_r_value-t_file[i])/(t_file[j]-t_file[i]);
    lum = (a*lum_file[i]+b*lum_file[j]);
    return lum;
}

double get_BC(double temp)
{
    double area, a, b;
    int i,j;
    i=get_index_BC(temp);
    j=i+1;
    a=(temp2_file[j]-temp)/(temp2_file[j]-temp2_file[i]);
    b=(temp-temp2_file[i])/(temp2_file[j]-temp2_file[i]);
    //printf("%i %i %lf %lf %lf %lf %lf \n", i,j,temp,temp2_file[i], temp2_file[j], a,b);
    
    BC_array[0]=a*u_file[i]+b*u_file[j];
    BC_array[1]=a*g_file[i]+b*g_file[j];
    BC_array[2]=a*r_file[i]+b*r_file[j];
    BC_array[3]=a*i_file[i]+b*i_file[j];
    BC_array[4]=a*z_file[i]+b*z_file[j];
    BC_array[5]=a*U_file[i]+b*U_file[j]+VEGA_U;
    BC_array[6]=a*B_file[i]+b*B_file[j]+VEGA_B;
    BC_array[7]=a*V_file[i]+b*V_file[j]+VEGA_V;
    BC_array[8]=a*R_file[i]+b*R_file[j]+VEGA_R;
    BC_array[9]=a*I_file[i]+b*I_file[j]+VEGA_I;
    //BC_array[10]=a*NUV_file[i]+b*NUV_file[j];
    //BC_array[11]=a*FUV_file[i]+b*FUV_file[j];
    //BC_array[12]=a*M2_file[i]+b*M2_file[j];
    //BC_array[13]=a*W1_file[i]+b*W1_file[j];
    //BC_array[14]=a*W2_file[i]+b*W2_file[j];
    
    //printf("%i %i %lf %lf %lf %lf %lf %lf \n",i,j,temp,temp2_file[i],temp2_file[j],a,b,BC_array[0]);
    //return BC_array;
}


#include "BC_mag_td.h"
