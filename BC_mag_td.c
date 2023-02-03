#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "odeint.h"

#define RAD "/Users/annahaynie/Desktop/Carnegie/CSM/SNEC-m14.5_K2.5e17_R1900/06/rad_photo.dat"
#define LUM "/Users/annahaynie/Desktop/Carnegie/CSM/SNEC-m14.5_K2.5e17_R1900/06/lum_observed.dat"
#define TEMP "/Users/annahaynie/Desktop/Carnegie/CSM/SNEC-m14.5_K2.5e17_R1900/06/T_eff.dat"
#define BCFILE "/Users/annahaynie/SNEC/BolCorr3.txt"
#define OUTPUT "/Users/annahaynie/Desktop/Carnegie/CSM/SNEC-m14.5_K2.5e17_R1900/06/BC_mag_td_z2.txt"

#define START 10223
#define LINESrad 44619
#define LINEStemp 44619
#define LINESbc 999

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
// Calculating magnitudes in all 15 bands using bolometric corrections for speed
// possibly not as accurate as we need it to be around the time of sbo where a black body is not
// a good approximation for the star, but good and fast elsewhere
// for RSG+CSM models, we expect this to work everywhere becasue they are diffusions dominated events
// not time delay dominated
//
// uses equation m_band = m_sun - 2.5log(L/Lsun) - BC_band where L is the corrected luminosity,
// which gets fully time delay corrected here the same was as in updated_lum4.c, and BC_band
// is a function of the temperature at a given time in the TEMP input file - this temp may or may not be
// time delay corrected based on the input file - corrected temps will come from fluxtemp.c
//
//  Created by Anna Haynie on 2/3/20.
//

Ode_Int ODE;
void derivs(double x, double y[], double dydx[]);

void readdata();

int get_index_rad(double tvalue);
int get_index_BC(double temp);
int get_index_lum_temp(double tvalue);
double t_r(double tvalue, double theta);
double get_lum(double t_r_value);
double get_BC(double temp);

double time_file[70000], lum_file[70000];
double time2_file[70000], temp_file[70000];
double time3_file[70000], rad_file[70000];
double BC_array[7000], m[7000],BC_T[7000], BC_Tz[7000];
double temp2_file[75000], u_file[75000],g_file[75000], r_file[75000], i_file[75000], z_file[75000], U_file[75000], B_file[75000], V_file[75000], R_file[75000], I_file[75000],  NUV_file[75000], FUV_file[75000], M2_file[75000], W1_file[75000], W2_file[75000];

double time;
double tvalue;
double temp;
double T;
double error=1.0E+005;

int main()
{
    readdata();
   
    for(int i = START; i < LINEStemp; i++) {
       
        tvalue = time_file[i];
        T = temp_file[i];
        ODE.init(1);
        ODE.set_bc(1,0.0);
        ODE.go(0.0,PI/2,PI*1.0e-2,5.0e-7,derivs); //do theta integral to find updated luminosity
        
        FILE *results;
        results = fopen(OUTPUT,"a");
        double L, Lz, BC;
        
        L= ODE.get_y(1,ODE.kount); // get updated luminosity from integral
        
    /* Instert Redshift ********* */
        time = tvalue*(1+z);
        temp = T/(1+z);
        Lz = L/((1+z)*(1+z));
    /* ************************** */
        
        
        get_BC(temp); //find the bolometric corrections associated with the temperature at this time step
        
        for (int j=0; j<15; j++) {
            m[j]=m_bol_sun-2.5*log10(Lz/L_sun)-BC_array[j]; // calculate mag in each band
        }
        
        ODE.tidy();
        fprintf(results,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", time, m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14]);

        fclose(results);
    }
}

void readdata()
{
    FILE *data, *data2, *data3, *data4;
    data = fopen(LUM,"r");
    data2 = fopen(TEMP,"r");
    data3 = fopen(RAD,"r");
    data4 = fopen(BCFILE,"r");
    for (int i=0; i<LINEStemp; i++) {
        fscanf(data, "%lf %lf", &time_file[i], &lum_file[i]);
        fscanf(data2, "%lf %lf", &time2_file[i], &temp_file[i]);
    }
    for(int i=0;i<LINESrad;i++) {
        fscanf(data3, "%lf %lf", &time3_file[i], &rad_file[i]);
    }
    for(int i=0;i<LINESbc;i++) {
        fscanf(data4, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &temp2_file[i], &u_file[i], &g_file[i], &r_file[i], &i_file[i], &z_file[i], &U_file[i], &B_file[i], &V_file[i], &R_file[i], &I_file[i],  &NUV_file[i], &FUV_file[i], &M2_file[i], &W1_file[i], &W2_file[i]);
    }
    fclose(data);
    fclose(data2);
    fclose(data3);
    fclose(data4);
}

void derivs(double x, double y[], double dydx[]) //details for theta integral
{
    double theta, t_r_value, lum;
    theta = x;
    t_r_value= t_r(tvalue, theta); //get time delay at each angle
    if (t_r_value>0.0) {
        lum=get_lum(t_r_value); // find lum at delayed time
    } else {
        lum = 0.0;
    }
    dydx[1] = sin(theta)*lum; // integrate over surface area of star
    // no 2PI here because our integral is really of the flux but 2PI*R*R*F = L
}

int get_index_lum_temp(double tvalue)
{
    //find index of a given time in the temperature file
    int j = 1;
    for(int i=0;i<LINEStemp;i++) {
        if(time2_file[i]<tvalue)
            j = i;
    }
    return j;
}

int get_index_rad(double tvalue)
{
    //find index of a given time in the radius file
    int j = 1;
    for(int i=0;i<LINESrad;i++) {
        if(time3_file[i]<tvalue)
            j = i;
    }
    return j;
}

int get_index_BC(double temp)
{
    //find index of a given temperature in the Bolometric Corrections table
    int j = 1;
    for(int i=0;i<LINESbc;i++) {
        if(temp2_file[i]<temp)
            j=i;
    }

    return j;
}

double t_r(double tvalue, double theta)
{
    //taking time delay effects into account
    // calculates retardation time given some theta on the star
    double t_r_value, r, r2, r3, d,t, d2;
    int i,j;
    i=get_index_rad(tvalue);
    r=rad_file[i];
    d=(r/c)*(1-cos(theta));
    t_r_value=tvalue-d;
    //iterate through this same process until the error between r(t) and r(tr) is sufficiently small
    do {
        int i=get_index_rad(t_r_value);
        double r2=rad_file[i];
        d2=(r2/c)*(1-cos(theta));
        t_r_value=tvalue-d2;
        error = (fabs(r-r2))/r;
        r=r2;
    } while (error>1.0E-03);
    
    return t_r_value;
}

double get_lum(double t_r_value)
{
    //linear interpolation of luminosity
    //find L at a delayed time calculated from a given theta called by the integral
    double lum, a, b;
    int i, j;
    i=get_index_lum_temp(t_r_value);
    j=i+1;
    a = (time_file[j]-t_r_value)/(time_file[j]-time_file[i]);
    b = (t_r_value-time_file[i])/(time_file[j]-time_file[i]);
    lum = (a*lum_file[i]+b*lum_file[j]);
    return lum;
}

double get_BC(double temp)
{
    //linear interpolation of bolometric corrections in each band
    //find the BC in each band for a given temperature
    double area, a, b;
    int i,j;
    i=get_index_BC(temp);
    j=i+1;
    a=(temp2_file[j]-temp)/(temp2_file[j]-temp2_file[i]);
    b=(temp-temp2_file[i])/(temp2_file[j]-temp2_file[i]);
    
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
    BC_array[10]=a*NUV_file[i]+b*NUV_file[j];
    BC_array[11]=a*FUV_file[i]+b*FUV_file[j];
    BC_array[12]=a*M2_file[i]+b*M2_file[j];
    BC_array[13]=a*W1_file[i]+b*W1_file[j];
    BC_array[14]=a*W2_file[i]+b*W2_file[j];
    
}


#include "BC_mag_td.h"
