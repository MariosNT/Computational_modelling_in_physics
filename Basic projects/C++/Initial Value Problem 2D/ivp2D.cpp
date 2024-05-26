#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

///////funkcje do liczenia równañ ró¿niczkowych 2 rzêdu //////////

void Euler_D2(double x,double y1_p,double y2_p,double & y1_k,double & y2_k,
           double h,double(*Deriv_1)(double,double,double),double(*Deriv_2)(double,double,double))
           {
                y1_k=y1_p+h*Deriv_1(x,y1_p,y2_p);
                y2_k=y2_p+h*Deriv_2(x,y1_p,y2_p);
           }

void Runge_Kutty_D2(double x,double y1_p,double y2_p,double & y1_k,double & y2_k,
           double h,double(*Deriv_1)(double,double,double),double(*Deriv_2)(double,double,double))
           {
                double k_11 = h*Deriv_1(x,y1_p,y2_p);
                double k_21 = h*Deriv_2(x,y1_p,y2_p);

                double k_12 = h*Deriv_1(x+h/2.0, y1_p+k_11/2.0, y2_p+k_21/2.0);
                double k_22 = h*Deriv_2(x+h/2.0, y1_p+k_11/2.0, y2_p+k_21/2.0);

                double k_13 = h*Deriv_1(x+h/2.0, y1_p+k_12/2.0, y2_p+k_22/2.0);
                double k_23 = h*Deriv_2(x+h/2.0, y1_p+k_12/2.0, y2_p+k_22/2.0);

                double k_14 = h*Deriv_1(x+h, y1_p+k_13, y2_p+k_23);
                double k_24 = h*Deriv_2(x+h, y1_p+k_13, y2_p+k_23);

                y1_k = y1_p + (k_11+2.0*k_12+2.0*k_13+k_14)/6.0;
                y2_k = y2_p + (k_21+2.0*k_22+2.0*k_23+k_24)/6.0;
           }

//////////////////// zmienne u¿ytkownika i funkcje ukladów równañ ///////////////////////

#define N 10001
#define T_p 0.0
#define T_k 5.0
#define Y_0 0.0
#define V_0 2.0
#define k 4.0

double Second_Deriv(double x,double y,double v) //pochodna podstawienia
{
    return -k*y;
}

double First_Deriv(double x,double y,double v) //podstawienie za dy/dx
{
    return v;
}

////// funkcja glowna //////

int main()
{
    double h=(T_k-T_p)/(N-1);
    double t[N];
    double y[N];
    double v[N];
    double E[N];
    t[0]=T_p;
    y[0]=Y_0;
    v[0]=V_0;
    for(int i=1;i<N;i++) t[i]=t[i-1]+h;
    for(int i=0;i<N-1;i++)
    {
        Runge_Kutty_D2(t[i],y[i],v[i],y[i+1],v[i+1],h,First_Deriv,Second_Deriv);
        //Euler_D2(t[i],y[i],v[i],y[i+1],v[i+1],h,First_Deriv,Second_Deriv);
        E[i]=v[i]*v[i]/2+k*y[i]*y[i]/2;
    }
    E[N-1]=v[N-1]*v[N-1]/2+k*y[N-1]*y[N-1]/2;
    fstream plik;
    plik.open("ivp2D.data", ios::out);
    plik.precision(6);
    for(int i=0;i<N;i++) plik<<t[i]<<" "<<y[i]<<" "<<v[i]<<" "<<E[i]<<endl;
    return 0;
}

/// program napisany przez : Aleksander Cz¹stkiewicz