// Kamil Wrzos, Piotr Tokarczyk (23.10.2023)

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

/////////////////////////////////////////////////////////////////////////////
// VARIABLES
#define Np        1000      // number of points
#define rE        100.0     // range of calculations
#define A         2.76E16   // star's parameter nr 1
#define B         2.98E-11  // star's parameter nr 2
#define T         10        // temperature
#define N         10        // number of iterations of the density for which
                            // we calculate star's pressure in the centre and its mass
double Rho_zero = 1000;     // initial density value in the centre

// CONSTANTS
#define G         6.67E-11  // gravitational constant
#define pi        3.141592  // pi value
#define R         8.31      // gas constant
#define mi        0.001     // molar mass of the gas

//////////////////////////////////////////////////////////////////////////////
double r[Np];               // distance from the centre
double Rho[Np];             // mass density in the function of the radius
double F[Np];               // function F defined as in (11.2.9)
double M[N];                // star's total mass
double P_zero[N];           // star's pressure in the centre

//////////////////////////////////////////////////////////////////////////////
double h = rE / (Np - 1);   // grid constant
fstream plik;               // destination file

//////////////////////////////////////////////////////////////////////////////

double P (double rho)       // function calculating star's pressure
{
    double P;

    P = A * ( exp(B * rho * rho) - 1) + R * T * rho / mi; // equation (11.2.11) with the component representing the ideal gas law
    return P;
}

double P_prim (double rho)  // function calculating the derivative of pressure with respect to density
{
    double P_prim;

    P_prim =  2 * A * B * rho * exp(B * rho * rho) + R * T / mi; // equation (11.2.12) with the component representing the ideal gas law
    return P_prim;
}

double mass()       // function calculating star's mass
{
    double m = 0.0;

    for (int i=0; i<Np-1; i++)
    {
        m += 0.5*h*( pow(r[i+1],2)*Rho[i+1] + pow(r[i],2)*Rho[i] ); // integration with the usage of trapezoidal rule
    }

    m *= 4*pi;

    return m;
}

void INIT() // function initialising tables and values
{
    Rho[0] = Rho_zero;              // we set the density value in the centre
    Rho[1] = Rho[0];                // and in the distance h from the centre
    F[0] = 0.0;                     // and also value of F function in the centre
    
    for (int i = 0; i < Np; i++)  
    {
        r[i] = i * h;              
    }
}

void Adams_Bashforth()      // function calculating F value and consequently Rho
{
    for (int n = 1; n < Np; n++)    // loop calculating 1,2,3,...,Np-1 element (not including zero_th one) 
    {
        // CALCULATING THE SUM IN EQUATION (11.2.9)
        F[n] = 0.0;     // we set the initial value to 0
        for (int i = 1; i <= n; i++)
        {
            F[n] += (Rho[i - 1] * pow(i - 1, 2) + Rho[i] * pow(i, 2)); // we add consecutive expressions
        }
        
        // CALCULATING THE F[n] VALUE (11.2.9)
        F[n] = -2.0 * pi * G * h * Rho[n] / pow(n, 2) / P_prim(Rho[n]) * F[n];

        // CALCULATING THE DENSITY Rho[n+1] (11.2.10)
        Rho[n + 1] = Rho[n] + 0.5 * h * (3.0 * F[n] - F[n - 1]);
    }
}

void Save_data()        // function saving data
{
    fstream plik ("star.data", ios::out);     // we open "star.data" file
        for (int i=0; i<N; i++)
        {
            plik
            << P_zero[i] << " "         // 1 - pressure in the centre
            << M[i]      << " "         // 2 - total mass
            << endl;
        }
        
    plik.close();                       // we close the file
}

int main()
{
    for (int i=0; i<N; i++)
    {
        INIT();                      // we initialise the table r and values
        Adams_Bashforth();           // calculate star's density
        M[i] = mass();               // calculate the total mass 
        P_zero[i] = P(Rho_zero);     // calculate the pressure in the centre
        Rho_zero += 1000;            // and change in the suitable way the density in the centre
    }

    Save_data();                     // we save the data
    return 0;
}