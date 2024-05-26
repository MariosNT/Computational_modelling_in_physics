// Kamil Wrzos, Piotr Tokarczyk (23.10.2023)

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

/////////////////////////////////////////////////////////////////////////////
// VARIABLES
#define Np       10000      // number of points
#define rE       1E9        // range of calculations
#define A        2.76E+16   // star's parameter nr 1
#define B        2.98E-11   // star's parameter nr 2
#define Rho_zero 1000       // density in the star's centre
#define T        1000       // temperature

// CONSTANTS
#define G        6.67E-11   // gravitational constant
#define pi       3.141592   // pi value
#define R        8.31       // gas constant
#define mi       0.001      // molar mass of the gas

//////////////////////////////////////////////////////////////////////////////
double r[Np];               // distance from the centre
double Rho[Np];             // mass density in the function of the radius
double Rho_an[Np];          // mass density in the function of the radius from the chapter nr 5 (5.2.4)
double F[Np];               // function F defined as in (11.2.9)
double a[Np];               // acceleration of gravity in the function of the radius

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

        // IF THE DENSITY TURNS OUT TO BE NEGATIVE, WE WRITE 0 (see: "Remark concerning the influence of grid constant on calculations")
        /*
            if (Rho[n+1] < 0)
                Rho[n+1] = 0;
        */
    }

    for (int i=0; i<Np; i++)    // calculating Rho_an from chapter nr 5 (5.2.4) (excercise nr 1)
        Rho_an[i] = 0.04 * exp(-1.0*i*h); 
}

void calculate_a()     // function calculating the acceleration of gravity value from the formula (11.1.5)
{
    for (int i=1; i<Np-1; i++)
    {
        a[i] = 0.0;

        for (int j=0; j<i+1; j++)
        {
            a[i] += 0.5*h*( pow(r[j+1],2)*Rho[j+1] + pow(r[j],2)*Rho[j] ); // integration with the usage of trapezoidal rule
        }

        a[i] *= -4.0*G*pi/pow(i*h,2);        
    }

    a[0] = a[1];
    a[Np-1] = a[Np-2];
}

double radius()     // function calculating the star's radius
{
  int i = 3;

  while (a[i-1] >= a[i] && i<Np-1)
  {
    i++;            // wee seek the extremum
  }
    
  return r[i-1];
}

void Save_data()        // function saving data
{
    fstream plik ("star.data", ios::out);     // we open "star.data" file

    for (int i = 0; i < Np; i++)    // we save the data
    {
        plik
        << r[i]      << " "         // 1 - distance from the centre of the star
        << Rho[i]    << " "         // 2 - star's density
        << a[i]      << " "         // 3 - acceleration of gravity
        << Rho_an[i] << " "         // 4 - star's density from the chapter nr 5 (excercise nr 1)
        << endl;
    }

    plik.close();                   // we close the file
}

int main()
{
    // functions calculating the density
    INIT();                 
    Adams_Bashforth();      

    // function calculating the acceleration of gravity
    calculate_a();

    // we save the data
    Save_data();      

    // we print the mass and the radius of the star
    cout << "The star's mass is: " << mass() << endl;
    cout << "The star's radius is: " << radius() << endl;

    return 0;
}