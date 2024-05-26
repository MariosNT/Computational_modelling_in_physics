#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
///////////////////////////////////////////////////////////////////////////////
// ZMIENNE UŻYTKOWNIKA
#define X_0 0.0		//poczatek liczonego zakresu
#define X_N 2	//koniec liczonego zakresu
#define Y_0 1.0		//y(0)-warunek poczatkowy row. rozniczkowego
#define N 10000		//ilosc liczonych pkt
///////////////////////////////////////////////////////////////////////////////
//ZMIENNE PROGRAMOWE
double x[N];                    //tablica argumentow
double y_Eulera[N];             //tablica wartosci funkcji otrzymanymi metoda eulera
double y_Runge_Kutty[N];        //tablica wartosci funkcji otrzymanymi metoda Runge Kutty
double y_Adamsa_Bashforta[N];   //tablica wartosci funkcji otrzymanymi metoda Adamsa Bashforta
double h = (X_N-X_0)/(N-1);     //szerokosc miedzy liczonymi argumentami
fstream plik;
///////////////////////////////////////////////////////////////////////////////
double Euler(double x, double y_0, double h, double (*F)(double, double) )
{
    double y_1 = y_0 + h * F(x, y_0);
    return y_1;
}

double Adam_Beshfort( double x, double y_1, double y_2, double y_3, double h, double (*F)(double, double))
{
    double y = y_3 + h * (5.0 * F(x - 2.0 * h, y_1) - 16.0 * F(x - h, y_2) + 23.0 * F(x, y_3)) / 12.0;
    return y;
}

double Runge_Kutte(double x, double y_0, double h, double (*F)(double, double) )
{
    double k_1 = h * F(x, y_0);
    double k_2 = h * F(x + h / 2.0, y_0 + k_1 / 2.0);
    double k_3 = h * F(x + h / 2.0, y_0 + k_2 / 2.0);
    double k_4 = h * F(x + h, y_0 + k_3);
    double y_1 = y_0 + (k_1 + 2.0 * k_2 + 2.0 * k_3 + k_4) / 6.0;
    return y_1;
}

double Right_side_function(double x, double y) //postac funkcji po praej stornie (dy/dx=f(x,y) ta f(x,y) jest funkcja ktorej piszemy tu jawna postac)
{
    double out = y;
    return out;
}
//////////////////////////////////////////////////////////////////////////////
int main(){
    for (int i = 0; i < N; i++) x[i] = X_0 + h * i;

    y_Eulera[0] = Y_0;
  for (int i = 0; i < N - 1; i++) y_Eulera[i + 1] = Euler(x[i], y_Eulera[i], h, Right_side_function);;

  y_Runge_Kutty[0] = Y_0;
  for (int i = 0; i < N - 1; i++) y_Runge_Kutty[i + 1] = Runge_Kutte(x[i], y_Runge_Kutty[i], h, Right_side_function);

  y_Adamsa_Bashforta[0] = Y_0;
  y_Adamsa_Bashforta[1] = y_Runge_Kutty[1];
  y_Adamsa_Bashforta[2] = y_Runge_Kutty[2];
  for (int i = 2; i < N - 1; i++) y_Adamsa_Bashforta[i + 1] = Adam_Beshfort(x[i], y_Adamsa_Bashforta[i - 2], y_Adamsa_Bashforta[i - 1], y_Adamsa_Bashforta[i], h, Right_side_function);

  plik.open("ivp.data", ios::out);
  plik.precision(10);
  plik.setf( ios::showpoint );


  for(int i=0; i<N; i++)
  {
      plik << x[i] << "   " << y_Eulera[i] << "   " << y_Adamsa_Bashforta[i] << "   " << y_Runge_Kutty[i] << "   " << endl; //Runge Kutty -najlepszy
  }

  plik.close();

  return 0;
}

///kod poprawiony przez : Aleksander Cząstkiewicz
