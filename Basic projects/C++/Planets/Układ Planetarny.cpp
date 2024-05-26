#include <iostream>
#include "Stellar_Body.h"
//#include "Stellar_Body.cpp" //zawrzeć tą linijke jeśli plik nie jest w projekcie z resztą

using namespace std;

int main()
{
	vector<double> time; //wektor z przedziałami czasowymi
	int proc = 0;
	int N = 10000;  // ilość próbek
	double h = 0.1; //przedział czasu
	Stellar_Body* Planety[1]; //tablica obiektów nie zafiksowanych

	Planety[0] = new Stellar_Body( 2, 100, 0, 0, 7);  //deklaracja obiektu (masa , x0 -pozycja początkowa , y0 -pozycja początkowa , vx0 - prędkośc początkowa , vy0  -prędkośc początkowa)
	//Planety[1]=new Stellar_Body(3,0,4,2,0); // aby dodac kolejne obiekty trzeba je zinicjalizowac w tablicy

	Planety[0]->lista_ob = Planety; //ważny el. w klasie ta tablice jest jako statyczny ** wskaznik wiec zeby metody miały dostep trzeba przypisać

	Planety[0]->Calculate_q1(h); //przed petla trzeba dla wszytskich policzyc 1 krok bo mamy algorytm 3 punktowy
	//Planety[1]->Calculate_q1(h);
	time.push_back(0);
	for (int i = 0; i < N; i++)
	{
		if (i % (N / 100) == 0)
		{
			cout << proc << "%" << " ";
			if (proc % 10 == 0) cout << endl;
			proc++;
		}
		Planety[0]->Calculate_Step(h); //dla kazdego obiektu liczymy jeden krok
		//Planety[1]->Calculate_Step(h);
		time.push_back(time.at(time.size()-1)+h);
	}
	Planety[0]->Write_Data(1); //wypisanie danych (dac inne wartosci do argumentów żeby pliki miały inną nazwę)
	//Planety[1]->Write_Data(2);

	Planety[0]->Calc_and_Write_E_ANGM(time); //wypisanie energi i momentu pędu dla całego układu od czasu 
											//(wywołac tą metode tylko raz bo ona liczy dla całego układu , więc potrzeba poprostu jednego wywołania żeby dostać się do tablicy planet)
}

///program napisany przez Aleksander Cząstkiewicz

