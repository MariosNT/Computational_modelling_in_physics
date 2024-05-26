#pragma once
#include <vector>
#include <fstream>
#include <string>
#include <math.h>
class Stellar_Body
{
public:
	static Stellar_Body** lista_ob; //tablica planet nie zafiksowanych
	static int ilosc_ob;			//liczba tych planet nie zafiksowanych
	Stellar_Body(double,double,double,double,double); //konstruktor
	void Calculate_q1(double);							//liczenie pierwszego kroku
	void Calculate_Step(double);						//liczenie dowolnego kroku
	void Write_Data(int);									//wypisanie pozycji
	void Calc_and_Write_E_ANGM(std::vector<double>);	//wypisanie czasu , energi i momentu pêdu ca³ego uk³adu
	~Stellar_Body();

private:
	double m; //masa planety
	std::vector<double> qx, qy; //wektory pozycji
	std::vector<double> vx,vy; //wektory prêdkoœci
	std::fstream plik; //plik do którego wypisuje
};

///program napisany przez Aleksander Cz¹stkiewicz

