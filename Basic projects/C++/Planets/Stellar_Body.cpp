#include "Stellar_Body.h"

#define M 10000 // masa zafiksowanego obiektu
#define X 0.0 //pozycja zafiksowanego obiektu
#define Y 0.0 // pozycja zafiksowanego obiektu
#define G 1 // stala grawitacyjna

int Stellar_Body::ilosc_ob = 0;
Stellar_Body** Stellar_Body::lista_ob = NULL; //inicjalizacja static'ów


Stellar_Body::Stellar_Body( double b, double c, double d, double e, double f)
{
	ilosc_ob++;
	m = b;
	qx.push_back(c); qy.push_back(d);
	vx.push_back(e); vy.push_back(f);

}

void Stellar_Body::Calculate_q1(double h)
{
	double r,rx,ry;
	rx = qx.at(0) - X;  // wsp x promienia miezy zafiksowanym obiektem a liczonym obiektem
	ry = qy.at(0) - Y; //wsp y promienia miezy zafiksowanym obiektem a liczonym obiektem
	r = sqrt(rx * rx + ry * ry); // dlugosc tego promienia
	double ax = -G * M / r / r / r * rx; //liczenie przyspieszenia
	double ay = -G * M / r / r / r * ry; //
	for (int i = 0; i < ilosc_ob; i++)
	{
		if (lista_ob[i] != NULL && this!=lista_ob[i])
		{
			rx = qx.at(0) - lista_ob[i]->qx.at(0); // wsp x promienia miêdzy innym obiektem i a liczonym obiektem
			ry = qy.at(0) - lista_ob[i]->qy.at(0); // wsp y promienia miêdzy innym obiektem i a liczonym obiektem
			r = sqrt(rx * rx + ry * ry);
			ax += -G * lista_ob[i]->m / r / r / r * rx; //liczenie przyspieszenia
			ay += -G * lista_ob[i]->m / r / r / r * ry; //
		}
	}
	double q1x = qx.at(0) + vx.at(0) * h + ax * h * h / 2;
	double q1y = qy.at(0) + vy.at(0) * h + ay * h * h / 2;
	qx.push_back(q1x); qy.push_back(q1y);
}

void Stellar_Body::Calculate_Step(double h)
{
	int step = qx.size()-1;
	double r, rx, ry;
	rx = qx.at(step) - X;
	ry = qy.at(step) - Y;
	r = sqrt(rx * rx + ry * ry);
	double ax = -G * M / r / r / r * rx;
	double ay = -G * M / r / r / r * ry;
	for (int i = 0; i < ilosc_ob; i++)
	{
		if (lista_ob[i] != NULL && this !=lista_ob[i])
		{
			rx = qx.at(step) - lista_ob[i]->qx.at(step);
			ry = qy.at(step) - lista_ob[i]->qy.at(step);
			r = sqrt(rx * rx + ry * ry);
			ax += -G * lista_ob[i]->m / r / r / r * rx;
			ay += -G * lista_ob[i]->m / r / r / r * ry;
		}
	}
	double qx_n_1 = 2 * qx.at(step) - qx.at(step - 1) + h * h * ax; // tea funkcja jest taka sama jak Calculate_q1(double h) tylko ze w wyliczeniu nastepnego punktu korzystamy z
	double qy_n_1 = 2 * qy.at(step) - qy.at(step - 1) + h * h * ay; // punktów poprzednich - formu³a 3 punktowa
	qx.push_back(qx_n_1); qy.push_back(qy_n_1);
	double vnx = (qx.at(step+1) - qx.at(step - 1)) / h / 2; //formu³a 3 punktowa dla prêdkoœæi
	double vny = (qy.at(step+1) - qy.at(step - 1)) / h / 2;
	vx.push_back(vnx); vy.push_back(vny);
}

void Stellar_Body::Write_Data(int id)
{
	std::string  pocz = "obiekt", kon = ".data";
	std::string fin;
	int a = id % 10,b=0;
	if (id > 9)
	{
		b = (id - a) / 10;
		fin = pocz + " " + char(b + 48) + char(a + 48) + kon;
	}
	else
	{
		fin = pocz + " " + char(a + 48) + kon;
	}
	 //u¿yæ tego kawa³ka kodu jeœli nastêpna czêœæ nie dzia³a przez starsz¹ wersje c++
	//////////
	//std::string fin = pocz + " " + std::to_string(id) + kon; //mo¿na u¿yæ przy nowszej wersji cpp dla wygody
	//////////
	plik.open(fin.c_str(),std::ios::out);
	for (unsigned int i = 0; i < qx.size(); i++)
	{
		plik << qx.at(i) << " " << qy.at(i);
		if (i < qx.size() - 1) plik << "\t\t\t" << vx.at(i) << " " << vy.at(i);
		plik << std::endl;
	}
	plik.close();
}

void Stellar_Body::Calc_and_Write_E_ANGM(std::vector<double> time)
{
	std::vector<double>  Energy, AngMom;
	double E = 0, AM = 0,rx,ry,vx,vy;
	for (int i = 0; i < qx.size() - 1; i++)
	{
		E = 0; AM = 0;
		for (int j = 0; j < ilosc_ob; j++)//sumujemy energie i moment pêdu wszystkich niezafiksowanych cia³
		{
			rx = lista_ob[j]->qx.at(i);
			ry = lista_ob[j]->qy.at(i);
			vx = lista_ob[j]->vx.at(i);
			vy = lista_ob[j]->vy.at(i);
			double _m_ = lista_ob[j]->m;
			E += _m_ * (vx * vx + vy * vy) / 2 - G * M * _m_ / sqrt(rx * rx + ry * ry); // energia kinetyczna plus potencjalna
			AM += _m_ * (rx * vy - ry * vx); //z-towa skladowa moemntu pêdu policzona z iloczynu wektorowego r x p
		}
		Energy.push_back(E); AngMom.push_back(AM);
	}
	plik.open("czas,Energia,Moment_pedu.data", std::ios::out);
	for (int i = 0; i < qx.size() - 1; i++)
	{
		plik << time.at(i) << " " << Energy.at(i) << " " << AngMom.at(i) << std::endl;
	}
	plik.close();
}

Stellar_Body::~Stellar_Body()
{
	ilosc_ob--;
	plik.close();
	qx.clear();
	qy.clear();
}


///program napisany przez Aleksander Cz¹stkiewicz