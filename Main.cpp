#include<iostream>
#include<cmath>
const int QUANTITY_Y = 6;
const int QUANTITY_X = 3;


double MatrixElem(unsigned int, unsigned int, double[QUANTITY_Y][QUANTITY_X]);
double CalcB(unsigned short, double[QUANTITY_Y], double[QUANTITY_Y][QUANTITY_X]);
double* GaussMethod(double[3][4]);
void AproksimFunction();

int main()
{
	AproksimFunction();
	system("pause");
}

double MatrixElem(unsigned int n,unsigned int m, double BazF[QUANTITY_Y][QUANTITY_X])
{
	setlocale(LC_ALL, "russian");
	double result = 0;
	for (unsigned short i(0); i < QUANTITY_Y; i++)
	{
		result += (BazF[i][n - 1] * BazF[i][m - 1]);
	}
	return result;
}
double CalcB(unsigned short n, double Y[QUANTITY_Y], double BazF[QUANTITY_Y][QUANTITY_X])
{
	double result = 0;
	for (unsigned short i(0); i < QUANTITY_Y; i++)
	{
		result += (Y[i] * BazF[i][n-1]);
	}
	return result;
}


double* GaussMethod(double Matrix[3][4])
{
	double* result = new double[3];
	for (unsigned short j(1); j < 3; j++)
	{
		double subdiv = Matrix[j][0] / Matrix[0][0];
		for (unsigned short i(0); i < 4; i++)
		{
			Matrix[j][i] = Matrix[j][i] - (Matrix[0][i] * subdiv);
		}
	}
	for (unsigned short j(2); j < 3; j++)
	{
		double subdiv = Matrix[j][1] / Matrix[1][1];
		for (unsigned short i(1); i < 4; i++)
		{
			Matrix[j][i] = Matrix[j][i] - (Matrix[1][i] * subdiv);
		}
	}

	result[2] = Matrix[2][3] / Matrix[2][2];
	result[1] = (Matrix[1][3] - Matrix[1][2] * result[2]) / Matrix[1][1];
	result[0] = (Matrix[0][3] - Matrix[0][1] * result[1] - Matrix[0][2] * result[2]) / Matrix[0][0];
	return result;
}

void AproksimFunction()
{
	double Y[QUANTITY_Y] = { 0.1f, 0.4f, 0.9f, 1.6f, 2.5f, 3.6f };
	double BazF[QUANTITY_Y][QUANTITY_X] = {
			{1, 1.1f, 0.33287108369f},
			{1, 2.1f, 0.12245642825f},
			{1, 3.1f, 0.04504920239f},
			{1, 4.1f, 0.0165726754f},
			{1, 5.1f, 0.00609674656f},
			{1, 6.1f, 0.00224286771f},
	};
	double Matrix[3][4];
	for (unsigned short i(0); i < 3; i++)
	{
		for (unsigned short j(0); j < 3; j++)
		{
			Matrix[i][j] = MatrixElem(i+1, j+1, BazF);
		}
		Matrix[i][3] = CalcB(i+1, Y, BazF);
	}

	//Вывод изначальной матрицы
	std::cout << "Изначальная матрица\n";
	for (unsigned short i(0); i < 3; i++)
	{
		for (unsigned short j(0); j < 4; j++)
		{
			std::cout << Matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
	double* answers = GaussMethod(Matrix);
	double f[QUANTITY_Y] = {};
	double sigma[QUANTITY_Y] = {};
	double J = 0;
	for (unsigned short i(0); i < QUANTITY_Y; i++)
	{
		f[i] = *answers * BazF[i][0] + *(answers + 1) * BazF[i][1] + *(answers + 2) * BazF[i][2];
		sigma[i] = abs(f[i] - Y[i]);
		J += pow(sigma[i], 2);
	}

	//Отсюда идет вывод данных
	std::cout << "\nМатрица после обработки\n";
	for (unsigned short i(0); i < 3; i++)
	{
		for (unsigned short j(0); j < 4; j++)
		{
			std::cout << Matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "\nX,Y и Z в результате выполнения метода Гаусса\n";
	for (unsigned short i(0); i < 3; i++)
	{
		std::cout << *(answers + i) << " ";
	}
	std::cout << "\n";
	std::cout << "\nЗначения апрокс. функций\n";
	for (unsigned short i(0); i < QUANTITY_Y; i++)
	{
		std::cout << f[i] << " ";
	}
	std::cout << "\n";
	std::cout << "\nЗначения сигмы\n";
	for (unsigned short i(0); i < QUANTITY_Y; i++)
	{
		std::cout << sigma[i] << " ";
	}
	std::cout << "\n";
	std::cout << "\nJ = " << J << "\n";
	delete[] answers;
}
