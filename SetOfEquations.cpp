#include <vector>
#include <math.h>
#include <string>

#include "exceptions.h"
#include "SetOfEquations.h"

//private functions for this file 
bool checkCriterion(Matrix& M);
bool isErrorGreaterThan(Matrix& x1, Matrix& x2,double error);


//operators 
std::ostream& operator<<(std::ostream& os, SetOfEquations& set)
{
	for (int i = 0; i < set.A.getRows(); i++)
	{
		for (int j = 0; j < set.A.getCol(); j++)
		{
			os << " x" << j << "*" << set.A[i][j] << " ";
			if (j != set.A.getRows()- 1)
				os << "+";
		}
		os << " = " << set.B[i][0] << "\n";
	}

	if (set.solution)
	{
		os << "\n";
		for (int i = 0; i < set.A.getRows(); i++)
			os << " x" << i << " = " << set.X[i][0]<<"\n";
	}

	return os;
}

//methods
void SetOfEquations::LU()
{

}

void SetOfEquations::Cramer()
{
	if (A.getRows() != A.getCol())
	{
		throw exception_setofequations(" Cant find solution by Cramer's method. Matrix isnt square.\n");
		return;
	}

	int size = A.getRows();
	double W = A.GiveDeterminantGaussMethod();
	if (!W)
	{
		throw exception_setofequations(" Cant find solution by Cramer's method. Determinant is 0.\n");
		return;
	}

	for (int i = 0; i < size; i++)
	{
		Matrix Wi = A;
		for (int j = 0; j < size; j++)
			Wi[j][i] = B[j][0];

		X[i][0] = Wi.GiveDeterminantGaussMethod() / W;
	}

	solution = true;
}

void SetOfEquations::GaussElimination()
{
	int size = A.getRows();

	//create extended matrix M = A|B -> convert to upper trianguklation -> alg backwards =
	Matrix M(size, size + 1);
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			M[i][j] = A[i][j];
	for (int i = 0; i < size; i++)
		M[i][size] = B[i][0];

	M.ConvertToUpperTriangulation();

	X[size - 1][0] = M[size - 1][size] / M[size - 1][size - 1];
	for (int i = size - 2; i >= 0; i--)
	{
		double E = 0;
		for (int s = i + 1; s < size; s++)
			E += M[i][s] * X[s][0];
		X[i][0] = (M[i][size] - E) / M[i][i];
	}

	solution = true;
}

void SetOfEquations::SimpleIteration(double error)
{
	//divide matrix A to D and R 
	int size = A.getRows();
	Matrix D(size);
	Matrix R(size);

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			if (i == j)
			{
				D[i][j] = -1 / A[i][j];
			}
			else
			{
				R[i][j] = A[i][j];
			}

	//W = -1/D*R
	Matrix W = D * R;

	for (int i = 0; i < size; i++)
				D[i][i] = -D[i][i];

	if (checkCriterion(W))
	{
		//Z=1/D*B
		Matrix Z = D * B;//vector 

		for (int i = 0; i < size; i++)
			D[i][i] = 1 / D[i][i];

		X.FillWithRandomNumbers(-2,2);

		int iter = 0;
		Matrix x_tmp = X;
		while (1)
		{
			//iterate 
			X = W * x_tmp + Z;

			if (!isErrorGreaterThan(X, x_tmp,error))
				break;

			x_tmp = X;
			iter++;
		}

		solution = true;
		return;
	}
	throw exception_setofequations(" Cant find solution by simple iteration method.\n");
}

void SetOfEquations::GaussaSeidelaIteration(double error)
{
	//divide matrix A to D,L(ower) and U(pper)  
	int size = A.getRows();
	Matrix D(size);
	Matrix L(size);
	Matrix U(size);

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			if (i == j)
			{
				D[i][j] = -1 / A[i][j];
			}
			else if (i > j)
			{
				L[i][j] = A[i][j];
			}
			else
			{
				U[i][j] = A[i][j];
			}

	Matrix Wl = D * L;
	Matrix Wu = D * U;

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			if (i == j)
				D[i][j] = -D[i][j];

	if (checkCriterion(Wl))
	{
		Matrix Z = D * B;
	
		X.FillWithRandomNumbers(-2, 2);

		Matrix X_next(size, 1);

		int iter = 1;
		while (1)
		{
			//iterate
			//X_next = Wl * X_next + Wu * X + Z;//idk 
			for (int i = 0; i < size; i++)
			{
				double s = 0;
				for (int k = 0; k < size; k++)
					s += Wl[i][k] * X_next[k][0] + Wu[i][k] * X[k][0];
				X_next[i][0] = s + Z[i][0];
			}

			/*std::cout << "iteration nr " << iter << "\n";
			for (int i = 0; i < size; i++)
				std::cout << "x" << i << " = " << X[i][0] << "\n";
			std::cout << std::endl;
			*/

			iter++;

			if (!isErrorGreaterThan(X_next, X,error))
			{
				X = X_next;
				break;
			}
			X = X_next;
		}

		solution = true;
	}
	throw exception_setofequations(" Cant find solution by Gauss-Sidel iteration method.\n");
}

void SetOfEquations::InsertSetFromFile(FILE* f)
{
	for (int i = 0; i < A.getRows(); i++)
	{
		for (int j = 0; j < A.getCol(); j++)
			if (!fscanf(f, "%lg", &A[i][j]))	return;
		if (!fscanf(f, "%lg", B[i]))	return;//if fscanf fails return 
	}
}

void SetOfEquations::PrintLayout()
{
	for (int i = 0; i < A.getRows(); i++)
	{
		for (int j = 0; j < A.getCol(); j++)
		{
			printf(" x%d*%g", j, A[i][j]);
			if (j != A.getRows() - 1)
				printf(" + ");
		}

		printf(" = %g\n", B[i][0]);
	}

	printf("\n");
	if (solution)
	{
		for (int i = 0; i < A.getRows(); i++)
			printf(" x%d = %f\n", i, X[i][0]);
		printf("\n");
	}
}

void SetOfEquations::PrintSolution()
{
	if (solution)
	{
		for (int i = 0; i < A.getRows(); i++)
			printf(" x%d = %f\n", i, X[i][0]);
	}
	else
		printf(" Solution wasnt found\n");
	printf("\n");
}


//private funcions for this file 
bool checkCriterion(Matrix& M)
{
	double w = 0;
	//first criterion, check max in row
	for (int i = 0; i < M.getRows(); i++)
	{
		double tmp = 0;
		for (int j = 0; j < M.getCol(); j++)
			tmp += fabs(M[i][j]);
		if (tmp > w)
			w = tmp;
	}
	if (w < 1)
		return true;

	//next criterion, check max in col
	w = 0;
	for (int j = 0; j < M.getCol(); j++)
	{
		double tmp = 0;
		for (int i = 0; i < M.getRows(); i++)
			tmp += fabs(M[j][i]);
		if (tmp > w)
			w = tmp;
	}
	if (w < 1)
		return true;

	//next criterion, check sqrt(x*x+...)
	w = 0;
	for (int i = 0; i < M.getRows(); i++)
	{
		for (int j = 0; j < M.getCol(); j++)
			w += M[i][j] * M[i][j];
	}
	if (sqrt(w) < 1)
		return true;

	return false;
}

bool isErrorGreaterThan(Matrix& x1, Matrix& x2, double error)
{
	for (int i = 0; i < x1.getRows(); i++)
		if (fabs(x1[i][0] - x2[i][0]) > error)
			return true;
	return false;
}
