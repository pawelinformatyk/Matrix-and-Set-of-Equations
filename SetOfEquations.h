#pragma once

#include <iostream>

#include "Matrix.h"

//to do 
	//aplication to manage set of equations 
	//LU 

class SetOfEquations
{
	//A*x=B 
		//n1*x1+n2*x2+...nn*xn=b1
		//n1*x1+n2*x2+...nn*xn=b2
		//.......................
		//.......................
		//.......................
		//n1*x1+n2*x2+...nn*xn=bn

	bool solution;

public:
	Matrix A;
	Matrix X;//solution / vector x
	Matrix B;//vector 


	SetOfEquations(const Matrix a, const Matrix b) :A(a), B(b), X(a.getRows(), 1)
	{
		//if size is dif 
		solution = false;
		if (A.getRows() != B.getRows())
		{
			//throw exception_setofequations("Matrix A has more rows than vector B\n");
			exit( 0 );
		}

	}

	//operators
	friend std::ostream& operator<<(std::ostream& os, SetOfEquations& set);

	//methods 
	void PrintLayout();
	void PrintSolution();
	void InsertSetFromFile(FILE*);
	//bool SolutionFound();

	//find solution to system of equations 
	void Cramer();
	void LU();
	void GaussElimination();
	void SimpleIteration(double error=0.1);
	void GaussaSeidelaIteration(double error=0.1);
};