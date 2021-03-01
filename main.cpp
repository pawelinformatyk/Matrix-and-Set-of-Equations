#include <Windows.h>
#include <iostream>
#include <time.h>
#include <array>
#include <vector>

#include "Matrix.h"
#include "SetOfEquations.h"


int main()
{
	Matrix m1( 5 );
	m1.FillWithRandomNumbers();
	Matrix m2( 5,1 );
	m2.FillWithRandomNumbers();

	std::cout << m1 << "\n\n" << m2 << "\n\n" << m1 * m2 << "\n\n";

	SetOfEquations set( m1, m2 );
	try
	{
		set.GaussaSeidelaIteration();

	}
	catch( std::exception& e)
	{
		std::cout << e.what() << std::endl;
		return 1;
	}

	set.PrintLayout();


	return 0;
}

