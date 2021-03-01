#pragma once

class Matrix
{
public:

	Matrix( int h , int w );
	Matrix( int size  );
	Matrix( const Matrix& M );
	~Matrix();
	
	//operators 
	friend Matrix operator+( const Matrix& M1 , const Matrix& M2 );
	friend Matrix operator-( const Matrix& M1 , const Matrix& M2 );
	friend Matrix operator*( const Matrix& M1 , const Matrix& M2 );
	friend Matrix operator*( double d , const Matrix& M2 );
	friend Matrix operator*( const Matrix& M1 , double d );
	friend Matrix operator^( const Matrix& M1 , int n );

	//operator += *= -= ^=
	//operator == != 
	friend std::ostream& operator<<( std::ostream& os , const Matrix M );

	Matrix& operator=( const Matrix& M );
	double& operator()( int x , int y );
	double* operator[]( int y );

	//methods 
	int getRows()const
	{
		return rows;
	};
	int getCol()const
	{
		return columns;
	}
	double GiveMaxValue();
	double GiveMinValue();
	double GiveDeterminantGaussMethod();//fast
	double GiveDeterminantLaplaceMethod();//slow
	int GiveIndexOfMaxValueInCol( int start );
	double GiveProductOfDiagonal();
	double GiveDeterminantLUMethod();

	Matrix& ConvertToUpperTriangulation();//eliminacja gaussa 
	Matrix& ConvertToLowerTriangulation();
	Matrix& ConvertToMatrixOfAlgebraicComplements();

	Matrix& SwapRows( int r1 , int r2 );
	Matrix& SwapColumns( int c1 , int c2 );

	std::array<Matrix,2> DecompositionLU();
	Matrix& Transpose();
	Matrix& Inverse();

	void CopyMatrixWithoutRowAndCol( const Matrix& M , int row , int col );

	Matrix& InsertNumbers();
	Matrix& InsertNumbersFromFile(FILE* f);
	Matrix& FillMatrixWith0();
	Matrix& FillWithRandomNumbers(int left=0,int right=10);
	Matrix& FillWithNumber( double d );

private:
	double** data;
	int rows;
	int columns;
};