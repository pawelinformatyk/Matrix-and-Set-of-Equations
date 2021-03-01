#include <iostream>
#include <array>
#include <time.h>
#include <iomanip>

#include "Matrix.h"

//konstruktory i destr
Matrix::Matrix( int h , int w )
{
	rows = h;
	columns = w;

	data = new double* [rows];
	for(int i = 0; i < rows; i++)
	{
		data[i] = new double[columns];
		for(int j = 0; j < columns; j++)
			data[i][j] = 0;
	}
}

Matrix::Matrix( int size )
{
	rows = size;
	columns = size;

	data = new double* [rows];
	for(int i = 0; i < rows; i++)
	{
		data[i] = new double[columns];
		for(int j = 0; j < columns; j++)
			data[i][j] = 0;
	}
}

Matrix::Matrix( const Matrix& M )
{
	columns = M.columns;
	rows = M.rows;

	data = new double* [rows];
	for(int i = 0; i < rows; i++)
	{
		data[i] = new double[columns];
		for(int j = 0; j < columns; j++)
			data[i][j] = M.data[i][j];
	}
}

Matrix::~Matrix()
{
	for(int i = 0; i < rows; i++)
		delete data[i] ;
	delete[]data;
}


//operatory 
Matrix operator+( const Matrix& M1 , const Matrix& M2 )
{
	if(M1.rows != M2.rows || M1.columns != M2.columns)
		printf( "aa" );//wyjatek

	Matrix tmp( M2.rows , M2.columns );

	for(int i = 0; i < M1.rows; i++)
		for(int j = 0; j < M1.columns; j++)
			tmp.data[i][j] = M1.data[i][j] + M2.data[i][j];

	return tmp;
}

Matrix operator-( const Matrix& M1 , const Matrix& M2 )
{
	if(M1.rows != M2.rows || M1.columns != M2.columns)
		printf( "aa" );//wyjatek

	Matrix tmp( M2.rows , M2.columns );

	for(int i = 0; i < M1.rows; i++)
		for(int j = 0; j < M1.columns; j++)
			tmp.data[i][j] = M1.data[i][j] - M2.data[i][j];

	return tmp;
}

Matrix operator*( const Matrix& M1 , const Matrix& M2 )
{
	if(M1.columns != M2.rows)
		printf( "aa" );//wyjatek

	Matrix tmp( M1.rows , M2.columns );

	for(int i = 0; i < M1.rows; i++)
		for(int j = 0; j < M2.columns; j++)
		{
			double s = 0;
			for(int k = 0; k < M1.columns; k++)
				s += M1.data[i][k] * M2.data[k][j];
			tmp.data[i][j] = s;
		}

	return tmp;
}

Matrix operator*( double d , const Matrix& M )
{
	Matrix tmp( M.rows , M.columns );

	for(int i = 0; i < M.rows; i++)
		for(int j = 0; j < M.columns; j++)
			tmp.data[i][j] = d * M.data[i][j];

	return tmp;
}

Matrix operator*( const Matrix& M , double d )
{
	Matrix tmp( M.rows , M.columns );

	for(int i = 0; i < M.rows; i++)
		for(int j = 0; j < M.columns; j++)
			tmp.data[i][j] = d * M.data[i][j];

	return tmp;

}

std::ostream& operator<<( std::ostream& os , const Matrix M )
{
	for(int i = 0; i < M.rows; i++)
	{
		for(int j = 0; j < M.columns; j++)
			os << " " << std::setprecision(4) << M.data[i][j]  ;
		os << std::endl;
	}

	return os;
}

Matrix& Matrix::operator=( const Matrix& M )
{
	if(this == &M)// pozwala na M=M
		return *this;

	for(int i = 0; i < rows; i++)
		delete data[i] ;
	delete[]data;

	rows = M.rows;
	columns = M.columns;

	data = new double* [rows];
	for(int i = 0; i < rows; i++)
		data[i] = new double[columns];

	for(int i = 0; i < M.rows; i++)
		for(int j = 0; j < M.columns; j++)
			data[i][j] = M.data[i][j];

	return *this;
}

double& Matrix::operator()( int x , int y )
{
	if(x >= columns || y >= rows)
		return data[0][0];

	return data[y][x];
}

double* Matrix::operator[]( int y )
{
	if(y >= rows)
		return NULL;

	return data[y];
}

Matrix operator^( const Matrix& M , int n )
{
	if(M.rows != M.columns)
		return M;

	if(n == 0)
	{
		Matrix tmp( M.rows , M.columns );

		for(int i = 0; i < M.columns; i++)
			tmp[i][i] = 1;

		return tmp;
	}
	else if(n == 1)
		return M;

	Matrix tmp = M;

	for(int i = 1; i < n/2 ; i++)
		tmp = tmp * M;

	tmp = tmp * tmp;

	return (n < 0) ? tmp.Inverse() : tmp;
}


//metody
Matrix& Matrix::ConvertToUpperTriangulation()
{
	for(int s = 0; s < rows - 1; s++)
		for(int i = s + 1; i < rows; i++)
		{
			if(fabs( data[s][s] ) <= 1e-12)
			{									  //h raize 0 na diagonalnej
				int r = GiveIndexOfMaxValueInCol( s );		  //znajduje max wartosc h kol
				SwapRows( s , r );				  //przestawiam s wiersz z tym
			}
			double m = data[i][s] / data[s][s];
			for(int j = s; j < columns; j++)
				data[i][j] -= data[s][j] * m;
		}

	return *this;
}

Matrix& Matrix::ConvertToMatrixOfAlgebraicComplements()
{
	if(rows != columns)
		return *this;

	Matrix tmp( rows - 1 );
	Matrix copy = *this;

	for(int i = 0; i < rows; i++)
	{
		int sign = (i % 2) ? -1 : 1;

		for(int j = 0; j < columns; j++)
		{
			tmp.CopyMatrixWithoutRowAndCol( copy , i , j );
			data[i][j] = sign * tmp.GiveDeterminantGaussMethod(); 
			sign = -sign;
		}
	}

	return *this;
}

double Matrix::GiveProductOfDiagonal()
{
	double d = 1;
	for(int i = 0; i < rows; i++)
		d *= data[i][i];

	return d;
}

double Matrix::GiveDeterminantGaussMethod()
{
	if(rows != columns)
		return 0;

	if(rows == 1)
		return data[0][0];
	else if(rows == 2)
		return data[0][0] * data[1][1] - data[0][1] * data[1][0];

	Matrix M = *this;
	int swap = 1;//if swapped h
	for(int s = 0; s < rows - 1; s++)
		for(int i = s + 1; i < rows; i++)
		{
			if(fabs( M[s][s] ) < 1e-12)
			{									  //h raize 0 na diagonalnej
				int r = M.GiveIndexOfMaxValueInCol( s );		  //znajduje max wartosc h kol
				M.SwapRows( s , r );				  //przestawiam s wiersz z tym
				swap = ++swap % 2;
			}
			double m = M[i][s] / M[s][s];
			for(int j = s; j < rows; j++)
				M[i][j] -= M[s][j] * m;
		}

	double d = (swap) ? 1 : -1;
	for(int i = 0; i < rows; i++)
		d *= M[i][i];

	return d;
}

double Matrix::GiveDeterminantLaplaceMethod()
{
	if(rows != columns)
		return 0;

	if(rows == 1)
		return data[0][0];
	else if(rows == 2)
		return data[0][0] * data[1][1] - data[0][1] * data[1][0];

	//0th h

	double d = 0;
	int sign = 1;

	Matrix tmp( rows - 1 );

	for(int i = 0; i < columns; i++ )
	{
		tmp.CopyMatrixWithoutRowAndCol( *this , 0 , i );
		d += sign * data[0][i] * tmp.GiveDeterminantLaplaceMethod();  // v+= (-1)^(i+j) * a(i,j) * det(minorAij)
		sign = -sign;
	}

	return d;
}

double Matrix::GiveDeterminantLUMethod()
{
	if(rows != columns)
		return 0;

	if(rows == 1)
		return data[0][0];
	else if(rows == 2)
		return data[0][0] * data[1][1] - data[0][1] * data[1][0];

	//double d = LU[0].GiveProductOfDiagonal() * LU[1].GiveProductOfDiagonal();

	return 0;
}			

int Matrix::GiveIndexOfMaxValueInCol( int start )
{
	int r = start + 1;
	for(int i = start + 2; i < rows; i++)
		if(data[r][start] < data[i][start])
			r = i;

	return r;
}

double Matrix::GiveMaxValue()
{
	return 0;
}

double Matrix::GiveMinValue()
{
	return 0;
}

Matrix& Matrix::SwapRows( int r1 , int r2 )
{
	if(r1 >= rows || r2 >= rows)
		return *this;//wyjatek

	double* tmp = data[r1];
	data[r1] = data[r2];
	data[r2] = tmp;

	return *this;
}

Matrix& Matrix::SwapColumns( int c1 , int c2 )
{
	if(c1 >= columns|| c2 >= columns)
		return *this;//wyjatek

	for(int i = 0; i < rows; i++)
	{
		double tmp = data[i][c1];
		data[i][c1] = data[i][c2];
		data[i][c2] = tmp;
	}

	return *this;
}

std::array<Matrix, 2> Matrix::DecompositionLU()
{
	//initz array 
	std::array<Matrix, 2> array = { Matrix(this->rows),Matrix(this->rows) };
	
	if(rows != columns)
	{
		std::cout << "rows != colums\n";
		return array;
	}

	//create matrix LU and insert values
	Matrix LU = *this;
	for(int k = 0; k < rows - 1; k++)
	{
		if(fabs( LU[k][k] ) < 1e-12)
			return array;

		for(int i = k + 1; i < rows; i++)
			LU[i][k] /= LU[k][k];

		for(int i = k + 1; i < rows; i++)
			for(int j = k + 1; j < rows; j++)
				LU[i][j] -= LU[i][k] * LU[k][j];
	}

	//insert values to array 
	for(int i = 0; i < rows; i++)
		for(int j = 0; j <= i; j++)
			(i != j) ? (array[0])[i][j] = LU[i][j] : (array[0])[i][j] = 1;

	for(int i = 0; i < rows; i++)
		for(int j = i; j < columns; j++)
			(array[1])[i][j] = LU[i][j];

	return array;
}

Matrix& Matrix::Transpose()
{
	Matrix M(columns,rows);

	for(int i = 0; i < rows; i++)
		for(int j = 0; j < columns; j++)
			M[j][i] = data[i][j];

	*this = M;

	return *this;
}

Matrix& Matrix::Inverse()
{
	if(rows != columns)
		return *this;

	double det = this->GiveDeterminantGaussMethod();
	if(fabs(det) < 1e-12)
	{
		std::cout << "Det equals 0 , no inverse \n";
		return *this;
	}

	det = 1;

	this->ConvertToMatrixOfAlgebraicComplements( );
	this->Transpose();

	for(int i = 0; i < rows; i++ )
	{
		for(int j = 0; j < columns; j++ )
			data[i][j] /= det; 
	}

	return *this;
}

void Matrix::CopyMatrixWithoutRowAndCol( const Matrix& M , int row , int col )
{
	//wycinanie wiersza i kolumny
	if(M.rows - 1 != this->rows || M.columns != this->columns || row >= M.rows || col >= M.columns)
		return;

	int k = 0;
	for(int i = 0; i < M.rows; i++)
	{
		if(i == row)  continue;

		int l = 0;
		for(int j = 0; j < M.columns; j++)
		{
			if(j == col)	continue;

			data[k][l] = M.data[i][j];
			l++;
		}
		k++;
	}
}

Matrix& Matrix::InsertNumbers()
{
	printf( "Insert data to %dx%d matrix\n" , rows , columns );
	for(int i = 0; i < rows; i++)
		for(int j = 0; j < columns; j++)
			std::cin >> data[i][j];

	return *this;
}

Matrix& Matrix::InsertNumbersFromFile(FILE* f)
{
	if (!f)
		return *this;

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			fscanf(f, "%lf ", &data[i][j]);
		}
		
		fscanf(f, "\n");
	}
}

Matrix& Matrix::FillMatrixWith0()
{
	for(int i = 0; i < rows; i++)
		for(int j = 0; j < columns; j++)
			data[i][j] = 0;

	return *this;
}

Matrix& Matrix::FillWithNumber(double d)
{
	for(int i = 0; i < rows; i++)
		for(int j = 0; j < columns; j++)
			data[i][j] = d;

	return *this;
}

Matrix& Matrix::FillWithRandomNumbers(int left, int right)
{
	srand( (unsigned int)time( NULL ) );

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < columns; j++)
			data[i][j] = rand() % (right - left + 1) + left;

	return *this;
}

