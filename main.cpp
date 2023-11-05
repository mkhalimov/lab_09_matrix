#include "matrix.hpp"
#include <iostream>

using std::cout;
using std::cin;
using std::endl;

int main()
{
    using my_type = double;
    Matrix<my_type> A;
    cin >> A;

    cout << "We got matrix" << endl;
    
    Matrix<my_type> M(A);
    M = M * 0;
    Matrix<my_type> L(M);
    Matrix<my_type> U(M);
    A.LU(&U, &L, &M);
    Matrix<my_type> temp = A.M();
    cout << "M:\n" << temp << endl 
         << "L:\n"     <<       L << endl
         << "U:\n"     <<       U << endl 
         << "M*A:\n"   <<     M*A << endl 
         << "M*L*U:\n" << M*(L*U) << endl 
         << "A:\n"     <<       A << endl;
    Matrix<my_type> b (2, 1);
    b(0,0) = 1;
    b(1,0) = 4;
    cout << "b is :\n" << b << endl;
    
    Matrix<my_type> answer(A.Solve(b));
    cout << "Ax = b..." << endl;
    cout << "x is :\n" << answer << endl;
    
    return 0;
}