#include <stdio.h>
#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/LU>

extern "C"{
  void calc_eigen(int n, double* mat, double* peigenval, double* neigenval, double* peigenvector, double* neigenvector){
    Eigen::MatrixXd A = Eigen::Map<Eigen::MatrixXd>(mat,n,n);
    Eigen::EigenSolver<Eigen::MatrixXd> s(A);
    for ( int i=0; i<n ; i++ ){
      std::complex<double> val = s.eigenvalues()(i);
      peigenval[i] = (double)val.real() ;
      neigenval[i] = (double)val.imag() ;
      for ( int j=0; j<n ; j++ ){
	peigenvector[i*n+j] = s.eigenvectors().col(i)(j).real() ;
	neigenvector[i*n+j] = s.eigenvectors().col(i)(j).imag() ;
      }
    }
  }
}

extern "C"{
  void calc_inverse_matrix(int n, double* mat, double* ret){
    Eigen::MatrixXd A = Eigen::Map<Eigen::MatrixXd>(mat,n,n);
    Eigen::Map<Eigen::MatrixXd>(ret,n,n) = A.inverse();
  }
}

extern "C"{
  double calc_determinant(int n, double* mat){
    Eigen::MatrixXd A = Eigen::Map<Eigen::MatrixXd>(mat,n,n);
    double ret ;
    ret = A.determinant();
    return ret;
  }
}

extern "C" {
  int debug_dump_float(float* v, int n){
    std::cout << v[0];
    for ( int i=1 ; i<n ; i++ ){
      std::cout << " " << v[i];
      v[i]++;
    }
    std::cout << std::endl;
    return 0;
  }
}

extern "C" {
  int debug_dump_double(double* v, int n){
    std::cout << v[0];
    for ( int i=1 ; i<n ; i++ ){
      std::cout << " " << v[i];
      v[i]++;
    }
    std::cout << std::endl;
    return 0;
  }
}

extern "C" {
  int debug_dump_int(int* v, int n){
    std::cout << v[0];
    for ( int i=1 ; i<n ; i++ ){
      std::cout << " " << v[i];
      v[i]++;
    }
    std::cout << std::endl;
    return 0;
  }
}

extern "C" {
  int print_args(double* ds, float f, double d, int i, long l, long* ls) {
    std::cout << __func__ << std::endl;
    printf("  -- double string=%2.16f\n", ds[0]);
    printf("  -- float        =%2.16f\n", f);
    printf("  -- double       =%2.16f\n", d);
    printf("  -- integer      =%2.16f\n", i);
    printf("  -- long         =%2.16f\n", l);
    printf("  -- long string  =%2.16f\n", ls[0]);
    ds[0] = 1.11111111111111111111111111111111;
    return 0;
  }
}
