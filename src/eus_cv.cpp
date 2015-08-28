#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/LU>

#include <opencv2/core/core.hpp>
#include <opencv2/flann/flann.hpp>

// void copy_array_to_cvmat(cv::Mat &mat, double *ar) {
// }

int cv_kmeans(int k, int sample_n, int sample_length,
	      float* _p, int* _id, float* _center,
	      int max_iter, double eps) {
  cv::Mat p(sample_n,sample_length,CV_32F,_p);
  cv::Mat id(sample_n,1,CV_32S,_id);
  cv::Mat center(k,sample_length,CV_32F,_center);
  //
  CV_Assert(center.refcount == NULL);
  CV_Assert(id.refcount == NULL);
  CV_Assert(p.refcount == NULL);
  //
  cv::kmeans(p, k, id, cvTermCriteria(CV_TERMCRIT_EPS|CV_TERMCRIT_ITER, max_iter, eps),
	     1, cv::KMEANS_PP_CENTERS, center);
  //
  return 0;
}

int main() {
  int k=2;
  int sample_n = 10;
  int sample_length = 3;
  float p[sample_n*sample_length];
  int id[sample_n];
  float center[k*sample_length];
  //
  for(int i=0 ; i<sample_n ; i++) {
    for ( int j=0; j<sample_length ; j++ ){
      Eigen::Vector2f pt = Eigen::Vector2f::Random();
      p[j+i*sample_length] = pt.coeff(0);
      if ( i > sample_n/2 ) {
	p[j+i*sample_length] += 100;
      }
    }
  }
  //
  cv_kmeans(k, sample_n, sample_length, p, id, center, 10, 0.01);
  //
  std::cout << "points:" << std::endl;
  for(int i=0 ; i<sample_n ; i++) {
    for ( int j=0; j<sample_length ; j++ ){
      std::cout << " " << p[j+i*sample_length];
    }
    std::cout << std::endl;
  }
  //
  std::cout << "id:" << std::endl;
  for(int i=0 ; i<sample_n ; i++) {
    std::cout << " " << id[i];
  }
  std::cout << std::endl;
  //
  std::cout << "centers:" << std::endl;
  for(int i=0 ; i<k ; i++) {
    for ( int j=0; j<sample_length ; j++ ){
      std::cout << " " << center[j+i*sample_length];
    }
    std::cout << std::endl;
  }
  //
  return 0;
}

