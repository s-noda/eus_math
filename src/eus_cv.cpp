#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/LU>

#include <math.h>
#include <string.h>

#include <opencv2/core/core.hpp>
#include <opencv2/flann/flann.hpp>

int kmeans_log(int k, int sample_n, int sample_length, float* _p, int* _id, float* _center){
  std::cout << "points:" << std::endl;
  for(int i=0 ; i<sample_n ; i++) {
    for ( int j=0; j<sample_length ; j++ ){
      std::cout << " " << _p[j+i*sample_length];
    }
    std::cout << std::endl;
  }
  //
  std::cout << "id:" << std::endl;
  for(int i=0 ; i<sample_n ; i++) {
    std::cout << " " << _id[i];
  }
  std::cout << std::endl;
  //
  std::cout << "centers:" << std::endl;
  for(int i=0 ; i<k ; i++) {
    for ( int j=0; j<sample_length ; j++ ){
      std::cout << " " << _center[j+i*sample_length];
    }
    std::cout << std::endl;
  }
  return 0;
}

extern "C" {
int cv_kmeans(int k, int sample_n, int sample_length,
	      float* _p, int* _id, float* _center,
	      int max_iter, double eps, int debug) {
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
  if ( debug ) kmeans_log(k, sample_n, sample_length, _p,_id,_center);
  //
  return 0;
}
}

// extern "C" {
//   float cv_calc_variance(int n, float* p, float center){
//     float ret = 0;
//     for ( int i = 0; i<n ; i++ ){
//       ret += (p[i]-center)*(p[i]-center);
//       ret /= n;
//       }
//     return sqrt(ret);
//   }
// }

extern "C" {
  int cv_calc_centroid(int n, int d, float* p, float* center){
    for ( int i=0 ; i<d ; i++ ) center[i] = 0;
    for ( int i = 0; i<n ; i++ ){
      for ( int j=0 ; j<d; j++ ) {
	center[j] += p[j+i*d]/n;
      }
    }
    return 0;
  }
}

extern "C" {
  int cv_calc_variance(int n, int d, float* p, float* center, float* var){
    for ( int i=0 ; i<d ; i++ ) var[i] = 0;
    for ( int i = 0; i<n ; i++ ){
      for ( int j=0 ; j<d; j++ ) {
	var[j] += (p[j+i*d]-center[j])*(p[j+i*d]-center[j]);
	var[j] /= n;
      }
    }
    for ( int i=0 ; i<d ; i++ ) {
      var[i] = sqrt(var[i]);
      if ( var[i] < 1e-12 ) var[i] = 1e-12;
    }
    return 0;
  }
}

// extern "C" {
//   float cv_calc_gauss(float p, float center, float var) {
//     return ( 1/(sqrt(2*M_PI)*var)
// 	     * exp( - (p-center)*(p-center) / (2*var*var) ) );
//   }
// }

extern "C" {
  double cv_calc_gauss(int d, float* p, float* center, float* var) {
    double ret = 1;
    for ( int i=0; i<d ; i++ ) {
      ret *= 1/(sqrt(2*M_PI)*var[i])
	* exp( - (p[i]-center[i])*(p[i]-center[i]) / (2*var[i]*var[i]) ) ;
    }
    return ret;
  }
}

extern "C" {
  double cv_calc_likelihood(int n, int d, float* p, float* center, float* var) {
    double ret = 1;
    // cv_calc_variance(n,d,p,center,var);
    for ( int i=0 ; i<n ; i++ ) {
      ret *= cv_calc_gauss(d,(p+i*d),center,var);
    }
    ret = log(ret);
    return ret;
  }
}

extern "C" {
  double cv_calc_bic(int n, int d, float* p, float* center, float* var) {
    double ret = cv_calc_likelihood(n,d,p,center,var);
    ret = -2 * ret + 2*d*log(n);
    return ret;
  }
}

extern "C" {
  int cv_xmeans(int (*stop)(int, float*,int*,float*),
		int sample_n, int sample_length,
		float* _p, int* _id, float* _center, float* var,
		int max_iter, double eps,
		int max_depth, int* depth,
		double bic_stop_thre, double c_dist_thre,
		int debug) {
    int ret = 0;
    if ( depth[0] >= max_depth ) return ret;
    if ( sample_n <= 1 ) return ret;
    // kmeans
    if ( debug ) {
      std::cout << "[xmeans] depth:" << depth[0] << std::endl;
    }
    cv_kmeans(2, sample_n, sample_length, _p, _id, _center, max_iter, eps, 0);
    // swap
    float buf[sample_length];
    int p1=0;
    //for ( int i=0; i<sample_n; i++ ) {
    while ( p1 < sample_n ){
      if ( _id[p1] == 0 ){
	p1++;
      } else {
	bool find = false;
	for ( int j=p1+1; j<sample_n; j++ ) {
	  if ( _id[j] == 0 ) {
	    memcpy(buf, (_p+j*sample_length), sizeof(float)*sample_length);
	    memcpy((_p+j*sample_length), (_p+p1*sample_length),
		   sizeof(float)*sample_length);
	    memcpy((_p+p1*sample_length), buf, sizeof(float)*sample_length);
	    //
	    int tmp = _id[j];
	    _id[j] = _id[p1];
	    _id[p1] = tmp;
	    //
	    find=true;
	    break;
	  }
	}
	if ( ! find ) break;
	p1++;
      }
    }
    for ( int j=0; j<p1; j++ ) _id[j]=depth[0];
    for ( int j=p1; j<sample_n; j++ ) _id[j]=depth[0]+1;
    // log
    if ( debug ) kmeans_log(2, sample_n, sample_length, _p,_id,_center);
    // check condition
    bool revert = false;
    if ( stop ){
      if ( (*stop)(sample_n,_p,_id,_center) ) revert = true;
    } else {
      cv_calc_centroid(sample_n, sample_length, _p, buf);
      cv_calc_variance(sample_n, sample_length, _p, buf, var);
      float bic_org = cv_calc_bic(sample_n, sample_length, _p, buf, var);
      //
      cv_calc_variance(p1, sample_length, _p, _center, var);
      cv_calc_variance((sample_n - p1), sample_length,
		       (_p + p1*sample_length),
		       (_center + sample_length),
		       (var + sample_length));
      float bic1 =cv_calc_likelihood(p1, sample_length, _p, _center, var);
      float bic2 =cv_calc_likelihood((sample_n - p1), sample_length,
				     (_p+p1*sample_length),
				     (_center + sample_length),
				     (var+sample_length));
      float bic = -2 * (bic1+bic2) + (2*sample_length + 2*sample_length)*log(sample_n);
      if ( bic - bic_org > bic_stop_thre ) revert = true;
      //
      float dist = 0;
      for ( int j=0 ; j<sample_length; j++ ) {
	float dif = _center[j+0*sample_length] - _center[j+1*sample_length];
	dist += dif*dif;
      }
      if ( dist < c_dist_thre ) revert = true;
      if ( debug ) {
	std::cout << " stop?:"  << std::endl;
	std::cout << "  -- bic: " << bic_org << " vs " << bic << std::endl;
	std::cout << "  -- dst: " << dist << " vs " << c_dist_thre << std::endl;
	std::cout << "  --  " << revert << std::endl;
      }
    }
    // revert
    if ( revert ) {
      cv_calc_centroid(sample_n, sample_length, _p, _center);
      cv_calc_variance(sample_n, sample_length, _p, _center, var);
      for ( int i=0 ; i<sample_n ; i++ ){
	_id[i] = depth[0];
      }
      return ret;
    }
    // next depth
    ret += cv_xmeans(stop,p1,sample_length,
		     _p,_id,_center,var,
		     max_iter,eps,max_depth,depth,
		     bic_stop_thre, c_dist_thre,
		     debug);
    depth[0]++;
    ret += cv_xmeans(stop,(sample_n - p1),sample_length,
		     (_p+p1*sample_length),
		     (_id+p1),
		     (_center+1*sample_length),
		     (var+1*sample_length),
		     max_iter,eps,max_depth,depth,
		     bic_stop_thre, c_dist_thre,
		     debug);
    return ret;
  }
}

extern "C" {
int cv_hcluster(int sample_n, int sample_length,
		float* _p, float* _center,
		int max_iter, int branch, double eps, int debug) {
  cvflann::KMeansIndexParams k_params(branch, max_iter,
				      cvflann::FLANN_CENTERS_KMEANSPP, eps);
  cv::Mat p(sample_n,sample_length,CV_32F,_p);
  cv::Mat center(sample_n,sample_length,CV_32F,_center);
  int cnt;
  //
  CV_Assert(center.refcount == NULL);
  CV_Assert(p.refcount == NULL);
  //
  cnt = cv::flann::hierarchicalClustering<cv::flann::L2<float> >(p,center,k_params);
  if ( debug ) {
    std::cout << "points:" << std::endl;
    for(int i=0 ; i<sample_n ; i++) {
      for ( int j=0; j<sample_length ; j++ ){
	std::cout << " " << _p[j+i*sample_length];
      }
      std::cout << std::endl;
    }
    //
    std::cout << "id: 0-" << cnt << std::endl;
    //
    std::cout << "centers:" << std::endl;
    for(int i=0 ; i<cnt ; i++) {
      for ( int j=0; j<sample_length ; j++ ){
	std::cout << " " << _center[j+i*sample_length];
      }
      std::cout << std::endl;
    }
  }
  return cnt;
}
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
  cv_kmeans(k, sample_n, sample_length, p, id, center, 10, 0.01, 1);
  //
  return 0;
}

