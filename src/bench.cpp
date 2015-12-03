#include <iostream>
#include <sys/time.h>

int main() {
  int cnt = 100000;
  double val = 1.0;
  double rate = 5.3;
  int i;

  struct timeval s, t;

  std::cout << "multiplication " << cnt << "times" << std::endl;
  val = 1.0;
  i = 0;
  gettimeofday(&s, NULL);
  // for ( int i=0 ; i<cnt ; i++ ){
  while ( i++ < cnt ){
    val = (val * rate);
  }
  gettimeofday(&t, NULL);
  std::cout << 1 << "-->" << val << std::endl;
  std::cout <<  ((t.tv_sec - s.tv_sec) * 1000.0 + (t.tv_usec - s.tv_usec) / 1000.0) << "ms" << std::endl;

  std::cout << "addition " << cnt << "times" << std::endl;
  val = 1.0;
  i = 0;
  gettimeofday(&s, NULL);
  while ( i++ < cnt ){
  // for ( int i=0 ; i<cnt ; i++ ){
    val = (val + rate);
  }
  gettimeofday(&t, NULL);
  std::cout << 1 << "-->" << val << std::endl;
  std::cout <<  ((t.tv_sec - s.tv_sec) * 1000.0 + (t.tv_usec - s.tv_usec) / 1000.0) << "ms" << std::endl;
}
