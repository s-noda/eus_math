#include <iostream>
#include <sys/time.h>

int main() {
  int cnt = 1000000;
  double val = 1.0;
  double rate = 5.3;

  struct timeval s, t;

  std::cout << "multiplication " << cnt << "times" << std::endl;
  gettimeofday(&s, NULL);
  for ( int i=0 ; i<cnt ; i++ ){
    val = (val * rate);
  }
  std::cout << 1 << "-->" << val << std::endl;
  gettimeofday(&t, NULL);
  std::cout <<  ((t.tv_sec - s.tv_sec) * 1000 + (t.tv_usec - s.tv_usec) / 1000) << "ms" << std::endl;

  std::cout << "addition " << cnt << "times" << std::endl;
  val = 1.0;
  gettimeofday(&s, NULL);
  for ( int i=0 ; i<cnt ; i++ ){
    val = (val + rate);
  }
  std::cout << 1 << "-->" << val << std::endl;
  gettimeofday(&t, NULL);
  std::cout <<  ((t.tv_sec - s.tv_sec) * 1000 + (t.tv_usec - s.tv_usec) / 1000) << "ms" << std::endl;
}
