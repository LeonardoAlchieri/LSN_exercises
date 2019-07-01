//
//
//  **** ESERCIZIO 1.2 ****
//
//  Corso di Simulazione Numerica - AA 2018/2019
//
//  Leonardo Alchieri, 886810
//
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132

using namespace std;

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

//
//
//  Method to create exponentially random variables in [0,∞).
//
//
double Random :: Exp(double rate){
    double x;       //exponentially generated random
    double y = Rannyu();    //random variable from uniform distribution [0,1)
    x = (-1/rate) * log(1-y);   //invers of CDF for y = λe^(-λx)
    
    return x;
}

//
//
//  Method to create random variables in a Cauchy-Lorentz distribution.
//
//
double Random :: Lorentz(double mean, double width){
    double x=11;      //random generated in Cauchy-Lorentz distribution
    double y;
    while(x>10 || x<-10){
        y = Rannyu();    //random variable from uniform distribution [0,1)
        x = mean + width * tan(PI * (y - 1/2));
    }
    return x;
}
