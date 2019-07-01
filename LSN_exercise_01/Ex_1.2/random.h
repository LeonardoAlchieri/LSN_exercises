//
//
//  **** ESERCIZIO 1.2 ****
//
//  Corso di Simulazione Numerica - AA 2018/2019
//
//  Leonardo Alchieri, 886810
//
//

#ifndef __Random__
#define __Random__

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
    // constructors
    Random();
    // destructor
    ~Random();
    // methods
    //
    //  The first 5 methods have been provided.
    //
    void SetRandom(int * , int, int);
    void SaveSeed();
    double Rannyu(void);
    double Rannyu(double min, double max);
    double Gauss(double mean, double sigma);
    //
    //
    //  The following 2 methods have been created as part of the exercise.
    //
    //
    double Exp(double rate);
    double Lorentz(double mean, double width);
};

#endif // __Random__


