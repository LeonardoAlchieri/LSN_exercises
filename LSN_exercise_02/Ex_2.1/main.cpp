//
//  Exercise 2.1
//
//  Leoanrdo Alchieri, 886810
//
//  Laboratorio di Simulazione Numerica
//
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

double error(double*, double*, int);        //function to calculate errors
 
int main (int argc, char *argv[]){

    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
      Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    
    rnd.SaveSeed();
    
    //
    //
    //  PARTE 1
    //
    //
    
    int M=1E4;      //number of random numbers to generate in each blocks
    int N=1E2;      //number of blocks on which I average
    
    //
    //
    //  I evaluate the integral ∫ π/2 cos(x π/2) dx
    //
    //
    double f=0;
    double x=0;
    
    double avg[N];      //array with the N averages for the integral
    for(int i=0; i<N; i++)
        avg[i] = 0;
    
    double sq_avg[N];      //as above, but with the squares of the averages
    for(int i=0; i<N; i++)
        sq_avg[i] = 0;
        
    for(int j=0; j<N; j++){
        for(int i=0; i<M; i++){
            x = rnd.Rannyu();
            f += cos(x*M_PI/2);
        }
        avg[j] = M_PI/2. * f / M;
        sq_avg[j] = avg[j] * avg[j];
        f=0;
    }
    
    ofstream out;
    out.open("part_1.txt");
    
    double progressive_sums[N];  //vector with progressive sums of the block averages
    for(int i=0; i<N; i++)
        progressive_sums[i]=0;
    
    double progressive_sq_sums[N];   //as above, but of the squared averages
    for(int i=0; i<N; i++)
        progressive_sq_sums[i]=0;
    
    double progressive_errors[N];   //the errors for each progressive sum
    for(int i=0; i<N; i++)
        progressive_errors[i]=0;
    
    out << "Throws" << " " << "Averages " << "Errors" << endl;
    for(int i=0; i<N; i++){
        for(int j=0; j<i+1; j++){
            progressive_sums[i] += avg[j];
            progressive_sq_sums[i] += sq_avg[j];
        }
        progressive_sums[i] = progressive_sums[i] / (i+1);        //cumulative averages
        progressive_sq_sums[i] = progressive_sq_sums[i] / (i+1);  //cumulative averages of the squares
        progressive_errors[i] = error(progressive_sums, progressive_sq_sums, i);
        out << i*M << " " << progressive_sums[i] << " " << progressive_errors[i] << endl;
    }
    
    //
    //
    //  PARTE 2
    //
    //
    double f2=0;
    double x2=0;
    
    double avg2[N];
    for(int i=0; i<N; i++)
        avg2[i] = 0;
    double sq_avg2[N];
    for(int i=0; i<N; i++)
        sq_avg2[i] = 0;
    
    for(int j=0; j<N; j++){
        
        for(int i=0; i<M; i++){
            x2 = 1 - sqrt(1-rnd.Rannyu());
            
            f2 += cos(x2*M_PI/2)/ ( 2 * ( 1 - x2));
        }
        avg2[j] = M_PI/2. * f2 / M;
        sq_avg2[j] = avg2[j] * avg2[j];
        f2=0;
    }
    
    ofstream out2;
    out2.open("part_2.txt");
    
    double progressive_sums2[N];  //vector with progressive sums of the block averages
    for(int i=0; i<N; i++)
        progressive_sums2[i]=0;
    
    double progressive_sq_sums2[N];   //as above, but of the squared averages
    for(int i=0; i<N; i++)
        progressive_sq_sums2[i]=0;
    
    double progressive_errors2[N];   //the errors for each progressive sum
    for(int i=0; i<N; i++)
        progressive_errors2[i]=0;
    
    out2 << "Averages " << "Errors" << endl;
    for(int i=0; i<N; i++){
        for(int j=0; j<i+1; j++){
            progressive_sums2[i] += avg2[j];
            progressive_sq_sums2[i] += sq_avg2[j];
        }
        progressive_sums2[i] = progressive_sums2[i] / (i+1);        //cumulative averages
        progressive_sq_sums2[i] = progressive_sq_sums2[i] / (i+1);  //cumulative averages of the squares
        progressive_errors2[i] = error(progressive_sums2, progressive_sq_sums2, i);
        out2 << progressive_sums2[i] << " " << progressive_errors2[i] << endl;
    }
    return 0;
}

double error(double* avg, double* sq_avg, int n){
    
    if(n==0)
        return 0;
    else
        return sqrt((sq_avg[n]-avg[n]*avg[n])/n);
}
