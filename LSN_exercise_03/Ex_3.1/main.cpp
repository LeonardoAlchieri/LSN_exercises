//
//
//  *** ES 3.1 ***
//
//  ©Leonardo Alchieri, 886810
//
//  Laboratorio di Simulazione Numerica
//
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "random.h"

using namespace std;

double error(double*, double*, int);        //function to calculate errors
 
int main (int argc, char *argv[]){

    //
    //
    //  Setting up the generator for random numbers
    //
    //
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
    //
    //
    //  *** PART 1 ***
    //
    //  Calculate the final price (at time T) directly using one sample of S(t).
    //
    //
    
    //
    //  __variable declaration__
    //
    double S_0 = 100.;      //asset price at t=0
    double T = 1.;          //delivery time
    double K = 100.;        //strike price
    double r = 0.1;        //risk-free interest rate
    double s = 0.25;       //volatility
    
    int M=10000;          //number of assets
    int N=100;          //number of blocks
    double S_T = 0.;
    
    double C[100] = {0.};        //call results
    double C_sq[100] = {0.};
    double P[100] = {0.};        //put results
    double P_sq[100] = {0.};
        
    double W=0;             //normally distribuited variable through N(0,T)
    
    
    
    int L = M/N;

    double prog_sum_C[100] = {0.};
    double prog_sq_sum_C[100] = {0.};
    double prog_sum_P[100] = {0.};
    double prog_sq_sum_P[100] = {0.};
    
    //
    //
    // __implementation__
    //
    //
    ofstream out;
    out.open("part2.txt");
    
    
    for(int j=0; j<N; j++){
        for(int i=0; i<L; i++){
            W = rnd.Gauss(0.,1.);
            S_T = S_0 * exp((r - 1./2. * s*s) + s * W);
            
            C[j] += exp(-r*T)*max(0., S_T - K);
            P[j] += exp(-r*T)*max(0., K - S_T);
        }
        C[j] /= double(L);
        P[j] /= double(L);
        
        C_sq[j] = C[j]*C[j];
        P_sq[j] = P[j]*P[j];
        
    }
    
    
    
    out << "OPTION PRICING: DIRECT EVALUATION" << endl;
    out << "Call " << "err_call " << "Put " << "err_put " << endl;
    for(int i=0; i<N; i++){
        for(int j=0; j<i+1; j++){
            prog_sum_C[i] += C[j];
            prog_sq_sum_C[i] += C[j] * C[j];
            
            prog_sum_P[i] += P[j];
            prog_sq_sum_P[i] += P[j] * P[j];
        }
        prog_sum_C[i] /= double(i+1);
        prog_sq_sum_C[i] /= double(i+1);
        
        prog_sum_P[i] /= double(i+1);
        prog_sq_sum_P[i] /= double(i+1);
        
        out << i << " ";
        out << prog_sum_C[i] << " " << error(prog_sum_C, prog_sq_sum_C, i);
        out << " " << prog_sum_P[i] << " " << error(prog_sum_P, prog_sq_sum_P, i) << endl;
    }
    
    
    //
    //
    //
    //  *** PART 2 ***
    //
    //
    //
    for(int i=0; i<N; i++){
        C[i] = C_sq[i] = P[i] = P_sq[i] = prog_sum_C[i] = prog_sum_P[i] = prog_sq_sum_C[i] = prog_sq_sum_P[i] = 0.;
    }
    
    double t=0.01;          //time interval
    
    //
    //
    // __implementation__
    //
    //
    out.close();
    out.open("part1.txt");
    
    
    
    for(int j=0; j<N; j++){
        for(int i=0; i<L; i++){
            
            S_T=S_0;
            
            for(int n=0; n<100; n++){
                
                W = rnd.Gauss(0.,1.);
                S_T = S_T * exp((r - 1./2. * s*s)*(t) + s * W * sqrt(t));
                
            }
            
            C[j] += exp(-r*T)*max(0., S_T - K);
            P[j] += exp(-r*T)*max(0., K - S_T);
        }
        C[j] /= double(L);
        P[j] /= double(L);
        
        C_sq[j] = C[j]*C[j];
        P_sq[j] = P[j]*P[j];
        
    }
    
    
    
    out << "OPTION PRICING: DESCRETE EVALUATION" << endl;
    out << "Call " << "err_call " << "Put " << "err_put " << endl;
    for(int i=0; i<N; i++){
        for(int j=0; j<i+1; j++){
            prog_sum_C[i] += C[j];
            prog_sq_sum_C[i] += C[j] * C[j];
            
            prog_sum_P[i] += P[j];
            prog_sq_sum_P[i] += P[j] * P[j];
        }
        prog_sum_C[i] /= double(i+1);
        prog_sq_sum_C[i] /= double(i+1);
        
        prog_sum_P[i] /= double(i+1);
        prog_sq_sum_P[i] /= double(i+1);
        
        out << i << " ";
        out << prog_sum_C[i] << " " << error(prog_sum_C, prog_sq_sum_C, i);
        out << " " << prog_sum_P[i] << " " << error(prog_sum_P, prog_sq_sum_P, i) << endl;
    }
    
    out.close();
    out.clear();
    rnd.SaveSeed();
    
    return 0;
}

//
//
//
//
double error(double* avg, double* sq_avg, int n){
    
    if(n==0)
        return 0;
    else
        return sqrt((sq_avg[n]-avg[n]*avg[n])/n);
}
