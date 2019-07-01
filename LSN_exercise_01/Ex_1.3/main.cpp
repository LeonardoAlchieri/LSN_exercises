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
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

double error(double*, double*, int);
 
int main (int argc, char *argv[]){
    //
    //  This set of the code has been provided as it is by the professor of the course.
    //
    //  My intention is to build on top of it as the exercise requires.
    //
    //

    Random rnd;  //element Random-type.
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
    //  The following code has been done as part of the exercise.
    //
    //
    
    int throws=1E4;       //number of throws
    //int d=1;                //length of lines. Using 1 simplifies the calculations
    double L=0.5;           //length of the needle.
    int N=10000;              //interval of the x-axis chosen
    //int n_lines = N/d;      // == N if d=1
    
    int hits=0;             //number of hits
    int misses=0;           //number of misses
    double x1;              //starting point of the segment - to be randomly generated ∈ [0,N)
    double theta;           //inclination of segment with respect to x axis - randomly generated ∈ [0,360º)
    double x2 ;             //endline point of segment on x-axis - generated as x1 + L*cos(theta)
    
    int blocks=100;              //number of blocks
    
    double sum_pi[blocks];        //array with the averages
    for(int i=0; i<blocks; i++)
        sum_pi[i] = 0;
    
    double sum_sq_pi[blocks];      //as above, but with the squares
    for(int i=0; i<blocks; i++)
        sum_sq_pi[i] = 0;
    
    for(int j=0; j<blocks; j++){
        hits=0;
        misses=0;
        for(int i=0; i<throws; i++){
            x1 = rnd.Rannyu(0,N);
            theta = rnd.Rannyu(0,360);
            x2 = x1 + L * cos(theta);
            if(int(x1) != int(x2))
                hits++;
            else
                misses++;
            
            
        }
        sum_pi[j] += (double)(2 * L * throws) / (hits);
        sum_sq_pi[j] += (double)(2 * L * throws) / (hits) * (double)(2 * L * throws) / (hits);
        
        //to check that everything is okay
        if(hits + misses != throws){
            cerr << "-- Error" << endl;
            return -1;
        }
        
        cout << "-- Progress: " << int((double)j/blocks*100) << "%\r";
    }
    cout << "-- Progress: 100%\n";
    
    ofstream out;
    out.open("pies.txt");
    
    
    double progressive_sums[blocks];
    for(int i=0; i<blocks; i++)
        progressive_sums[i] =0;
    
    double progressive_sq_sums[blocks];
    for(int i=0; i<blocks; i++)
        progressive_sq_sums[i] =0;
    
    double progressive_errors[blocks];
    for(int i=0; i<blocks; i++)
        progressive_errors[i] =0;
    
    for(int i=0; i<blocks; i++){
        for(int j=0; j<i+1; j++){
            progressive_sums[i] += sum_pi[j];
            progressive_sq_sums[i] += sum_sq_pi[j];
        }
        progressive_sums[i] = progressive_sums[i] / (i+1);        //cumulative averages
        progressive_sq_sums[i] = progressive_sq_sums[i] / (i+1);  //cumulative averages of the squares
        progressive_errors[i] = error(progressive_sums, progressive_sq_sums, i);
        out << progressive_sums[i] << " " << progressive_errors[i] << endl;
    }
    
    out.close();
    
    rnd.SaveSeed();
    return 0;
}

double error(double* avg, double* sq_avg, int n){
    
    if(n==0)
        return 0;
    else
        return sqrt((sq_avg[n]-avg[n]*avg[n])/n);
}
