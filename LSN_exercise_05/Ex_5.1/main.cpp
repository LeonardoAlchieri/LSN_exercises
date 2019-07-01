/*
 ***********************
 ––––– EXERCISE 5 ––––––
 ***********************
 
 Leonardo Alchieri, 886810
 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Laboratorio di Simulazione Numerica
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

double evaluate_H_1s(double, double);
double evaluate_H_2p(double, double, double, double);
double error(double, double, int);
 
int main (int argc, char *argv[]){

    //
    //  ** Setting up random number generator **
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

    /*
     
        ~~~~~ PART 1 ~~~~~
     
     */
    
    ofstream out;
    ofstream threeDmap;
    out.open("results/1s_distances.txt");
    threeDmap.open("results/3Dmap_1s.txt");
    //
    //  ** Declarations **
    //
    int M=1E6;      //number of throws
    int N=1E2;      //number of blocks
    
    //
    //  ## NOTE: the coordinate system is sat so that Bohr's radius is a_0=1
    //
    double step=1.;
    
    //
    //  Starting positions
    //
    double x=0.;
    double y=0.;
    double z=0.;
    //
    //  I gave a "bad" position, i.e. easily visible when plotted.
    //
    
    double r, x_new, y_new, z_new;
    double r_prime;
    double alpha;
    double acceptor;
    
    
    double prog_sum=0.;
    double prog_sq_sum=0.;
    
    double sum[N], sq_sum[N];
    for(int i=0; i<N; i++)
        sum[i] = sq_sum[i] = 0;
    
    
    r = sqrt(x*x + y*y + z*z);
    
    int ncount=0;
    
    threeDmap << x << " " << y << " " << z << endl;
    
    
    for(int j=0; j<N; j++){
        
        for(int i=1; i<int(M/N); i++){

            x_new = rnd.Rannyu(x-step, x+step);
            y_new = rnd.Rannyu(y-step, y+step);
            z_new = rnd.Rannyu(z-step, z+step);
            
            
            r_prime = sqrt(x_new * x_new + y_new * y_new + z_new * z_new);
            
            alpha = evaluate_H_1s(r, r_prime);
            
            acceptor = rnd.Rannyu();
            
            if(acceptor <= alpha){
                x = x_new;
                y = y_new;
                z = z_new;
                
                ncount++;
                if(ncount%10==0)
                    threeDmap << x << " " << y << " " << z << endl;
            }
            
            r = sqrt(x*x + y*y + z*z);
            sum[j] += r;
        }

        sum[j] /= int(M/N);
        sq_sum[j] = sum[j]*sum[j];
        
    }
    
    
    out << "averages errors" << endl;
    for(int i=0; i<N; i++){
        prog_sum=0.;
        prog_sq_sum=0.;
        for(int j=0; j<i+1; j++){
            prog_sum += sum[j];
            prog_sq_sum += sq_sum[j];
        }
        prog_sum /= (i+1);
        prog_sq_sum /= (i+1);
        
        out << prog_sum << " " << error(prog_sum, prog_sq_sum, i) << endl;
    }
    
    out.close();
    threeDmap.close();
    /*
     
     ~~~~~ PART 2 ~~~~~
     
     */
    out.open("results/2p_distances.txt");
    
    threeDmap.open("results/3Dmap_2p.txt");
    
    
    for(int i=0; i<N; i++)
        sum[i] = sq_sum[i] = 0;
    
    prog_sum=0.;
    prog_sq_sum=0.;
    r=0.;
    x=y=0.;
    z=0.5;
    
    ncount=0;
    
    threeDmap << x << " " << y << " " << z << endl;
    
    for(int j=0; j<N; j++){
        for(int i=1; i<int(M/N); i++){
            
            x_new = rnd.Rannyu(x-step, x+step);
            y_new = rnd.Rannyu(y-step, y+step);
            z_new = rnd.Rannyu(z-step, z+step);
            
            r_prime = sqrt(x_new * x_new + y_new * y_new + z_new * z_new);
            
            alpha = evaluate_H_2p(r, r_prime, z, z_new);
            
            acceptor = rnd.Rannyu();
            
            if(acceptor <= alpha){
                x = x_new;
                y = y_new;
                z = z_new;
                
                ncount++;
                if(ncount%10==0)
                    threeDmap << x << " " << y << " " << z << endl;
                
            }
            
            r = sqrt(x*x + y*y + z*z);
            
            sum[j] += r;
            
        }
        sum[j] /= int(M/N);
        sq_sum[j] = sum[j]*sum[j];
        
    }
    
    
    out << "averages errors" << endl;
    for(int i=0; i<N; i++){
        prog_sum=0.;
        prog_sq_sum=0.;
        for(int j=0; j<i+1; j++){
            prog_sum += sum[j];
            prog_sq_sum += sq_sum[j];
        }
        prog_sum /= (i+1);
        prog_sq_sum /= (i+1);
        
        out << prog_sum << " " << error(prog_sum, prog_sq_sum, i) << endl;
    }
    
    threeDmap.close();
    out.close();
    
    threeDmap.clear();
    out.clear();
    
    rnd.SaveSeed();
    return 0;
}
//
//
//
//
//

double evaluate_H_1s(double x, double x_prime){
    double alpha;
    double move;
    
    move = exp(-2 * x_prime) / exp(-2 * x);
    
    alpha = min(1., move);
    
    return alpha;
}

double evaluate_H_2p(double r, double r_prime, double z, double z_prime){
    double alpha;
    double move;
    
    move = ( z_prime * z_prime * exp(- r_prime) ) / ( z * z * exp(-r));
    
    alpha = min(1., move);
    
    return alpha;
}



double error(double avg, double sq_avg, int n){
    
    if(n==0)
        return 0;
    else{
        return sqrt((sq_avg-avg*avg)/n);
        
    }
}


/*
 ***********************
 ––––– EXERCISE 5 ––––––
 ***********************
 
 Leonardo Alchieri, 886810
 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Laboratorio di Simulazione Numerica
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
