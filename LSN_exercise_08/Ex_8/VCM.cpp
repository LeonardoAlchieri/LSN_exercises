/*
 ***********************
 ––––– EXERCISE 8 ––––––
 ***********************
 
 Leonardo Alchieri, 886810
 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Laboratorio di Simulazione Numerica
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "VCM.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    // Initialize the system
    
    
    Input(argc, argv);
    
    
    
    for(int i_blk = 1; i_blk <= n_blk; i_blk++)
    {
        Reset(i_blk);   //Reset block averages
        for(int i_step = 1; i_step <= n_step; i_step++)
        {
            Move();
            Measure();
            Accumulate();
        }
        
        
        Averages(i_blk);
    }
    cout << "Acceptance rate (for final block): " << (double)accepted/(double)attempted << endl;
    
    return 0;
}

void Input(int argc,char ** argv)
{
    ifstream ReadInput,ReadConf;
    
    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    
    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
    
    if(argc == 1)
    {
        /*
        cout << "-- Input from: input.dat" << endl;
        ReadInput.open("input.dat");
         */
    }
    else
    {
        /*
        cout << "-- Input from: " << argv[1] << endl;
        ReadInput.open(argv[1]);
         */
        sigma = atof(argv[1]);
        mu = atof(argv[2]);
    }
    cout << "sigma: " << sigma << " mu: " << mu << endl;
    /*
     
     read N_part, sigma, mu, x, delta;
     
    */
    
    bin_length = (sup - inf)/(double)n_bins;
    
    return;
}

void Move()
{
    // square of the wave function |ψ|² before and after the move
    double p_old, p_new;
    // new proposed position
    double x_new;
    
    double test;
    
    for(int i=0; i<N_part; i++)
    {
        p_old = wave_function(x) * wave_function(x);
        
        x_new = x + delta*(rnd.Rannyu(-1.,1.));
        p_new = wave_function(x_new) * wave_function(x_new);
        
        // Metropolis test
        test = min( p_new/p_old, 1. );
        
        if(rnd.Rannyu() <= test)
        {
            x = x_new;
            
            accepted++;
        }
        attempted++;
    }
    return;
}

double wave_function(double x)
{
    double psi;
    
    psi = exp( - (x - mu)*(x - mu) / (2. * sigma*sigma) ) + exp( - (x + mu)*(x + mu) / (2. * sigma*sigma) );
    
    return psi;
}

void Measure()
{
    // |ψ|² / ∫|ψ|²dx
    double p;
    double H=0.;
    int position =-1;
    
    p = wave_function(x) * wave_function(x);
        
    H = (-0.5 * second_derivative(x)  + potential(x) * wave_function(x)) / (wave_function(x)) ;
    if( x > inf && x < sup)
    {
        position = int( (x - inf)/bin_length );
        histogram[position]++;
    }
    else
    {
        cout << "problemo" << endl;
    }
    walker = H;
    return;
}

double second_derivative(double x)
{
    double h;
    double first_term;
    double second_term;
    
    
    first_term = exp( - (x - mu)*(x - mu) / (2. * sigma*sigma) ) * ( (x-mu)*(x-mu) / (sigma*sigma*sigma*sigma) - 1./(sigma*sigma));
    
    second_term = exp( - (x + mu)*(x + mu) / (2. * sigma*sigma) ) * ( (x+mu)*(x+mu) / (sigma*sigma*sigma*sigma) - 1./(sigma*sigma));
    
    h = first_term + second_term;
    
    return h;
}

double potential(double x)
{
    double V;
    
    V = pow(x,4.)-2.5*pow(x,2.);
    
    return V;
}

void Accumulate()
{
    blk_avg += walker;
    blk_norm++;
    
    for(int i=0; i<n_bins; i++)
    {
        blk_histo[i] += histogram[i];
        tot += histogram[i];
    }
        
    return;
}

void Averages(int i)
{
    double estimate_H=0.;
    double x_position;
    double estimate_histo=0.;
    double err_histo;
    ofstream Hamiltonian, histogram;
    
    Hamiltonian.open("output.Hamiltonian.0", ios::app);
    
    estimate_H = blk_avg/((double)blk_norm);
    glob_avg += estimate_H;
    glob_avg_sq += estimate_H*estimate_H;
    
    error_H = Error(glob_avg, glob_avg_sq, i);
    
    Hamiltonian << " " << i << " " << estimate_H << " " << glob_avg/(double)i << " " << error_H << endl;
    
    Hamiltonian.close();
    
    histogram.open("output.histo.0");
    
    for(int k=0; k<n_bins; k++)
    {
        x_position = k * bin_length + inf;
        estimate_histo = blk_histo[k]/(bin_length * tot);
        glob_histo[k] += estimate_histo;
        glob_histo_sq[k] += estimate_histo*estimate_histo;
        
        err_histo = Error(glob_histo[k], glob_histo_sq[k], i);
        
        histogram << " " << x_position << " " << estimate_histo << " " << glob_histo[k]/(double)i << " " << err_histo << endl;
    }
    
    return;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Reset(int iblk) //Reset block averages
{
    blk_avg = 0.;
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
    tot =0;
    for(int i=0; i<n_bins; i++)
    {   
        histogram[i] = 0;
        blk_histo[i] = 0;
    }
}
