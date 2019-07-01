

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
//#include "walker.h"

using namespace std;

double distance(double, double, double);

double error(double* , double* , int );
 
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
    //
    //
    //
    //
    
    ofstream out;
    out.open("dati.txt");
    
    int T=1E2;      //number of steps.
    int N=1E4;      //blocks over which I average
    int dir=-1E2;   //the direction I take. To be chosen randomly between forward of backward in all 3D directions.
    
    double walk[N][4];      //matrix with the N averages and the 3 directions. discerte.
    for(int j=0; j<N; j++){
        for(int i=0; i<4; i++)
            walk[j][i] =0;
    }
    
    double walk_cont[N][4];      //matrix with the N averages and the 3 directions. continous.
    for(int j=0; j<N; j++){
        for(int i=0; i<4; i++)
            walk_cont[j][i] =0;
    }
    //walker jimmy;
    
    double avg_dist[T];     //array with the averages distances at a certain step.
    for(int i=0; i<T; i++)
        avg_dist[i] =0;
    
    double sq_avg_dist[T];     //array with the averages distances at a certain step.
    for(int i=0; i<T; i++)
        sq_avg_dist[i] =0;
    
    double cont_avg_dist[T];    //analogous as above, but for the continous path.
    for(int i=0; i<T; i++)
        cont_avg_dist[i] =0;
    
    double sq_cont_avg_dist[T];    //analogous as above, but for the continous path.
    for(int i=0; i<T; i++)
        sq_cont_avg_dist[i] =0;
    
    // random angles for continous random walk. θ ∈ [0,π) and φ ∈ [0,2π)
    double theta=0;
    double phi=0;
    
    double discr_err=0;
    double cont_err=0;
    
    
    out << "Step" << " " << "Average discrete" << " " << "Error avg Discrete" << " " << "Average cont." << " " << "Error avg cont." << endl;
    for(int i=0; i<T; i++){
        for(int j=0; j<N; j++){
            //
            //  *** PART 1 ***
            //
            //  Random Descrete Walk
            //
            //  Firstly, I generate a random variable ∈[0,6). It decides in which direction the step will be taken.
            //
            dir = int(rnd.Rannyu(0,6));
            //
            //  Series of conditions to figure out where to move
            //
            if(dir == 0)
                walk[j][0]++;
            if(dir == 1)
                walk[j][0]--;
            if(dir == 2)
                walk[j][1]++;
            if(dir == 3)
                walk[j][1]--;
            if(dir == 4)
                walk[j][2]++;
            if(dir == 5)
                walk[j][2]--;
            
            avg_dist[i] += distance(walk[j][0], walk[j][1], walk[j][2]);    //I sum all of the distances
            sq_avg_dist[i] += distance(walk[j][0], walk[j][1], walk[j][2]) * distance(walk[j][0], walk[j][1], walk[j][2]);
            //
            //  *** PART 2 ***
            //
            //  Random Continuous Walk
            //
            //  I generate two angle, and upgrade the coordinates accordingly.
            //
            theta = rnd.Rannyu(0,M_PI);
            phi = rnd.Rannyu(0, 2*M_PI);
            //
            //  I upgrade using spherical coordinates.
            //
            walk_cont[j][0] += sin(theta) * cos(phi);        // x
            walk_cont[j][1] += sin(theta) * sin(phi);        // y
            walk_cont[j][2] += cos(theta);                   // z
            
            cont_avg_dist[i] += distance(walk_cont[j][0], walk_cont[j][1], walk_cont[j][2]);
            sq_cont_avg_dist[i] += distance(walk_cont[j][0], walk_cont[j][1], walk_cont[j][2]) * distance(walk_cont[j][0], walk_cont[j][1], walk_cont[j][2]);
        }
        //
        //  I average
        //
        cont_avg_dist[i] /= N;
        sq_cont_avg_dist[i] /= N;
        
        avg_dist[i] /= N;
        sq_avg_dist[i] /= N;
        
        //
        //  Evaluate the statistical uncertainties
        //
        cont_err = sqrt((sq_cont_avg_dist[i]-cont_avg_dist[i]*cont_avg_dist[i]));
        discr_err = sqrt((sq_avg_dist[i]-avg_dist[i]*avg_dist[i]));
        //
        //  Output to external file
        //
        out << i << " " <<  avg_dist[i] << " " << discr_err << " " <<  cont_avg_dist[i] << " " << cont_err << endl;
        
        cout << "-- Progress: " << (double)(i/T *100.) << "% \r";
    }
    cout << "-- Progress: 100% \n";
    
    //
    //  I try to figure out a way to evaluate the errors.
    //
    
    out.close();
    return 0;
}

double distance(double x, double y, double z){
    return sqrt(x*x + y*y + z*z);
}

double error(double* avg, double* sq_avg, int n){
    
    if(n==0)
        return 0;
    else
        return sqrt((sq_avg[n]-avg[n]*avg[n])/n);
}
