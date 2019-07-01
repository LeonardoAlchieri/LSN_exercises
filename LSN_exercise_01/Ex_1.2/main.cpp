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
#include <ctgmath>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){
    cout << endl;
    //
    //  This set of the code has been provided as it is by the professor of the course.
    //
    //  My intention is to build on top of it as the exercise requires.
    //
    //

    Random rnd;     //element Random-type.
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
   
    ofstream out;
    out.open("data.txt");
    
    int n=1E4;      //number of variables to be generated
    
    double Uniform[n][4];   //vector with uniformally distributed varibles in [0,1)
    double Exp[n][4];       //exponentially distributed variables in [0,∞)
    double Lorentz[n][4];   //Lorentzly distributed variables in (-∞,∞)
    for(int i=0; i<4; i++){
        for(int j=0; j<n; j++){
            Uniform[j][i] =0;
            Exp[j][i] =0;
            Lorentz[j][i] =0;
        }
    }
    
    int count =0;
    int N=1;
    
    
    for(int i=0; i<n; i++){
        while(count < 4){
            for(int j=0; j<N; j++){
                Uniform[i][count] += rnd.Rannyu();
                Exp[i][count] += rnd.Exp(1);
                Lorentz[i][count] += rnd.Lorentz(0,1);
            }
            
            out << Uniform[i][count]/N << " " << Exp[i][count]/N << " " << Lorentz[i][count]/N << " ";
            
            if( (Uniform[i][count]/N) != (Uniform[i][count]/N)){
                cerr << "-- Errore: i=" << i << " N="<< N << endl;
                return -1;
            }
            if(N==10){
                N=100;
            }
            if(N==2){
                N=10;
            }
            if(N==1){
                N=2;
            }
            
            
            count++;
        }
        out << endl;
        cout << "-- Progress: " << int((double)i/n*100) << "%\r";
        
        count=0;
        N=1;
    }
    cout << "-- Progress: 100%\n";
    cout << " *** Data saved on external file 'data.txt'. ***" << endl;
    cout << endl;
    
    
    rnd.SaveSeed();
    out.close();
    return 0;
}
