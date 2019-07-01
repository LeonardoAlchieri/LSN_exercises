#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

double error(double*, double*, int);        //function to calculate errors

int main (int argc, char *argv[]){
    
    ofstream out;
    out.open("averages.txt");
    
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
    // ### PRIMA PARTE ###
    //
    //
    
    int M = 1E4;      //number of random numbers to generate
    int N = 100;        //number of blocks
    int L = int(M/N);    //number of throws in each block
    
    cout << "M: " << M << endl;
    cout << "N: " << N << endl;
    //generating M random numbers in [0,1), using the provided library.
    double numbers[M];    //vector with the M random numbers
    for(int i=0; i<M; i++){
        numbers[i] = rnd.Rannyu();
    }
    
    
    double averages[N];     //vector with the averages of L throws in each block
    for(int i=0; i<N; i++)
        averages[i]=0;
        
    double sq_averages[N];  //vector with the squared averages (to calculate error)
    
    double sum = 0;     //auxiliary variable
    for(int i=0; i<N; i++){
        sum = 0;
        for(int j=0; j<L; j++){
            sum += numbers[j+i*L];
        }
        averages[i] = sum / L ;
        sq_averages[i] = averages[i] * averages[i] ;
    }
    
    double progressive_sums[N];  //vector with progressive sums of the block averages
    for(int i=0; i<N; i++)
        progressive_sums[i]=0;
    
    double progressive_sq_sums[N];   //as above, but of the squared averages
    for(int i=0; i<N; i++)
        progressive_sq_sums[i]=0;
    
    double progressive_errors[N];   //the errors for each progressive sum
    for(int i=0; i<N; i++)
        progressive_errors[i]=0;
    
    out << "Averages " << "Errors" << endl;
    for(int i=0; i<N; i++){
        for(int j=0; j<i+1; j++){
            progressive_sums[i] += averages[j];
            progressive_sq_sums[i] += sq_averages[j];
        }
        progressive_sums[i] = progressive_sums[i] / (i+1);        //cumulative averages
        progressive_sq_sums[i] = progressive_sq_sums[i] / (i+1);  //cumulative averages of the squares
        progressive_errors[i] = error(progressive_sums, progressive_sq_sums, i);
        out << progressive_sums[i] << " " << progressive_errors[i] << endl;
    }
    
    //
    // ### SECONDA PARTE ###
    //
    //
    //
    // Esegueo praticamente uguale a prima. Si risparmia il commento assiduo
    //
    ofstream out2;
    out2.open("errors.txt");
    
    
    double averages_err[N];
    double averages_sq_err[N];
    
    for(int i=0; i<N; i++){
        sum = 0;
        for(int j=0; j<L; j++){
            sum += (numbers[j+i*L] - 0.5)*(numbers[j+i*L] - 0.5);
        }
        averages_err[i] = sum/L;
        averages_sq_err[i] = averages_err[i]*averages_err[i];
    }
    
    double progressive_sums_err[N];
    for(int i=0; i<N; i++)
        progressive_sums_err[i]=0;
    
    double progressive_sq_sums_err[N];
    for(int i=0; i<N; i++)
        progressive_sq_sums_err[i]=0;
    
    double progressive_errors_err[N];
    for(int i=0; i<N; i++)
        progressive_errors_err[i]=0;
    
    out2 << "Averages " << "Errors" << endl;
    for(int i=0; i<N; i++){
        for(int j=0; j<i+1; j++){
            progressive_sums_err[i]+=averages_err[j];
            progressive_sq_sums_err[i]+=averages_sq_err[j];
        }
        progressive_sums_err[i] = progressive_sums_err[i]/(i+1);
        progressive_sq_sums_err[i] = progressive_sq_sums_err[i]/(i+1);
        progressive_errors_err[i] = error(progressive_sums_err, progressive_sq_sums_err, i);
        out2 << progressive_sums_err[i] << " " << progressive_errors_err[i] << endl;
    }
    
    //
    //
    // ### TERZA PARTE ###
    //
    //
    
    ofstream out3;
    out3.open("ki.txt");
    
    int m=100;                                  //numero di intervalli
    int n=0;                                    //variabile delle frequenze
    double number = -1;                         //variabile dove imagazzinare i numeri random generati
    //
    //  N_TAL non mi serve salvare i numeri dentro un array: devo solametne andare a controllare la sequenza
    //
    double avg_ki2[100];
    for(int i=0; i<100; i++)
        avg_ki2[i] =0;
    
    double progressive_ki2[100];
    for(int i=0; i<100; i++)
        progressive_ki2[i] =0;
    
    double progressive_sq_ki2[100];
    for(int i=0; i<100; i++)
        progressive_sq_ki2[i] =0;
    
    double progressive_err_ki2[100];
    for(int i=0; i<100; i++)
        progressive_err_ki2[i] =0;
    
    
    double ki2 =0;
    int N_tot = 1E4;                            //numero di casuali calcolati ogni ciclo
    double Expected = double(N_tot/m);
    for(int nn=0; nn<100; nn++){
        for(int k=0; k<m; k++){
            ki2 =0;
            n=0;                                    //riazzero la frequenza a ogni ciclo
            for(int i=0; i<N_tot; i++){
                number = int(rnd.Rannyu(0,m));
                if(number == k){                    //aumento la frequenza solamente se i numeri coincidono
                    n++;
                }
                if(number == -1){
                    cerr << "error" << endl;
                }
            }
        
            ki2 = double((n-Expected)*(n-Expected)/(Expected));       //esegueo il test del ki^2
            
            progressive_ki2[nn] += ki2;
            //progressive_sq_ki2[nn] += ki2 * ki2;
        }
        
        for(int i=0; i<nn+1; i++){
            avg_ki2[nn] += progressive_ki2[i];
            progressive_sq_ki2[nn] += progressive_ki2[i]*progressive_ki2[i];
        }
        progressive_sq_ki2[nn] = progressive_sq_ki2[nn] / (nn+1);
        avg_ki2[nn] = avg_ki2[nn] / (nn+1);
        progressive_err_ki2[nn] = error(avg_ki2, progressive_sq_ki2, nn);
        cout << "Progress: " << nn << "%\r" << endl;
        out3 << progressive_ki2[nn] << " " << avg_ki2[nn] << " " << progressive_err_ki2[nn] << endl;
    }
    
    out.close();
    out2.close();
    out3.close();
    
    rnd.SaveSeed();
    return 0;
}

double error(double* avg, double* sq_avg, int n){
    
    if(n==0)
        return 0;
    else
        return sqrt((sq_avg[n]-avg[n]*avg[n])/n);
}
