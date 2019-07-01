/*
 ***********************
 ––––– EXERCISE 9 ––––––
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
#include <vector>
#include "classes.h"
#include "metropolis.h"
#include "random.h"

using namespace std;




//
//
//
//  ***FUNCTIONS***
//
//
//

void Input(int, char**);
int check(vector<salesman> );
void fitness(vector<salesman> & , cities);

//void mutate(vector<salesman> &, double, double, double, double, int);
void mutate(vector<salesman> &, double, double, double, double, int, int);
void mutate(vector<salesman> &, int);
void sex(vector<salesman> &, int, int, int);

int partition(vector<salesman> &, int, int, cities);
void quickSort(vector<salesman> &, int , int , cities);

int selection(double);

void scambia(vector<salesman> &, int, int);

void print(vector<salesman>, cities, int);

//
//
//
//
//
//
//

//
// 1 for squares, 0 for circ
//


int main(int argc, char ** argv)
{
    
    Input(argc, argv);
    
    cities citta(N_cities, type_input);
    //
    //  Per formattazione delle funzioni, devo definire un vettore lungo 1.
    //
    
    vector<salesman> uomini(1, salesman(N_cities));
    
    uomini[0].shuffle();
    
    
    if (check(uomini) != 0)
    {
        return -1;
    }
    
    ofstream out;
    out.open("results.txt");
    double temp=1.;
    int d = 15;
    
    out << setw(d) << "Temperature" << setw(d) << "Paths" << endl;
    for(int i=1; i<=1000; i++)
    {
        attempted=0;
        accepted=0;
        temp = 1./(double)i;
        
        for(int step=0; step<1000; step++)
        {
            
        
            Move(uomini, temp, citta, N_cities);
        }
        
        out << setw(d) << temp << setw(d) << uomini[0].loss(citta) << endl;
    }
    cout << endl<< setw(50) << "--- FINAL RESULTS ---" << setw(50) << endl;
    cout << setw(d) << "Temperature" << setw(d) << "Paths" << endl;
    cout << setw(d) << temp << setw(d) << uomini[0].loss(citta) << endl;
    out.close();
    
    out.open("best_conf.txt");
    
    out << setw(d) << "x" << setw(d) << "y" << endl;
    for(int i=0; i<N_cities; i++)
        out << setw(d) << citta.get_x(uomini[0].get_position(i)) << setw(d) << citta.get_y(uomini[0].get_position(i)) << endl;
    out << setw(d) << citta.get_x(uomini[0].get_position(0)) << setw(d) << citta.get_y(uomini[0].get_position(0)) << endl;
    out.close();
    
    
    cout << endl << "*** PROGRAM COMPLETED *** " << endl << endl;
    return 0;
}

//
//
//

void Input(int argc, char ** argv)
{
    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    
    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
    
    if(argc != 1)
    {
        type_input = atoi(argv[1]);
    }

    
    return;
}

//
//
//

int check(vector<salesman> individui)
{
    double sum=0.;
    double expected_sum = ((double)N_cities-1.)*((double)N_cities)/2.;
    for(int ind=0; ind<n_individuals; ind++)
    {
        sum =0.;
        for(int i=0; i<N_cities-1; i++)
        {
            sum += individui[ind].get_position(i);
            for(int j=i+1; j<N_cities; j++)
            {
                if( individui[ind].get_position(i) == individui[ind].get_position(j) )
                {
                    cerr << "ERRORE: un individuo ha due citta uguali. " << endl;
                    cerr << individui[ind].get_position(i) << " " << individui[ind].get_position(j) << endl;
                    return -1;
                }
                
                    
            }
        }
        sum += individui[ind].get_position(N_cities-1);
        
        if (sum != expected_sum)
        {
            cerr << "ERRORE: un individuo ha qualcosa che non va." << endl;
            for(int i=0; i<N_cities; i++)
            {
                cerr << individui[ind].get_position(i) << endl;
            }
            return -1;
        }
        
    }
    return 0;
}

//
//
//


int partition (vector<salesman> & arr, int low, int high, cities citta)
{
    
    double pivot = arr[high].loss(citta);    // pivot
    int i = (low - 1);  // Index of smaller element
    
    for (int j = low; j <= high- 1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if (arr[j].loss(citta) <= pivot)
        {
            i++;    // increment index of smaller element
            scambia(arr, i, j);
        }
    }
    scambia(arr, i+1, high);
    return (i + 1);
}

/* The main function that implements QuickSort
 arr[] --> Array to be sorted,
 low  --> Starting index,
 high  --> Ending index */
void quickSort(vector<salesman> & arr, int low, int high, cities citta)
{
    
    if (low < high)
    {
        
        /* pi is partitioning index, arr[p] is now
         at right place */
        
        int pi = partition(arr, low, high, citta);
        
        // Separately sort elements before
        // partition and after partition
        
        quickSort(arr, low, pi - 1, citta);
        quickSort(arr, pi + 1, high, citta);
        
    }
}

void scambia(vector<salesman> & individui, int i, int j)
{
    salesman scambio(N_cities);
    
    scambio = individui[i];
    individui[i] = individui[j];
    individui[j] = scambio;
    
    return;
}

void print(vector<salesman> individuo, cities citta, int index)
{
    cout << "Best: \n";
    for(int i=0; i<N_cities; i++)
        cout << individuo[index].get_position(i) << endl;
    
    cout << "Loss: " << individuo[index].loss(citta) << endl;
    cout << endl;
    return;
}

//
//
//
//
//
//

