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

using namespace std;

int type_input=0;

int N_cities = 30;
int n_individuals = 900;

int SEX_LENGTH = 4;
// cycles in each mutation
int X_RAYS = 3;
int BIG_SHIFT = 10;

double p_perm = 0.05;
double p_shift = 0.01;
double p_shift_cont = 0.04;
double p_sex = 0.5;
double p_big = 0.04;


int N_for_sex = 100;

// exponent for the selection process
double p = 5;

int N_generations = 100;

void Input(int, char**);
int check(vector<salesman> );
void fitness(vector<salesman> & , cities);

void mutate(vector<salesman> &, double, double, double, double, int);
void mutate(vector<salesman> &, double, double, double, double, int, int);
void mutate(vector<salesman> &, int);
void sex(vector<salesman> &, int, int, int);

int partition(vector<salesman> &, int, int, cities);
void quickSort(vector<salesman> &, int , int , cities);

int selection(double);

void scambia(vector<salesman> &, int, int);

void print(vector<salesman>, cities, int);


int main(int argc, char ** argv)
{
    
    Input(argc, argv);
    
    //
    // 1 for squares, 0 for circ
    //
    
    cities citta(N_cities, type_input);
    
    vector<salesman> uomini(n_individuals, salesman(N_cities));
    
    //
    //  Inzializzando con la classe <vector>, gli individui (<uomini>) sono tutti uguali. Richiamo una volta su tutti uno shuffle.
    //
    for(int i=0; i<n_individuals; i++)
    {
        uomini[i].shuffle();
    }
    
    if (check(uomini) != 0)
    {
        return -1;
    }
    
    //
    //  Order the data based on the loss value.
    //
    
    quickSort(uomini, 0, n_individuals-1, citta);
        
    //
    //
    //
    
    //
    //  Controlla che funzioni bene la copiatura nel quoickSort
    //
    
    
    vector<double> best_loss;
    vector<int> mamma(N_for_sex, 0);
    vector<int> papa(N_for_sex, 0);
    vector<salesman> new_son_1(N_for_sex, salesman(N_cities));
    vector<salesman> new_son_2(N_for_sex, salesman(N_cities));
    vector<salesman> old_best(4, salesman(N_cities));
    vector<double> avgs;
    
    int papa_safe=0;
    int mamma_safe=0;
    int n=0;
    
    cout << endl << " --- Best initial loss: " << uomini[0].loss(citta) << endl;
    best_loss.push_back(uomini[0].loss(citta));
    
    for(int i=0; i<N_generations; i++)
    {
        
        
        
        for(int j=0; j<n_individuals; j++)
        {
            mutate(uomini, p_perm, p_shift, p_shift_cont, p_big, j, X_RAYS);
        }
        quickSort(uomini, 0, n_individuals-1, citta);
        for(int j=0; j<1; j++)
        {
            old_best[j] = uomini[j];

        }
        
        for(int j=0; j<N_for_sex; j++)
        {
            n=0;
            
            while(n == 0)
            {
                mamma_safe = selection(p);
                papa_safe = selection(p);
                
                n=1;
                //
                //  I know it's stupid, but without this cicle, for some reason,
                //  the code gives worse results (but not my much).
                //  It obviously runs faster though.
                //
                for(int count=0; count<j; count++)
                {
                    if(mamma_safe == papa_safe)
                        n=0;        
                    
                }
            }
            
            mamma[j] = mamma_safe;
            papa[j] = papa_safe;
            
        }
        
        
       
        
        for(int j=0; j<N_for_sex; j++)
        {
            if (rnd.Rannyu() < p_sex)
            {
                sex(uomini, papa[j], mamma[j], SEX_LENGTH);
            }
            new_son_1[j] = uomini[papa[j]];
            new_son_2[j] = uomini[mamma[j]];
        }
        
        /*
        for(int j=0; j<n_individuals; j++)
        {
            uomini[j].shuffle();
        }
        */
        for(int j=0; j<N_for_sex; j++)
        {
            uomini[j] = new_son_1[j];
            uomini[j+N_for_sex] = new_son_2[j];
        }
        
        for(int j=0; j<1; j++)
            uomini[j] = old_best[j];
        
        quickSort(uomini, 0, n_individuals-1, citta);
        
        
        
        if (check(uomini) != 0)
        {
            return -1;
        }
        double avg=0;
        for(int j=0; j<n_individuals/2; j++)
        {
            avg += uomini[i].loss(citta);
        }
        avg = avg/(double)((double)n_individuals/2.);
        best_loss.push_back(uomini[0].loss(citta));
        avgs.push_back(avg);
    }
   
    
    cout << " --- Best final loss: " << uomini[0].loss(citta) << endl;
    
    ofstream out;
    out.open("losses.txt");
    
    for(int i=0; i<N_generations; i++)
        out << best_loss[i] << endl;
    
    out.close();
    

    out.open("medie.txt");
    
    for(int i=0; i<N_generations; i++)
        out << avgs[i] << endl;
    
    
    out.close();
    out.clear();
    
    ofstream last_conf;
    last_conf.open("last_configuration.txt");
    
    for(int i=0; i<N_cities; i++)
        last_conf << citta.get_x(uomini[0].get_position(i)) << " " << citta.get_y(uomini[0].get_position(i)) << endl;
    
    last_conf << citta.get_x(uomini[0].get_position(0)) << " " << citta.get_y(uomini[0].get_position(0)) << endl;
    
    last_conf.close();
    last_conf.clear();
    
    
    
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
    
    ifstream readInput;
    readInput.open("input.txt");
    
    readInput >> type_input;
    cout << type_input << endl;
    if(argc != 1)
    {
        type_input = atoi(argv[1]);
    }
    
    readInput >> N_cities;
    readInput >> n_individuals;
    readInput >> N_generations;
    readInput >> N_for_sex;
    readInput >> p_perm;
    readInput >> p_shift;
    readInput >> p_shift_cont;
    readInput >> p_sex;
    readInput >> p_big;
    readInput >> p;
    readInput >> SEX_LENGTH;
    readInput >> X_RAYS;
    
    
    cout << " --- N_cities: " << N_cities << endl;
    cout << " --- n_individuals: " << n_individuals << endl;
    cout << " --- N_generations: " << N_generations << endl;
    cout << " --- N_for_sex: " << N_for_sex << endl;
    cout << " --- p_perm: " << p_perm << endl;
    cout << " --- p_shift: " << p_shift << endl;
    cout << " --- p_shift_cont: " << p_shift_cont << endl;
    cout << " --- p_sex: " << p_sex << endl;
    cout << " --- p_big: " << p_big << endl;
    cout << " --- p: " << p << endl;
    cout << " --- SEX_LENGTH: " << SEX_LENGTH << endl;
    cout << " --- X_RAYS: " << X_RAYS << endl << endl;
    readInput.close();
    readInput.clear();
    
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

void fitness(vector<salesman> & individui, cities le_mie_citta)
{
    salesman scambio(N_cities);
    for(int i=0; i<n_individuals-1; i++)
    {
        for(int j=i+1; j<n_individuals; j++)
        {
            if( individui[j].loss(le_mie_citta) < individui[i].loss(le_mie_citta) ){
                
                scambio = individui[i];
                individui[i] = individui[j];
                individui[j] = scambio;
                
            }
        }
    }
    
    return;
}

//
//
//
void mutate(vector<salesman> & individui, double prob_perm, double prob_shift, double prob_shift_cont, double prob_big, int i, int count)
{
    
    double start=0, finish=1;
    
    if( prob_perm > 1. || prob_shift > 1. || prob_shift_cont >1.)
    {
        cerr << "Invalid probabilities given. The code will not mutate." << endl;
        return;
    }
    
    if( rnd.Rannyu() < prob_perm )
    {
        individui[i].permutation(count);
    }
    if( rnd.Rannyu() < prob_shift )
    {
        for(int j=0; j<count; j++)
            individui[i].shift();
    }
    if( rnd.Rannyu() < prob_shift_cont )
    {
        for(int j=0; j<count; j++)
        {
            start = int(rnd.Rannyu(0,N_cities - 5));
            finish = start + int(rnd.Rannyu(0,5));
            if( start > N_cities || finish > N_cities)
            {
                cerr << "-- ! -- Errore nello shift_perm: sono stati generati dei numeri casuali maligni. Abort Mutation." << endl;
            }
            individui[i].shift_cont(start, finish);
        }
    }
    if(rnd.Rannyu() < prob_big)
    {
        individui[i].big_shift(BIG_SHIFT);
    }
    
    
    return;
}

void mutate(vector<salesman> & individui, double prob_perm, double prob_shift, double prob_shift_cont, double prob_big, int i)
{
    
    double start=0, finish=1;
    
    if( prob_perm > 1. || prob_shift > 1. || prob_shift_cont >1.)
    {
        cerr << "Invalid probabilities given. The code will not mutate." << endl;
        return;
    }
    
    if( rnd.Rannyu() < prob_perm )
    {
        individui[i].permutation(10);
    }
    if( rnd.Rannyu() < prob_shift )
    {
        individui[i].shift(3);
    }
    if( rnd.Rannyu() < prob_shift_cont )
    {
        start = int(rnd.Rannyu(0,N_cities - 20));
        finish = start + int(rnd.Rannyu(0,20));
        if( start > N_cities || finish > N_cities)
        {
            cerr << "-- ! -- Errore nello shift_perm: sono stati generati dei numeri casuali maligni. Abort Mutation." << endl;
        }
        individui[i].shift_cont(start, finish);
    }
    if(rnd.Rannyu() < prob_big)
    {
        individui[i].big_shift(BIG_SHIFT);
    }

    
    return;
}

void mutate(vector<salesman> & individui, int i)
{
    double prob_perm=0.05, prob_shift=0.05, prob_shift_cont=0.05;
    
    double start=0, finish=1;
    

    if( rnd.Rannyu() < prob_perm )
    {
        individui[i].permutation();
    }
    if( rnd.Rannyu() < prob_shift )
    {
        individui[i].shift();
    }
    if( rnd.Rannyu() < prob_shift_cont )
    {
        start = int(rnd.Rannyu(0,N_cities - 5));
        finish = start + int(rnd.Rannyu(0,5));
        if( start > N_cities || finish > N_cities)
        {
            cerr << "-- ! -- Errore nello shift_perm: sono stati generati dei numeri casuali maligni. Abort Mutation." << endl;
        }
        individui[i].shift_cont(start, finish);
    }

    
    
    return;
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


int selection(double p)
{
    int s=0;
    double random=0;
    
    random = pow(rnd.Rannyu(),p);
    
    s = int(n_individuals * random);
    
    return s;
}

void sex(vector<salesman> & persona, int dad, int mum, int length)
{
    // pick a part of the array that changes with sex. I pick it randomly.
    //int length = 4;
    int start = int(rnd.Rannyu(0,N_cities-length));
    int count1=0;
    int count2=0;
    int index=0;
    int index_2=0;
    
    int pos_mum[10] = {0};
    int pos_dad[10] = {0};
    
    for(int i=start; i<start + length; i++)
    {
        
        count2=0;
        count1=0;
        for(int j=0; j<N_cities; j++)
        {
            if( persona[dad].get_position(i) == persona[mum].get_position(j))
            {
                pos_dad[index] = j;
                count1++;
            }
            if( persona[dad].get_position(j) == persona[mum].get_position(i))
            {
                pos_mum[index] = j;
                count2++;
                
            }
            if(count1 > 1 || count2 > 1)
            {
                cerr << "-- ! -- ERRORE: c'è qualcosa che non va nella riproduzione. >ABORT" << endl;
                return;
            }
        }
        index++;
    }
    index=0;
    
    // Just do a "BAD SORT".
    
    double aux;
    for(int i=start; i<(start+length); i++)
    {
        index_2=index+1;
        for(int j=i+1; j<start+length; j++)
        {
            
            if(pos_dad[index_2] < pos_dad[index])
            {
                persona[dad].swap(i,j);
                
                aux = pos_dad[index];
                pos_dad[index] = pos_dad[index_2];
                pos_dad[index_2] = aux;
                
            }
            if(pos_mum[index] < pos_mum[index_2])
            {
                persona[mum].swap(i,j);
                aux = pos_mum[index];
                pos_mum[index] = pos_mum[index_2];
                pos_mum[index_2] = aux;
                
            }
            index_2++;
        }
        index++;
    }

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

