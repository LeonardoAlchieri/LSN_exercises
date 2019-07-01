#ifndef metropolis_h
#define metropolis_h

#include <stdio.h>
#include <algorithm>    // std::random_shuffle
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include "random.h"
#include "classes.h"


int N_cities = 30;
int type_input=0;
int n_individuals=1;

//
//
//
//
//

int attempted=0;
int accepted=0;

void mutate(vector<salesman> & individui, int i)
{
    
    double start=0, finish=1;
    int rnd_shift;
    double prob=0.2;
    //
    //  *** PERMUTATION ***
    //
    if(rnd.Rannyu() < prob)
        individui[i].permutation();
    //
    //  *** SHIFT ***
    //
    rnd_shift = int(rnd.Rannyu(0,29));
    
    if(rnd.Rannyu() < prob)
        individui[i].shift(rnd_shift);
    //
    //  *** SHIFT CONT ***
    //
    start = int(rnd.Rannyu(0,N_cities - 20));
    finish = start + int(rnd.Rannyu(0,20));
    if( start > N_cities || finish > N_cities)
    {
        cerr << "-- ! -- Errore nello shift_perm: sono stati generati dei numeri casuali maligni. Abort Mutation." << endl;
    }
    if(rnd.Rannyu() < prob)
        individui[i].shift_cont(start, finish);
    //
    //  *** SHIFT BIG ***
    //
    if(rnd.Rannyu() < prob)
        individui[i].big_shift(10);


    return;
    }

//
//
//
/*          METROPOLIS          */
//
//
//

void Move(vector<salesman> & walker, double temperature, cities citta, int N_cities)
{
    attempted++;
    
    double beta = 1/temperature;
    double energy_old, energy_new;
    double r;
    double A;
    
    vector<salesman> new_walker(1, salesman(N_cities));
    
    new_walker = walker;
    
    energy_old = walker[0].loss(citta);
    
    //
    //
    //
    mutate(new_walker, 0);
    
    energy_new = new_walker[0].loss(citta);;
    
    if( (energy_new - energy_old) < 0 )
    {
        accepted++;
        walker = new_walker;
    }
    else
    {
        r = rnd.Rannyu();
        A = exp(-beta * (energy_new - energy_old));
        
        if(r < min(1.,A))
        {
            accepted++;
            walker = new_walker;
            
        }
    }
    
    
    
}


#endif /* metropolis_h */
