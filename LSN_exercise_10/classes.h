#ifndef altro_h
#define altro_h

#include <stdio.h>
#include <algorithm>    // std::random_shuffle
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include "random.h"
int seed[4];
Random rnd;

using namespace std;

//
//
//

class cities
{
public:
    // constructor
    cities(int N_cities, int type)
    {
        double theta;
        m_N_cities = N_cities;
        if(type == 0)
        {
            cout << "--- Circumference ---" << endl;
            for(int i=0; i<N_cities; i++)
            {
                
                theta = rnd.Rannyu(0,2*M_PI);
                x.push_back(cos(theta));
                y.push_back(sin(theta));
                
            }
        }
        else if(type == 1)
        {
            cout << "--- Square ---" << endl;
            for(int i=0; i<N_cities; i++)
            {
                
                x.push_back(rnd.Rannyu());
                y.push_back(rnd.Rannyu());
            }
        }
    }
    //~cities();
    
    // functions
    double get_x(int city)
    {
        return x[city];
    }
    double get_y(int city)
    {
        return y[city];
    }
    
private:
    int m_N_cities;
    std::vector<double> x;
    std::vector<double> y;
};

//
//
//

class salesman
{
public:
    salesman(int N_cities)
    {
        //
        //  I save in this one only the random order of the cities
        //
        m_N_cities = N_cities;
        for(int i=0; i<m_N_cities; i++)
        {
            order.push_back(i);
        }
    }
    //
    //
    //
    
    void shuffle()
    {
        random_shuffle(order.begin(), order.end());
        return;
    }
    
    
    //
    //
    //

    int get_N_cities()
    {
        return m_N_cities;
    }
    
    //
    //
    //
    
    int get_position(int i)
    {
        return order[i];
    }
    
    //
    //
    //
    
    void swap(int i, int j)
    {
        double ausiliario;
        ausiliario = order[i];
        order[i] = order[j];
        order[j] = ausiliario;
        
        return;
    }
    
    //
    //
    //
    
    double loss(cities cities)
    {
        double loss=0;
        double x_1, y_1, x_2, y_2;
        double last_x, last_y, first_x, first_y;
        
        for(int i=0; i <m_N_cities-1; i++)
        {
            x_1 = cities.get_x(order[i]);
            y_1 = cities.get_y(order[i]);
            
            x_2 = cities.get_x(order[i+1]);
            y_2 = cities.get_y(order[i+1]);
            
            loss += sqrt( (x_1 - x_2)*( x_1 - x_2) + (y_1 - y_2)*(y_1 - y_2) );
        }
        last_x = cities.get_x(order[m_N_cities-1]);
        last_y = cities.get_y(order[m_N_cities-1]);
        
        first_x = cities.get_x(order[0]);
        first_y = cities.get_y(order[0]);
        
        
        loss += sqrt( (first_x - last_x)*( first_x - last_x) + (first_y - last_y)*(first_y - last_y) );
        
        m_loss = loss;
        return loss;
    }
    
    //
    //
    //
    
    void permutation()
    {
        int rand_position;
        double scambio;

        rand_position = int(rnd.Rannyu(0,m_N_cities));
        
        if( rand_position == m_N_cities-1)
        {
            scambio = order[rand_position];
            order[rand_position] = order[0];
            order[0] = scambio;
        }
        else
        {
            scambio = order[rand_position];
            order[rand_position] = order[rand_position + 1];
            order[rand_position+1] = scambio;
        }
    
        return;
    }
    
    //
    //
    //
    
    void shift(int n_shift)
    {
        rotate(order.begin(), order.begin() + n_shift, order.end());
        
        return;
    }
    //
    //
    //
    
    // I have to specify the indexes of the contiguous cities I am shifting
    void shift_cont(int n_shift, int start, int finish)
    {
        vector<double> ausiliario(30);
        double excess;
        
        if(finish >= m_N_cities || start >= m_N_cities)
        {
            cerr << "-- ! -- Bad index given. No shift will happen." << endl;
            return;
        }
        for(int i=start; i <= finish; i++)
        {
            ausiliario[i] = order[i];
            
            if( i + n_shift >= m_N_cities)
            {
                excess = (i+n_shift) - (m_N_cities);
                order[i] = order[excess];
            }
            order[i] = order[i+n_shift];
        }
        for(int i=start; i <= finish; i++)
        {
            if( i + n_shift >= m_N_cities)
            {
                excess = (i+n_shift) - (m_N_cities);
                order[excess] = ausiliario[i];
            }
            order[i+n_shift] = ausiliario[i];
        }
        
        return;
    }
    
    void shift_cont(int start, int finish)
    {
        int n_shift=1;
        int ausiliario;
        double excess;
        
        if(finish >= m_N_cities || start >= m_N_cities)
        {
            cerr << "-- ! -- Bad index given. No shift will happen." << endl;
            return;
        }
        for(int i=start; i <= finish; i++)
        {
            ausiliario = order[i];
            if( i + n_shift >= m_N_cities)
            {
                excess = (i+n_shift) - (m_N_cities);
                order[i] = order[excess];
                
                order[excess] = ausiliario;
            }
            order[i] = order[i+n_shift];
            order[i+n_shift] = ausiliario;
            
        }
        return;
    }
    
    // I move the first N cities to the end, and vice versa
    void big_shift(int N)
    {
        double aux_swap;
        int step;
        step = (m_N_cities-1) - N;
        
        for(int i=0; i<N; i++)
        {
            aux_swap = order[i];
            order[i] = order[i+step];
            order[i+step] = aux_swap;
            
        }
        
        return;
    }
    
    
    
    //
    //
    //
    
    //~salesman();
    
private:
    int m_N_cities;
    double m_loss;
    std::vector<int> order;
    
};
#endif /* altro_h */
