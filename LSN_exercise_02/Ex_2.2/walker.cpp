#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "walker.h"

using namespace std;

walker :: walker(double step_lenght, double x_0, double y_0, double z_0){
    step = step_lenght;
    x=x_0;
    y=y_0;
    z=z_0;
    
}

walker :: walker(){
    cerr << "-- No step length specified: set to 1. " << endl;
    cerr << "-- No origin specified: set to (0,0,0)." << endl;
    step = 1.;
    x=0.;
    y=0.;
    z=0.;
}
walker :: walker(double step_lenght){
    step = step_lenght;
    cerr << "-- No origin specified: set to (0,0,0)." << endl;
    x=0.;
    y=0.;
    z=0.;
}
walker :: ~walker(){}

double walker :: get_step(){
    return step;
}

void walker :: movement_discrete(int direction){
    if(direction == 0)
        x += step;
    if(direction == 1)
        y += step;
    if(direction == 2)
        z += step;
    if(direction > 2){
        cerr << "ERROR: wrong direction." << endl;
    }
        
}

void walker :: movement_continuum(double theta, double phi){
    if( theta > M_PI ){
        cerr << "** ERROR: theta > π." << endl;
    }
    if( phi > 2*M_PI ){
        cerr << "** ERROR: phi > 2π." << endl;
    }
    
        
    x = step * sin(theta) * cos(phi);
    y = step * sin(theta) * sin(phi);
    z = step * cos(theta);
}

double walker :: get_distance(){
    double distance=0.;
    distance = sqrt(x*x + y*y + z*z);
    return distance;
}

double walker :: get_x(){
    return x;
}

double walker :: get_y(){
    return y;
}

double walker :: get_z(){
    return z;
}
