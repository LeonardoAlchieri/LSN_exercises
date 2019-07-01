//
//  VCM.hpp
//  
//
//  Created by Leonardo Alchieri on 19/05/2019.
//

#ifndef VCM_hpp
#define VCM_hpp

#include <stdio.h>
#include "random.h"
int seed[4];
Random rnd;

// FUcNTIONS
void Input(int, char **);
double wave_function(double);
void Move(void);
void Measure(void);
double potential(double);
double second_derivative(double);
void Accumulate(void);
void Averages(int);
void Reset(int);

double Error(double, double, int);

// VARIABLES
double mu = 0.76;
double sigma = 0.6;
int N_part = 1;
double delta = 2.2;

int accepted=0;
int attempted=0;

int blk_norm=0;
double blk_avg=0.;

int n_blk = 100;
int n_step = 10000;

double walker = 0.;

double glob_avg=0.;
double glob_avg_sq=0.;
double error_H=0.;

int n_bins = 100;

int histogram[100];
int blk_histo[100];
double glob_histo[100];
double glob_histo_sq[100];
// limits for the x in the histogram
double inf=-3;
double sup=3;
double bin_length;

//number of total points in the histogram
double tot=0;

// starting position
double x=0.;
#endif /* VCM_hpp */
