//
//  ODEsolver.cpp
//  CFDstochasticGrowth
//
//  Created by Victor Vidal Suriel on 4/2/16.
//  Copyright © 2016 Victor Vidal Suriel. All rights reserved.
//
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
using namespace std;

/*
 ODE solver solves the forward order difference equation and returns a double* y which points
an array of function values y[i] = solution at ih. The inputs are
 
 N    = number of intervals
 T    = interval length
 sig  = variance of k_n
 A    = correlation length
 y0   = initial conditions
 flag = 1 for forward euler, 4 for 4th order
 */
double *ODEsolver(int N, double T, double mu, double sig, double A, double y0, int flag)
{
    double h = T/N;      //timestep
    double C;            //correlation between k_n and k_n+1
    
    if( A == 0 )
        C = 0;           //uncorrelated
    else if( A == -1)
        C = 1;           //perfect correlation
    else
        C = exp( -h/A ); //partial correlation
    
    //Create random number generator and RNG seed
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    
    //Creates normal distribution
    normal_distribution<double> distribution (0.0,1.0);
    
    //Create N normal random variables with mean mu and variance sig
    double * k;
    k = new double[N];
    
    k[0] = mu + sig*distribution(generator);
    
    for( int i=1; i<N; i++ ) 
        k[i] = C*(k[i-1]- mu) + sig * sqrt(1-C*C) * distribution(generator) + mu;
    

    //Creates N+1 array y which stores solutions
    double * y;
    y = new double[N+1];
    
    y[0] = y0;
    
    if( flag == 1 ) //forward Euler 1st order
    {
        for ( int i=1; i<N+1; i++)
            y[i] = y[i-1] + h * k[i-1]*y[i-1];
    }
    else if( flag == 4 ) //Adams-Bashforth 4th order
    {
        y[1] = y[0] + h * k[0]*y[0];
        
        y[2] = y[1] + h * ( (3/2)*k[1]*y[1] - (1/2)*k[0]*y[0] );
        
        y[3] = y[2] + h * ( (23/12)*k[2]*y[2] - (4/3)*k[1]*y[1] + (5/12)*k[0]*y[0] );
        
        for( int i=4; i<N+1; i++ )
            y[i] = y[i-1] + h * ( (55/24)*k[i-1]*y[i-1] - (59/24)*k[i-2]*y[i-2] +
                                 (37/24)*k[i-3]*y[i-3] - (3/8)*k[i-4]*y[i-4] );
    }
    else // if flag is anything besides 1 or 4
    {
        cout << "Error\n";
        exit(-1);
    }
    delete [] k;
    return y; //pointer which points to solution vector
}
