//
//  main.cpp
//  CFDstochasticGrowth
//
//  Created by Victor Vidal Suriel on 4/2/16.
//  Copyright © 2016 Victor Vidal Suriel. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double * ODEsolver(int N, double T, double mu, double sig, double A, double y0, int flag);

int main(int argc, const char * argv[]) {
    int N      = atoi(argv[1]); //number of sub intervals
    double T   = atof(argv[2]); //total length of interval
    double mu  = 0;             //mean
    double sig = 1;             //variance
    double A   = atof(argv[3]);
    double y0  = 1;             //initial conditions
    int trials = atoi(argv[4]);
    int flag   = atoi(argv[5]); //flag = 1 for forward Euler, 4 for Adam-Bashforth 
    
    double h = T/N;      // time step
    double C, Gamma;
    
    if( A == 0 ) {       //uncorrelated
        C = 0;
        Gamma = h*T;
    }
    else if( A == -1 ) { //perfect correlation
        C = 1;
        Gamma = T*T;
    }
    else {               //partial correlation
        C = exp( -h/A );
        Gamma = h*( T*(1+C)/(1-C) - 2*C*h*(1-pow(C,N))/pow(1-C,2) );
    }
    
    //opens file that saves all data for plotting in MATLAB
    ofstream outFile;
    outFile.open("ODEsolution.txt");
    
    outFile << N << " " << trials << " ";
    
    for( int n=0; n<N+1; n++ )
        outFile << n*h << " ";
    
    double * y = nullptr;   //vector of solutions. For n = 0,...,N, y[n] = y(nh)
    double y_final[trials]; //y_final[i] stores y_N at trial i for Monte Carlo estimates.
    
    for( int i=0; i<trials; i++ ) {
        delete [] y;
        
        y = ODEsolver(N, T, mu, sig, A, y0, flag); // main program, see function file
        
        for( int n=0; n<N+1; n++ )
            outFile << y[n] << " ";
        
        y_final[i] = y[N];
    }
     //compute Monte Carlo mean
    double mean = 0;
    for( int i=0; i<trials; i++)
        mean += y_final[i];
    
    mean /= trials;
    
    //compute Monte Carlo variance
    double var = 0;
    for( int i=0; i<trials; i++)
        var += (y_final[i] - mean) * (y_final[i] - mean);
    
    var /= trials - 1;
    
    outFile << mean << " " << var;
    
    outFile.close();
    
    //display the theoretical and computed values on command line
    cout << "Computed mean = " << mean;
    cout << "\nTheoretical mean: " << y0*exp(mu*T + sig*sig*Gamma/2);
    cout << "\nComputed variance: " << var;
    cout << "\nTheoretical variance: " << (exp(sig*sig*Gamma) - 1) * pow(y0*exp(mu*T + sig*sig*Gamma/2),2) << endl;
    return 0;
}
