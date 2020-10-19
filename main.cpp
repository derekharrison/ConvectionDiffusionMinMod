/*
 * main.cpp
 *
 *  Created on: Oct 19, 2020
 *      Author: d-w-h
 */

#include <math.h>
#include <stdio.h>

double min(double a, double b) {
    double r;
    if(a > b) {
        r = b;
    }
    else {
        r = a;
    }

    return r;
}

double psi(double r) {
    double result;
    if(r > 0) {
        result = min(r, 1);
    }
    else {
        result = 0;
    }

    return result;
}

int main(int argc, char* argv[]) {
    double rho, L, U, gamma, del_x, phi_a, phi_b;
    int N;

    /* Parameters */
    N = 50;          //Number of nodes
    L = 1.0;         //Length of domain
    rho = 1.0;       //Density
    U = 1.0;         //Fluid Velocity
    gamma = 1.0;     //Diffusion coefficient
    phi_a = 1.0;     //Value of transport property at inlet
    phi_b = 0.5;     //Value of transport property at outlet

    /* Start simulation */
    del_x = L / N;
    double* phi = new double[N];
    double* x_c = new double[N];

    /* Initialize */
    for(int i = 0; i < N; ++i) {
        phi[i] = 0.0;
        x_c[i] = i*del_x + 0.5*del_x;
    }

    /* Start Gauss-Seidel iterations */
    int max_iter = 2000;
    int it = 0;
    while(it < max_iter) {
        double re, rw;
        //Left most node
        re = 2*(phi[0] - phi_a)/(phi[1] - phi[0] + 1e-20);
        phi[0] = (gamma/(0.5*del_x)*phi_a + gamma/del_x*phi[1] + rho*U*phi_a - rho*U*0.5*psi(re)*(phi[1] - phi[0]))/(gamma/(0.5*del_x) + gamma/del_x + rho*U);

        //Second node
        rw = (2*phi[0] - 2*phi_a)/(phi[1] - phi[0] + 1e-20);
        re = (phi[1] - phi[0])/(phi[2] - phi[1] + 1e-20);
        phi[1] = (gamma/del_x*phi[0] + gamma/del_x*phi[2] + rho*U*phi[0] + rho*U*0.5*psi(rw)*(phi[1] - phi[0]) - rho*U*0.5*psi(re)*(phi[2] - phi[1]))/(gamma/del_x + gamma/del_x + rho*U);

        //Central nodes
        for(int j = 2; j < N - 1; ++j) {
            rw = (phi[j-1] - phi[j-2])/(phi[j] - phi[j-1] + 1e-20);
            re = (phi[j] - phi[j-1])/(phi[j+1] - phi[j] + 1e-20);
            phi[j] = (gamma/del_x*phi[j-1] + gamma/del_x*phi[j+1] + rho*U*phi[j-1] + rho*U*0.5*psi(rw)*(phi[j] - phi[j-1]) - rho*U*0.5*psi(re)*(phi[j+1] - phi[j]))/(gamma/del_x + gamma/del_x + rho*U);
        }

        //Last node
        rw = (phi[N-2] - phi[N-3])/(phi[N-1] - phi[N-2] + 1e-20);
        phi[N-1] = (gamma/del_x*phi[N-2] + gamma/(0.5*del_x)*phi_b + rho*U*phi[N-2] + rho*U*0.5*psi(rw)*(phi[N-1] - phi[N-2]) - rho*U*phi_b)/(gamma/del_x + gamma/(0.5*del_x));

        it++;
    }

    /* Print results */
    for(int i = 0; i < N; ++i) {
        printf("i: %i, x_c: %f, phi: %f\n", i, x_c[i], phi[i]);
    }

    return 0;
}
