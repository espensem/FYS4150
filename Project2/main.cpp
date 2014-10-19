#include <iostream>
#include <string>
#include <armadillo>
#include <sstream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <time.h>
#include <cmath>

using namespace std;
using namespace arma;


double maxoffdiag(double **A, int * k, int * l, int n){
    // Function to find the maximum matrix element

    double max = 0.0;

    // we search area only over the diagonal, remembering that A is symmetric
    for (int i = 0; i < n; i++){
        for (int j = i+1; j < n; j++){
            if (fabs(A[i][j]) > max){
                max = fabs(A[i][j]);
                *l = i;
                *k = j;
            }
        }
    }
    return max;
}


void rotate(double **A, double **R, int k, int l, int n){
    // Function to find the values of cos and sin

    double s, c;

    if (A[k][l] != 0.0){
        double t, tau;
        tau = (A[l][l] - A[k][k]) / (2 * A[k][l]);
        if (tau > 0){
            t = 1.0 / (tau + sqrt(1.0 + tau*tau));
        }
        else{
            t = -1.0 / (-tau + sqrt(1.0 + tau*tau));
        }

        c = 1.0 / sqrt(1 + t*t);
        s = c*t;
    }
    else{
        c = 1.0;
        s = 0.0;
    }

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A[k][k];
    a_ll = A[l][l];

    // changing the matrix elements with indices k and l
    A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
    A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
    A[k][l] = 0.0; // hard-coding of the zeros
    A[l][k] = 0.0;

    // changing the remaining elements
    for (int i = 0; i < n; i++){
        if (i != k && i != l){
            a_ik = A[i][k];
            a_il = A[i][l];
            A[i][k] = c*a_ik - s*a_il;
            A[k][i] = A[i][k];
            A[i][l] = c*a_il + s*a_ik;
            A[l][i] = A[i][l];
        }

        // computing the new eigenvectors
        r_ik = R[i][k];
        r_il = R[i][l];
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;
    }
    return;
}


void jacobi_method(double ** A, double **R, int n, int n_new){
    /*
     * Jacobi's method for finding eigenvalues and
     * eigenvectors of the symmetric matrix A.
     * The eigenvalues of A will be on the diagonal of A,
     * with eigenvalue i being A[i][i].
     * The j-th component of the i-th eigenvector is
     * stored in R[i][j].
     *
     * @param double ** A: input matrix (n x n)
     * @param double ** R: empty matrix for eigenvectors (n x n)
     * @param int n: dimention of matrices
     */



    // If n != 0 write out the number of iterations in terminal
    if (n != 0){

        // setting up the eigenvector matrix
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (i == j){
                    R[i][j] = 1.0;
                }
                else{
                    R[i][j] = 0.0;

                }
            }
        }

        int k, l;
        double epsilon = 1.0e-10;
        double max_number_iterations = pow(n,3);
        int iterations = 0;
        double max_offdiag = maxoffdiag(A, &k, &l, n);

        while(fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations){
            max_offdiag = maxoffdiag(A, &k, &l, n);
            rotate(A, R, k, l, n);
            iterations++;
        }
        cout << "Number of iterations: " << iterations << "\n" << endl;
    }

    // If n = 0 write number of iterations and n to file, using n_new for each run
    else if (n==0){

        // setting up the eigenvector matrix
        for (int i = 0; i < n_new; i++){
            for (int j = 0; j < n_new; j++){
                if (i == j){
                    R[i][j] = 1.0;
                }
                else{
                    R[i][j] = 0.0;

                }
            }
        }

        int k, l;
        double epsilon = 1.0e-10;
        double max_number_iterations = pow(n_new,3);
        int iterations = 0;
        double max_offdiag = maxoffdiag(A, &k, &l, n_new);

        while(fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations){
            max_offdiag = maxoffdiag(A, &k, &l, n_new);
            rotate(A, R, k, l, n_new);
            iterations++;
        }

        stringstream s;
        string filename;

        s << "n_vs_iterations" << ".dat";

        filename = s.str();

        ofstream outfile;
        outfile.open(filename.c_str(), ofstream::out|ofstream::app);

        outfile << n_new << " " << iterations << endl;

        cout << "File written" << endl;

        outfile.close();
    }

    return;
}


int main(int argc, char *argv[]){

    /*
     * Command line arguments:
     * Input 1: n - nember of steps (integer).
     * Input 2: rho_final - end value of rho (float).
     * Input 3: omega_r - parameter reflecting the strength of the oscillator potential (float).
     * Input 4: n_electrons - number of electrons in potential (1 or 2 for this exercise) (integer).
     * Input to run ex. a: program n rho_final 1.0 1
     * Input to run ex. b: program n rho_final 1.0 1
     */

    // missing command line arguments
    if (argc < 5){
        cout << "Command line arguments not correct!\nRun as: program integer float float integer." << endl;
        cout << "Argument 1: n - nember of steps." << "\n"
             << "Argument 2: rho_final - end value of rho." << "\n"
             << "Argument 3: omega_r - parameter reflecting the strength of the oscillator potential." << "\n"
             << "Argument 4: n_electrons - number of electrons in potential (1 or 2 for this exercise)." << "\n" << endl;
        return 1;
      }

    // convert command line input strings to integer or float
    int n = atoi(argv[1]);
    if (n==0){
        // Running for multiple n if n==0
        cout << "Warning: Setting n = 0 will run for multiple n (exercise. b)" << endl;
    }

    double rho_final = atof(argv[2]);
    if (rho_final==0.0){
        // unsuccessfull run for rho_final = 0
        cout << "rho_final (second argument) must be greater than 0" << endl;
        return 1;
    }

    double omega_r = atof(argv[3]);
    if (omega_r==0.0){
        cout << "Warning: Setting omega_r = 0 (third argument) will result in a zero potential" << endl;
    }

    int n_electrons = atoi(argv[4]);
    if (n_electrons < 1 || n_electrons > 2){
        // Unsuccessfull run if n_electrons is not 1 or 2
        cout << "n_electrons (fourth argument) need to be 1 or 2 for successfull run" << endl;
        return 1;
    }

    int n_new;
    n_new = n;
    int counter;
    counter = 0;    

    // Running program for given n
    // ----------------------------------------------------
    if (n!=0){

        // Declaring lower boundary and step size
        double rho0 = 0.0;
        double h = (rho_final - rho0) / n;

        // Declaring the second derivative matrix A
        double **A = new double* [n];
        for (int i = 0; i < n; i++){
            A[i] = new double[n];
        }

        for (int i = 0; i < n; i++){
            double rho_i = rho0 + (i+1)*h;
            for (int j = 0; j < n; j++){
                if (i == j){
                    if (n_electrons==2){
                        A[i][j] = (2.0 / pow(h,2)) + (pow(omega_r, 2) * pow(rho_i, 2)) + 1./rho_i;
                    }
                    else if (n_electrons==1){
                        A[i][j] = (2.0 / pow(h,2)) + pow(rho_i, 2);
                    }
                }
                else if (i - j == 1){
                    A[i][j] = -1.0 / (h*h);
                }
                else if (j - i == 1){
                    A[i][j] = -1.0 / (h*h);

                }
            }
        }

        // Initializing eigenvector matrix R
        double **R = new double* [n];
        for (int i = 0; i < n; i++){
            R[i] = new double[n];
        }

        jacobi_method(A, R, n, 0);

        // Collecting diagonal values of final A (eigenvalues) in an eigenvalue vector
        vec eigenvalues;
        eigenvalues.zeros(n);
        for (int i = 0; i < n; i++){
            eigenvalues[i] = A[i][i];
        }

        // Finding the eigenvector corresponding to the ground state (lowest eigenvalue)
        int ground_state_number = 0;
        double smallest = 100;          // We know that the lowest eigenvalue is 3 so set to a number higher than this
        for (int i = 0; i < n; i++){
            if (eigenvalues[i] < smallest){
                smallest = eigenvalues[i];
                ground_state_number = i;

            }
        }

        cout << smallest << " " << ground_state_number << endl;

        vec R_ground_state;
        R_ground_state.zeros(n);
        for (int i = 0; i < n; i++){
            R_ground_state[i] = R[i][ground_state_number];

        }

        // Sort array eigenvalues for print out
        eigenvalues = sort(eigenvalues);
        cout << eigenvalues[0] << " " << eigenvalues[1] << " " << eigenvalues[2] << endl;

        // Writing eigenvectors to file
        stringstream s;
        string filename;

        int run = 0;
        if (omega_r == 0.01){
            run = 1;
        }
        else if (omega_r == 0.5){
            run = 2;
        }
        else if (omega_r == 1.0){
            run = 3;
        }
        else if (omega_r == 5.0){
            run = 4;
        }

        s << "Wavefunc_omega_r" << run << ".dat";

        filename = s.str();

        ofstream outfile;
        outfile.open(filename.c_str());

        for (int i = 0; i < n; i++){
            double rho_i = rho0 + (i+1)*h;
            outfile << rho_i << " " << R_ground_state[i] << endl;
        }

        cout << "File written" << endl;
        outfile.close();



    }
    // ---------------------------------------------------

    // Running program for multiple n
    // ---------------------------------------------------
    else if (n==0){

        // Clear content of file before doing a new run
        stringstream s;
        string filename;

        s << "n_vs_iterations" << ".dat";

        filename = s.str();

        ofstream outfile;
        outfile.open(filename.c_str(), ofstream::trunc);

        outfile.close();

        // Doing the Jacobi method for different n's
        for (int i = 20; i <= 300; i = i + 20){
            int n_new = i;
            // Declaring lower boundary and step size
            double rho0 = 0.0;
            double h = (rho_final - rho0) / n_new;

            // Declaring the second derivative matrix A
            double **A = new double* [n_new];
            for (int i = 0; i < n_new; i++){
                A[i] = new double[n_new];
            }

            for (int i = 0; i < n_new; i++){
                double rho_i = rho0 + (i+1)*h;
                for (int j = 0; j < n_new; j++){
                    if (i == j){
                        if (n_electrons==2){
                            A[i][j] = (2.0 / pow(h,2)) + (pow(omega_r, 2) * pow(rho_i, 2)) + 1./rho_i;
                        }
                        else if (n_electrons==1){
                            A[i][j] = (2.0 / pow(h,2)) + pow(rho_i, 2);
                        }
                    }
                    else if (i - j == 1){
                        A[i][j] = -1.0 / (h*h);
                    }
                    else if (j - i == 1){
                        A[i][j] = -1.0 / (h*h);

                    }
                }
            }

            // Initializing eigenvector matrix R
            double **R = new double* [n_new];
            for (int i = 0; i < n_new; i++){
                R[i] = new double[n_new];
            }

            jacobi_method(A, R, n, n_new);

            counter ++;
        }

        cout << "To run plotting_b.py give number of runs (number of n_values) as input." << "\n"
                "Number of n_values = " << counter << endl;

    }

    return 0;
}
