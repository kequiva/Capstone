#include <iostream>
#include <fstream>
#include "cosmo.h"

using namespace std;

int main(){
    // Open a file
    string filename = "redshifts.txt";
    ifstream file;
    file.open(filename);

    // Frist line of file contains the three parameters
    double C1,C2,C3;
    file >> C1 >> C2 >> C3;

    // Second line of file contains the number of redshifts
    int N;
    file >> N;

    cout << N << endl;

    // Create an array for the redshifts
    double* redshifts = new double[N]; 

    // Then we read the redshifts
    for (int i=0;i<N;i++){
        file >> redshifts[i];
    }

    // Close the file
    file.close();

    // Open the output file
    ofstream outfile;
    outfile.open("results.csv");

    // Write a header
    outfile << "Angular Diameter Distance (Mpc), Luminosity Distance (Mpc), Comoving Radial Distance (Mpc), Comoving Transverse Distance (Mpc)" << endl;

    // Now create the cosmo stuff
    Cosmo* c = new Cosmo(C1,C2,C3);

    // For all the redshift values print out the necessary constants on the file
    for (int i=0;i<N;i++){
        c->setRedshift(redshifts[i]);

        // Write the data to the file
        outfile << c->dA() << "," << c->dL() << "," << c->dC() << "," << c->dM() << endl;
        cout << redshifts[i] << endl;
    }

    // Close the output file
    outfile.close();

    return 0;
}

