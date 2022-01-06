/**
 * @file calculation.cpp
 * @author Charalampos Mavidis (iphpxs@gmail.com)
 * @brief Example of transfer matrix implementation
 * @version 1.0
 * @date January 5, 2022
 * 
 * @copyright MIT license (c) 2022
 * 
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "tmm.h" 

using namespace std;
using namespace tmm;


int main() {
    const double pi = 3.1415926534;

    long double startFreq = 1.0; // THz
    long double endFreq = 10.0; // THz
    double step = 0.1; // Hz

    vector<complex<double>> rs; 
    vector<complex<double>> ts;
    vector<double> freqList;

    string filename = "diel_slabSi_normalInc.txt";
    ofstream calcFile; 
    calcFile.open(filename);


    Material  air =  Material();
    Material Si = Material(12.1);
    Layer  superstrate =  Layer(1e-6, air,"superstrate");
    Layer  dielectricLayer =  Layer(1e-6, Si,"Si");
    Layer substrate = Layer(1e-6, air);
    vector<Layer> structure{superstrate,dielectricLayer,substrate};

    
    for (double frequency = startFreq; frequency <= endFreq; frequency+=step){
        long double w = 2.0*pi*frequency*1e12;
        TransferMatrix * system = new  TransferMatrix(w,0.0,structure);

        system->calculate();

        freqList.push_back(frequency);
        rs.push_back(system->getRs());
        ts.push_back(system->getTs());

        delete system;

    }

    cout << "Calculation completed for " << freqList.size() << " frequencies.";

    // Write to file
    for (int i=0; i != freqList.size(); ++i) {
        calcFile << freqList[i] << ", " << rs[i].real() << ", " << rs[i].imag() << ", ";
        calcFile << ts[i].real() << ", " << ts[i].imag() << "\n";

    }
    calcFile.close();

    return 0; 
}