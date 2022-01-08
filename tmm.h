/**
 * @file tmm.h
 * @author Charalampos Mavidis (iphpxs@gmail.com)
 * @brief Transfer Matrix method implementation
 * @version 1.0
 * @date January 5, 2022
 * 
 * @copyright MIT license (c) 2022
 * 
 */
#ifndef TMM_H
#define TMM_H 

#include <vector>
#include <cmath>
#include <complex>
#include <string>

using namespace std;

namespace tmm {
// frequency should be in rad/sec// frequency should be in rad/sec
double const c_const = 299792458.0; // m/s
double const eps0 = 8.85418782e-12; 
double const mu0 = 1.25663706e-6;
double const pi = 3.1415926534;
    
class Material { 
public:
    Material() : epsilon(1.0), mu(1.0) { n = sqrt(epsilon*mu); }; // vacuum constructor
    Material(complex<double> epsilon) : epsilon(epsilon), mu(1.0) { n = sqrt(epsilon*mu);};
    Material(complex<double> epsilon, complex<double> mu, string name="") 
    : epsilon(epsilon), mu(mu), name(name) { n = sqrt(epsilon*mu); };

    complex<double> getEpsilon() { return epsilon; };
    complex<double> getMu() { return mu; };
    complex<double> getN() {return sqrt(epsilon*mu); };

private:
    complex<double> epsilon; 
    complex<double> mu;
    complex<double> n;
    string name = ""; 
};

class Layer { 
public:
    Layer(double thickness,Material const &  material, string name="") : thickness(thickness), material(material), name(name) {};
    double getThickness() { return thickness; }; 
    Material getMaterial() { return material; };
    complex<double> getKz() { return kz; };
    string getName() { return name; }; 
    void setKz(complex<double> kz_) { kz = kz_; }; 

private:
    complex<double> kz; 
    double thickness; 
    Material material; 
    string name;


};

class TransferMatrix {
public:
    TransferMatrix(double frequency, double angle, vector<Layer> & structure_ ) : 
    frequency(frequency), angle(angle), structure(structure_) { runSetup = false;  };; 
    TransferMatrix(double frequency, double angle) : frequency(frequency), angle(angle) { runSetup = false;  };

    void calculate(); // calculate transfer matrix for the given structure

    // get S-parameters
    complex<double> getRs() { return rs; };
    complex<double> getRp() { return rp; };
    complex<double> getTs() { return ts; };
    complex<double> getTp() { return tp; };
    vector<complex<double>> getSparams() { return {rs,ts,rp,tp}; }; // get all the above

    // Composite matrices
    vector<vector<complex<double>>> Ms; // s polarization 2x2 matrix
    vector<vector<complex<double>>> Mp; // p polarization 2x2 matrix

    // Elementary interface matrices
    vector<vector<complex<double>>> interfaceMatrix_s(Layer  &  layer1, Layer  &  layer2);
    vector<vector<complex<double>>> interfaceMatrix_p(Layer  &  layer1, Layer  &  layer2);
    vector<vector<complex<double>>> propagate(Layer  &  layer, double distance); // propagation matrix

private:
    double frequency; 
    double angle;
    double k0;
    complex<double> kx;
    bool runSetup;
    vector<Layer> structure;
    // S parameters
    complex<double> rs = 0; 
    complex<double> rp = 0;
    complex<double> ts = 0;
    complex<double> tp = 0;   


};

}

#endif // define TMM_H

