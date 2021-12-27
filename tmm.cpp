#include <vector>
#include <cmath>
#include <complex>
#include <string>


#include "tmm.h"

using namespace std;
// frequency should be in rad/sec// frequency should be in rad/sec
double const c_const = 299792458.0; // m/s
double const eps0 = 8.85418782e-12; 
double const mu0 = 1.25663706e-6;


// matrix multiplication @to-do
template <typename M>
M multiply(M & matrix1, M & matrix2) { 
    M u;
  
    return u;
}



void TransferMatrix::calculate() {

    // Step 1: for a given frequency and angle of incidence 
    // set kz for all layers
    double k0 = frequency/c_const; 

    complex<double> kx = k0*sin(angle);

   

    for (Layer  l : structure ) { 
        complex<double> kz = sqrt(l.getMaterial().getN()*k0*k0 - kx*kx);
  
        l.setKz(kz);
    }

    // Step 2: calculate transfer matrix
    Ms = interfaceMatrix_s(structure[0], structure[1]);
    Mp = interfaceMatrix_p(structure[0], structure[1]);

    //for (int i=1; i != structure.size()-1; ++i) {
    //    Ms = multiply(Ms,propagate(structure[i], structure[i].getThickness()));
    //    Mp = multiply(Mp,propagate(structure[i], structure[i].getThickness()));

    //    Ms = multiply(Ms, interfaceMatrix_s(structure[i].getMaterial(), structure[i+1].getMaterial()));
    //    Mp = multiply(Mp, interfaceMatrix_p(structure[i].getMaterial(), structure[i+1].getMaterial()));
//
  //  }

    // step 3: populate S params 
    //rs = Ms[1][0]/Ms[0][0];
    //rp = Mp[1][0]/Mp[0][0];
    //ts = 1.0/Ms[0][0];
    //tp = 1.0/Mp[0][0];
    runSetup = true;

}

vector<vector<complex<double>>> TransferMatrix::interfaceMatrix_s(Layer const &  l1, Layer const & l2) {

    complex<double> eta = 0;

    complex<double> M11 = 0.5 + 0.5*eta; 
    complex<double> M12 = 0.5 - 0.5*eta; 
    complex<double> M21 = 0.5 - 0.5*eta;
    complex<double> M22 = 0.5 + 0.5*eta;

    vector<vector<complex<double>>> Ms; 
    Ms.push_back({M11,M12});
    Ms.push_back({M21,M22});

    return Ms;   
}

vector<vector<complex<double>>> TransferMatrix::interfaceMatrix_p(Layer  const & l1, Layer  const & l2) {

    complex<double> eta; 

    complex<double> M11 = 0.5 + 0.5*eta; 
    complex<double> M12 = 0.5 - 0.5*eta; 
    complex<double> M21 = 0.5 - 0.5*eta;
    complex<double> M22 = 0.5 + 0.5*eta;

    vector<vector<complex<double>>> Mp; 
    Mp.push_back({M11,M12});
    Mp.push_back({M21,M22});

    return Mp;   
}

vector<vector<complex<double>>> propagate(Layer  & layer, double distance) {
    complex<double> kz = layer.getKz();
    complex<double> imagI = complex<double>{0.0,1.0};
    vector<vector<complex<double>>> P = {{exp(-imagI*kz*distance),0},{0,exp(+imagI*kz*distance)}};
    return P;

}