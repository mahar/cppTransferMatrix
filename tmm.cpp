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
template <typename Matrix>
Matrix matmul(Matrix & m1, Matrix & m2) {
    /**
     * Matrix multiplication utility function
     * 
     */
    int n = m1.size();
    int m = m1[0].size();
    int p = m2[1].size();

    Matrix newMatrix(n, vector<complex<double>>(p));

    for (int i=0; i < n; ++i) {
        for (int j=0; j < m; ++j) {
            complex<double> sum(0.0,0.0);
            for (int k=0; k < m; ++k) {
                sum += m1[i][k]*m2[k][j];
            }
            newMatrix[i][j] = sum;
        }
    }
    return newMatrix;
} 



void tmm::TransferMatrix::calculate() {

    // Step 1: for a given frequency and angle of incidence 
    // set kz for all layers
    if (!structure.empty() && frequency > 0) {
        k0 = frequency/c_const;

        // Initialize matrices
        Ms = {{1.0,0.0}, {0.0,1.0}};
        Mp = {{1.0,0.0}, {0.0,1.0}};
    } else {
        return;
    }

    kx = k0*sin(angle);

   
    
    for (Layer  l : structure ) { 
        complex<double> kz = sqrt(l.getMaterial().getN()*k0*k0 - kx*kx);

        l.setKz(kz);
    }
  

    // Step 2: calculate transfer matrix
    Ms = interfaceMatrix_s(structure[0], structure[1]);
    Mp = interfaceMatrix_p(structure[0], structure[1]);

    for (int i=1; i < structure.size()-1; ++i) {

        auto propagMatrix = propagate(structure[i], structure[i].getThickness());
        auto matrixS = interfaceMatrix_s(structure[i], structure[i+1]);
        auto matrixP = interfaceMatrix_p(structure[i], structure[i+1]);
        Ms = matmul(Ms,propagMatrix);
        Mp = matmul(Mp,propagMatrix);

        Ms = matmul(Ms, matrixS);
        Mp = matmul(Mp, matrixP);

     

    }

    // step 3: populate S params 
    rs = Ms[1][0]/Ms[0][0];
    rp = Mp[1][0]/Mp[0][0];
    ts = 1.0/Ms[0][0];
    tp = 1.0/Mp[0][0];
    runSetup = true;

}

vector<vector<complex<double>>> tmm::TransferMatrix::interfaceMatrix_s(Layer  &  layer1, Layer  & layer2) {
    // Material parameters
    complex<double> eps1 = layer1.getMaterial().getEpsilon();
    complex<double> mu1 = layer1.getMaterial().getMu();
    complex<double> k1 = sqrt(eps1*mu1)*k0;

    complex<double> k1z = sqrt(norm(k1) - kx*kx); // norm(z) = ||z||^2

    complex<double> eps2 = layer2.getMaterial().getEpsilon();
    complex<double> mu2 = layer2.getMaterial().getMu();
    complex<double> k2 = sqrt(eps1*mu1)*k0;
    complex<double> k2z = sqrt(norm(k2) - kx*kx); // norm(z) = ||z||^2

    complex<double> eta = k2z*mu1/(k1z*mu2);

    complex<double> M11 = 0.5 + 0.5*eta; 
    complex<double> M12 = 0.5 - 0.5*eta; 
    complex<double> M21 = 0.5 - 0.5*eta;
    complex<double> M22 = 0.5 + 0.5*eta;

    vector<vector<complex<double>>> Ms; 
    Ms.push_back({M11,M12});
    Ms.push_back({M21,M22});

    return Ms;   
}

vector<vector<complex<double>>> tmm::TransferMatrix::interfaceMatrix_p(Layer  & layer1, Layer  & layer2) {
    // Material parameters
    complex<double> eps1 = layer1.getMaterial().getEpsilon();
    complex<double> mu1 = layer1.getMaterial().getMu();
    complex<double> k1 = sqrt(eps1*mu1)*k0;

    complex<double> k1z = sqrt(norm(k1) - kx*kx); // norm(z) = ||z||^2

    complex<double> eps2 = layer2.getMaterial().getEpsilon();
    complex<double> mu2 = layer2.getMaterial().getMu();
    complex<double> k2 = sqrt(eps1*mu1)*k0;
    complex<double> k2z = sqrt(norm(k2) - kx*kx); // norm(z) = ||z||^2

    complex<double> eta = k2z*eps1/(k1z*eps2);

    complex<double> M11 = 0.5 + 0.5*eta; 
    complex<double> M12 = 0.5 - 0.5*eta; 
    complex<double> M21 = 0.5 - 0.5*eta;
    complex<double> M22 = 0.5 + 0.5*eta;

    vector<vector<complex<double>>> Mp; 
    Mp.push_back({M11,M12});
    Mp.push_back({M21,M22});

    return Mp;   
}

vector<vector<complex<double>>> tmm::TransferMatrix::propagate(Layer  & layer, double distance) {
    complex<double> kz = layer.getKz();
    complex<double> imagI = complex<double>{0.0,1.0};
    vector<vector<complex<double>>> P = {{exp(-imagI*kz*distance),0},{0,exp(+imagI*kz*distance)}};
    return P;

}