/**
 * @file tmm.cpp
 * @author Charalampos Mavidis (iphpxs@gmail.com)
 * @brief Transfer Matrix method implementation
 * @version 1.0
 * @date January 5, 2022
 * 
 * @copyright MIT license (c) 2022
 * 
 */
#include <vector>
#include <cmath>
#include <complex>
#include <string>
#include <iostream>
#include "tmm.h"

using namespace std;


/**
 * @brief Print 2D matrix.
 * 
 * @tparam matrix 
 * @param m 
 */
template <typename matrix>
void printMatrix(matrix & m) {
    for (auto v : m) { 
        for (auto el : v) { 
            std::cout << el << " - ";
        }
        std::cout << endl;
    }
}

/**
 * @brief Matrix multiplication function for 2x2 matrices
 * 
 * @tparam matrix 
 * @param m1 
 * @param m2 
 * @return matrix 
 */
template <typename matrix>
matrix matmul( matrix & m1,  matrix & m2) {
    /**
     * Matrix multiplication utility function
     * 
     */
    int N = m1.size(); // m1 rows
    int M = m1[0].size(); // m1 cols
    int P = m2[0].size(); // m2 cols

    matrix  newMatrix = {{0.0,0}, {0,0.0}};

    newMatrix[0][0] = m1[0][0]*m2[0][0] + m1[0][1]*m2[1][0];
    newMatrix[1][1] = m1[1][0]*m2[0][1] + m1[1][1]*m2[1][1];
    newMatrix[0][1] = m1[0][0]*m2[0][1] + m1[0][1]*m2[1][1];
    newMatrix[1][0] = m1[1][0]*m2[0][0] + m1[1][1]*m2[1][0];
    return newMatrix;
} 

/**
 * @brief Calculation function.
 * 
 */
void tmm::TransferMatrix::calculate() {

    // Step 1: for a given frequency and angle of incidence 
    // set kz for all layers
    if (!structure.empty() && frequency > 0) {
        k0 = frequency/c_const;

        // Initialize matrices
        Ms = {{1.0,0.0}, {0.0,1.0}};
        Mp = {{1.0,0.0}, {0.0,1.0}};
    } else {
        throw "Frequency must be positive and structure array must not be empty.";
    }

    kx = k0*sin(angle);
    
    for (Layer &  l : structure ) { 
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

/**
 * @brief Interface Matrix for s polarization
 * 
 * @param layer1 
 * @param layer2 
 * @return vector<vector<complex<double>>> 
 */
vector<vector<complex<double>>> tmm::TransferMatrix::interfaceMatrix_s(Layer  &  layer1, Layer  & layer2) {
    // Material parameters
    complex<double> eps1 = layer1.getMaterial().getEpsilon();
    complex<double> mu1 = layer1.getMaterial().getMu();
    complex<double> k1 = sqrt(eps1*mu1)*k0;

    complex<double> k1z = sqrt(norm(k1) - kx*kx); // norm(z) = ||z||^2

    complex<double> eps2 = layer2.getMaterial().getEpsilon();
    complex<double> mu2 = layer2.getMaterial().getMu();
    complex<double> k2 = sqrt(eps2*mu2)*k0;
    complex<double> k2z = sqrt(norm(k2) - kx*kx); // norm(z) = ||z||^2

    complex<double> normImpedance = k2z*mu1/(k1z*mu2); // normalized impedance

    complex<double> M11 = 0.5 + 0.5*normImpedance; 
    complex<double> M12 = 0.5 - 0.5*normImpedance; 
    complex<double> M21 = 0.5 - 0.5*normImpedance;
    complex<double> M22 = 0.5 + 0.5*normImpedance;

    vector<vector<complex<double>>> Ms; 
    Ms.push_back({M11,M12});
    Ms.push_back({M21,M22});

    return Ms;   
}

/**
 * @brief Interface matrix for p polarization.
 * 
 * @param layer1 
 * @param layer2 
 * @return vector<vector<complex<double>>> 
 */
vector<vector<complex<double>>> tmm::TransferMatrix::interfaceMatrix_p(Layer  & layer1, Layer  & layer2) {
    // Material parameters
    complex<double> eps1 = layer1.getMaterial().getEpsilon();
    complex<double> mu1 = layer1.getMaterial().getMu();
    complex<double> k1 = sqrt(eps1*mu1)*k0;

    complex<double> k1z = sqrt(norm(k1) - kx*kx); // norm(z) = ||z||^2

    complex<double> eps2 = layer2.getMaterial().getEpsilon();
    complex<double> mu2 = layer2.getMaterial().getMu();
    complex<double> k2 = sqrt(eps2*mu2)*k0;
    complex<double> k2z = sqrt(norm(k2) - kx*kx); // norm(z) = ||z||^2

    complex<double> normImpedance = k2z*eps1/(k1z*eps2);

    complex<double> M11 = 0.5 + 0.5*normImpedance; 
    complex<double> M12 = 0.5 - 0.5*normImpedance; 
    complex<double> M21 = 0.5 - 0.5*normImpedance;
    complex<double> M22 = 0.5 + 0.5*normImpedance;

    vector<vector<complex<double>>> Mp; 
    Mp.push_back({M11,M12});
    Mp.push_back({M21,M22});

    return Mp;   
}

/**
 * @brief Propagation matrix.
 * 
 * @param layer 
 * @param distance 
 * @return vector<vector<complex<double>>> 
 */
vector<vector<complex<double>>> tmm::TransferMatrix::propagate(Layer  & layer, double distance) {
    //complex<double> kz = layer.getKz();

    complex<double> eps1 = layer.getMaterial().getEpsilon();
    complex<double> mu1 = layer.getMaterial().getMu();
    complex<double> k1 = sqrt(eps1*mu1)*k0;

    complex<double> kz = sqrt(norm(k1) - kx*kx); // norm(z) = ||z||^2
    layer.setKz(kz);

    complex<double> imagI = complex<double>{0.0,1.0};
    vector<vector<complex<double>>> P = {{exp(-imagI*kz*distance),0},{0,exp(+imagI*kz*distance)}};
    return P;

}

