#include <iostream>
#include <vector> 
#include "tmm.h" 

using namespace std;
using namespace tmm;


int main() {
    const double pi = 3.1415926534;

    Material  air =  Material();
    Material Si = Material(12.0);
    Layer  l1 =  Layer(1e-6, air,"haha");
    Layer  l2 =  Layer(1e-6, Si,"Si");
    Layer l3 = Layer(1e-6, air);
    vector<Layer> structure{l1,l2,l3};

    long double frequency = 2.0*pi*1e12;
    TransferMatrix * system = new  TransferMatrix(frequency,0.0,structure);

    system->calculate();

    cout.precision(14);
    cout << "Frequency = " << frequency << endl;
    cout << " rs = " << system->getRs() << endl;
    cout << " rp = " << system->getRp() << endl;
    cout << " ts = " << system->getTs() << endl;
    cout << " tp = " << system->getTp() << endl;

    cout << "interface matrix";
    // Elementary interface matrices
    auto vs = system->interfaceMatrix_s(l1,  l2);
    auto vp = system->interfaceMatrix_p(l1, l2);

    for (auto v : vp) { 
        for (auto el : v) { 
            std::cout << el << " - ";
        }
        std::cout << endl;
    }
   


    delete system;


    return 0; 
}