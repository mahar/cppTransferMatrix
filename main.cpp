#include <iostream>
#include <vector> 
#include "tmm.h" 

using namespace std;
using namespace tmm;


int main() {
    Material  air =  Material();
    Material Si = Material(12.0);
    Layer  l1 =  Layer(1e-6, air,"haha");
    Layer  l2 =  Layer(1e-6, Si,"Si");
    Layer l3 = Layer(1e-6, air);
    vector<Layer> structure{l1,l2,l3};

    double frequency = 2.0*3.14*1e12;
    TransferMatrix * system = new  TransferMatrix(frequency,0.0,structure);


    cout << " rs = " << system->getRs() << endl;

    delete system;


    return 0; 
}