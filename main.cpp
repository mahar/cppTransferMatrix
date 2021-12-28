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
    vector<Layer> st{l1,l2};


 
    TransferMatrix * system = new  TransferMatrix(1.0,0.0,st);
    system->calculate();

    cout << "Layer name : " << l1.getName() << endl;
    cout << "material eps : " << air.getEpsilon() << endl;

    for (auto v : system->Ms) {
        for (auto u : v) {
            cout << u << " - "; 
        }
        cout << endl;
    }

    cout << " rs = " << system->getTs() << endl;

    delete system;

  

   




    return 0; 
}