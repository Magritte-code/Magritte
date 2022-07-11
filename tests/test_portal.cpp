#include <iostream>
using std::cout;
using std::endl;

#include "model/model.hpp"
#include "solver/solver.hpp"



int main (int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Provide model file!" << endl;
        return (-1);
    }

    const string model_file = argv[1];

    cout << "Model file:" << endl;
    cout <<  model_file   << endl;

    Model model = Model(model_file);



    return (0);
}
