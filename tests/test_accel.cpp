#include <iostream>
using std::cout;
using std::endl;

#include "model/model.hpp"
#include "tools/timer.hpp"


int main (int argc, char **argv)
{
    Model model;

    model.set();
    model.add();

    cout << "Done." << endl;

    return (0);
}
