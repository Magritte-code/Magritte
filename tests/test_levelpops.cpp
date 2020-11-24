#include <iostream>
using std::cout;
using std::endl;

#include "model/model.hpp"
#include "tools/timer.hpp"


int main (int argc, char **argv)
{
    const string modelName = argv[1];

    cout << "Running test_levelpops..."                              << endl;
    cout << "-------------------------"                              << endl;
    cout << "Model name: " << modelName                              << endl;
    cout << "n threads = " << pc::multi_threading::n_threads_avail() << endl;

    Model model (modelName);
    model.compute_spectral_discretisation ();
    model.compute_LTE_level_populations   ();
    model.compute_inverse_line_widths     ();

    Timer timer("compute level populations");
    timer.start();
    model.compute_level_populations (true, 100);
    timer.stop();
    timer.print();

    cout << "Done." << endl;

    return (0);
}
