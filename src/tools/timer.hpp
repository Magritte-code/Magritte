#pragma once

#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <chrono>

/// TIMER: class for precise process timing
///////////////////////////////////////////
class singleTimer {

  private:
    std::chrono::high_resolution_clock::time_point start_ =
        std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point stop_ =
        std::chrono::high_resolution_clock::now();

  public:
    inline void start() { start_ = std::chrono::high_resolution_clock::now(); }
    inline void stop() { stop_ = std::chrono::high_resolution_clock::now(); }

    inline double get_interval() const {
        std::chrono::duration<double> interval = stop_ - start_;

        return interval.count();
    }
};

/// TIMER: class for precise process timing
///////////////////////////////////////////

class Timer {

  private:
    string name;

    vector<std::chrono::high_resolution_clock::time_point> starts;
    vector<std::chrono::high_resolution_clock::time_point> stops;

    vector<std::chrono::duration<double>> intervals;

    std::chrono::duration<double> total;

  public:
    ///  Constructor for TIMER
    //////////////////////////

    Timer(const string timer_name) { name = timer_name; }

    ///  start: start timer i.e. set initial time stamp
    ///////////////////////////////////////////////////

    void start() { starts.push_back(std::chrono::high_resolution_clock::now()); }

    ///  stop: stop timer and calculate interval for every process
    //////////////////////////////////////////////////////////////

    void stop() {
        stops.push_back(std::chrono::high_resolution_clock::now());

        const std::chrono::duration<double> interval = stops.back() - starts.back();

        intervals.push_back(interval);

        total += interval;
    }

    ///  print_to_file: print time interval to file
    ///////////////////////////////////////////////

    // void print_to_file ()
    //{
    //  string file_name = output_folder + "timer_" + name + ".txt";

    //	ofstream outFile (file_name, ios_base::app);

    //  outFile << interval.count() << endl;

    //	outFile.close();
    //}

    ///  print: print time interval to screen
    /////////////////////////////////////////

    string get_print_string() {
        return ("T   | " + name + " : " + to_string(intervals.back().count()) + " seconds");
    }

    string get_print_total_string() {
        return ("Tot | " + name + " : " + to_string(total.count()) + " seconds");
    }

    void print() { cout << get_print_string() << endl; }

    void print_total() { cout << get_print_total_string() << endl; }
};

#if (MAGRITTE_MPI_PARALLEL)

/// MPI_TIMER: class for precise process timing when using MPI
//////////////////////////////////////////////////////////////
//
// class MPI_TIMER
//{
//  private:
//
//    string name;
//
//    int world_size;
//    int world_rank;
//
//    chrono::duration <double> interval;
//    chrono::high_resolution_clock::time_point initial;
//
//	  vector<double> times;
//
//
//  public:
//
//
//	  ///  Constructor for TIMER
//	  //////////////////////////
//
//	  MPI_TIMER(string timer_name)
//	  {
//      name = timer_name;
//
//	  	MPI_Comm_size (MPI_COMM_WORLD, &world_size);
//	  	MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
//
//	  	times.resize (world_size);
//	  }
//
//
//	  ///  start: start timer i.e. set initial time stamp
//	  ///////////////////////////////////////////////////
//
//    void start ()
//    {
//      initial = chrono::high_resolution_clock::now();
//    }
//
//
//	  ///  stop: stop timer and calculate interval for every process
//	  //////////////////////////////////////////////////////////////
//
//    void stop ()
//    {
//      interval = chrono::high_resolution_clock::now() - initial;
//
//	  	times[world_rank] = interval.count();
//
//	  	MPI_Allgather (MPI_IN_PLACE,
//	  			           0,
//	  								 MPI_DATATYPE_NULL,
//	  								 times.data(),
//	  								 1,
//	  								 MPI_DOUBLE,
//	  								 MPI_COMM_WORLD);
//    }
//
//
//	  ///  print_to_file: let rank 0 process print times for every rank to file
//	  ///    @param[in] file_name: name of the file to print to
//	  /////////////////////////////////////////////////////////////////////////
//
//
//	  void print_to_file ()
//	  {
//	  	if (world_rank == 0)
//	  	{
//        string file_name = output_folder + "timer_" + name + ".txt";
//
//	  		ofstream outFile (file_name, ios_base::app);
//
//	  	  for (int w = 0; w < world_size; w++)
//	  	  {
//          outFile << world_size << "\t" << w << "\t" << times[w] << endl;
//	  	  }
//
//	  		outFile.close();
//	  	}
//	  }
//
//
//	  ///  print: let rank 0 process print times for every rank
//	  /////////////////////////////////////////////////////////
//
//    void print ()
//    {
//	  	if (world_rank == 0)
//	  	{
//        cout << "Timer " << name << ":" << endl;
//
//	  	  for (int w = 0; w < world_size; w++)
//	  	  {
//          cout << "rank[" << w << "]: time = " << times[w] << " seconds" << endl;
//	  	  }
//	  	}
//    }
//
//};

#endif
