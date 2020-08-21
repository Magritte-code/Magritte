#include <fstream>
#include <sys/stat.h>
#include <iomanip>

#include "io_cpp_text.hpp"


///  Constructor for IoText
///    @param[in] io_file : file to read from and write to
//////////////////////////////////////////////////////////

IoText :: IoText (const string &io_file) : Io (io_file)
{

}   // END OF CONSTRUCTOR




///  Check if the given path exists
///    @param[in] path : path to check
//////////////////////////////////////

bool pathExist (const string &path)
{
  struct stat buffer;

  return (stat (path.c_str(), &buffer) == 0);
}




///  Reader for the length of a text file
///    @param[in]  file_name : path to file containing the data
///    @param[out] length    : length to be read
///////////////////////////////////////////////////////////////

int IoText :: read_length (const string file_name, Size &length) const
{
  string fname = io_file + file_name;

  length = 0;

  if (pathExist (fname + ".txt"))
  {
    std::ifstream file (fname + ".txt");

    string line;

    while (std::getline (file, line))
    {
      length++;
    }

    file.close();
  }

  else
  {
    while (pathExist (fname + std::to_string (length)))
    {
      length++;
    }
  }

  return (0);
}




///  Reader for the number of columns (width) of a text file
///  or the number of files with a similar file name
///    @param[in]  file_name : path to file containing the data
///    @param[out] width     : width to be read
///////////////////////////////////////////////////////////////

int IoText :: read_width (const string file_name, Size &width) const
{
  string fname = io_file + file_name;

  width = 0;

  if (pathExist (fname + ".txt"))
  {
    // Count the number of columns in the file

    std::ifstream file (io_file + file_name + ".txt");

    string line, elem;

    std::getline (file, line);

    std::stringstream ss (line);

    while (std::getline (ss, elem, '\t'))
    {
      width++;
    }

    file.close();
  }

  else
  {
    // Count the number of files with a similar name

    while (pathExist (fname + std::to_string (width)))
    {
      width++;
    }
  }

  return (0);
}




///  Reader for a single (long integer) number from a text file
///    @param[in]  file_name : path to the file containing the number
///    @param[out] number    : number to be read
/////////////////////////////////////////////////////////////////////

int IoText :: read_number (const string file_name, size_t &number) const
{
    std::ifstream file (io_file + file_name + ".txt");

    if (file.is_open())
    {
        file >> number;

        return ( 0);
    }
    else
    {
        return (-1);
    }
}




///  Writer for a single (long integer) number to a text file
///    @param[in]  file_name : path to the file to be written
///    @param[out] number    : number to be written
/////////////////////////////////////////////////////////////

int IoText :: write_number (const string file_name, const size_t &number) const
{
    std::ofstream file (io_file + file_name + ".txt");

    if (file.is_open())
    {
        file << std::scientific << std::setprecision(16);
        file << number;

        return ( 0);
    }
    else
    {
        return (-1);
    }
}




///  Reader for a single (long integer) number from a text file
///    @param[in]  file_name : path to the file containing the number
///    @param[out] number    : number to be read
/////////////////////////////////////////////////////////////////////

int IoText :: read_number (const string file_name, long &number) const
{
    std::ifstream file (io_file + file_name + ".txt");

    if (file.is_open())
    {
        file >> number;

        return ( 0);
    }
    else
    {
        return (-1);
    }
}




///  Writer for a single (long integer) number to a text file
///    @param[in]  file_name : path to the file to be written
///    @param[out] number    : number to be written
/////////////////////////////////////////////////////////////

int IoText :: write_number (const string file_name, const long &number) const
{
    std::ofstream file (io_file + file_name + ".txt");

    if (file.is_open())
    {
        file << std::scientific << std::setprecision (16);
        file << number;

        return ( 0);
    }
    else
    {
        return (-1);
    }
}




///  Reader for a single (double) number from a text file
///    @param[in]  file_name : path to the file containing the number
///    @param[out] number    : number to be read
/////////////////////////////////////////////////////////////////////

int IoText :: read_number (const string file_name, double &number) const
{
    std::ifstream file (io_file + file_name + ".txt");

    if (file.is_open())
    {
        file >> number;

        return ( 0);
    }
    else
    {
        return (-1);
    }
}




///  Writer for a single (double) number to a text file
///    @param[in]  file_name : path to the file to be written
///    @param[out] number    : number to be written
/////////////////////////////////////////////////////////////

int IoText :: write_number (const string file_name, const double &number) const
{
    std::ofstream file (io_file + file_name + ".txt");

    if (file.is_open())
    {
        file << std::scientific << std::setprecision (16);
        file << number;

        return ( 0);
    }
    else
    {
        return (-1);
    }
}




///  Reader for a single string from a text file
///    @param[in]  file_name : path to the file containing the string
///    @param[out] word      : string to be read
/////////////////////////////////////////////////////////////////////

int IoText :: read_word (const string file_name, string &word) const
{
    std::ifstream file (io_file + file_name + ".txt");

    if (file.is_open())
    {
        file >> word;

        return ( 0);
    }
    else
    {
        return (-1);
    }
}




///  Writer for a single string to a text file
///    @param[in]  file_name : path to the file to be written
///    @param[out] word      : string to be written
/////////////////////////////////////////////////////////////

int IoText :: write_word (const string file_name, const string &word) const
{
    std::ofstream file (io_file + file_name + ".txt");

    if (file.is_open())
    {
        file << word;

        return ( 0);
    }
    else
    {
        return (-1);
    }
}




///  Reader for a single boolean from a text file
///    @param[in]  file_name : path to the file containing the boolean
///    @param[out] value     : value to be read
//////////////////////////////////////////////////////////////////////

int IoText :: read_bool (const string file_name, bool &value) const
{
    // Treat booleans as text in io
    string word;

    int err = read_word (file_name, word);

    if      (word.compare("true" ) == 0) {value = true; }
    else if (word.compare("false") == 0) {value = false;}
    else                                 {  err = -1;   }

    return err;
}




///  Writer for a single boolean to a text file
///    @param[in]  file_name : path to the file to be written
///    @param[out] value     : value to be written
/////////////////////////////////////////////////////////////

int IoText :: write_bool (const string file_name, const bool &value) const
{
  // Treat booleans as text in io
  string word = "false";

  if (value) {word = "true";}

  return write_word (file_name, word);
}




///  Reader for a list of long integers from a text file
///     @param[in] file_name : path to file containing the list
///     @param[in] list      : list to be read
///////////////////////////////////////////////////////////////

int IoText :: read_list (const string file_name, Long1 &list) const
{
  std::ifstream file (io_file + file_name + ".txt");

  long n = 0;

  while (file >> list[n])
  {
    n++;
  }

  file.close();

  return (0);
}




///  Writer for a list of long integers to a text file
///    @param[in] file_name : path to file to be written
///    @param[in] list      : list to be written
////////////////////////////////////////////////////////

int IoText :: write_list (const string file_name, const Long1 &list) const
{
  std::ofstream file (io_file + file_name + ".txt");

  file << std::scientific << std::setprecision (16);

  for (long n = 0; n < list.size(); n++)
  {
    file << list[n] << endl;
  }

  file.close();

  return (0);
}




///  Reader for a list of doubles from a text file
///    @param[in] file_name : path to file containing the list
///    @param[in] list      : list to be read
//////////////////////////////////////////////////////////////

int IoText :: read_list (const string file_name, Double1 &list) const
{
  std::ifstream file (io_file + file_name + ".txt");

  long   n = 0;

  while (file >> list[n])
  {
    n++;
  }

  file.close();

  return (0);
}




///  Writer for a list of strings to a text file
///    @param[in] file_name : path to file to be written
///    @param[in] list      : list to be written
////////////////////////////////////////////////////////

int IoText :: write_list (const string file_name, const Double1 &list) const
{
  std::ofstream file (io_file + file_name + ".txt");

  file << std::scientific << std::setprecision (16);

  for (long n = 0; n < list.size(); n++)
  {
    file << list[n] << endl;
  }

  file.close();

  return (0);
}




///  Reader for a list of doubles from a text file
///    @param[in] file_name : path to file containing the list
///    @param[in] list      : list to be read
//////////////////////////////////////////////////////////////

int IoText :: read_list (const string file_name, Size_t1 &list) const
{
    std::ifstream file (io_file + file_name + ".txt");

    long   n = 0;

    while (file >> list[n])
    {
        n++;
    }

    file.close();

    return (0);
}




///  Writer for a list of strings to a text file
///    @param[in] file_name : path to file to be written
///    @param[in] list      : list to be written
////////////////////////////////////////////////////////

int IoText :: write_list (const string file_name, const Size_t1 &list) const
{
    std::ofstream file (io_file + file_name + ".txt");

    file << std::scientific << std::setprecision (16);

    for (long n = 0; n < list.size(); n++)
    {
        file << list[n] << endl;
    }

    file.close();

    return (0);
}




///  Reader for a list of doubles from a text file
///    @param[in] file_name : path to file containing the list
///    @param[in] list      : list to be read
//////////////////////////////////////////////////////////////

int IoText :: read_list (const string file_name, Size1 &list) const
{
    std::ifstream file (io_file + file_name + ".txt");

    long   n = 0;

    while (file >> list[n])
    {
        n++;
    }

    file.close();

    return (0);
}




///  Writer for a list of strings to a text file
///    @param[in] file_name : path to file to be written
///    @param[in] list      : list to be written
////////////////////////////////////////////////////////

int IoText :: write_list (const string file_name, const Size1 &list) const
{
    std::ofstream file (io_file + file_name + ".txt");

    file << std::scientific << std::setprecision (16);

    for (long n = 0; n < list.size(); n++)
    {
        file << list[n] << endl;
    }

    file.close();

    return (0);
}




///  Reader for a list of doubles from a text file
///    @param[in] file_name : path to file containing the list
///    @param[in] list      : list to be read
//////////////////////////////////////////////////////////////

int IoText :: read_list (const string file_name, Real1 &list) const
{
    std::ifstream file (io_file + file_name + ".txt");

    long   n = 0;

    while (file >> list[n])
    {
        n++;
    }

    file.close();

    return (0);
}




///  Writer for a list of strings to a text file
///    @param[in] file_name : path to file to be written
///    @param[in] list      : list to be written
////////////////////////////////////////////////////////

int IoText :: write_list (const string file_name, const Real1 &list) const
{
    std::ofstream file (io_file + file_name + ".txt");

    file << std::scientific << std::setprecision (16);

    for (long n = 0; n < list.size(); n++)
    {
        file << list[n] << endl;
    }

    file.close();

    return (0);
}




///  Reader for a list of strings from a text file
///    @param[in] file_name : path to file containing the list
///    @param[in] list      : list to be read
//////////////////////////////////////////////////////////////

int IoText :: read_list (const string file_name, String1 &list) const
{
  std::ifstream file (io_file + file_name + ".txt");

  long   n = 0;

  while (file >> list[n])
  {
    n++;
  }

  file.close();

  return (0);
}




///  Writer for a list of strings to a text file
///    @param[in] file_name : path to file to be written
///    @param[in] list      : list to be written
////////////////////////////////////////////////////////

int IoText :: write_list (const string file_name, const String1 &list) const
{
  std::ofstream file (io_file + file_name + ".txt");

  file << std::scientific << std::setprecision (16);

  for (long n = 0; n < list.size(); n++)
  {
    file << list[n] << endl;
  }

  file.close();

  return (0);
}




///  Reader for an array of long integers from a text file
///    @param[in] file_name : path to file containing the array
///    @param[in] array     : array to be read
///////////////////////////////////////////////////////////////

int IoText :: read_array (const string file_name, Long2 &array) const
{
  std::ifstream file (io_file + file_name + ".txt");

  string line;

  for (long n1 = 0; n1 < array.size(); n1++)
  {
    std::getline (file, line);

    std::stringstream ss (line);

    for (long n2 = 0; n2 < array[n1].size(); n2++)
    {
      ss >> array[n1][n2];
    }
  }

  file.close();

  return (0);
}




///  Writer for an array of long integers to a text file
///    @param[in] file_name : path to file to be written
///    @param[in] array     : array to be written
////////////////////////////////////////////////////////

int IoText :: write_array (const string file_name, const Long2 &array) const
{
  std::ofstream file (io_file + file_name + ".txt");

  file << std::scientific << std::setprecision (16);

  for (long n1 = 0; n1 < array.size(); n1++)
  {
    for (long n2 = 0; n2 < array[n1].size(); n2++)
    {
      file << array[n1][n2] << "\t";
    }

    file << endl;
  }

  file.close();

  return (0);
}




///  Reader for an array of doubles from a text file
///    @param[in] file_name : path to file containing the array
///    @param[in] array     : array to be read
///////////////////////////////////////////////////////////////

int IoText :: read_array (const string file_name, Double2 &array) const
{
  std::ifstream file (io_file + file_name + ".txt");

  string line;

  for (long n1 = 0; n1 < array.size(); n1++)
  {
    std::getline (file, line);

    std::stringstream ss (line);

    for (long n2 = 0; n2 < array[n1].size(); n2++)
    {
      ss >> array[n1][n2];
    }
  }

  file.close();

  return (0);
}




///  Writer for an array of doubles from a text file
///    @param[in] file_name : path to file to be written
///    @param[in] array     : array to be written
////////////////////////////////////////////////////////

int IoText :: write_array (const string file_name, const Double2 &array) const
{
  std::ofstream file (io_file + file_name + ".txt");

  file << std::scientific << std::setprecision (16);

  for (long n1 = 0; n1 < array.size(); n1++)
  {
    for (long n2 = 0; n2 < array[n1].size(); n2++)
    {
      file << array[n1][n2] << "\t";
    }

    file << endl;
  }

  file.close();

  return (0);
}




///  Reader for a list of 3-vectors of doubles from a text file
///    @param[in] file_name : path to file containing the vectors
///    @param[in] x         : x component of the vector to be read
///    @param[in] y         : y component of the vector to be read
///    @param[in] z         : z component of the vector to be read
//////////////////////////////////////////////////////////////////

int IoText :: read_3_vector (
        const string   file_name,
              Double1 &x,
              Double1 &y,
              Double1 &z         ) const
{
  std::ifstream file (io_file + file_name + ".txt");

  long   n = 0;

  while (file >> x[n] >> y[n] >> z[n])
  {
    n++;
  }

  file.close();

  return (0);
}




///  Writer for a list of 3-vectors of doubles to a text file
///    @param[in] file_name : path to file containing the vectors
///    @param[in] x         : x component of the vector to be written
///    @param[in] y         : y component of the vector to be written
///    @param[in] z         : z component of the vector to be written
/////////////////////////////////////////////////////////////////////

int IoText :: write_3_vector (
        const string   file_name,
        const Double1 &x,
        const Double1 &y,
        const Double1 &z         ) const
{
  const long length = x.size();

  // Check if all 3 vectors are the same size
  if (   (length != y.size())
      || (length != z.size()) )
  {
    return (-1);
  }

  std::ofstream file (io_file + file_name + ".txt");

  file << std::scientific << std::setprecision (16);

  for (long n = 0; n < length; n++)
  {
    file << x[n] << "\t" << y[n] << "\t" << z[n] << endl;
  }

  file.close();

  return (0);
}
