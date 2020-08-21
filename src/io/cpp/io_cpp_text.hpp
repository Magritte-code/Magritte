#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"


///  IoText: io specified by text files
///////////////////////////////////////

struct IoText : public Io
{
    // Constructor
    IoText (const string &io_file);

    int  read_length   (const string fname,       Size    &length) const override;

    int  read_width    (const string fname,       Size    &width ) const override;

    int  read_number   (const string fname,       size_t  &number) const override;
    int write_number   (const string fname, const size_t  &number) const override;

    int  read_number   (const string fname,       long    &number) const override;
    int write_number   (const string fname, const long    &number) const override;

    int  read_number   (const string fname,       double  &number) const override;
    int write_number   (const string fname, const double  &number) const override;

    int  read_word     (const string fname,       string  &word  ) const override;
    int write_word     (const string fname, const string  &word  ) const override;

    int  read_bool     (const string fname,       bool    &value ) const override;
    int write_bool     (const string fname, const bool    &value ) const override;

    int  read_list     (const string fname,       Long1   &list  ) const override;
    int write_list     (const string fname, const Long1   &list  ) const override;

    int  read_list     (const string fname,       Double1 &list  ) const override;
    int write_list     (const string fname, const Double1 &list  ) const override;

    int  read_list     (const string fname,       Size_t1 &list  ) const override;
    int write_list     (const string fname, const Size_t1 &list  ) const override;

    int  read_list     (const string fname,       Real1   &list  ) const override;
    int write_list     (const string fname, const Real1   &list  ) const override;

    int  read_list     (const string fname,       Size1   &list  ) const override;
    int write_list     (const string fname, const Size1   &list  ) const override;

    int  read_list     (const string fname,       String1 &list  ) const override;
    int write_list     (const string fname, const String1 &list  ) const override;

    int  read_array    (const string fname,       Long2   &array ) const override;
    int write_array    (const string fname, const Long2   &array ) const override;

    int  read_array    (const string fname,       Double2 &array ) const override;
    int write_array    (const string fname, const Double2 &array ) const override;

    int  read_3_vector (const string fname,       Double1 &x,
                                                  Double1 &y,
                                                  Double1 &z     ) const override;
    int write_3_vector (const string fname, const Double1 &x,
                                            const Double1 &y,
                                            const Double1 &z     ) const override;
};