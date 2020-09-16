#pragma once


#include "tools/types.hpp"


///  Io: Magritte's input/output interface
//////////////////////////////////////////

struct Io
{
    const string io_file;

    // Constructor
    Io (const string &io_file): io_file (io_file) {};

    virtual int  read_length   (const string fname,       Size    &length) const = 0;
    virtual Size  get_length   (const string fname                       ) const = 0;

    virtual int  read_width    (const string fname,       Size    &width ) const = 0;
    virtual Size  get_width    (const string fname                       ) const = 0;

    virtual int  read_number   (const string fname,       long    &number) const = 0;
    virtual int write_number   (const string fname, const long    &number) const = 0;

    virtual int  read_number   (const string fname,       Size    &number) const = 0;
    virtual int write_number   (const string fname, const Size    &number) const = 0;

    virtual int  read_number   (const string fname,       Real    &number) const = 0;
    virtual int write_number   (const string fname, const Real    &number) const = 0;

    virtual int  read_word     (const string fname,       string  &word  ) const = 0;
    virtual int write_word     (const string fname, const string  &word  ) const = 0;

    virtual int  read_bool     (const string fname,       bool    &value ) const = 0;
    virtual int write_bool     (const string fname, const bool    &value ) const = 0;

    virtual int  read_list     (const string fname,       Long1   &list  ) const = 0;
    virtual int write_list     (const string fname, const Long1   &list  ) const = 0;

    virtual int  read_list     (const string fname,       Double1 &list  ) const = 0;
    virtual int write_list     (const string fname, const Double1 &list  ) const = 0;

    virtual int  read_list     (const string fname,       Size_t1 &list  ) const = 0;
    virtual int write_list     (const string fname, const Size_t1 &list  ) const = 0;

    virtual int  read_list     (const string fname,       Real1   &list  ) const = 0;
    virtual int write_list     (const string fname, const Real1   &list  ) const = 0;

    virtual int  read_list     (const string fname,       Size1   &list  ) const = 0;
    virtual int write_list     (const string fname, const Size1   &list  ) const = 0;

    virtual int  read_list     (const string fname,       String1 &list  ) const = 0;
    virtual int write_list     (const string fname, const String1 &list  ) const = 0;

    virtual int  read_array    (const string fname,       Long2   &array ) const = 0;
    virtual int write_array    (const string fname, const Long2   &array ) const = 0;

    virtual int  read_array    (const string fname,       Double2 &array ) const = 0;
    virtual int write_array    (const string fname, const Double2 &array ) const = 0;

    virtual int  read_array    (const string fname,       Real2   &array ) const = 0;
    virtual int write_array    (const string fname, const Real2   &array ) const = 0;

    virtual int  read_3_vector (const string fname,       Double1 &x,
                                                          Double1 &y,
                                                          Double1 &z     ) const = 0;
    virtual int write_3_vector (const string fname, const Double1 &x,
                                                    const Double1 &y,
                                                    const Double1 &z     ) const = 0;


    int read_list (const string fname, Vector<Size>& v) const
    {
        int err = read_list (fname, v.vec);
        v.set_dat ();
        return err;
    }


    int write_list (const string fname, const Vector<Size>& v) const
    {
        return write_list (fname, v.vec);
    }


    int read_list (const string fname, Vector<Real>& v) const
    {
        int err = read_list (fname, v.vec);
        v.set_dat ();
        return err;
    }


    int write_list (const string fname, const Vector<Real>& v) const
    {
        return write_list (fname, v.vec);
    }
};