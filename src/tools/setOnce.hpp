#pragma once


#include <exception>
#include <iostream>

#include "paracabs.hpp"


struct DoubleSetException : public std::exception
{
	const char* what () const throw ()
    {
	    return "Tried to overwrite SetOnce object with another value.";
    }
};


struct GetBeforeSetException : public std::exception
{
	const char* what () const throw ()
    {
    	return "Tried to get SetOnce object before setting it.";
    }
};


template <typename type>
class SetOnce
{
    private:
        bool already_set = false;
        type value;

    public:
        inline SetOnce () {}
        inline SetOnce (const type new_value): value (new_value) {}
        inline SetOnce (const SetOnce& s):
            already_set (s.already_set),
            value       (s.value) {};

        inline void set (const type new_value)
        {
            if (already_set)
            {
                if (value != new_value)
                {
                    std::cout << "ERROR value = " << value << " new value = " << new_value << std::endl;
                    throw DoubleSetException ();
                }
            }
            else
            {
                already_set = true;
                value       = new_value;
            }
        }


        accel inline type get () const
        {
            return value;
//            if (already_set) {return value;                  }
//            else             {throw GetBeforeSetException ();}
        }


};
