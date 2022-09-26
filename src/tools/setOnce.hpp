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
        bool default_set = false;
        type value;
				string var_name;//should be set when setting the parameter

    public:
        inline SetOnce () {}
        inline SetOnce (const type new_value, const string name): value (new_value), var_name(name) {}
        inline SetOnce (const SetOnce& s):
            already_set (s.already_set),
            value       (s.value),
						var_name(s.var_name) {};

        inline void set (const type new_value, const string name)
        {
            if (already_set)
            {
                if (value != new_value)
                {
                    std::cout << "variable: " << var_name << " ERROR value = " << value << " new value = " << new_value << std::endl;
                    throw DoubleSetException ();
                }
            }
            else
            {
								var_name = name;
                already_set = true;
                value       = new_value;
            }
        }


        inline void set_default (const type new_value, const string name)
        {
            if (already_set)
            {
                if (value != new_value)
                {
                    std::cout << "variable: " << var_name << " ERROR value = " << value << " new value = " << new_value << std::endl;
                    throw DoubleSetException ();
                }
            }
            else
            {
								var_name = name;
                default_set = true;
                value       = new_value;
            }
        }


        accel inline type get () const
        {
            if (already_set || default_set) {return value;                  }
            else                            {throw GetBeforeSetException ();}
        }


};
