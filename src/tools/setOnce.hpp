#pragma once

#include "paracabs.hpp"

#include <exception>
#include <iostream>

struct DoubleSetException : public std::exception {
    const char* what() const throw() { return "Tried to overwrite SetOnce object with another value."; }
};

struct GetBeforeSetException : public std::exception {
    const char* what() const throw() { return "Tried to get SetOnce object before setting it."; }
};

// essential parameters will not be derived, so must be set by the user
// themselves when creating a model. Note: currently does not actually do much
template <typename type, bool essential> class SetOnce {
  private:
    bool already_set = false;
    bool default_set = false;
    type value;
    const string var_name; // should be set when setting the parameter

  public:
    inline SetOnce() {}
    inline SetOnce(const string name) : var_name(name) {}
    inline SetOnce(const type new_value, const string name) : value(new_value), var_name(name) {}
    inline SetOnce(const SetOnce& s) : already_set(s.already_set), value(s.value), var_name(s.var_name){};

    inline void set(const type new_value) {
        if (already_set) {
            if (value != new_value) {
                std::cout << "variable: " << var_name << " ERROR value = " << value << " new value = " << new_value
                          << std::endl;
                throw DoubleSetException();
            }
        } else {
            already_set = true;
            value       = new_value;
        }
    }

    inline void set_default(const type new_value) {
        if (already_set) {
            if (value != new_value) {
                std::cout << "variable: " << var_name << " ERROR value = " << value << " new value = " << new_value
                          << std::endl;
                throw DoubleSetException();
            }
        } else {
            default_set = true;
            value       = new_value;
        }
    }

    accel inline type get() const {
        if (already_set || default_set) {
            return value;
        } else {
            throw GetBeforeSetException();
        }
    }
};
