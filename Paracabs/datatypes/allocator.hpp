#pragma once


#include <memory>


// References
// https://medium.com/@vgasparyan1995/what-is-an-allocator-c8df15a93ed


namespace paracabs
{
    template <class value_type, class memory_type>
    class allocator : public memory_type, public std::allocator<value_type>
    {
        public:
            allocator() = default;

            template <typename U>
            allocator (const allocator<U, memory_type>& other) : std::allocator<value_type>(other) {};

            template<typename U>
            struct rebind
            {
                using other = allocator<U, memory_type>;
            };

            value_type* allocate (size_t num)
            {
                return static_cast<value_type*> (this->malloc (num*sizeof(value_type)));
            }

            void deallocate (value_type* ptr, size_t num) noexcept
            {
                this->free (ptr);
            }
    };


//    template <class T, class U>
//    bool operator== (allocator<T> const& x, allocator<U> const& y) noexcept
//    {
//        return true;
//    }


//    template <class T, class U>
//    bool operator!= (allocator<T> const& x, allocator<U> const& y) noexcept
//    {
//        return !(x == y);
//    }
}
