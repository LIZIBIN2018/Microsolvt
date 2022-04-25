#pragma once

#include <iostream>

template<int dim, class T>
class Array {
    T* arr;
public:
    size_t len[dim];
    size_t size;
    Array(size_t *len)
    : arr(nullptr),size(1)
    {
        for(int d = 0;d < dim; d++ ) // current length
        {
            this->len[d] = len[d];
            size        *= len[d];
            if (len[d] == 0)
            {
                size = 0;
                return;
            }
        }
        arr = new T[size];
    }
    
    template<class ... Args>
    constexpr T& at(Args ... args) const
    {
        return T[index(args)];
    }

    template<class ... Args>
    int index(Args ... args) const
    {
        return index(args, 0, 0, 1);
    }

    template<class ... Args>
    int index(const size_t &i, Args ... args, size_t ind_cur, size_t d_cur, size_t multiple) const
    {
        if(d == dim)
            return ind_cur;
        ind_cur += i*multiple;
        multiple *= len[d_cur];
        d_cur++;
        return index(args, inc_cur, d_cur, multiple);
    }

    ~Array() {
        delete arr;
        arr = nullptr;
    }
};


template<class T>
using Array1D = Array<1,T>;

template<class T>
using Array2D = Array<2,T>;

template<class T>
using Array3D = Array<3,T>;