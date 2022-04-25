#pragma once

template<class T>
class Array1D {
    T* arr;
public:
    size_t len;
    Array1D(size_t len = 0)
    : arr(nullptr)
    {
        if(len == 0) return;
        this->len = len;
        arr = new T[len];
    }
    
    constexpr T& at(size_t i) const {
        return arr[i];
    }
    ~Array1D() {
        delete[] arr;
        arr = nullptr;
    }
};
