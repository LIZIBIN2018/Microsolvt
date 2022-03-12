#pragma once

template<class T>
class Array2D {
    T** arr;
public:
    size_t lenX, lenY;
    Array2D(size_t lenX = 0, size_t lenY = 0)
    : arr(nullptr)
    {
        if(lenX == 0 || lenY == 0) return;
        this->lenX = lenX;
        this->lenY = lenY;
        arr = new T*[lenX];
        for (size_t i = 0; i < lenX; i++) {
            arr[i] = new T[lenY];
        }
    }
    
    constexpr T& at(size_t i, size_t j) const {
        return arr[i][j];
    }
    ~Array2D() {
        for (size_t i = 0; i < lenX; i++) {
            delete[] arr[i];
        }
        delete[] arr;
        arr = nullptr;
    }
};
