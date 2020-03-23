#ifndef DifferentiableFunction_h
#define DifferentiableFunction_h

#include <functional>

struct DifferentiableFunction {
    
    //Function (must be differentiable)
    std::function<double(double)> func;

    //Derivative
    std::function<double(double)> d_func;
};

#endif