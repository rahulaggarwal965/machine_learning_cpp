#ifndef DeepNeuralNetwork_hpp
#define DeepNeuralNetwork_hpp

#include <vector>
#include <string>
#include "Matrix.hpp"
#include "DifferentiableFunction.h"


class DeepNeuralNetwork {
    private:
        //Layer
        size_t nLayers;
        std::vector<MatrixD> weights;
        std::vector<MatrixD> biases;

       //Activation Function
       DifferentiableFunction activation_func;

       //Learning Rate
       double learning_rate;

    public: 
        
        DeepNeuralNetwork(const std::vector<int>& shape, const DifferentiableFunction& func, double learning_rate);
        DeepNeuralNetwork(const std::vector<int>& shape, const DifferentiableFunction& func);

        MatrixD predict(const MatrixD& input) const;
        void learn(const MatrixD& input, const MatrixD& output, std::vector<MatrixD>& dWeights, std::vector<MatrixD>& dBiases);
};

#endif