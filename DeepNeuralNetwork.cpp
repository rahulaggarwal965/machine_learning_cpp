#include "DeepNeuralNetwork.hpp"

DeepNeuralNetwork::DeepNeuralNetwork(const std::vector<int>& shape, const DifferentiableFunction& activation_func, double learning_rate) :
    nLayers(shape.size() - 1),
    activation_func(activation_func) ,
    learning_rate(learning_rate)
{
    weights.reserve(nLayers);
    biases.reserve(nLayers);
    for (size_t i = 0; i < this->nLayers; i++) {
        weights.emplace_back(shape[i + 1], shape[i]);
        biases.emplace_back(shape[i+1], 1);
    }
}

DeepNeuralNetwork::DeepNeuralNetwork(const std::vector<int>& shape, const DifferentiableFunction& activation_func) :
    DeepNeuralNetwork(shape, activation_func, 0.1) {}

MatrixD DeepNeuralNetwork::predict(const MatrixD& input) const {
    if(input.nCols() != 1 || input.nRows() != this->weights[0].nCols()) 
        throw std::invalid_argument("Invalid Input Shape: Entered = " + std::to_string(input.nRows()) + "x" + std::to_string(input.nCols()) + "Expected = " + std::to_string(this->weights[0].nCols()) + "x1");
    MatrixD outputs = input;
    for (int i = 0; i < this->nLayers; i++) {
        outputs = this->weights[i] * outputs + this->biases[i];
        outputs.map(this->activation_func.func);
    }
    return outputs; //CAN I GET RVO (fuck me if i get an extra copy)
}

void DeepNeuralNetwork::learn(const MatrixD& input, const MatrixD& output, std::vector<MatrixD>& dWeights, std::vector<MatrixD>& dBiases) {
    if(input.nCols() != 1 || input.nRows() != this->weights[0].nCols()) 
        throw std::invalid_argument("Invalid Input Shape: Entered = " + std::to_string(input.nRows()) + "x" + std::to_string(input.nCols()) + "Expected = " + std::to_string(this->weights[0].nCols()) + "x1");
    if(output.nCols() != 1 || input.nRows() != this->weights[this->nLayers - 1].nCols())
         throw std::invalid_argument("Invalid Output Shape: Entered = " + std::to_string(output.nRows()) + "x" + std::to_string(output.nCols()) + "Expected = " + std::to_string(this->weights[this->nLayers - 1].nCols()) + "x1");

    std::vector<MatrixD> intermediates, activations;
    intermediates.reserve(this->nLayers);
    activations.reserve(this->nLayers + 1);

    activations.emplace_back(input); //Inefficient
    for (int i = 0; i < this->nLayers; i++) {
        const MatrixD& t = intermediates.emplace_back(this->weights[i] * activations[i] + this->biases[i]);
        activations.emplace_back(t.apply(this->activation_func.func));
    }

    Matrix delta = ((output - activations[this->nLayers]) * 2).hadamard(intermediates[this->nLayers - 1].apply(this->activation_func.d_func));
    dWeights.emplace_back(delta * activations[nLayers - 1].transpose());
    dBiases.emplace_back(delta);

    for (size_t i = this->nLayers - 1; i >= 0; i--) {
        delta = this->weights[i].transpose() * delta;
        dWeights.emplace_back(delta * activations[i - 1].transpose());
        dBiases.emplace_back(delta);
    }


    //BackProp


}

    