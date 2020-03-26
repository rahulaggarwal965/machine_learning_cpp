#ifndef Matrix_hpp
#define Matrix_hpp

#include <vector>
#include <string>
#include <stdexcept>
#include <ostream>
#include <iomanip>
#include <functional>
#include <type_traits>

template <typename T>
class Matrix {
    private: 
        //FIELDS
        size_t rows;
        size_t cols;
        std::vector<T> data;

    public: 
        //CONSTRUCTORS
        Matrix() {
            rows = cols = 0;
        }

        Matrix(size_t dimension) {
            Matrix(dimension, dimension);
        }

        Matrix(size_t rows, size_t cols) :
            rows(rows),
            cols(cols), 
            data(rows * cols) {
        } 

	//Rule of Five
	~Matrix() = default;
	Matrix(const Matrix& m) = default;
	Matrix& operator=(const Matrix& m) = default;
	Matrix(Matrix&& m)  = default;
	Matrix& operator=(Matrix&& m) = default;

        size_t nRows() const {
            return this->rows;
        }

        size_t nCols() const {
            return this->cols;
        }

        Matrix(size_t rows, size_t cols, const std::vector<T>& data) :
            rows(rows),
            cols(cols) {
            if(data.size() != rows * cols)
                throw std::invalid_argument("Dimensions do not match arguments");
            this->data = data;
        }

        //INDEXING OPERATORS
        inline T operator()(size_t i, size_t j) const {
            return this->data[j * this->cols + i];
        }

        inline T& operator()(size_t i, size_t j) {
            return this->data[j * this->cols + i];
        }

        //ADDITION OPERATOR
        Matrix& operator+=(T value) {
            for (size_t i = 0; i < this->data.size(); i++) {
                this->data[i] += value;
            }
            return *this;
        }
        friend Matrix operator+(const Matrix& m, T value) {
            Matrix r(m.rows, m.cols);
            for (size_t i = 0; i < m.data.size(); i++) {
                r.data[i] = m.data[i] + value;
            }
            return r;
        }
        friend Matrix operator+(T value, const Matrix& m) {
            return m + value;
        }
        friend Matrix operator+(Matrix&& m, T value) {
            return m += value;
        }
        friend Matrix operator+(T value, Matrix&& m) {
            return m += value;
        }
        Matrix& operator+=(const Matrix &m) {
             if(this->rows != m.rows || this->cols != m.cols)
                throw std::invalid_argument("Dimensions do not match: A = " + std::to_string(this->rows) + "x" + std::to_string(this->cols) + 
                    ", B = " + std::to_string(m.rows) + "x" + std::to_string(m.cols));
            for (size_t i = 0; i < this->data.size(); i++) {
                this->data[i] += m.data[i];
            }
            return *this;
        }
        friend Matrix operator+(const Matrix& m0, const Matrix&m1) {
            return Matrix(m0) += m1;
        }
        friend Matrix operator+(const Matrix& m0, Matrix&& m1) {
            return m1 += m0;
        }
        friend Matrix operator+(Matrix&& m0, const Matrix& m1) {
            return m0 += m1;
        }
        friend Matrix operator+(Matrix&& m0, Matrix&& m1) {
            return m0 += m1;
        }

        //SUBBTRACTION OPERATORS
        Matrix& operator-=(T value) {
            for (size_t i = 0; i < this->data.size(); i++) {
                this->data[i] -= value;
            }
            return *this;
        }
        friend Matrix operator-(const Matrix& m, T value) {
            Matrix r(m.rows, m.cols);
            for (size_t i = 0; i < m.data.size(); i++) {
                r.data[i] = m.data[i] - value;
            }
            return r;
        }
        friend Matrix operator-(T value, const Matrix& m) {
            Matrix r(m.rows, m.cols);
            for (size_t i = 0; i < m.data.size(); i++) {
                r.data[i] = value - m.data[i];
            }
            return r;
        }
        friend Matrix operator-(Matrix&& m, T value) {
            return m -= value;
        }
        Matrix& operator-=(const Matrix &m) {
             if(this->rows != m.rows || this->cols != m.cols)
                throw std::invalid_argument("Dimensions do not match: A = " + std::to_string(this->rows) + "x" + std::to_string(this->cols) + 
                    ", B = " + std::to_string(m.rows) + "x" + std::to_string(m.cols));
            for (size_t i = 0; i < this->data.size(); i++) {
                this->data[i] -= m.data[i];
            }
            return *this;
        }
        friend Matrix operator-(const Matrix& m0, const Matrix&m1) {
            return Matrix(m0) -= m1;
        }
        friend Matrix operator-(Matrix&& m0, const Matrix& m1) {
            return m0 -= m1;
        }
        friend Matrix operator-(Matrix&& m0, Matrix&& m1) {
            return m0 -= m1;
        }

        //MULTIPLICATION OPERATIONS
        Matrix& operator*=(T value) {
            for (size_t i = 0; i < this->data.size(); i++) {
                this->data[i] *= value;
            }
            return *this;
        }
        friend Matrix operator*(const Matrix& m, T value) {
            Matrix r(m.rows, m.cols);
            for (size_t i = 0; i < m.data.size(); i++) {
                r.data[i] = m.data[i] * value;
            }
            return r;
        }
        friend Matrix operator*(T value, const Matrix& m) {
            return m * value;
        }
        friend Matrix operator*(Matrix&& m, T value) {
            return m *= value;
        }
        friend Matrix operator*(T value, Matrix&& m) {
            return m *= value;
        }
        Matrix& operator*=(const Matrix &m) {
             if(this->rows != m.rows || this->cols != m.cols)
                throw std::invalid_argument("Dimensions do not match: A = " + std::to_string(this->rows) + "x" + std::to_string(this->cols) + 
                    ", B = " + std::to_string(m.rows) + "x" + std::to_string(m.cols));
            for (size_t i = 0; i < this->data.size(); i++) {
                this->data[i] *= m.data[i];
            }
            return *this;
        }
        friend Matrix operator*(const Matrix& m0, const Matrix& m1) {
            if(m0.cols != m1.rows)
              throw std::invalid_argument("Rows of A do not match Cols of B: A = " + std::to_string(m0.rows) + "x" + std::to_string(m0.cols) + 
                    ", B = " + std::to_string(m1.rows) + "x" + std::to_string(m1.cols));
            Matrix r(m0.rows, m1.cols);
            for (size_t i = 0; i < m0.rows; i++) {
                for (size_t j = 0; j < m1.cols; j++) {
                   T sum = 0;
                   for (int k = 0; k < m1.rows; k++) {
                       sum += m0(k, i) * m1(j, k);
                   }
                   r(i, j) = sum;
                }
            }
            return r;
        }
        //COME UP WITH BETTER SOLUTIONS FOR THIS
        /*
        friend Matrix operator*(const Matrix& m0, Matrix&& m1) {
            return m0 * m1;
        }
        friend Matrix operator*(Matrix&& m0, const Matrix& m1) {
            return m0 -= m1;
        }
        friend Matrix operator*(Matrix&& m0, Matrix&& m1) {
            return m0 -= m1;
        }
        */

        //Element-Wise multiplication
        Matrix hadamard(const Matrix &m) const {
            return Matrix(*this) *= m;
        }

        Matrix hadamard(Matrix&& m) const {
            return m *= (*this);
        }
        
        //DIVISION OPERATORS
        Matrix& operator/=(T value) {
            for (size_t i = 0; i < this->data.size(); i++) {
                this->data[i] /= value;
            }
            return *this;
        }
        friend Matrix operator/(const Matrix& m, T value) {
            Matrix r(m.rows, m.cols);
            for (size_t i = 0; i < m.data.size(); i++) {
                r.data[i] = m.data[i] / value;
            }
            return r;
        }
        friend Matrix operator/(T value, const Matrix& m) {
            Matrix r(m.rows, m.cols);
            for (size_t i = 0; i < m.data.size(); i++) {
                r.data[i] = value / m.data[i];
            }
            return r;
        }
        friend Matrix operator/(Matrix&& m, T value) {
            return m /= value;
        }
        Matrix& operator/=(const Matrix &m) {
             if(this->rows != m.rows || this->cols != m.cols)
                throw std::invalid_argument("Dimensions do not match: A = " + std::to_string(this->rows) + "x" + std::to_string(this->cols) + 
                    ", B = " + std::to_string(m.rows) + "x" + std::to_string(m.cols));
            for (size_t i = 0; i < this->data.size(); i++) {
                this->data[i] /= m.data[i];
            }
            return *this;
        }
        friend Matrix operator/(const Matrix& m0, const Matrix&m1) {
            return Matrix(m0) /= m1;
        }
        friend Matrix operator/(Matrix&& m0, const Matrix& m1) {
            return m0 /= m1;
        }
        friend Matrix operator/(Matrix&& m0, Matrix&& m1) {
            return m0 /= m1;
        }

        Matrix operator-() const {
            Matrix r(this->rows, this->cols);
            for (size_t i = 0; i < this->data.size(); i++) {
                r.data[i] = -this->data[i];
            }
            return r;
        }

        bool operator==(const Matrix &m) const {
            if(this->rows != m.rows || this->cols != m.cols) return false;

            for (size_t i = 0; i < this->data.size(); i++) {
                if(this->data[i] != m.data[i]) return false;
            }
            return true;
        }

        bool operator!=(const Matrix &m) const {
            return !(*this == m);
        }

        Matrix transpose() const {
            Matrix r(this->cols, this->rows);
            for (size_t i = 0; i < this->cols; i++) {
                for (size_t j = 0; j < this->rows; j++) {
                    r(j, i) = operator()(i, j);
                }
            }
            return r;
        }

        void map(const std::function<T(T)>& f) {
            std::transform(this->data.begin(), this->data.end(), this->data.begin(), f);
        }

        Matrix apply(const std::function<T(T)>& f) const {
            Matrix result(this->rows, this->cols);
            std::transform(this->data.begin(), this->data.end(), result.data.begin(), f);
            return result;
        }
 
        friend std::ostream& operator<<(std::ostream &os, const Matrix &m) {
            const int numWidth = 13;
            char fill = ' ';

            for (size_t i = 0; i < m.rows; i++) {
                for (size_t j = 0; j < m.cols; j++) {
                    os << std::left << std::setw(numWidth) << std::setfill(fill) << to_string(m(i, j));
                }
                os << std::endl;
            }

            return os;
        }
        

};

typedef Matrix<double> MatrixD;
typedef Matrix<int> MatrixI;

#endif
