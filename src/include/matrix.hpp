#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <cassert>
#include <stdexcept>
#include <vector>

namespace mat {

template <typename T>
class row_view {
private:
    T* data_begin_;
    int size_;

public:
    row_view(T* data_begin, int size)
        : data_begin_(data_begin), size_(size) {};

    [[nodiscard]] int get_cnt() const {
        return size_;
    }

    double& operator()(int inx) {
        if(inx >= size_ || inx < 0)
            throw std::out_of_range("index out of range");

        return data_begin_[inx];
    }
};

template <typename T>
class col_view {
private:
    T* data_begin_;
    int delta_;
    int size_;

public:
    col_view(T* data_begin, int size, int delta)
        : data_begin_(data_begin), delta_(delta), size_(size) {};

    [[nodiscard]] int get_cnt() const {
        return size_;
    }

    double& operator()(int inx) {
        if(inx >= size_ || inx < 0)
            throw std::out_of_range("index out of range");

        return data_begin_[inx * delta_];
    }

    const double& operator()(int inx) const {
        if(inx >= size_ || inx < 0)
            throw std::out_of_range("index out of range");

        return data_begin_[inx * delta_];
    }
};

template <typename T>
class h_mat_view;

template <typename T>
class vector;

template <typename T>
class matrix {
private:
    std::vector<std::remove_const_t<T>> mat_;
    int rows_cnt_;
    int cols_cnt_;

public:
    //constructors
    explicit matrix() {}

    explicit matrix(int row, int col) {
        rows_cnt_ = row;
        cols_cnt_ = col;
        mat_.resize(row * col);
    }

    explicit matrix(int side) {
        rows_cnt_ = side;
        cols_cnt_ = side;
        mat_.resize(side * side);
    }

    explicit matrix(const h_mat_view<T>& other) {
        rows_cnt_ = other.get_row_size();
        cols_cnt_ = other.get_col_size();
        mat_.resize(rows_cnt_ * cols_cnt_);

        for(int i = 0; i < rows_cnt_; i++) {
            for(int j = 0; j < cols_cnt_; j++) {
                mat_[i * cols_cnt_ + j] = other(i, j);
            }
        }
    }

    explicit matrix(vector<T>& diag) {
        rows_cnt_ = diag.get_size();
        cols_cnt_ = diag.get_size();
        mat_.resize(diag.get_size() * diag.get_size());

        for(int i = 0; i < rows_cnt_; i++) {
            (*this)(i, i) = diag(i);
        }
    }

    matrix(const matrix &other) {
        rows_cnt_ = other.rows_cnt_;
        cols_cnt_ = other.cols_cnt_;
        mat_ = other.mat_;
    }

    matrix(matrix &&other) noexcept {
        rows_cnt_ = other.rows_cnt_;
        cols_cnt_ = other.cols_cnt_;
        mat_ = std::move(other.mat_);
    }

    matrix &operator=(const matrix &other) = default;
    matrix &operator=(matrix &&other) = default;

    //access operations
    [[nodiscard]] int get_row_size() const {
        return rows_cnt_;
    }

    [[nodiscard]] int get_col_size() const {
        return cols_cnt_;
    }

    double& operator()(int row, int col) {
        if(row < 0 || row >= rows_cnt_) {
            throw std::out_of_range("row index out of range");
        }

        if(col < 0 || col >= cols_cnt_) {
            throw std::out_of_range("col index out of range");
        }

        return mat_[row * cols_cnt_ + col];
    }

    const double& operator()(int row, int col) const {
        if(row < 0 || row >= rows_cnt_) {
            throw std::out_of_range("row index out of range");
        }

        if(col < 0 || col >= cols_cnt_) {
            throw std::out_of_range("col index out of range");
        }

        return mat_[row * cols_cnt_ + col];
    }

    row_view<T> rows(int row_inx) {
        if(row_inx >= rows_cnt_ || row_inx < 0)
            throw std::out_of_range("row index out of range");

        return {&mat_[row_inx * cols_cnt_], cols_cnt_};
    }

    row_view<const T> rows(int row_inx) const {
        if(row_inx >= rows_cnt_ || row_inx < 0)
            throw std::out_of_range("row index out of range");

        return {&mat_[row_inx * cols_cnt_], cols_cnt_};
    }

    col_view<T> cols(int col_inx) {
        if(col_inx >= cols_cnt_ || col_inx < 0)
            throw std::out_of_range("col index out of range");

        return {&mat_[col_inx], rows_cnt_, cols_cnt_};
    }

    col_view<const T> cols(int col_inx) const {
        if(col_inx >= cols_cnt_ || col_inx < 0)
            throw std::out_of_range("col index out of range");

        return {&mat_[col_inx], rows_cnt_, cols_cnt_};
    }

    h_mat_view<T> extract_cols(std::vector<int> const& col_indices) {
        h_mat_view<T> result(rows_cnt_);

        for(int const& col_indice : col_indices) {
            result.append(cols(col_indice));
        }

        return result;
    }

    h_mat_view<T> extract_all_cols() {
        h_mat_view<T> result(rows_cnt_);

        for(int i = 0; i < cols_cnt_; i++) {
            result.append(cols(i));
        }

        return result;
    }

    //BREAKS VIEWS
    void swap_rows(int i, int j) {
        for(int k = 0; k < cols_cnt_; k++) {
            std::swap(mat_[i * cols_cnt_ + k], mat_[j * cols_cnt_ + k]);
        }
    }

    //multiplication
    template <typename Matrix>
    matrix<T> operator*(Matrix& other) {
        if(get_col_size() != other.get_row_size()) {
            throw std::out_of_range("invalid dimentions");
        }

        matrix<T> result(get_row_size(), other.get_col_size());

        for(int i = 0; i < get_row_size(); i++) {
            for(int j = 0; j < other.get_col_size(); j++) {
                for(int k = 0; k < cols_cnt_; k++) {
                    result(i, j) += (*this)(i, k) * other(k, j);
                }
            }
        }

        return result;
    }

    vector<T> operator*(col_view<T> other) {
        if(cols_cnt_ != other.get_cnt()) {
            throw std::out_of_range("invalid dimentions");
        }

        vector<T> result(get_row_size());

        for(int i = 0; i < get_row_size(); i++) {
            for(int j = 0; j < other.get_cnt(); j++) {
                result(i) += (*this)(i, j) * other(j);
            }
        }

        return result;
    }

    //utils
    static matrix<T> E(int size) {
        matrix<double> result(size);

        for(int i = 0; i < size; i++) {
            result(i, i) = 1;
        }

        return result;
    }
};

template <typename T>
class h_mat_view {
private:
    std::vector<col_view<std::remove_const_t<T>>> data_;
    int rows_cnt_;

public:
    explicit h_mat_view(int rows_cnt) : rows_cnt_(rows_cnt) {};
    explicit h_mat_view() : rows_cnt_(0) {};

    h_mat_view(h_mat_view const& other) : rows_cnt_(other.rows_cnt_), data_(other.data_) {};
    h_mat_view(h_mat_view&& other) : rows_cnt_(other.rows_cnt_), data_(std::move(other.data_)) {};
    h_mat_view& operator=(const h_mat_view& other) = default;
    h_mat_view& operator=(h_mat_view&& other) = default;

    void append(col_view<T> data) {
        if(data.get_cnt() != rows_cnt_)
            throw std::out_of_range("h_mat_view size mismatch");

        data_.emplace_back(data);
    }

    void remove(int col_inx) {
        if(col_inx >= data_.size()) {
            data_.erase(data_.begin() + col_inx);
        }
    }

    [[nodiscard]] int get_row_size() const {
        return rows_cnt_;
    }

    [[nodiscard]] int get_col_size() const {
        return data_.size();
    }

    double& operator()(int row, int col) {
        if(row < 0 || row >= rows_cnt_) {
            throw std::out_of_range("row index out of range");
        }

        if(col < 0 || col >= data_.size()) {
            throw std::out_of_range("col index out of range");
        }

        return data_[col](row);
    }

    const double& operator()(int row, int col) const {
        if(row < 0 || row >= rows_cnt_) {
            throw std::out_of_range("row index out of range");
        }

        if(col < 0 || col >= data_.size()) {
            throw std::out_of_range("col index out of range");
        }

        return data_[col](row);
    }

    col_view<T> cols(int col_inx) {
        if(col_inx >= get_col_size() || col_inx < 0)
            throw std::out_of_range("col index out of range");

        return data_[col_inx];
    }

    h_mat_view<T> extract_all_cols() {
        return {this};
    }

    void remove_at(int col) {
        if(col < 0 || col >= get_col_size() || col >= data_.size())
            throw std::out_of_range("col index out of range");

        data_.erase(data_.begin() + col);
    }

    //inserts before col index
    void insert_at(int col, col_view<T> other) {
        if(col < 0  || col > data_.size())
            throw std::out_of_range("col index out of range");

        data_.insert(data_.begin() + col, other);
    }

    // col_view<const T> cols(int col_inx) const {
    //     if(col_inx >= get_col_size() || col_inx < 0)
    //         throw std::out_of_range("col index out of range");
    //
    //     return data_[col_inx];
    // }

    template <typename Matrix>
    matrix<T> operator*(Matrix& other) {
        if(get_col_size() != other.get_row_size()) {
            throw std::out_of_range("invalid dimentions");
        }

        matrix<T> result(get_row_size(), other.get_col_size());

        for(int i = 0; i < get_row_size(); i++) {
            for(int j = 0; j < other.get_col_size(); j++) {
                for(int k = 0; k < get_col_size(); k++) {
                    result(i, j) += (*this)(i, k) * other(k, j);
                }
            }
        }

        return result;
    }
};

template <typename T>
class vector {
private:
    matrix<T> mat_;

public:
    template<typename U>
    friend void print_matrix(vector<U> const& vec);

    // constructors
    explicit vector() {}
    explicit vector(int size) : mat_(size, 1) {}
    explicit vector(int size, T value) : mat_(size, 1) {
        for(int i = 0; i < size; ++i) {
            mat_(i, 0) = value;
        }
    }

    vector(const vector &other) : mat_(other.mat_) {}
    vector(vector &&other) : mat_(std::move(other.mat_)) {}
    vector &operator=(const vector &other) = default;
    vector &operator=(vector &&other) = default;

    // access
    double &operator()(int inx) {
        if (inx >= mat_.get_row_size() || inx < 0)
            throw std::out_of_range("index out of range");
        return mat_(inx, 0);
    }

    const double &operator()(int inx) const {
        if (inx >= mat_.get_row_size() || inx < 0)
            throw std::out_of_range("index out of range");
        return mat_(inx, 0);
    }

    [[nodiscard]] int get_size() const { return mat_.get_row_size(); }

    [[nodiscard]] vector extract(std::vector<int> const& indices) {
        vector result(indices.size());

        for(int i = 0; i < indices.size(); i++) {
            result(i) = (*this)(indices[i]);
        }

        return result;
    }

    [[nodiscard]] vector splice(int begin, int end) {
        vector result(end - begin + 1);

        for(int i = 0; i <= end - begin; i++) {
            result(i) = (*this)(i + begin);
        }

        return result;
    }

    //operations
    [[nodiscard]] double operator*(vector &other) {
        if (other.get_size() != get_size())
            throw std::out_of_range("vector size mismatch");

        double sum = 0;
        for (int i = 0; i < get_size(); i++) {
            sum += (*this)(i) *other(i);
        }

        return sum;
    }

    [[nodiscard]] vector operator+(vector const&other) const {
        if (other.get_size() != get_size())
            throw std::out_of_range("vector size mismatch");

        vector sum(get_size());

        for (int i = 0; i < get_size(); i++) {
            sum(i) = (*this)(i) + other(i);
        }

        return sum;
    }

    [[nodiscard]] vector operator-(vector const&other) const {
        if (other.get_size() != get_size())
            throw std::out_of_range("vector size mismatch");

        vector sum(get_size());

        for (int i = 0; i < get_size(); i++) {
            sum(i) = (*this)(i) - other(i);
        }

        return sum;
    }

    template<typename U>
    [[nodiscard]] double operator*(col_view<U> const& other) const {
        if (this->get_size() != other.get_cnt())
            throw std::out_of_range("vector size mismatch");

        double sum = 0;
        for (int i = 0; i < get_size(); i++) {
            sum += (*this)(i) *other(i);
        }

        return sum;
    }

    template<typename Matrix>
    [[nodiscard]] vector operator*(Matrix& other) {
        if (this->get_size() != other.get_row_size())
            throw std::out_of_range("matrix size mismatch");

        vector result(other.get_col_size());

        for (int i = 0; i < other.get_col_size(); i++) {
            result(i) = *this * other.cols(i);
        }

        return result;
    }
};

template <typename T>
vector<T> operator*(double scalar, vector<T> const& v) {
    vector<T> result(v.get_size());

    for (int i = 0; i < v.get_size(); i++) {
        result(i) = scalar * v(i);
    }

    return result;
}

template <typename T>
vector<T> operator*(const vector<T>& v, double scalar) {
    vector result(v.get_size());

    for (int i = 0; i < v.get_size(); i++) {
        result(i) = scalar * v(i);
    }

    return result;
}

template<typename Matrix>
void print_matrix(Matrix const& mat) {
    for(int i = 0; i < mat.get_row_size(); i++) {
        for(int j = 0; j < mat.get_col_size(); j++) {
            std::cout << mat(i, j) << "\t";
        }
        std::cout << std::endl;
    }
}

template <typename T>
void print_matrix(vector<T> const& vec) {
    print_matrix(vec.mat_);
}

template<
    template<typename> class MatrixA,
    template<typename> class MatrixB,
    typename T>
h_mat_view<T> h_concat(MatrixA<T>& mat_a, MatrixB<T>& mat_b) {
    if(mat_a.get_row_size() != mat_b.get_row_size())
        throw std::out_of_range("matrix size mismatch");

    h_mat_view<T> result = mat_a.extract_all_cols();
    h_mat_view<T> result2 = mat_b.extract_all_cols();

    for(int i = 0; i < mat_b.get_col_size(); i++) {
        result.append(result2.cols(i));
    }

    return result;
}

template<typename T>
vector<T> concat(vector<T> const& v1, vector<T> const& v2) {
    vector<T> result(v1.get_size() + v2.get_size());

    for(int i = 0; i < v1.get_size(); i++) {
        result(i) = v1(i);
    }
    for(int i = 0; i < v2.get_size(); i++) {
        result(v1.get_size() + i) = v2(i);
    }

    return result;
}
}

#endif //MATRIX_HPP
