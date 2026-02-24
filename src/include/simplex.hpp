#ifndef SIMPLEX_HPP
#define SIMPLEX_HPP
#include <cmath>
#include <limits>
#include <algorithm>


#include "matrix.hpp"

namespace smp {

constexpr double EPS = 1e-10;

struct task {
    mat::matrix<double> A;
    mat::vector<double> b;
    mat::vector<double> c;
};

enum status_codes {
    CONTINUE,
    FINISH,
    UNBOUNDED,
    EMPTY
};

struct status {
    mat::vector<double> x;
    status_codes code;
};

//TODO: Optimize
// Function to compute the rank of a matrix using Gaussian elimination
template <typename ClonableMatrix>
int rank(ClonableMatrix const& A) {
    mat::matrix<double> A_buff(A);

    int m = A_buff.get_row_size();
    int n = A_buff.get_col_size();
    int rank = 0;

    for (int col = 0, row = 0; col < n && row < m; ++col) {
        int pivot = row;
        for (int i = row; i < m; ++i)
            if (fabs(A_buff(i, col)) > fabs(A_buff(pivot, col)))
                pivot = i;

        if (fabs(A_buff(pivot, col)) < EPS) continue;

        A_buff.swap_rows(row, pivot);

        for (int i = row + 1; i < m; ++i) {
            double factor = A_buff(i, col) / A_buff(row, col);
            for (int j = col; j < n; ++j)
                A_buff(i, j) -= factor * A_buff(row, j);
        }

        ++row;
        ++rank;
    }

    return rank;
}

template<typename ClonableMatrix>
mat::matrix<double> invert_matrix(ClonableMatrix &A) {
    int n = A.get_row_size();

    mat::matrix<double> A_buff(A);
    auto inverse = mat::matrix<double>::E(A.get_row_size());

    for (int i = 0; i < n; i++) {
        int pivot_row = i;
        for (int r = i + 1; r < n; ++r) {
            if (fabs(A_buff(r, i)) > fabs(A_buff(pivot_row, i)))
                pivot_row = r;
        }

        if (fabs(A_buff(pivot_row, i)) < EPS)
            throw std::runtime_error("Matrix is not invertible (zero pivot).");

        if (pivot_row != i) {
            A_buff.swap_rows(pivot_row, i);
            inverse.swap_rows(pivot_row, i);
        }

        auto pivot = A_buff(i, i);

        for (int j = 0; j < n; j++) {
            A_buff(i, j) /= pivot;
            inverse(i, j) /= pivot;
        }

        for (int k = 0; k < n; k++) {
            if (k == i)
                continue;

            double factor = A_buff(k, i);
            for (int j = 0; j < n; j++) {
                A_buff(k, j) -= factor * A_buff(i, j);
                inverse(k, j) -= factor * inverse(i, j);
            }
        }
    }

    return inverse;
}

inline mat::vector<double> calculate_move_direction(int size,
    mat::matrix<double>& B, mat::matrix<double>& A,
    std::vector<int> const& n_ext, std::vector<int> const& l, int neg_index) {
    mat::vector<double> result(size);

    auto BA = B * A.cols(neg_index);

    for(int i = 0; i < n_ext.size(); ++i) {
        result(n_ext[i]) = BA(i);
    }

    for(int i = 0; i < l.size(); ++i) {
        result(l[i]) = 0;
    }

    result(neg_index) = -1;

    return result;
}

//returns n_ext, l
std::pair<std::vector<int>, std::vector<int>> fill_till_MxM(mat::h_mat_view<double>& A_np, std::vector<int>& np,
    mat::h_mat_view<double>& A_nn, std::vector<int>& nn) {
    int prev_rk = A_np.get_col_size();

    std::vector<int> n_ext(np);
    std::vector<int> l;

    for(int i = 0; i < A_nn.get_col_size(); ++i) {
        A_np.append(A_nn.cols(i));
        auto new_rk = rank(A_np);

        if(new_rk == prev_rk) {
            A_np.remove_at(A_np.get_col_size() - 1);
            l.emplace_back(nn[i]);
        } else {
            n_ext.emplace_back(nn[i]);
            prev_rk = new_rk;
        }
    }

    return {n_ext, l};
}

inline int find_first_neg_index(mat::vector<double>& d, std::vector<int>& l) {
    for(int i = 0; i < d.get_size(); ++i) {
        if(d(i) < 0) {
            return l[i];
        }
    }

    return -1;
}

inline double calculate_move_coefficient(
    mat::vector<double> const& x_prev, mat::vector<double>& direction, std::vector<int> const& n_ext) {

    double theta = std::numeric_limits<double>::infinity();

    for(int i = 0; i < n_ext.size(); ++i) {
        auto new_theta = x_prev(n_ext[i]) / direction(n_ext[i]);
        if(direction(n_ext[i]) > EPS && new_theta < theta) {
            theta = new_theta;
        }
    }

    return theta;
}

inline bool all_non_positive(mat::vector<double>& direction, std::vector<int>& n_ext, std::vector<int>& np) {
    for(int i = 0; i < n_ext.size(); ++i) {
        if(std::find(np.begin(), np.end(), n_ext[i]) != np.end()) {
            continue;
        }

        if(direction(n_ext[i]) >= 0)
            return false;
    }

    return true;
}

inline void find_new_basis(mat::matrix<double>& A, mat::h_mat_view<double>& A_np,
    std::vector<int>& np, std::vector<int>& l_used, std::vector<int>& l, std::vector<int>& n_ext) {
    auto current_rk = rank(A_np);
    bool have_changed_basis = false;

    for(int i = 0; i < n_ext.size(); ++i) {
        if(std::find(np.begin(), np.end(), n_ext[i]) != np.end()) {
            continue;
        }

        for(int j = 0; j < l.size(); ++j) {
            if(std::find(n_ext.begin(), n_ext.end(), l[j]) != n_ext.end()) {
                continue;
            }
            if(std::find(l_used.begin(), l_used.end(), l[j]) != l_used.end()) {
                continue;
            }
            //j-ый вектор из l ещё не положили в n_ext, попробуём поменять с i-ым добавленным вектором из n_ext

            A_np.remove_at(i);
            A_np.insert_at(i, A.cols(l[j]));

            auto new_rk = rank(A_np);

            if(new_rk < current_rk) {
                A_np.remove_at(i);
                A_np.insert_at(i, A.cols(n_ext[i]));
            } else {
                auto tmp1 = n_ext[i];
                auto tmp2 = l[j];

                n_ext.erase(n_ext.begin() + i);
                n_ext.insert(n_ext.begin() + i, l[j]);
                l.erase(l.begin() + j);
                l.insert(l.begin() + j, tmp1);

                l_used.push_back(tmp1);

                have_changed_basis = true;
            }

            if(have_changed_basis) {
                break;
            }
        }
        if(have_changed_basis) {
            break;
        }
    }
}

inline bool is_unbounded(mat::vector<double>& direction, std::vector<int>& n_ext) {
    for(int i = 0; i < n_ext.size(); ++i) {
        if(direction(n_ext[i]) > EPS)
            return false;
    }
    return true;
}

inline mat::vector<double> remove_neg_error(mat::vector<double> const& x) {
    mat::vector<double> res(x);

    for(int i = 0; i < x.get_size(); ++i) {
        if(x(i) < 0 && x(i) > -EPS) {
            res(i) = 0;
        }
    }

    return res;
}

inline status simplex(mat::vector<double>& x_prev, mat::matrix<double> &A,
                                    mat::vector<double> &c) {
    const auto N = x_prev.get_size();
    const auto M = A.get_row_size();

    std::vector<int> np;
    std::vector<int> nn;
    std::vector<int> n_ext;
    std::vector<int> l;
    std::vector<int> l_used;


    for(int i = 0; i < x_prev.get_size(); i++) {
        if(x_prev(i) > EPS) {
            np.push_back(i);

            if(rank(A.extract_cols(np)) != np.size()) {
                nn.push_back(i);
                np.pop_back();
            }
        } else {
            nn.push_back(i);
        }
    }

    auto A_np = A.extract_cols(np);
    auto A_nn = A.extract_cols(nn);

    if(A_np.get_col_size() != M) { //дополняем до квадратной
        std::tie(n_ext, l) = fill_till_MxM(A_np, np, A_nn, nn);
    } else {
        n_ext = np;
        l = nn;
    }

    A_np = A.extract_cols(n_ext);

    while(true) {
        auto B = invert_matrix(A_np);
        auto tmp = A.extract_cols(l);
        auto d_l = remove_neg_error(c.extract(l) - c.extract(n_ext) * B * tmp);
        auto j = find_first_neg_index(d_l, l);

        if(j == -1) { //нет отрицательных коэфф — оптимальное решение
            return {x_prev, FINISH};
        }

        auto direction = calculate_move_direction(N, B, A, n_ext, l, j);

        if(is_unbounded(direction, n_ext)) { //все вектора отрицательны — тета растёт бесконечно
            return {x_prev, UNBOUNDED};
        }

        if(np.size() < M && !all_non_positive(direction, n_ext, np)) { //надо менять базис
            find_new_basis(A, A_np, np, l_used, l, n_ext);
        } else { //тэта существует
            auto theta = calculate_move_coefficient(x_prev, direction, n_ext);

            return {remove_neg_error(x_prev - theta * direction), CONTINUE};
        }
    }
}

inline bool any_greater_than_zero(mat::vector<double> const& y) {
    for(int i = 0; i < y.get_size(); ++i) {
        if(y(i) > EPS) {
            return true;
        }
    }
    return false;
}

inline status find_vector(task& task) {
    const auto M = task.A.get_row_size();
    const auto N = task.A.get_col_size();

    mat::matrix<double> A_ext = mat::matrix<double>::E(M);
    mat::matrix new_A(mat::h_concat(task.A, A_ext));

    mat::vector<double> x_base(N);
    mat::vector<double> new_x = mat::concat(x_base, task.b);

    mat::vector<double> c_base(N);
    mat::vector<double> c_ext(M, 1);
    mat::vector<double> new_c = mat::concat(c_base, c_ext);

    status_codes code = CONTINUE;

    do {
        auto status = simplex(new_x, new_A, new_c);
        new_x = status.x;
        code = status.code;
    } while(code == CONTINUE);

    if(any_greater_than_zero(new_x.splice(N - 1, N + M - 1))) {
        return {new_x, EMPTY};
    }

    return {new_x.splice(0, N - 1), code};
}

inline status simplex_main(task& task) {
    const auto N = task.A.get_col_size();

    auto pivot_status = find_vector(task);

    if(pivot_status.code != FINISH) {
        return pivot_status;
    }

    auto x = pivot_status.x;
    status_codes code = CONTINUE;

    do {
        auto status = simplex(x, task.A, task.c);
        x = status.x;
        code = status.code;
    } while(code == CONTINUE);

    return {x, code};
}

}

#endif //SIMPLEX_HPP
