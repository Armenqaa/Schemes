#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

double H = 0.1, TAU = 0.0005, LEFT_X = -1., RIGHT_X = 1., LEFT_T = 0., RIGHT_T = 1., OMEGA = 0.1, A = 1. / 2.;
int X_SIZE = 0;

double find_norm(const std::vector<double> &vec) {
    /*double res = vec[0];
    for (double elem : vec) {
        res = std::max(elem, res);
    }
    return res;*/
    double res = 0.;
    for (int i = 0; i < vec.size(); ++i) {
        res += vec[i] * vec[i];
    }
    return res;
}

std::vector<double> plus(const std::vector<double> &v1, const std::vector<double> &v2) {
    std::vector<double> res;
    // Размер дельта х меньше размера v
    for (int i = 1; i < v1.size() - 1; ++i) {
        res.push_back(v1[i] + v2[i - 1]);
    }
    return res;
}

// Аналитическое решение задачи для линейного случая
double linear_analytic_solution(double x, double t) {
    return (x <= A * t ? 0 : 1);
}

// Аналитическое решение задачи для нелинейного случая
double non_linear_analytic_solution(double x, double t) {
    return (x <= 0 ? 0 : (x <= t ? x / t : 1));
}

void calculate_G_linear(std::vector<double> &G, const std::vector<double> &v_next, const std::vector<double> &v) {
    for (int i = 1; i < X_SIZE - 1; ++i) {
        double right_side_of_G = -TAU * A / (2 * H) * v_next[i - 1]  + v_next[i] - v[i] +
                                 TAU * A / (2 * H) * v_next[i + 1] - OMEGA * H * H * (v[i + 1] - 2 * v[i] + v[i - 1]);
        G[i - 1] = right_side_of_G;
    }
}

void calculate_G_non_linear(std::vector<double> &G, const std::vector<double> &v_next, const std::vector<double> &v) {
    for (int i = 1; i < X_SIZE - 1; ++i) {
        double right_side_of_G = -TAU * A / (2 * H) * v_next[i - 1] * v_next[i - 1] + v_next[i] - v[i] +
                TAU * A / (2 * H) * v_next[i + 1] * v_next[i + 1] - OMEGA * H * H * (v[i + 1] - 2 * v[i] + v[i - 1]);
        G[i - 1] = right_side_of_G;
    }
}

void calculate_J_non_linear(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, const std::vector<double> &v_next) {
    for (int i = 1; i < X_SIZE - 1; ++i) { // b size == X_SIZE - 2
        b.push_back(1.);
    }

    for (int i = 1; i < X_SIZE - 2; ++i) { // a size == c size == X_SIZE - 3
        a.push_back(-TAU * A / H * v_next[i - 1]);
        c.push_back(TAU * A/ H * v_next[i + 1]);
    }
}

// a - поддиагональ, b - диагональ, c - наддиагональ
void calculate_J_linear(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, const std::vector<double> &v_next) {
    for (int i = 1; i < X_SIZE - 1; ++i) { // b size == X_SIZE - 2
        b.push_back(1.);
    }

    for (int i = 1; i < X_SIZE - 2; ++i) { // a size == c size == X_SIZE - 3
        a.push_back(-TAU * A / (2 * H));
        c.push_back(TAU * A/ (2 * H));
    }
}

// a - поддиагональ, b - диагональ, c - наддиагональ
std::vector<double> solve_equations(std::vector<double> &v_next,
                                    const std::vector<double> &v, std::vector<double> &delta_x, bool linear) {
    std::vector<double> a, b, c, G(X_SIZE - 2);
    if (!linear) {
        calculate_J_non_linear(a, b, c, v_next);
        calculate_G_non_linear(G, v_next, v);
    } else {
        calculate_J_linear(a, b, c, v_next);
        calculate_G_linear(G, v_next, v);
    }
    /*
    std::vector<double> a, b, c, G_minus;
    c = {1, 1, 1};
    b = {1, 1, 1, 1};
    a = {0, 0, 0};
    G_minus = {1, 1, 1, 1};
    X_SIZE = b.size() + 2;
    */
    G[0] = -G[0];
    for (int i = 1; i < X_SIZE - 2; ++i) {
        b[i] -= a[i - 1] * c[i - 1] / b[i - 1];
        G[i] = -G[i] - a[i - 1] * G[i - 1] / b[i - 1];
    }

    delta_x.resize(G.size());
    delta_x.back() = G.back() / b.back();
    // Решение СЛАУ методом прогонки
    for (int i = X_SIZE - 4; i >= 0; --i) {
        delta_x[i] = (G[i] - c[i] * delta_x[i + 1]) / b[i];
    }

    std::vector<double> v_next_for_isaev_sonin = v_next;
    std::vector<double> G_isaev_sonin = G;
    for (int k = 0; k < 20; ++k) {
        for (int i = 1; i < X_SIZE - 1; ++i) {
            v_next_for_isaev_sonin[i] = delta_x[i - 1] * pow(2, -k) + v[i];
        }
        if (linear) {
            calculate_G_linear(G_isaev_sonin, v_next_for_isaev_sonin, v);
        } else {
            calculate_G_non_linear(G_isaev_sonin, v_next_for_isaev_sonin, v);
        }
        if (find_norm(G_isaev_sonin) < find_norm(G)) {
            break;
        }
    }
    v_next = v_next_for_isaev_sonin;
    /*
    std::cout << "delta_x:";
    for (int i = 0; i < delta_x.size(); ++i) {
        std::cout << delta_x[i] << " ";
    }
    std::cout << std::endl;
    */
    return G_isaev_sonin;
}

void prepare_to_next_iteration_of_newton(double t, std::vector<double> &v_next, std::vector<double> &v,
                                         std::vector<double> &delta_x, bool linear) {
    std::vector<double> tmp = v_next;
    v_next.clear();
    if (!linear) {
        v_next.push_back(non_linear_analytic_solution(LEFT_X, t + TAU));
    } else {
        v_next.push_back(linear_analytic_solution(LEFT_X, t + TAU));
    }
    std::vector<double> plus_res = plus(v, delta_x);
    for (int i = 0; i < plus_res.size(); ++i) {
        v_next.push_back(plus_res[i]);
    }

    if (!linear) {
        v_next.push_back(non_linear_analytic_solution(RIGHT_X, t + TAU));
    } else {
        v_next.push_back(linear_analytic_solution(RIGHT_X, t + TAU));
    }
    v = tmp;
}

void print_norm_and_current_layers(double t, std::vector<double> &G, std::vector<double> &v_next, std::vector<double> &v) {
    std::cout << t << ":" << find_norm(G) << std::endl;
    std::cout << "v:";
    for (int i = 0; i < v.size(); ++i) {
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "v_next:";
    for (int i = 0; i < v_next.size(); ++i) {
        std::cout << v_next[i] << " ";
    }
    std::cout << std::endl;
}

void start_newton(std::vector<double> &v, std::vector<double> &v_next, double t) {
    std::vector<double> delta_x, G;
    bool linear = false;
    do {
        G = solve_equations(v_next, v, delta_x, linear);
        prepare_to_next_iteration_of_newton(t, v_next, v, delta_x, linear);
        // print_norm_and_current_layers(t, G, v_next, v);
    } while (find_norm(G) > 1.e-3);
}

void initialize_for_zero_iteration(std::vector<double>& v_next, std::vector<double>& v) {
    for (double x = LEFT_X; x <= RIGHT_X; x += H) {
        // std::cout << x << " " << non_linear_analytic_solution(x, 0) << std::endl;
        v.push_back(non_linear_analytic_solution(x, 0));
        v_next.push_back(non_linear_analytic_solution(x, 0));
        X_SIZE++;
    }
}

void print_result(std::vector<double> &v) {
    std::ofstream out("res.txt");
    double x = LEFT_X;
    for (int i = 0; i < v.size(); ++i, x += H) {
        out << x << " " << v[i] << " " << std::endl;
    }
}

void print_analytic_solutions() {
    std::ofstream analytic_linear("analytic_linear.txt"), analytic_non_linear("analytic_non_linear.txt");
    for (double x = LEFT_X; x < RIGHT_X; x += H) {
        analytic_linear << x << " " << linear_analytic_solution(x, RIGHT_T) << std::endl;
    }

    for (double x = LEFT_X; x < RIGHT_T; x += H) {
        analytic_non_linear << x << " " << non_linear_analytic_solution(x, RIGHT_T) << std::endl;
    }
}

double Ch_norm(const std::vector<double>& grid) {
    if (grid.size() == 0) {
        return -1;
    }

    double res = fabs(grid[0]);
    for (int i = 1; i < grid.size(); ++i) {
        res = std::max(res, fabs(grid[i]));
    }

    return res;
}

double L1_norm (const std::vector<double>& grid) {
    double res = 0.;
    for (int i = 0; i < grid.size(); ++i) {
        res += abs(grid[i]);
    }

    return res * H;
}

std::vector<double> find_difference(const std::vector<double>& grid, bool linear = true) {
    std::vector<double> difference;
    double x = LEFT_X;
    for (int i = 0; i < grid.size(); ++i, x += H) {
        if (linear) {
            difference.push_back(grid[i] - linear_analytic_solution(x, 1));
        } else {
            difference.push_back(grid[i] - non_linear_analytic_solution(x, 1));
        }
    }
    return difference;
}

void find_norms(const std::vector<double>& grid, bool linear = true) {
    if (linear) {
        std::cout << "linear" << std::endl;
    } else {
        std::cout << "non linear" << std::endl;
    }

    { // Delta_Ch
        std::cout << Ch_norm(find_difference(grid, linear)) << " & ";
    }

    { // Delta_L1
        std::cout << L1_norm(find_difference(grid, linear)) << " & ";
    }

    { // delta_Ch
        std::cout << Ch_norm(find_difference(grid, linear)) / find_norm(grid) << " & ";
    }

    { // delta_L1
        std::cout << L1_norm(find_difference(grid, linear)) / find_norm(grid) << std::endl;
    }


}

int main() {
    std::vector<double> v, v_next;
    initialize_for_zero_iteration(v_next, v);
    for (double t = LEFT_T; t <= RIGHT_T; t += TAU) {
        start_newton(v, v_next, t);
    }

    print_result(v);
    print_analytic_solutions();
    /*for (int i = 0; i < v.size(); ++i) {
        std::cout << v_next[i] << " ";
    }*/
    find_norms(v_next, false);
    return 0;
    //std::cout << X_SIZE << std::endl;
    //std::vector<double> v1, v2, delta_x;
    //solve_equations(v1, v2, delta_x, true);
}
