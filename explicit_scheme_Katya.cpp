#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

double H = 0.1, TAU = 0.001, LEFT_X = -1., RIGHT_X = 1., LEFT_T = 0., RIGHT_T = 1., A = 1. / 2.;
size_t X_SIZE = 0;

// Аналитическое решение задачи для линейного случая
double linear_analytic_solution(double x, double t) {
    return (x <= A * t ? 0 : 1);
}

// Аналитическое решение задачи для нелинейного случая
double non_linear_analytic_solution(double x, double t) {
    return (x <= 0 ? 0 : (x <= t ? x / t : 1));
}

// Просчет нулевого слоя сетки с помощью начального условия
void calculate_initial_condition_linear(std::vector<double> &grid) {
    for (double x = LEFT_X; x <= RIGHT_X; x += H) {
        double val = linear_analytic_solution(x, 0.);
        // if (val < 1.e-10) val = 0.;
        grid.push_back(val);
    }
    X_SIZE = grid.size();
}

// Просчет нулевого слоя сетки с помощью начального условия
void calculate_initial_condition_non_linear(std::vector<double> &grid) {
    for (double x = LEFT_X; x <= RIGHT_X; x += H) {
        double val = non_linear_analytic_solution(x, 0.);
        // if (val < 1.e-10) val = 0.;
        grid.push_back(val);
    }
    X_SIZE = grid.size();
}

// Просчет слоев после 1го в линейном случае
void calculate_interior_layers_linear(std::vector<double> &grid) {
    std::vector<double> next_layer(X_SIZE, 0.);
    for (double t = LEFT_T + TAU; t <= RIGHT_T; t += TAU) {
        for (size_t i = 0; i < X_SIZE; ++i) {
            if (i == 0 || i == X_SIZE - 1) {
                double val = (i == 0 ? linear_analytic_solution(LEFT_X, t) : linear_analytic_solution(RIGHT_X, t));
                next_layer[i] = val;
                continue;
            }
            next_layer[i] = -A * TAU / (2 * H) * (grid[i + 1] - grid[i - 1]) +
                    grid[i] + A * A * TAU * TAU / (2 * H * H) * (grid[i + 1] - 2 * grid[i] + grid[i - 1]);
        }
        grid = next_layer;
    }
}

// Просчет слоев после 1го в нелинейном случае
void calculate_interior_layers_non_linear(std::vector<double> &grid) {
    std::vector<double> next_layer(X_SIZE, 0.);
    for (double t = LEFT_T + TAU; t <= RIGHT_T; t += TAU) {
        for (size_t i = 0; i < X_SIZE; ++i) {
            if (i == 0 || i == X_SIZE - 1) {
                double val = (i == 0 ? non_linear_analytic_solution(LEFT_X, t) : non_linear_analytic_solution(RIGHT_X, t));
                next_layer[i] = val;
                continue;
            }
            next_layer[i] =
                    grid[i] - A * TAU / (2 * H) * (grid[i + 1] * grid[i + 1] - grid[i - 1] * grid[i - 1]) +
                    A * A * TAU * TAU / (H * H) * (grid[i] * grid[i + 1] - grid[i - 1] * grid[i] - grid[i] * grid[i] + grid[i - 1] * grid[i - 1]);
            /*std::string str = std::to_string(next_layer[i]);
            if (!strchr(str.c_str(), '.') && str != "0") {
                std::cout << "not double: " << t << " " << i << std::endl;
            }*/
        }
        grid = next_layer;
    }
}

// Результат в файл
void print_result(std::vector<double> &grid, bool linear = true) {
    size_t i = 0;
    double x = LEFT_X;
    std::ofstream output;
    if (linear) {
        output.open("res_linear.txt");
    } else {
        output.open("res_non_linear.txt");
    }
    for (; i < X_SIZE; ++i, x += H) {
        output << x << " " << grid[i] << std::endl;
    }
}

// Аналитическое решение в файл
void print_analytic_solution() {
    {
        size_t i = 0;
        double x = LEFT_X;
        std::ofstream output("analytic_linear.txt");
        for (; i < X_SIZE; ++i, x += H) {
            output << x << " " << linear_analytic_solution(x, RIGHT_T) << std::endl;
        }
    }

    {
        size_t i = 0;
        double x = LEFT_X;
        std::ofstream output("analytic_non_linear.txt");
        for (; i < X_SIZE; ++i, x += H) {
            output << x << " " << non_linear_analytic_solution(x, RIGHT_T) << std::endl;
        }
    }
}

std::vector<double> calculate_linear_case() {
    std::vector<double> grid;
    calculate_initial_condition_linear(grid);
    calculate_interior_layers_linear(grid);
    print_result(grid);
    return grid;
}

std::vector<double> calculate_non_linear_case() {
    std::vector<double> grid;
    calculate_initial_condition_non_linear(grid);
    /*for (const auto row : grid) {
        for (const auto elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }*/
    calculate_interior_layers_non_linear(grid);
    print_result(grid, false);
    return grid;
}

double find_norm(const std::vector<double>& grid) {
    double res = 0;
    for (int i = 0; i < grid.size(); ++i) {
        res += grid[i] * grid[i];
    }

    return res;
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


}/*
double find_norm(const std::vector<double>& grid) {
    double res = 0;
    for (int i = 0; i < grid.size(); ++i) {
        res += grid[i] * grid[i];
    }

    return res;
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

std::vector<double> find_difference(const std::vector<double>& grid, const std::vector<double>& grid_, bool linear = true) {
    std::vector<double> difference;
    double x = LEFT_X;
    for (int i = 0; i < grid.size(); ++i, x += H) {
        difference.push_back(grid[i] - grid_[i]);
    }
    return difference;
}

void find_norms(const std::vector<double>& grid, const std::vector<double>& grid_, bool linear = true) {
    if (linear) {
        std::cout << "linear" << std::endl;
    } else {
        std::cout << "non linear" << std::endl;
    }

    { // Delta_Ch
        std::cout << Ch_norm(find_difference(grid, grid_, linear)) << " & ";
    }

    { // Delta_L1
        std::cout << L1_norm(find_difference(grid, grid_, linear)) << " & ";
    }

    { // delta_Ch
        std::cout << Ch_norm(find_difference(grid, grid_, linear)) / find_norm(grid) << " & ";
    }

    { // delta_L1
        std::cout << L1_norm(find_difference(grid, grid_, linear)) / find_norm(grid) << std::endl;
    }


}*/

int main() {
    std::vector<double> grid_linear = calculate_linear_case();
    std::vector<double> grid_non_linear = calculate_non_linear_case();
    find_norms(grid_non_linear, false);
    print_analytic_solution();
}
