#include <fstream>
// #include <cmath>
#include <vector>
#include <string>
#include <iostream>

double H = 1.e-2, TAU = 1.e-3, LEFT_X = -1., RIGHT_X = 1., LEFT_T = 0., RIGHT_T = 1., A = -1. / 2.;
size_t X_SIZE = 0;

// Аналитическое решение задачи для линейного случая
double linear_analytic_solution(double x, double t) {
    return (x + t / 2. <= 0. ? 0. : (x + t / 2. <= 1. / 4. ? 4. * x + 2. * t : 1.));
}

// Аналитическое решение задачи для нелинейного случая
double non_linear_analytic_solution(double x, double t) {
    return (x <= 0. ? 0. : (x <= 1. / 4. + t ? 4. * x / (1. + 4. * t) : 1.));
}

// Просчет нулевого слоя сетки с помощью начального условия
void calculate_initial_condition(std::vector<std::vector<double>> &grid) {
    for (double x = LEFT_X; x <= RIGHT_X; x += H) {
        double val = (x < 0. ? 0. : (x < 0.25 ? 4. * x : 1.));
        // if (val < 1.e-10) val = 0.;
        grid[0].push_back(val);
    }
    grid[1].resize(grid[0].size(), 0);
    X_SIZE = grid[0].size();
}

// Просчет первого слоя сетки в линейном случае (в моей задаче это нужно было делать другой схемой, так как начальная схема использует 3 слоя)
void calculate_first_layer_linear(std::vector<std::vector<double>> &grid) {
    for (size_t i = 0; i < X_SIZE; ++i) {
        if (i == 0 || i == X_SIZE - 1) {
            double val = (i == 0 ? linear_analytic_solution(LEFT_X, TAU) : linear_analytic_solution(RIGHT_X, TAU));
            grid[1][i] = val;
            continue;
        }
        grid[1][i] = grid[0][i] + (grid[0][i + 1] - grid[0][i - 1]) * TAU * -A / (2 * H);
    }
}

// Просчет первого слоя сетки в нелинейном случае (в моей задаче это нужно было делать другой схемой, так как начальная схема использует 3 слоя)
void calculate_first_layer_non_linear(std::vector<std::vector<double>> &grid) {
    for (size_t i = 0; i < X_SIZE; ++i) {
        if (i == 0 || i == X_SIZE - 1) {
            double val = (i == 0 ? non_linear_analytic_solution(LEFT_X, TAU) : non_linear_analytic_solution(RIGHT_X, TAU));
            grid[1][i] = val;
            continue;
        }
        // grid[1][i] = grid[0][i] + TAU / H * (grid[0][i] * grid[0][i] - grid[0][i - 1] * grid[0][i - 1]); // сапог Самохина
        grid[1][i] = grid[0][i] + A * TAU / H * (grid[0][i + 1] * grid[0][i + 1] - grid[0][i - 1] * grid[0][i - 1]); // Шляпа без центра
        // std::cout << grid[1][i] << " ";
    }
}

// Просчет слоев после 1го в линейном случае
void calculate_interior_layers_linear(std::vector<std::vector<double>> &grid) {
    std::vector<double> next_layer(X_SIZE, 0.);
    for (double t = LEFT_T + 2 * TAU; t <= RIGHT_T; t += TAU) {
        for (size_t i = 0; i < X_SIZE; ++i) {
            if (i == 0 || i == X_SIZE - 1) {
                double val = (i == 0 ? linear_analytic_solution(LEFT_X, t) : linear_analytic_solution(RIGHT_X, t));
                next_layer[i] = val;
                continue;
            }
            next_layer[i] = (grid[1][i + 1] - grid[1][i - 1]) * TAU / (2 * H) + grid[0][i];
        }
        grid[0] = grid[1];
        grid[1] = next_layer;
    }
}

// Просчет слоев после 1го в нелинейном случае
void calculate_interior_layers_non_linear(std::vector<std::vector<double>> &grid) {
    std::vector<double> next_layer(X_SIZE, 0.);
    for (double t = LEFT_T + 2 * TAU; t <= RIGHT_T; t += TAU) {
        for (size_t i = 0; i < X_SIZE; ++i) {
            if (i == 0 || i == X_SIZE - 1) {
                double val = (i == 0 ? non_linear_analytic_solution(LEFT_X, t) : non_linear_analytic_solution(RIGHT_X, t));
                next_layer[i] = val;
                continue;
            }
            next_layer[i] =
                    (grid[1][i + 1] * grid[1][i + 1] - grid[1][i - 1] * grid[1][i - 1]) * TAU / (2 * H) + grid[0][i];
            std::string str = std::to_string(next_layer[i]);
            if (!strchr(str.c_str(), '.') && str != "0") {
                std::cout << "not double: " << t << " " << i << std::endl;
            }
        }
        grid[0] = grid[1];
        grid[1] = next_layer;
    }
}

// Результат в файл
void print_result(std::vector<std::vector<double>> &grid, bool linear = true) {
    size_t i = 0;
    double x = LEFT_X;
    std::ofstream output;
    if (linear) {
        output.open("res_linear.txt");
    } else {
        output.open("res_non_linear.txt");
    }
    for (; i < X_SIZE; ++i, x += H) {
        output << x << " " << grid[1][i] << std::endl;
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

void calculate_linear_case() {
    std::vector<std::vector<double>> grid(2, std::vector<double>());
    calculate_initial_condition(grid);
    calculate_first_layer_linear(grid);
    calculate_interior_layers_linear(grid);
    print_result(grid);
}

void calculate_non_linear_case() {
    std::vector<std::vector<double>> grid(2, std::vector<double>());
    calculate_initial_condition(grid);
    calculate_first_layer_non_linear(grid);
    /*for (const auto row : grid) {
        for (const auto elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }*/
    calculate_interior_layers_non_linear(grid);
    print_result(grid, false);
}

void find_linear_system_solution() {
    double A_i, B_i, C_i;

}

int main() {
    calculate_linear_case();
    calculate_non_linear_case();
    print_analytic_solution();
}