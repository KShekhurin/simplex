#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <functional>
#include <climits>

using namespace std;

namespace szlab {

    // 1) Структура ячейки транспортной таблицы
    // Хранит тариф, текущий план перевозки и признак базисности
    struct Cell {
        double cost;      // Стоимость перевозки (C_ij)
        double product;      // Количество перевезенного груза (X_ij)
        bool isBase;    // Является ли ячейка базисной (true) или свободной (false)

        // Конструктор по умолчанию
        Cell() : cost(0.0), product(0.0), isBase(false) {}

        // Конструктор с инициализацией стоимости
        Cell(double c) : cost(c), product(0.0), isBase(false) {}
    };

    struct Dot {
        int x;
        int y;

        Dot(int setx, int sety) : x(setx), y(sety) {}
    };

    class TransportProblem {
    public:
        int m; // Количество поставщиков (строки)
        int n; // Количество потребителей (столбцы)

        std::vector<double> supplies; // Вектор запасов (A)
        std::vector<double> demands;  // Вектор потребностей (B)

        std::vector<std::vector<Cell>> table;

        // 2) Векторы для метода потенциалов
        std::vector<double> u; // Потенциалы поставщиков (строки)
        std::vector<double> v; // Потенциалы потребителей (столбцы)
        vector<int> u_indexes; // Вектор индексов базисных элементов u
        vector<int> v_indexes; // Вектор индексов базисных элементов v

        // Вектор хранения базисных ячеек
        std::vector<Cell> basis_vec;

        // Конструктор класса
        TransportProblem(int rows, int cols) : m(rows), n(cols) {
            // Инициализация векторов запасов и потребностей
            supplies.resize(m, 0.0);
            demands.resize(n, 0.0);

            // Инициализация таблицы ячеек
            table.resize(m, std::vector<Cell>(n));

            // Инициализация векторов потенциалов
            u.resize(m, 0.0);
            v.resize(n, 0.0);
        }

        void balanceTable() {
            double sumDemands = accumulate(demands.begin(), demands.end(), 0.0);
            double sumSupplies = accumulate(supplies.begin(), supplies.end(), 0.0);
            double diff = sumSupplies - sumDemands;
            cout << "Начало балансировки..." << endl;
            if (diff > 0.0) {
                demands.push_back(diff);
                v.push_back(0.0);
                for (int i = 0; i < m; ++i) {\
                    table[i].push_back(Cell());
                }
                n++;
                cout << "Добавлен фиктивный потребитель с потребностью " << diff << endl;
            }
            if (diff < 0.0) {
                diff *= -1;
                supplies.push_back(diff);
                u.push_back(0.0);
                std::vector<Cell> newRow(n, Cell());
                table.push_back(newRow);
                m++;
                cout << "Добавлен фиктивный поставщик с количеством товара " << diff << endl;
            }
            cout << "Конец балансировки" << endl;
        }

        void NorthwestAngle() {
            vector<double> supplies_dynamic = supplies;
            vector<double> demands_dynamic = demands;
            int i = 0;
            int j = 0;
            cout << "Выполняю метод северо-западного угла..." << endl;
            while (i != m-1 || j != n-1) {
                double diff = supplies_dynamic[i] - demands_dynamic[j];
                table[i][j].product = (supplies_dynamic[i] <= demands_dynamic[j] ? supplies_dynamic[i] : demands_dynamic[j]);
                if (diff > 0.0 && demands_dynamic[j] != 0) {
                    supplies_dynamic[i] -= demands_dynamic[j];
                    table[i][j].isBase = true;
                    // для метода потенциалов нужны индексы базисных элементов
                    u_indexes.push_back(i);
                    v_indexes.push_back(j);
                    j++;
                };
                if (diff < 0.0 && supplies_dynamic[i] != 0) {
                    demands_dynamic[j] -= supplies_dynamic[i];
                    table[i][j].isBase = true;
                    u_indexes.push_back(i);
                    v_indexes.push_back(j);
                    i++;
                };
                if (diff == 0.0) {
                    table[i][j].isBase = true;
                    table[i][j+1].isBase = true;
                    table[i][j+1].product = 0.0;
                    u_indexes.push_back(i);
                    v_indexes.push_back(j);
                    u_indexes.push_back(i);
                    v_indexes.push_back(j+1);
                    j++;
                    i++;
                }
            }
            table[i][j].product = (supplies_dynamic[i] <= demands_dynamic[j] ? supplies_dynamic[i] : demands_dynamic[j]);
            table[i][j].isBase = true;
            u_indexes.push_back(i);
            v_indexes.push_back(j);
            cout << "Метод северо-западного угла выполнен" << endl;
        }

        void CalculatePotentials() {
            // нам доступны векторы индексов потенциалов
            // из-за того что они по порядку то можно по ним проходить
            // и вычислять все потенциалы и так далее
            u[0] = 0;
            v[0] = table[0][0].cost;
            cout << "Считаю потенциалы..." << endl;
            for (int current_potential = 0; current_potential < (n + m - 1); current_potential++) {
                for (int previous_potential = 0; previous_potential < current_potential; ++previous_potential) {
                    int u_cur = u_indexes[current_potential];
                    int v_cur = v_indexes[current_potential];
                    int u_prev = u_indexes[previous_potential];
                    int v_prev = v_indexes[previous_potential];
                    bool break_flag = false;
                    // сравним индексы каждой компоненты u и v
                    // они совпали - не бывает
                    // оба индекса разные - пропускаем
                    // иначе - вычисляем новую компоненту и выходим из этого цикла
                    if (u_cur != u_prev && v_cur == v_prev) {
                        u[u_cur] = table[u_cur][v_cur].cost - v[v_prev];
                        break_flag = true;
                    }
                    else if (u_cur == u_prev && v_cur != v_prev) {
                        //u[current_potential] = v[previous_potential] - table[u[current_potential]][v[current_potential]].cost;
                        v[v_cur] = table[u_cur][v_cur].cost - u[u_prev];
                        break_flag = true;
                    }
                    if (break_flag) break;
                }
            }
        }

        Dot findUnopt() {
                double min_delta = 0.0;  // Ищем только отрицательные значения
                int min_i = -1, min_j = -1;

                // Перебираем все свободные (небазисные) ячейки
                for (int i = 0; i < m; ++i) {
                    for (int j = 0; j < n; ++j) {
                        if (!table[i][j].isBase) {
                            // Вычисляем оценку ячейки: Δ = c_ij - u_i - v_j
                            double delta = table[i][j].cost - u[i] - v[j];

                            // Если оценка отрицательная и меньше текущей минимальной
                            if (delta < min_delta) {
                                min_delta = delta;
                                min_i = i;
                                min_j = j;
                            }
                        }
                    }
                }

                // Если отрицательных оценок нет — решение оптимально
                if (min_i == -1 && min_j == -1) {
                    return Dot(-1, -1);  // Специальное значение
                }

                // Возвращаем координаты ячейки с наибольшей возможностью улучшения
                return Dot(min_i, min_j);
            }

        void optiCycle(Dot matrix_pos) {
            // Поиск цикла пересчёта методом DFS
            // Цикл должен начинаться и заканчиваться в matrix_pos,
            // чередовать горизонтальные и вертикальные ходы,
            // проходить только через базисные ячейки (кроме стартовой)

            vector<pair<int, int>> cycle;
            cycle.push_back({matrix_pos.x, matrix_pos.y});

            // Рекурсивная функция DFS для поиска цикла
            function<bool(int, int, bool)> dfs = [&](int i, int j, bool move_in_row) -> bool {
                // Условие завершения: вернулись в стартовую ячейку с циклом из ≥4 ячеек
                if (cycle.size() > 1 && i == matrix_pos.x && j == matrix_pos.y) {
                    return true;
                }

                if (move_in_row) {
                    // Горизонтальный ход: та же строка, другой столбец
                    for (int nj = 0; nj < n; ++nj) {
                        if (nj == j) continue;

                        // Можно использовать только стартовую ячейку или базисные
                        if (!((i == matrix_pos.x && nj == matrix_pos.y) || table[i][nj].isBase)) {
                            continue;
                        }

                        // Проверка, что ячейка ещё не в цикле (кроме стартовой)
                        bool in_cycle = false;
                        for (size_t k = 1; k < cycle.size(); ++k) {
                            if (cycle[k].first == i && cycle[k].second == nj) {
                                in_cycle = true;
                                break;
                            }
                        }
                        if (in_cycle) continue;

                        cycle.push_back({i, nj});
                        if (dfs(i, nj, false)) return true;
                        cycle.pop_back();
                    }
                } else {
                    // Вертикальный ход: тот же столбец, другая строка
                    for (int ni = 0; ni < m; ++ni) {
                        if (ni == i) continue;

                        if (!((ni == matrix_pos.x && j == matrix_pos.y) || table[ni][j].isBase)) {
                            continue;
                        }

                        bool in_cycle = false;
                        for (size_t k = 1; k < cycle.size(); ++k) {
                            if (cycle[k].first == ni && cycle[k].second == j) {
                                in_cycle = true;
                                break;
                            }
                        }
                        if (in_cycle) continue;

                        cycle.push_back({ni, j});
                        if (dfs(ni, j, true)) return true;
                        cycle.pop_back();
                    }
                }
                return false;
            };

            // Попытка найти цикл: сначала с горизонтального хода, потом с вертикального
            if (!dfs(matrix_pos.x, matrix_pos.y, true)) {
                cycle.clear();
                cycle.push_back({matrix_pos.x, matrix_pos.y});
                dfs(matrix_pos.x, matrix_pos.y, false);
            }

            // Если цикл не найден или неполный — выход
            if (cycle.size() < 4 || cycle.back() != make_pair(matrix_pos.x, matrix_pos.y)) {
                return;
            }

            // Убираем дубликат конечной ячейки
            cycle.pop_back();

            // Находим θ = min{продукции} на позициях со знаком "-" (нечётные индексы: 1, 3, 5...)
            double theta = numeric_limits<double>::max();
            for (size_t k = 1; k < cycle.size(); k += 2) {
                int i = cycle[k].first;
                int j = cycle[k].second;
                theta = min(theta, table[i][j].product);
            }

            // Корректируем значения product вдоль цикла:
            // чётные индексы (0,2,4...): +θ, нечётные (1,3,5...): -θ
            for (size_t k = 0; k < cycle.size(); ++k) {
                int i = cycle[k].first;
                int j = cycle[k].second;
                if (k % 2 == 0) {
                    table[i][j].product += theta;
                } else {
                    table[i][j].product -= theta;
                }
            }

            // Обновляем базис: входящая ячейка становится базисной
            table[matrix_pos.x][matrix_pos.y].isBase = true;

            // Находим и помечаем выходящую ячейку (базисная, ставшая нулевой)
            for (size_t k = 1; k < cycle.size(); k += 2) {
                int i = cycle[k].first;
                int j = cycle[k].second;
                if (table[i][j].isBase && abs(table[i][j].product) < 1e-9) {
                    table[i][j].isBase = false;
                    table[i][j].product = 0.0;
                    break;
                }
            }

            // Обновляем векторы индексов для пересчёта потенциалов
            u_indexes.push_back(matrix_pos.x);
            v_indexes.push_back(matrix_pos.y);
        }

        void PotentialsMethod() {
            cout << "\n=== Запуск метода потенциалов ===" << endl;

            int iteration = 0;
            const int MAX_ITERATIONS = 100;  // Защита от зацикливания

            while (iteration < MAX_ITERATIONS) {
                cout << "\n--- Итерация " << (iteration + 1) << " ---" << endl;

                // Шаг 1: Вычисляем потенциалы для текущего базиса
                CalculatePotentials();

                // Шаг 2: Ищем неоптимальную ячейку
                Dot unopt = findUnopt();

                // Проверка на оптимальность: если отрицательных оценок нет
                if (unopt.x == -1 && unopt.y == -1) {
                    cout << "✓ Достигнуто оптимальное решение!" << endl;
                    break;
                }

                double delta = table[unopt.x][unopt.y].cost - u[unopt.x] - v[unopt.y];
                cout << "→ Найдена неоптимальная ячейка: (" << unopt.x << ", " << unopt.y
                     << "), Δ = " << delta << endl;

                // Шаг 3: Выполняем пересчёт по циклу
                optiCycle(unopt);

                // Вывод текущего состояния
                printState();
                functionResult();

                iteration++;
            }

            // Финальный вывод
            cout << "\n=== ИТОГОВЫЙ РЕЗУЛЬТАТ ===" << endl;
            functionResult();
            printState();
        }

        double functionResult() {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < m; i++) {
                    if (table[i][j].isBase) {
                        sum += table[i][j].cost * table[i][j].product;
                    }
                }
            }
            cout << "Сумма целевой функции текущей конфигурации : " << sum << endl;
            return sum;
        }


        void printState() const {
            std::cout << "=== Транспортная Таблица ===" << std::endl;

            // Заголовки столбцов (Потребители)
            std::cout << "      ";
            for (int j = 0; j < n; ++j) {
                std::cout << "B" << j << "(" << demands[j] << ") ";
            }
            std::cout << std::endl;

            // Строки (Поставщики) и ячейки
            for (int i = 0; i < m; ++i) {
                std::cout << std::fixed << std::setprecision(1) << "A" << i << "(" << supplies[i] << ") ";
                for (int j = 0; j < n; ++j) {
                    // Формат: [Цена:Поток*] где * означает базисность
                    std::cout << "["
                              << table[i][j].cost << ":"
                              << table[i][j].product
                              << (table[i][j].isBase ? "*" : " ") << "] ";
                }
                std::cout << std::endl;
            }
            std::cout << "============================" << std::endl;
            }
        };
    }

using namespace szlab;

int main(void) {
    bool balanced_table = true;
    TransportProblem problem3 = TransportProblem(4, 5);
    problem3.table = {
        {1, 6, 4, 8, 7},
        {2, 4, 2, 6, 3},
        {3, 1, 1, 3, 2},
        {2, 6, 8, 9, 6}
    };
    problem3.supplies = {8, 6, 7, 8};
    // сбалансированное условие:
    if (balanced_table) problem3.demands = {9, 5, 7, 5, 3};
    if (!balanced_table) problem3.demands = {10, 6, 8, 6, 4};
    problem3.balanceTable();
    // несбалансированное условие:
    problem3.printState();
    problem3.NorthwestAngle();
    problem3.printState();
    problem3.PotentialsMethod();
    return 0;
}
