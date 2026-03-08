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

        Dot findUnopt() { // Найдем первую свободную ячейку, нарушающую условие оптимальности.
            int entering_row = -1, entering_col = -1;
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < m; i++) {
                    if (!table[i][j].isBase) {
                        if (u[i] + v[j] > table[i][j].cost) {
                            entering_row = i;
                            entering_col = j;
                            break;
                        }
                    }
                }
                if (entering_row != -1) break;
            }
            if (entering_row == -1) return {-1,-1}; // решение оптимально
            Dot value(entering_row, entering_col);
            return value;
        }

        void optiCycle(Dot matrix_pos) {
            int i_opt = matrix_pos.x;
            int j_opt = matrix_pos.y;
            // счетчик для одинокости
            int counter = 0;
            // угловой вектор
            vector<Dot> angleVec = {};
            vector<vector<double>> basisplus(m, vector<double>(n, 0.0)); // вектор базисов + точки неоптимизации, хранящий cost таблицы

            basisplus[i_opt][j_opt] = table[i_opt][j_opt].cost; // добавили точку неоптимизации
            for (int basis_index = 0; basis_index < (n * m - 1); basis_index++) {
                basisplus[u_indexes[basis_index]][v_indexes[basis_index]] = table[u_indexes[basis_index]][v_indexes[basis_index]].cost;
            }
            // угловой называем полную ячейку участвующую в цикле пересчёта
            // у угловых ячеек есть замечательное свойство - они являются "одинокими"
            // это значит что кроме них нет никаких ячеек или в строке, или в столбце.
            // найдем все одинокие базисные ячейки, включая точку не-оптимума
            // проверяем столбцы и строки по базисным векторам

            // столбцы
            for (int i = 0; i < m; i++) {
                counter = 0;
                Dot buffer = {-1,-1};
                for (int j = 0; j < n; j++) {
                    if (basisplus[i][j] != 0) {
                        counter++;
                        buffer.x = i;
                        buffer.y = j;
                    }
                }
                if (counter == 1) {
                    basisplus[buffer.x][buffer.y] = 0.0;

                }
            }
            // строки
            for (int j = 0; j < n; j++) {
                counter = 0;
                Dot buffer = {-1,-1};
                for (int i = 0; i < m; i++) {
                    if (basisplus[i][j] != 0) {
                        counter++;
                        buffer.x = i;
                        buffer.y = j;
                    }
                }
                if (counter == 1) {
                    basisplus[buffer.x][buffer.y] = 0.0;
                }
            }
            // теперь работаем только с угловыми клетками в basisplus
            // по ним можно однозначно придти к нужной точке
            // создаем угловой вектор и вычисляем
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    if(basisplus[i][j]) angleVec.push_back({i,j});
                }
            }
        }

        void PotentialsMethod() {
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < m; i++) {
                    if (!table[i][j].isBase) {
                        if (u[i]+v[j]>table[i][j].cost) {
                            //RecalcCycle();
                            return;
                        }
                    }
                }
            }
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
    return 0;
}
