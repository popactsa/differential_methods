#include "Solver.h"
#include "Writer.h"

int main()
{
    double CFL = 0.5;
    double l = 1;
    int N = 120;
    double t_end = 0.1;
    Boundary bound = {"FREE", 0.0};
    Boundary* bounds = new Boundary[2];
    for(int i = 0; i < 2; i++)
        bounds[i] = bound;

    Parameters parameters = {1, 120, 0.5, 0.1, "LW", 1, 2};
    // parameters = Parameters({l, N, CFL, 0.1*N/CFL});         // параметры сетки и расчета
    Manager manager(parameters);
    manager.set_boundary_conditions(bounds);    // граничные условия
    manager.create_mesh();                      // создание сетки
    manager.create_fields();                    // выделение памяти и обьявление массивов
    manager.set_initial_conditions();           // задание начальных условий


    Solver solver(manager);

    while(solver.t < solver.manager.parameters.t_end)
    {
        solver.get_time_step();
        solver.solve_problem();
        solver.apply_boundary_conditions();     // граничные условия
        write_data(solver, solver.M);           // сохраниение данных если надо
    }

    int fict = solver.manager.fict;
    double* c = solver.manager.c;
    for(int i = fict; i < solver.manager.N_all-fict; i++)
    {
        std::cout << c[i] << std::endl;
    }
}
