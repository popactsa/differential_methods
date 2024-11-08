#include "Solver.h"
#include "Writer.h"

int main()
{
    double CFL = 0.5;
    double l = 1;
    int N = 100;
    double t_end = 0.1;
    Boundary bound1 = {"IN",  0.5, 0.1, 0, 0, 0, 1.0, 0};
    Boundary bound2 = {"OUT", 1, 0, 0, 0, 0, 0, 8.31};
    // Boundary bound3 = {"WALL", 1, 0, 0, 0, 0, 0, 8.31};
    Boundary* bounds = new Boundary[2];
    bounds[0] = bound1;
    bounds[1] = bound2;

    Parameters parameters = {0, 1.0, 200, 0.1, 10, "WENO3", 1, 10};// параметры сетки и расчета
    Manager manager(parameters);
    manager.set_boundary_conditions(bounds);    // граничные условия
    manager.create_mesh();                      // создание сетки
    manager.create_fields();                    // выделение памяти и обьявление массивов
    manager.set_initial_conditions();           // задание начальных условий


    Solver solver(manager);

    while(solver.manager.t < solver.manager.parameters.t_end && (solver.manager.t_done < 100))
    {
        solver.get_time_step();
        // std::cout << solver.manager.dt << std::endl;
        solver.apply_boundary_conditions();     // граничные условия
        solver.solve_problem();
        write_data(solver, solver.manager.t_done);           // сохраниение данных если надо
    }

    int fict = solver.manager.fict;
    // double* c = solver.manager.c;
    /*for(int i = fict; i < solver.manager.N_all-fict; i++)
    {
        std::cout << c[i] << std::endl;
    }*/
}
