#include "WENO3.h"
#include "Writer.h"

int main()
{
    Boundary bound1 = {IN,  0.5, 1, 0, 0, 0, 0.1, 0};
    Boundary bound2 = {OUT, 1, 0, 0, 0, 0, 0, 0.1*8.31};
    Boundary bound3 = {WALL, 1, 0, 0, 0, 0, 0, 8.31};
    Boundary* bounds = new Boundary[2];
    bounds[0] = bound1;
    bounds[1] = bound2;

    // Parameters parameters = {0, 1.0, 200, 0.1, 10, "WENO3", 1, 10};// параметры сетки и расчета
    std::string file_name = "Parameters_file.txt";
	Parameters parameters(file_name);

    Manager manager(parameters);
    manager.set_boundary_conditions(bounds);    // граничные условия
    manager.create_mesh();                      // создание сетки
    manager.create_fields();                    // выделение памяти и обьявление массивов
    manager.set_initial_conditions();           // задание начальных условий


    // Solver* solver = new WENO3(manager);
    WENO3 solver(manager);

    while(solver.manager.t < solver.manager.parameters.t_end && (solver.manager.t_done < 25))
    {
        solver.get_time_step();
        // std::cout << solver.manager.dt << std::endl;
        solver.apply_boundary_conditions();     // граничные условия
        solver.solve_problem();
        write_data(solver, solver.manager.t_done);           // сохраниение данных если надо
    }

    // int fict = solver.manager.fict;
    // double* c = solver.manager.c;
    /*for(int i = fict; i < solver.manager.N_all-fict; i++)
    {
        std::cout << c[i] << std::endl;
    }*/
    // delete(solver);
}
