#include "Writer.h"

void write_data(Solver& solver, int N_time)
{
    if(N_time%(solver.manager.parameters.n_write) == 0)
    {
        std::string file_name = std::to_string(N_time) + ".csv";
        std::ofstream file(file_name);
        int fict = solver.manager.fict;
        for(int i = fict; i < solver.manager.N_all - fict; i++)
        {
            file << solver.manager.xc[i] << " " << solver.manager.c[i] << "\n";
        }
    }
}
