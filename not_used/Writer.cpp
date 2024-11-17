#include "Writer.h"

void write_data(WENO3& solver, int N_time)
{
    if(N_time%(solver.manager.parameters.nt_write) == 0)
    {
        std::string file_name = "data/" + std::to_string(N_time) + ".csv";
        std::ofstream file(file_name);
        int fict = solver.manager.fict;
        double* c = solver.manager.rho;
        for(int i = fict; i < solver.manager.N_all - fict; i++)
        {
            file << solver.manager.xc[i] << " " << c[i] << "\n";
        }
    }
}
