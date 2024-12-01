#include "Solver_Godunov1D.h"
#include <cmath>
//#include <iostream>
//#include <string>



Solver_Godunov1D::Solver_Godunov1D(const Parameters& _par): par(_par)
{
    // Выделение памяти
    
    // Переменные
    p = new double[par.nx_all];
    rho = new double[par.nx_all];
    rho_u = new double[par.nx_all];
    rho_e = new double[par.nx_all];
    //U = new double[par.nx_all];
    //m = new double[par.nx_all];
    //v = new double[par.nx_all + 1];
    
    x = new double[par.nx_all + 1];
    
    //omega = new double[par.nx_all];   

    // Потоки и сетка
    F_m = new double[par.nx + 1];
    F_imp = new double[par.nx + 1];
    F_e = new double[par.nx + 1];
    
    //char ch;
    
    // Начальные условия
    set_initial_conditions();
    

    t = 0.0;

    for (step = 1; step <= par.nt; step++)
    {
        // Решение на текущем шаге
        solve_step();

        // Запись в файл
        if (step % par.nt_write == 0)
        {
            /*std::cout << "====" << step << "====" << std::endl;
            std::cout << "dt\t" << dt << std::endl;
            std::cout << "p\t" << p[par.nx / 2] << std::endl;
            std::cout << "rho\t" << rho[par.nx / 2] << std::endl;
            std::cout << "U\t" << U[par.nx / 2] << std::endl;
            std::cout << "m\t" << m[par.nx / 2] << std::endl;
            std::cout << "v\t" << v[par.nx / 2] << std::endl;
            std::cout << "x\t" << x[par.nx / 2] << std::endl;
            std::cout << std::endl;*/
            write_data();
        };
//		std::cin >> ch;
    };

    // Конец расчета
    std::cout << "Done" << std::endl;
};



Solver_Godunov1D::~Solver_Godunov1D()
{
    delete[] p;
    delete[] rho;
    delete[] rho_u;
    delete[] rho_e;

    delete[] x;
    
    delete[] F_m;
    delete[] F_e;
    delete[] F_imp;
    //omega;
};


// TODO переписать под выбор типа граничного условия
void Solver_Godunov1D::apply_boundary_conditions()
{
    // Свободный поток
    for (int i = 0; i < par.nx_fict; ++i)
    {
        rho[i] = rho[par.nx_fict];
        rho[par.nx_all - i - 1] = rho[par.nx_all - par.nx_fict - 1];

        rho_u[i] = rho_u[par.nx_fict];
        rho_u[par.nx_all - i - 1] = rho_u[par.nx_all - par.nx_fict - 1];
        
        rho_e[i] = rho_e[par.nx_fict];
        rho_e[par.nx_all - i - 1] = rho_e[par.nx_all - par.nx_fict - 1];

        p[i] = p[par.nx_fict];
        p[par.nx_all - i - 1] = p[par.nx_all - par.nx_fict - 1];
    }

    
    F_m[0] = rho_u[0];
    F_m[par.nx] = rho_u[par.nx_all - 1];
        
    F_imp[0] = pow(rho_u[0], 2) / rho[0] + p[0];
    F_imp[par.nx] = pow(rho_u[par.nx_all - 1], 2) / rho[par.nx_all - 1] + p[par.nx_all - 1];
     
    F_e[0] = (rho_e[0] + p[0]) * rho_u[0] / rho[0];
    F_e[par.nx] = (rho_e[par.nx_all - 1] + p[par.nx_all - 1]) * rho_u[par.nx_all - 1] / rho[par.nx_all - 1];
    
    // Скорость
    // for (unsigned int i = 0; i < 2; ++i)   // 1D -> 2 walls
    // {
    //     // Стенка
    //     if (par.walls[i].type == WALL)
    //     {
    //         if(par.ic_preset == IC_TEST1 || par.ic_preset == IC_TEST3  || par.ic_preset == IC_TEST4 )
    //         {
    //             v[par.nx_fict] = 0.0;
    //             v[par.nx_all - par.nx_fict] = 0.0;
    //         }
    //         else if(par.ic_preset == IC_TEST2)
    //         {
    //             v[par.nx_fict] = -2.0;
    //             v[par.nx_all - par.nx_fict] = 2.0;
    //         }
    //         else
    //         {
    //             std::cerr << "For nx_fict = " << par.nx_fict << " bc is not specified" << std::endl;
    //             exit(1);
    //         };
    //         int i_temp = 0;
    //         if(i)
    //             i_temp = par.nx_all;
    //         v[i_temp] = -v[i_temp+(1-2*i)*2] + 2.0 * v[i_temp+(1-2*i)*1];
    //         // v[0] = -v[2] + 2.0 * v[1];
    //         // v[par.nx_all] = -v[par.nx_all - 2] + 2.0 * v[par.nx_all - 1];
    //     }
    //     
    //
    //     // Свободный поток
    //     else if (par.walls[i].type == FLUX)
    //     {
    //         
    //
    //         int i_temp = 0;
    //         if(i)
    //             i_temp = par.nx_all;
    //         v[i_temp] = v[i_temp+(1-2*i)*2];
    //         v[i_temp+(1-2*i)*1] = v[i_temp+(1-2*i)*2];
    //     }
    //}

    // Плотность и внутренняя энергия
    // for (unsigned int i = 0; i < par.nx_fict; ++i)
    // {
    //     rho[i] = rho[par.nx_fict];
    //     rho[par.nx_all - par.nx_fict - i] = rho[par.nx_all - par.nx_fict - 1];
    //     U[i] = U[par.nx_fict];
    //     U[par.nx_all - par.nx_fict - i] = U[par.nx_all - par.nx_fict - 1];
    // };
};



void Solver_Godunov1D::solve_step()
{
    // Пересчитать шаг по времени
    get_time_step();
   

    // Общее время
    t += dt;


    // Применить граничные условия
    apply_boundary_conditions();
        

    // Параметры слева и справа на границе ячеек
    double *rho_left = new double[par.nx + 1];
    double *rho_u_left = new double[par.nx + 1];
    double *rho_e_left = new double[par.nx + 1];
    double *p_left = new double[par.nx + 1];
    
    double *rho_right = new double[par.nx + 1];
    double *rho_u_right = new double[par.nx + 1];
    double *rho_e_right = new double[par.nx + 1];
    double *p_right = new double[par.nx + 1];


    for (int i = 0; i < par.nx; ++i)
    {
        int i_temp = i + par.nx_fict;

        rho_left[i + 1] = rho[i_temp];
        rho_u_left[i + 1] = rho_u[i_temp];
        rho_e_left[i + 1] = rho_e[i_temp];
        p_left[i + 1] = p[i_temp];
        
        rho_right[i] = rho[i_temp];
        rho_u_right[i] = rho_u[i_temp];
        rho_e_right[i] = rho_e[i_temp];
        p_right[i] = p[i_temp];
    }

    
    // Решатель Римана - получим параметры на грани ячейки
    for (int i = 0; i < par.nx; ++i)
    {
        int N_iter = 20; // как в Toro, критерия по точности нет

        double p_l = p_left[i];
        double p_r = p_right[i];

        double u_l = rho_u_left[i] / rho_left[i];
        double u_r = rho_u_right[i] / rho_right[i];
        
        double rho_l = rho_left[i];
        double rho_r = rho_right[i];

        double c_l = (rho_l != 0.0) ? sqrt(par.gamma * p_l / rho_l) : 0.0;
        double c_r = (rho_r != 0.0) ? sqrt(par.gamma * p_r / rho_r) : 0.0;
        

        // Условие из Toro, 175 (153). Вывод, 165 (143)  
        if (2.0 / (par.gamma - 1) * (c_l + c_r) <= (u_r - u_l))
        {
            //std::cout << u_r << " " << u_l << " " << 2.0 / (par.gamma - 1) * (c_l + c_r) << std::endl;
            //std::cout << "On time " << t << std::endl;
            std::cerr << "Vaccum!" << std::endl;
        }
    

        double p_s = 0.5 * (p_l + p_r);
 

        // Метод Ньютона - получим давление на контактном разрыве
        double fl, fr;
        double dfl, dfr;
        for (int k = 0; k < N_iter; ++k)
        {
            // TODO переписать в виде вызова функции 'newton_iteration'
            if (p_s <= p_l) // Волна разрежения
            {
                fl = 2 * c_l / (par.gamma - 1) * (pow(p_s / p_l, ((par.gamma - 1) / 2 / par.gamma)) - 1);
                dfl = c_l / par.gamma * pow(p_s / p_l, (par.gamma - 1) / 2 / par.gamma - 1) / p_l;
            }


            else // Ударная волна
            {
                fl = (p_s - p_l) / sqrt(rho_l / 2 * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_l));
                dfl = rho_l / 4 * ((par.gamma + 1) * p_s + (3 * par.gamma + 1) * p_l) / pow(rho_l / 2 * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_l), 1.5);
            }


            if (p_s <= p_r) // Волна разрежения
            {
                fr = 2 * c_r / (par.gamma - 1) * (pow(p_s / p_r, (par.gamma - 1) / 2 / par.gamma) - 1);
                dfr = c_r / par.gamma * pow(p_s / p_r, (par.gamma - 1) / 2 / par.gamma - 1) / p_r;
            }


            else // Ударная волна
            {
                fr = (p_s - p_r) / sqrt(rho_r / 2 * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_r));
                dfr = rho_r / 4 * ((par.gamma + 1) * p_s + (3 * par.gamma + 1) * p_r) / pow(rho_r / 2 * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_r), 1.5);
            }


            p_s = p_s - (fl + fr - (u_l - u_r)) / (dfl + dfr);

            
            if (p_s <= 0.0)
            {
                p_s = 1e-6;
                break;
            }
        }

        // Волна разрежения
        if (p_s <= p_l) fl = 2 * c_l / (par.gamma - 1) * (pow(p_s / p_l, (par.gamma - 1) / 2 / par.gamma) - 1);

        // Ударная волна
        else fl = (p_s - p_l) / sqrt(rho_l / 2 * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_l));

        // Волна разрежения
        if (p_s <= p_r) fr = 2 * c_r / (par.gamma - 1) * (pow(p_s / p_r, (par.gamma - 1) / 2 / par.gamma) - 1);

        // Ударная волна
        else fr = (p_s - p_r) / sqrt(rho_r / 2 * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_r));


        double u_s = 0.5 * (u_l + u_r + fr - fl);

        
        // TODO вынести в отдельную функцию 'choose_solution'

        // Отбор решения
        double u_temp = 0.0;
        double p_temp = 0.0;
        double rho_temp = 0.0;
        
        double S = 0.0;
        
        double D;

        u_temp = u_s;
        p_temp = p_s;


        if (S <= u_s) // Левое
        {
            if (p_s > p_l) // Ударная волна
            {
                D = u_l - c_l * sqrt((par.gamma + 1) / 2 / par.gamma * p_s / p_l + (par.gamma - 1) / 2 / par.gamma); // Toro, 182 (161)

                if (S < D) // Невозмущенное
                {
                    rho_temp = rho_l;
                    u_temp = u_l;
                    p_temp = p_l;
                }


                else // Перед контактным разрывом с ад. Гюгонию
                {
                    rho_temp = rho_l * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_l) / ((par.gamma - 1) * p_s + (par.gamma + 1) * p_l);
                    u_temp = u_s;
                    p_temp = p_s;
                }
            }


            else // Волна разрежения
            {
                double c_s = c_l * pow(p_s / p_l, (par.gamma - 1) / 2 / par.gamma);
                

                // Левая граница волны разрежения
                double s_h = u_l - c_l;

                // Правая граница волны разрежения
                double s_t = u_s - c_s;


                if (S < s_h) // Невозмущенное
                {
                    rho_temp = rho_l;
                    u_temp = u_l;
                    p_temp = p_l;
                }


                else
                {
                    if (S > s_t) // Снаружи волны разрежения
                    {
                        rho_temp = rho_l * pow(p_s / p_l, 1 / par.gamma);
                        u_temp = u_s;
                        p_temp = p_s;
                    }


                    else // Внутри волны разрежения, Toro, 181 (160)
                    {
                        c_s = 2 / (par.gamma + 1) * (c_l + (par.gamma - 1) * (u_l - S) / 2);
                        rho_temp = rho_l * pow(c_s / c_l, 2 / (par.gamma - 1));
                        u_temp = 2 / (par.gamma + 1) * (c_l + (par.gamma - 1) * u_l / 2 + S);
                        p_temp = p_l * pow(c_s / c_l, 2 * par.gamma / (par.gamma - 1));
                    }
                }
            }
        }


        else // Правое
        {
            if (p_s > p_r) // Ударная волна
            {
                D = u_r + c_r * sqrt((par.gamma + 1) / 2 / par.gamma * p_s / p_r + (par.gamma - 1) / 2 / par.gamma); // 1-е соотношение Гюгонио, Toro, 182 (161)
                
                if (S > D) // Невозмущенное
                {
                    rho_temp = rho_r;
                    u_temp = u_r;
                    p_temp = p_r;
                } 


                else // Перед контактным разрывом с адиаб. Гюгонио
                {
                    rho_temp = rho_r * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_r) / ((par.gamma - 1) * p_s + (par.gamma + 1) * p_r);
                    u_temp = u_s;
                    p_temp = p_s;
                }
            }


            else // Волна разрежения
            {
                double c_s = c_r * pow(p_s / p_r, (par.gamma - 1) / 2 / par.gamma);
                double s_h = u_r + c_r; // Правая граница волны разрежения
                double s_t = u_s + c_s; // Левая граница волны разрежения
                

                if (S > s_h) // Невозмущенное
                {
                    rho_temp = rho_r;
                    u_temp = u_r;
                    p_temp = p_r;
                }


                else
                {
                    if (S < s_t) // Снаружи волны разрежения
                    {
                        rho_temp = rho_r * pow(p_s / p_r, 1 / par.gamma);
                        u_temp = u_s;
                        p_temp = p_s;
                    }


                    else // Внутри волны разрежения, Toro, 183 (161)
                    {
                        c_s = 2 / (par.gamma + 1) * (c_r - (par.gamma - 1) * (u_r - S) / 2);
                        rho_temp = rho_r * pow(c_s / c_r, 2 / (par.gamma - 1));
                        u_temp = 2 / (par.gamma + 1) * (-c_r + (par.gamma - 1) * u_r / 2 + S);
                        p_temp = p_r * pow(c_s / c_r, 2 * par.gamma / (par.gamma - 1));
                    }
                }
            }
        }


        // Поток Годунова
        F_m[i] = rho_temp * u_temp;
        F_imp[i] = rho_temp * pow(u_temp, 2) + p_temp;
        F_e[i] = (p_temp + rho_temp * pow(u_temp, 2) / 2.0 + p_temp / (par.gamma - 1)) * u_temp;
    }



    // Пересчет значений
    for (int i = par.nx_fict; i < par.nx + par.nx_fict; ++i)
    {
        rho[i] = rho[i] - (F_m[i] - F_m[i - 1]) * dt / par.dx;
        rho_u[i] = rho_u[i] - (F_imp[i] - F_imp[i - 1]) * dt / par.dx;
        rho_e[i] = rho_e[i] - (F_e[i] - F_e[i - 1]) * dt / par.dx;
        p[i] = (par.gamma - 1) * (rho_e[i] - pow(rho_u[i], 2) / rho[i] / 2.0);
    }


    // double v_last[par.nx_all + 1];
    // for(unsigned int i = 0; i < par.nx_all + 1; i++)
    // {
    //     v_last[i] = v[i];
    // }
    // 
   //  // Вязкость
   //  if (par.viscosity == VISC_NONE)
   //  {
   //      for (unsigned int i = 0; i < par.nx_all; ++i)
   //          omega[i] = 0.0;
   //  }
   //  else if (par.viscosity == VISC_NEUMAN)
   //  {
   //      for (unsigned int i = 0; i < par.nx_all; ++i)
   //          omega[i] = -par.mu0*rho[i]*fabs(v[i+1] - v[i])*(v[i+1] - v[i]);
   //  }
   //  else if (par.viscosity == VISC_LATTER)
   //  {
   //      for (unsigned int i = 0; i < par.nx_all; ++i)
   //          omega[i] = ((v[i+1] - v[i]) < 0.0) ? par.mu0*rho[i]*(v[i+1] - v[i])*(v[i+1] - v[i]) : 0.0;
   //  }
   //  else if (par.viscosity == VISC_LINEAR)
   //  {
   //      for (unsigned int i = 0; i < par.nx_all; ++i)
	  // {
   //          omega[i] = par.mu0 * rho[i] * (v[i + 1] - v[i]) * m[i];
	  // }
   //  }
   //  else if (par.viscosity == VISC_SUM)
   //  {
   //      for (unsigned int i = 0; i < par.nx_all; ++i)
	  // {
   //          omega[i] = par.mu0 * rho[i] * (v[i + 1] - v[i]) * m[i] - par.mu0 * rho[i]*fabs(v[i+1] - v[i])*(v[i+1] - v[i]);
	  // }
   //  }
   //  else
   //  {
   //      std::cerr << "Viscosity is not set: " << par.viscosity << std::endl;
   //      exit(1);
   //  };
    // velocity
    // for (unsigned int i = 1 + par.nx_fict; i < par.nx_all - par.nx_fict; ++i)
    // {
    //     double dx = 0.5 * (m[i] + m[i - 1]);
    //     v[i] = v[i] - ((p[i] + omega[i]) - (p[i - 1] + omega[i - 1])) * dt / dx;
    // };
    // // coordinate
    // for (unsigned int i = 0; i < par.nx_all + 1; ++i)
    // {
    //     x[i] += v[i] * dt;
    // };
    // ////// pressure at the boundaries for conservative law
    // double pc[par.nx_all];
    // for(unsigned int i = 1; i < par.nx_all; i++)
    // {
    //     pc[i] = 0.5*(p[i] + omega[i] + p[i-1] + omega[i-1]);
    // }
    // // density and internal energy
    // for (unsigned int i = par.nx_fict; i < par.nx_all - par.nx_fict; ++i)
    // {
    //     rho[i] = rho[i] / (rho[i] * (v[i + 1] - v[i]) * dt / m[i] + 1.0);
    //     if (par.is_conservative == false)
    //         U[i] = U[i] / (rho[i] * (v[i + 1] - v[i]) * (par.gamma - 1.0) * dt / m[i] + 1.0);
    //     else
    //     {
    //         double U_temp = U[i];
    //         U[i] = U[i] - (v[i+1]*pc[i+1] - v[i]*pc[i])*dt/m[i] + (v_last[i+1] + v_last[i])*(v_last[i+1] + v_last[i])/8.0 - (v[i+1] + v[i])*(v[i+1] + v[i])/8.0;
    //         if (U[i] < 0.0)
    //         {
    //             U[i] = U_temp / (rho[i] * (v[i + 1] - v[i]) * (par.gamma - 1.0) * dt / m[i] + 1.0);
    //         }
    //     }
    // };
    // // pressure
    // for (unsigned int i = 0; i < par.nx_all; ++i)
    // {
    //     p[i] = rho[i] * (par.gamma - 1.0) * U[i];
    //     // std::cout << "p " << p[i] << std::endl;
    // };
};



// TODO переписать под передачу массивов в качестве аргументов
void Solver_Godunov1D::set_initial_conditions()
{
    if (par.ic_preset == IC_TEST1)
    {
        for (unsigned int i = 0; i < par.nx_all; ++i)
        {
            if (i * par.dx <= 0.5 * (par.x_end - par.x_start) + par.nx_fict * par.dx)
            {
                rho[i] = 1.0;
                p[i] = 1.0;
            }
            else
            {
                rho[i] = 0.125;
                p[i] = 0.1;
            };
            
            rho_e[i] = rho[i] * ((pow(0.0, 2)) / 2.0  + p[i] / (par.gamma - 1.0) / rho[i]);
            rho_u[i] = rho[i] * 0.0;
            x[i] = (i + 0.5 - par.nx_fict) * par.dx;
        };
    }


    else if (par.ic_preset == IC_TEST2)
    {
        for (unsigned int i = 0; i < par.nx_all; ++i)
        {
            rho[i] = 1.0;


            if (i * par.dx <= 0.5 * (par.x_end - par.x_start) + par.nx_fict * par.dx) rho_u[i] = rho[i] * (-2.0);

            else rho_u[i] = rho[i] * 2.0;
            

            p[i] = 0.4;
            rho_e[i] = rho[i] * pow(2.0, 2) / 2.0 + p[i] / (par.gamma - 1.0);
            x[i] = (i + 0.5 - par.nx_fict) * par.dx;
        };
    }
    

    else if (par.ic_preset == IC_TEST3)
    {
        for (unsigned int i = 0; i < par.nx_all; ++i)
        {
            if (i * par.dx <= 0.5 * (par.x_end - par.x_start) + par.nx_fict * par.dx) p[i] = 1000.0;

            else p[i] = 0.01;
            

            rho[i] = 1.0;
            rho_e[i] = rho[i] * ((pow(0.0, 2)) / 2.0  + p[i] / (par.gamma - 1.0) / rho[i]);
            rho_u[i] = rho[i] * 0.0;
            x[i] = (i + 0.5 - par.nx_fict) * par.dx;
        };
    }
    

    // else if (par.ic_preset == IC_TEST4)
    // {
    //     for (unsigned int i = 0; i < par.nx_all; ++i)
    //     {
    //         if (i * par.dx <= par.x_start + 0.5 * (par.x_end - par.x_start) + par.nx_fict * par.dx)
    //         {
    //             p[i] = 0.01;
    //             rho[i] = 1.0;
    //         }
    //         else
    //         {
    //             p[i] =  100.0;
    //             rho[i] = 1.0;
    //         }
    //         U[i] = p[i] / (par.gamma - 1.0) / rho[i];
    //     };
    //     for (unsigned int i = 0; i < par.nx_all + 1; ++i)
    //     {
    //         if (i * par.dx <= par.x_start + 0.5 * (par.x_end - par.x_start) + par.nx_fict * par.dx)
    //         {
    //             v[i] = 0.0;
    //         }
    //         else
    //         {
    //             v[i] = 0.0;
    //         };
    //         x[i] = par.x_start - par.nx_fict * par.dx + i * par.dx;
    //     };
    // }
    

    else
    {
        std::cerr << "Initial condition is not set: " << par.ic_preset << std::endl;
        exit(1);
    };
    //for (unsigned int i = 0; i < par.nx_all; ++i) m[i] = rho[i] * (x[i + 1] - x[i]);
};



void Solver_Godunov1D::get_time_step()
{
    double min_dt = 1.0e6;
    for (unsigned int i = 1; i < par.nx_all; ++i)
    {
        double u = rho_u[i] / rho[i];
        double c = (rho[i] != 0.0) ? c = std::sqrt(par.gamma * p[i] / rho[i]) : 0.0;
        //double dx = x[i+1] - x[i];
        //double V = 0.5*(v[i+1] + v[i]);
        double dt_temp = par.CFL * par.dx / (c + fabs(u));
        if (dt_temp < min_dt) min_dt = dt_temp;
    }
    dt = min_dt;
    // std::cout << "dt = " << dt << std::endl;
}



void Solver_Godunov1D::write_data()
{
    std::string file_name = "data/" + std::to_string(step) + ".csv";
    std::ofstream file(file_name);
    

    //double x_grid[par.nx_all];
    double u[par.nx_all];
    

    for(unsigned int i = 0; i < par.nx_all; i++)
    {
         
         u[i] = rho_u[i] / rho[i];
    }
     

    //double* c = rho;
    

    for(int i = par.nx_fict; i < par.nx_all - par.nx_fict; i++)
    {
        file << x[i] << " " << rho[i]  << " " << p[i] << " " << u[i] << "\n";
    }
}
