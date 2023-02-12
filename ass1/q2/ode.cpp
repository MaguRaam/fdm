#include <library.h>

// ode integrator:
enum Scheme
{
    euler,
    rk4,
    ab2
};

int main()
{

    // read simulation parameters from input file:
    auto parameters = read_input("input");

    Scheme scheme = static_cast<Scheme>(parameters[0]);
    scalar y0 = parameters[1];
    scalar t0 = parameters[2];
    scalar tf = parameters[3];
    scalar dt = parameters[4];

    // compute no of time steps:
    int nt = static_cast<int>((tf - t0) / dt);

    // define f(y, t) system:
    auto f = [](scalar y, scalar t)
    { return -y; };

    // set up solvers:
    // euler:
    auto solve_euler = [=]()
    {
        write_data write("euler_" + std::to_string(dt) + ".dat");

        // set up initial state:
        scalar y = y0;
        scalar t = t0;

        // update:
        for (int n = 0; n < nt; ++n)
        {
            write(t, y);
            euler_update(f, y, t, dt);
            t += dt;
        }
    };

    // rk4:
    auto solve_rk4 = [=]()
    {
        write_data write("rk4_" + std::to_string(dt) + ".dat");

        // set up initial state:
        scalar y = y0;
        scalar t = t0;

        // update:
        for (int n = 0; n < nt; ++n)
        {
            write(t, y);
            rk4_update(f, y, t, dt);
            t += dt;
        }
    };

    // ab2:
    auto solve_ab2 = [=]()
    {
        write_data write("ab2_" + std::to_string(dt) + ".dat");

        // set up initial state:
        scalar y_old = y0;
        scalar t = t0;
        write(t, y_old);
        
        // use rk2 for first time step:
        scalar y_current = y_old;
        rk2_update(f, y_current, t, dt);
        t += dt;

        // use ab2 scheme for future steps:
        for (int n = 1; n < nt; ++n)
        {
            write(t, y_current);
            scalar y_new = ab2_update(f, y_current, y_old, t, dt);
            t += dt;

            y_old = y_current;
            y_current = y_new;

        }

    };


    // choose solvers:
    switch (scheme)
    {
    case euler:
        solve_euler();
        break;

    case rk4:
        solve_rk4();
        break;

    case ab2:
        solve_ab2();

    default:
        break;
    }

    return 0;
}
