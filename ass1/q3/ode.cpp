// Solve spring mass damper system using the RK4 method:

#include <library.h>

int main()
{
    // read simulation parameters:
    auto parameters = read_input("input");
    scalar m = parameters[0];
    scalar c = parameters[1];
    scalar k = parameters[2];
    scalar t0 = parameters[3];
    scalar tf = parameters[4];
    scalar dt = parameters[5];
    scalar x0 = parameters[6];
    scalar v0 = parameters[7];

    // create write object: TODO change filename here acoording to python script:
    write_data write(std::to_string(dt) + ".dat");

    // compute beta and omega sqr:
    scalar beta = c / m;
    scalar omega_sqr = k / m;

    // compute no of time steps:
    int nt = static_cast<int>((tf - t0) / dt);

    // define f(y, t) for spring mass damper system:
    auto f = [beta, omega_sqr](const state<scalar> &y, scalar t)
    {
        return state<scalar>{y.v, -beta * y.v - omega_sqr * y.x};
    };

    // set up initial state:
    state<scalar> y{x0, v0};

    // evolve ode:
    for (int n = 0; n < nt; ++n)
    {

        // compute t
        scalar t = n * dt;

        // write t and x to file
        write(t, y.x);

        // update state using rk4 method
        rk4_update(f, y, t, dt);
    }

    return 0;
}
