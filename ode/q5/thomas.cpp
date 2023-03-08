#include <library.h>

int main()
{
    // read simulation parameters:
    auto parameters = read_input("input");
    scalar S = parameters[0];

    // create grid pts:
    int N = static_cast<int>(parameters[1]);
    scalar dr = 1.0 / static_cast<scalar>(N - 1);
    std::vector<scalar> r = linspace(0.0, 1.0, N);

    // initialize matrix, rhs vector and solution vector:
    std::vector<scalar> rhs(N, -S*dr*dr), sol(N, 0.0), b(N, 0.0), d(N, 0.0), a(N, 0.0);

    // set rhs for boundary nodes: 
    rhs[0] = 0.0;            //left boundary
    rhs[N - 1] = 1.0;        //right boundary

    // assemble linear system for interior nodes:
    for (int i = 1; i < N -1; ++i){
        d[i] = -2.0;
		b[i] = 1.0 - (0.5*dr)/r[i];
		a[i] = 1.0 + (0.5*dr)/r[i];
    }

    //assemble for left boundary:
    d[0] = -1.0;
    a[0] = 1.0;

    //assemble for right boundary:
    d[N-1] = 1.0;

    thomas(b, d, a, rhs, sol, N);

     // write x vs T for x:
    write_data write(std::to_string(S) + ".dat");
    for (int i = 0; i < N; ++i)
        write(r[i], sol[i]);


    return 0;
}
