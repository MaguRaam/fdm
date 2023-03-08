// solve the given ode using implicit euler scheme:
#include <library.h>

inline void implicit_euler_update(scalar& y, const scalar x, const scalar dx){

    const double c = 200000.0; 
    const double xnew = x + dx;

    y = (y + (c*dx - 1)*exp(-xnew))/(1.0 + c*dx);

}


int main()
{
    // read simulation parameters from input file:
    auto parameters = read_input("input");
    const scalar y0  = parameters[0];
    const scalar x0 = parameters[1];
    const scalar xf = parameters[2];
    const scalar dx = parameters[3];


    //set initial condition:
    scalar x = x0;
    scalar y = y0;

    // create file object:
    write_data write("sol.dat");

    // evolve ode:
    while (x < xf){
        write(x, y);

        // explicit euler:
        //euler_update([](scalar y, scalar x){return -200000*y + 200000*exp(-x) - exp(-x);}, y, x, dx);
        
        // impicit euler:
        implicit_euler_update(y, x, dx);

        x += dx;

    }


    return 0;
}


