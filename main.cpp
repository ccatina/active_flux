#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

static double N = 10;
static double M = 20;

double fi(double x)
{
    if (x < 0 || x > 2)
        return 0;
    if ( x >= 0 && x <= 1 )
        return x;
    if ( x > 1 && x <= 2)
        return -x + 2;
//    if ((x <= 0) || (x >= 2))
//        return 0;
//    if (fabs(x - 0.1) > 0 && x <= 1.5)
//        return -x + 2;
//    if ((x > 1.5) && (fabs(x-2.9) > 0))
//        return x;

}



void left_angle(double h, double c, double tau)
{
    double gamma = c*tau/h;
    std::ofstream fout;
    fout.open("left_angle.txt");
    std::vector<double> x(M/h);
    x[0] = 0;
    for (int i = 1; i < M/h; i++)
        x[i] = x[i-1] + h;

    std::vector<double> y0(M/h);
    std::vector<double> y(M/h);

    for (int i = 0; i < M/h; i++)
        y0[i] = fi(x[i]);
    for (int i = 0; i < N/tau; i++) {
        y[0] = 0;
        for (int j = 1; j < M / h; j++)
            y[j] = (1 - gamma) * y0[j] + gamma * y0[j - 1];

        for(int j = 0; j != y.size(); j++ )
        {
            fout << y[j] << '\t';
        }
        fout << '\n';
        y0 = y;
    }
}

void linear_active_flux(double h, double c, double tau)
{
    std::vector<double> y0(M/h);
    std::vector<double> y(M/h);
    std::vector<double> x(M/h);
    std::vector<double> y_2(M/h + 1);
    std::vector<double> y_22(M/h +1);
    x[0] = 0;
    for (int i = 1; i < M/h; i++)
            x[i] = x[i - 1] + h;
    for (int i = 1; i < M/h; i++)
        y0[i] = fi(x[i]);
    for (int i = 0; i < (M/h-1); i++)
        y_2[i] = fi(x[i] + h/2);
    double v = c*tau/h;
    std::ofstream fout;
    fout.open("/home/irina/Documents/MetVich/Active_flux/left_angle.txt");
    for(int i = 0; i < N/tau; i++)
    {
        for(int j = 0; j != y.size(); j++ )
            fout << y0[j] << '\t';
        fout << '\n';
        std::cout << "New time level " << i<<"\n";
        // перезапись промежуточных ячеек
        for(int k = 0; k < (M/h-1); k++)
        {
            y_22[k] = 6*v*(1-v)*y0[k];
            if (k > 0)
                y_22[k] += v*(3*v -2)*y_2[k-1];

                y_22[k] += (1-v)*(1-3*v)*y_2[k+1];
        }
        //обновление вершин сетки
        for(int j = 0; j < M/h; j++)
        {
            y[j] = (1-v)*(1-v)*(1 + 2 * v) * y0[j];
                if (j > 0)
                    y[j] +=  v*v*(3-2*v) * y0[j-1] + v*(1-v)*y_2[j-1];;
                if (j >= 2)
                    y[j] += v*v*(v-1)*y_2[j-2];
                if (j !=  (M/h-1))
                    y[j] -=  v*(1-v)*(1-v)*y_2[j];
        }
        y0 = y;
        y_2 = y_22;
    }
    fout.close();
}

int main() {
    double h = 0.1;
    double c = 0.5;
    double tau = 0.1;
    //left_angle(h,c,tau);
    linear_active_flux(h,c,tau);

    return 0;
}
