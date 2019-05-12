#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
enum test
{
    LOL,
    LOL1
};
int main()
{
    double alpha;
    int j;
    double step_alpha = 0.5;
    FILE * output = fopen("Test.data", "w");
    for(alpha = 0.; alpha <= 25.; alpha += step_alpha)
    {
        fprintf(output, "\\subsection{$\\alpha = %e$} \n", alpha);
        fprintf(output, "Метод Ньютона сошелся с большого количества начальных точек к результатам, указанным ниже. \\\\ \n");
        fprintf(output, "\\begin{tabular}{ | c | c | c | c |} \n");
        fprintf(output, "\\hline \n");
        fprintf(output, "Initial point  & Newton method result & Right point & Functional \n \\\\ $x_1(0), \\; p_2(0)$ & $x_1(0), \\; p_2(0)$ & $x_1\\l(\\frac{\\pi}{2}\\r) - 0.8, \\; p_2\\l(\\frac{\\pi}{2}\\r)$ & $\\int_{0}^{\\frac{\\pi}{2}}u^2dt$  \\\\ \\hline \n");

        for(j = 0; j < 40; j++)
        {
            fprintf(output, "%e; %e & %e; %e & %e; %e & %e \\\\ \\hline \n", alpha, alpha, alpha, alpha, alpha, alpha, alpha);
        }
        fprintf(output, "\\end{tabular} \n");
    }
    fclose(output);
}