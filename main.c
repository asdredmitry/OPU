#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
const double EPS = 1e-50;
const double H_MIN =  __DBL_MAX__;

const unsigned int FULL_NORM = 1;
const unsigned int LOL = (1 << 1);
const unsigned int LAST_NORM = 1 << 2;
const unsigned int WRITE_IN_FILE = 1 << 3;
const unsigned int GET_LAST = 1 << 5;
const unsigned int PRINT_H = 1 << 6;
const unsigned int PRINT_GLOBAL_ERROR = 1 << 7;

const unsigned int WRITE_PERIODS = 1;
const unsigned int WRITE_ZEROES = 1 << 2;
const int MAX_STEPS_SHOOTING = 100000;


const double fac = 0.8;
const double facmin = 0.0;
const double facmax = 1.5;

const double alpha = 0.;

double c[13] = {0., 1./18., 1./12., 1./8., 5./16., 3./8., 59./400., 93./200., 5490023248./9719169821., 13./20., 1201146811./1299019798., 1., 1.};
double a[13][12] = {{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},\
                    {1./18., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},\
                    {1./48., 1./16., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},\
                    {1./32., 0., 3./32., 0., 0., 0., 0., 0., 0., 0., 0., 0.},\
                    {5./16., 0., -75./64, 75./64., 0., 0., 0., 0., 0., 0., 0., 0.},\
                    {3./80., 0., 0., 3./16., 3./20., 0., 0., 0., 0., 0., 0., 0.},\
                    {29443841./614563906., 0., 0., 77736538./692538347., -28693883./1125000000., 23124283./1800000000., 0., 0., 0., 0., 0., 0.},\
                    {16016141./946692911., 0., 0., 61564180./158732637., 22789713./633445777., 545815736./2771057229., -180193667./1043307555., 0., 0., 0., 0., 0.},\
                    {39632708./573591083., 0., 0., -433636366./683701615., -421739975./2616292301., 100302831./723423059., 790204164./839813087., 800635310./3783071287., 0., 0., 0., 0.},\
                    {246121993./1340847787., 0., 0., -37695042795./15268766246., -309121744./1061227803., -12992083./490766935., 6005943493./2108947869., 393006217./1396673457., 123872331./1001029789., 0., 0., 0.},\
                    {-1028468189./846180014., 0., 0., 8478235783./508512852., 1311729495./1432422823., -10304129995./1701304382., -48777925059./3047939560., 15336726248./1032824649., -45442868181./3398467696., 3065993473./597172653., 0., 0.},\
                    {185892177./718116043., 0., 0., -3185094517./667107341., -477755414./1098053517., -703635378./230739211., 5731566787./1027545527, 5232866602./850066563., -4093664535./808688257., 3962137247./1805957418., 65686358./487910083., 0},\
                    {403863854./491063109., 0., 0., -5068492393./434740067., -411421997./543043805., 652783627./914296604., 11173962825./925320556., -13158990841./6184727034., 3936647629./1978049680., -160528059./685178525., 248638103./1413531060., 0.}};
double b[13] = {14005451./335480064., 0., 0., 0., 0., -59238493./1068277825., 181606767./758867731., 561292985./797845732., -1041891430./1371343529., 760417239./1151165299., 118820643./751138087., -528747749./2220607170., 1./4.};
double b_[13] = {13451932./455176623., 0., 0., 0., 0., -808719846./976000145., 1757004468./5645159321., 656045339./265891186., -3867574721./1518517206., 465885868./322736535., 53011238./667516719., 2./45., 0.};

typedef struct 
{
    double* data;
    int alloc;
    int used;
}vector;

void init(vector* v, int n);
void reall(vector* v);
void push_back(vector* v, double val);
void free_vec(vector* v);

double max(double x1, double x2);
double min(double x1, double x2);
double fmax(double x1, double x2);
double fmin(double x1, double x2);

double f1(double t, double x1, double x2, double p1, double p2, double integr);
double f2(double t, double x1, double x2, double p1, double p2, double integr);
double f3(double t, double x1, double x2, double p1, double p2, double integr);
double f4(double t, double x1, double x2, double p1, double p2, double integr);
double f5(double t, double x1, double x2, double p1, double p2, double integr);

double check_sol(double t);

double norm(double x1, double x2, double p1, double p2, double integr);
double norm2(double x1, double x2);
double get_h(double h, double err, double tol);
int inverse_matrix(double * data, double * inversed);


double dorman_prince(double h, double t,  double x1_ic, double x2_ic, double p1_ic, double p2_ic, double integr_ic, double* x1, double* x2, double* p1, double* p2, double* integr);
void solve_dp(double l, double r, double tol, double x1_ic, double x2_ic, double p1_ic, double p2_ic, double integr_ic, double* x1, double* x2, double* p1, double* p2, double* integr, unsigned int attr, const char* name);
//int shooting_method(double l, double r, double tol, double shooting_method_accuracy, double x1_ic_l, double p1_ic_l, double x1_ic_r, double x2_ic_r, double* x1_ic_answer_l, double *x2_ic_answer_l, double* p1_ic_answer_l, double* p2_ic_answer_l, double* integral_ic_answer_l);
int shooting_method(double l, double r, double tol, double shooting_method_accuracy, double x1_l, double x2_l, double p1_l, double p2_l, double x1_r, double x2_r, double* x1_ic_answer_l, double *x2_ic_answer_l, double* p1_ic_answer_l, double* p2_ic_answer_l, double* integral_ic_answer_l);
int check_shooting_method(double l, double r, double tol, double shooting_method_accuracy, double x1_l, double x2_l, double x1_r, double x2_r, double* x1_l_answer, double* x2_l_answer);



int main(void)
{
    double x1, x2, p1, p2, integr;
    double x1_answer, x2_answer;
    solve_dp(0, 10, 1e-10, 1, 0, 0, 0, 0, &x1, &x2, &p1, &p2, &integr, WRITE_IN_FILE, "data.dat");
    check_shooting_method(0, M_PI, 1e-14, 1e-8, 100, -1344, 1, 0, &x1_answer, &x2_answer);
    printf("%lf %lf \n", x1_answer, x2_answer);
    //solve_dp(0, 10, 1e-10, 1, 0, -1, 0, 0, &x1, &x2, &p1, &p2, &integr, WRITE_IN_FILE, "data.dat");
    return 0;
}







double max(double x1, double x2)
{
    return (x1 > x2) ? x1 : x2;
}
double min(double x1, double x2)
{
    return (x1 > x2) ? x2 : x1;
}
double fmax(double x1, double x2)
{
    return max(fabs(x1), fabs(x2));
}
double fmin(double x1, double x2)
{
    return min(fabs(x1), fabs(x2));
}
double f1(double t, double x1, double x2, double p1, double p2, double integr)
{
    return x2;
    return x2;
    return x2;
}
double f2(double t, double x1, double x2, double p1, double p2, double integr)
{
    return cos(t);
    return cos(t);
    return p1;
    return p2/2. - x1*exp(-alpha*x1);
}
double f3(double t, double x1, double x2, double p1, double p2, double integr)
{
    return 0;
    return p2;
    return -p2*(-exp(-alpha*x1) + alpha*x1*exp(-alpha*x1));
}
double f4(double t, double x1, double x2, double p1, double p2, double integr)
{
    return 0;
    return sin(t);
    return -p1;
}
double f5(double t, double x1, double x2, double p1, double p2, double integr)
{
    return 0;
    return p2*p2/4.;
}
double check_sol(double t)
{
    return cos(t);
}
double norm(double x1, double x2, double p1, double p2, double integr)
{
    return fmax(x1, fmax(x2, fmax(p1, fmax(p2, integr))));
}
double norm2(double x1, double x2)
{
    return fmax(x1, x2);
}
double get_h(double h, double err, double tol)
{
    if(err < EPS)
        err = 2*EPS;
    h *= min(facmax, max(facmin, fac*pow(tol/err, 1./7.)));
    if(h < EPS)
        h = 2*EPS;
    return h;
}
int inverse_matrix(double* data, double* res)
{
    double det = 0.;
    det = data[0]*data[3] - data[2]*data[1];
    if(fabs(det) < EPS)
        return 0;
    res[0] = data[3]/det;
    res[1] = -data[2]/det;
    res[2] = -data[1]/det;
    res[3] = data[0]/det;
}
double dorman_prince(double h, double t,  double x1_ic, double x2_ic, double p1_ic, double p2_ic, double integr_ic, double* x1, double* x2, double* p1, double* p2, double* integr)
{
    double tmpt, tmpx1, tmpx2, tmpp1, tmpp2, tmpintegr;
    double kx1[13];
    double kx2[13];
    double kp1[13];
    double kp2[13];
    double kintegr[13];
    double x1_8, x2_8, p1_8, p2_8, integr_8;
    double x1_7, x2_7, p1_7, p2_7, integr_7;
    int i, j;
    kx1[0] = f1(t, x1_ic, x2_ic, p1_ic, p2_ic, integr_ic);
    kx2[0] = f2(t, x1_ic, x2_ic, p1_ic, p2_ic, integr_ic);
    kp1[0] = f3(t, x1_ic, x2_ic, p1_ic, p2_ic, integr_ic);
    kp2[0] = f4(t, x1_ic, x2_ic, p1_ic, p2_ic, integr_ic);
    kintegr[0] = f5(t, x1_ic, x2_ic, p1_ic, p2_ic, integr_ic);
    for(i = 0; i < 13; i++)
    {
        tmpx1 = x1_ic;
        tmpx2 = x2_ic;
        tmpp1 = p1_ic;
        tmpp2 = p2_ic;
        tmpintegr = integr_ic;
        tmpt = t + c[i]*h;
        for(j = 0; j < i; j++)
        {
            tmpx1 += a[i][j]*h*kx1[j];
            tmpx2 += a[i][j]*h*kx2[j];
            tmpp1 += a[i][j]*h*kp1[j];
            tmpp2 += a[i][j]*h*kp2[j];
            tmpintegr += a[i][j]*h*kintegr[j];
        }
        kx1[i] = f1(tmpt, tmpx1, tmpx2, tmpp1, tmpp2, tmpintegr);
        kx2[i] = f2(tmpt, tmpx1, tmpx2, tmpp1, tmpp2, tmpintegr);
        kp1[i] = f3(tmpt, tmpx1, tmpx2, tmpp1, tmpp2, tmpintegr);
        kp2[i] = f4(tmpt, tmpx1, tmpx2, tmpp1, tmpp2, tmpintegr);
        kintegr[i] = f5(tmpt, tmpx1, tmpx2, tmpp1, tmpp2, tmpintegr);
    }
    x1_8 = x2_8 = p1_8 = p2_8 = integr_8 = 0;
    for(i = 0; i < 13; i++)
    {
        x1_8 += b[i]*kx1[i];
        x2_8 += b[i]*kx2[i];
        p1_8 += b[i]*kp1[i];
        p2_8 += b[i]*kp2[i];
        integr_8 += b[i]*kintegr[i];
    }
    x1_8 *= h;
    x2_8 *= h;
    p1_8 *= h;
    p2_8 *= h;
    integr_8 *= h;
    x1_8 += x1_ic;
    x2_8 += x2_ic;
    p1_8 += p1_ic;
    p2_8 += p2_ic;
    integr_8 += integr_ic;
    x1_7 = x2_7 = p1_7 = p2_7 = integr_7 = 0;
    for(i = 0; i < 13; i++)
    {
        x1_7 += b_[i]*kx1[i];
        x2_7 += b_[i]*kx2[i];
        p1_7 += b_[i]*kp1[i];
        p2_7 += b_[i]*kp2[i];
        integr_7 += b_[i]*kintegr[i];
    }
    x1_7 *= h;
    x2_7 *= h;
    p1_7 *= h;
    p2_7 *= h;
    integr_7 *= h;tmpx1 = tmpx2 = tmpp1 = tmpp2 = tmpintegr = 0;
    x1_7 += x1_ic;
    x2_7 += x2_ic;
    p1_7 += p1_ic;
    p2_7 += p2_ic;
    integr_7 += integr_ic;
    *x1 = x1_7;
    *x2 = x2_7;
    *p1 = p1_7;
    *p2 = p2_7;
    *integr = integr_7;
    return norm(x1_7 - x1_8, x2_7 - x2_8, p1_7 - p1_8, p2_7 - p2_8, integr_7 - integr_8);  
}
void solve_dp(double l, double r, double tol, double x1_ic, double x2_ic, double p1_ic, double p2_ic, double integr_ic, double* x1, double* x2, double* p1, double* p2, double* integr, unsigned int attr, const char* name)
{
    double h, err;
    double glob_err = 0.f;
    FILE * output;
    double x1_tmp, x2_tmp, p1_tmp, p2_tmp, integr_tmp;
    h = 0.1;
    if(attr & WRITE_IN_FILE)
    {
        output = fopen(name, "w");
        if(output == NULL)
        {
            printf("Cannot open file to write \n");
            exit(EXIT_FAILURE);
        }
        fprintf(output, "%e %e %e %e %e %e \n", l, x1_ic, x1_ic - check_sol(l), p1_ic, p2_ic, integr_ic);
    }
    while(l < r)
    {
        if(l + h > r)
            h = r - l;
        err = dorman_prince(h, l, x1_ic, x2_ic, p1_ic, p2_ic, integr_ic, &x1_tmp, &x2_tmp, &p1_tmp, &p2_tmp, &integr_tmp);
        if(err < tol)
        {
            l += h;
            if(attr&WRITE_IN_FILE)
            {
                fprintf(output, "%e %e %e %e %e %e \n", l, x1_tmp, x1_tmp - check_sol(l), p1_tmp, p2_tmp, integr_tmp);
            }
            h = get_h(h, err, tol);
            x1_ic = x1_tmp; x2_ic = x2_tmp; p1_ic = p1_tmp; p2_ic = p2_tmp; integr_ic = integr_tmp;
        }
        else 
            h = get_h(h, err, tol);
    }
    if(attr&WRITE_IN_FILE)
        fclose(output);
    *x1 = x1_ic; *x2 = x2_ic; *p1 = p1_ic; *p2 = p2_ic; *integr = integr_ic;
}
int shooting_method(double l, double r, double tol, double shooting_method_accuracy, double x1_l, double x2_l, double p1_l, double p2_l, double x1_r, double x2_r, double* x1_ic_answer_l, double *x2_ic_answer_l, double* p1_ic_answer_l, double* p2_ic_answer_l, double* integral_ic_answer_l)
{
    double x1_r_, x2_r_;
    double trash;
    double jack[4];
    double eps = 1e-10;
    do
    {
        solve_dp(l, r, tol, x1_l, x2_l, p1_l, p2_l, 0, &x1_r_, &x2_r_, &trash, &trash, &trash, 0, NULL);
        double left_value_x1, right_value_x1;
        double left_value_x2, right_value_x2;
        solve_dp(l, r, tol, x1_l, x2_l - eps, p1_l, p2_l, 0, &left_value_x1, &left_value_x2, &trash, &trash, &trash, 0, NULL);
        solve_dp(l, r, tol, x1_l, x2_l + eps, p1_l, p2_l, 0, &right_value_x1, &right_value_x2, &trash, &trash, &trash, 0, NULL);
        jack[0] = (right_value_x1 - left_value_x1)/(2.*eps);
        jack[1] = (right_value_x2 - left_value_x2)/(2.*eps);
        solve_dp(l, r, tol, x1_l, x2_l, p1_l - eps, p2_l, 0, &left_value_x1, &left_value_x2, &trash, &trash, &trash, 0, NULL);
        solve_dp(l, r, tol, x1_l, x2_l, p1_l + eps, p2_l, 0, &right_value_x1, &right_value_x2, &trash, &trash, &trash, 0, NULL);
        jack[2] = (right_value_x1 - left_value_x1)/(2.*eps);
        jack[3] = (right_value_x2 - left_value_x2)/(2.*eps);
        x2_l = x2_l - (jack[0]*(x1_r_ - x1_r) + jack[1]*(x2_r_ - x2_r));
        p1_l = p1_l - (jack[2]*(x1_r_ - x1_r) + jack[3]*(x2_r_ - x2_r));
    } while (norm2(x1_r - x1_r_, x2_r - x2_r_) > shooting_method_accuracy);
    *x1_ic_answer_l = x1_l;
    *x2_ic_answer_l = x2_l;
    *p1_ic_answer_l = p1_l;
    *p2_ic_answer_l = p2_l;
}
int check_shooting_method(double l, double r, double tol, double shooting_method_accuracy, double x1_l, double x2_l, double x1_r, double x2_r, double* x1_l_answer, double* x2_l_answer)
{
    double x1_r_;
    double x2_r_;
    double trash;
    double jack[4];
    double eps = 1e-10;
    do
    {
        solve_dp(l, r, tol, x1_l, x2_l, 0, 0, 0, &x1_r_, &x2_r_, &trash, &trash, &trash, 0, NULL);
        double x1_left_value, x1_right_value;
        double x2_left_value, x2_right_value;
        solve_dp(l, r, tol, x1_l - eps, x2_l, 0, 0, 0, &x1_left_value, &x2_left_value, &trash, &trash, &trash, 0, NULL);
        solve_dp(l, r, tol, x1_l + eps, x2_l, 0, 0, 0, &x1_right_value, &x2_right_value, &trash, &trash, &trash, 0, NULL);
        jack[0] = (x1_right_value - x1_left_value)/(2.*eps);
        jack[1] = (x2_right_value - x2_left_value)/(2.*eps);
        solve_dp(l, r, tol, x1_l, x2_l - eps, 0, 0, 0, &x1_left_value, &x2_left_value, &trash, &trash, &trash, 0, NULL);
        solve_dp(l, r, tol, x1_l, x2_l + eps, 0, 0, 0, &x1_right_value, &x2_right_value, &trash, &trash, &trash, 0, NULL);
        jack[2] = (x1_right_value - x1_left_value)/(2.*eps);
        jack[3] = (x2_right_value - x2_left_value)/(2.*eps);
        //double det = jack[0]*jack[3] - jack[1]*jack[2];
        double inverse[4];
        inverse_matrix(jack, inverse);
        //inverse[0] = jack[3]/det;
        //inverse[1] = -jack[2]/det;
        //inverse[2] = -jack[1]/det;
        //inverse[3] = jack[0]/det;
        printf("%lf %lf \n", x1_l, x2_l);
        x1_l = x1_l - (inverse[0])*(x1_r_ - x1_r) - (inverse[1])*(x2_r_ - x2_r);
        x2_l = x2_l - (inverse[2])*(x1_r_ - x1_r) - inverse[3]*(x2_r_ - x2_r);
    } while (fabs(x1_r - x1_r_) > shooting_method_accuracy);
    *x1_l_answer = x1_l;    
    *x2_l_answer = x2_l;
}