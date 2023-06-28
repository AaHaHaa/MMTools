#include "gaussianElimination.cu"

/* ------------------------------------------------------------------------
* ---------------------------- mySpline_helper ----------------------------
* -------------------------------------------------------------------------
* This function is the helper function of mySpline().
* It uses Gaussian elimination to solve the coefficients of the polynomial 
* for each spline interval.
* Up to 6th order is implemented.
* -----------------------------------------------------------------------*/
__device__ double2 call_order1(const double target_x, const double upper_x,
                               const double2 lower_y, const double2 upper_y,
                               const double2 lower_dy, const double2 upper_dy) {
    const unsigned int num_order = 1;
    const unsigned int poly_order = 2*(num_order+1); // this is actually the order of polynomial plus one

    // real part
    double cx[poly_order];
    double AMx[poly_order*(poly_order+1)] = \
                                            {3*pow(upper_x,2),      2*upper_x,       1, 0, upper_dy.x,\
                                               pow(upper_x,3), pow(upper_x,2), upper_x, 1, upper_y.x,\
                                                            0,              0,       1, 0, lower_dy.x,\
                                                            0,              0,       0, 1, lower_y.x};
    gaussianElimination(AMx, poly_order, cx);

    // imaginary part
    double cy[poly_order];
    double AMy[poly_order*(poly_order+1)] = \
                                            {3*pow(upper_x,2),      2*upper_x,       1, 0, upper_dy.y,\
                                               pow(upper_x,3), pow(upper_x,2), upper_x, 1, upper_y.y,\
                                                            0,              0,       1, 0, lower_dy.y,\
                                                            0,              0,       0, 1, lower_y.y};
    gaussianElimination(AMy, poly_order, cy);

    double2 target_y;
    target_y.x = 0; target_y.y = 0;
    for (unsigned int i=0; i<poly_order; i++) {
        target_y.x += cx[i]*pow(target_x,double(poly_order-1-i));
        target_y.y += cy[i]*pow(target_x,double(poly_order-1-i));
    }

    return target_y;
}

__device__ double2 call_order2(const double target_x, const double upper_x,
                               const double2 lower_y, const double2 upper_y,
                               const double2 lower_dy, const double2 upper_dy,
                               const double2 lower_d2y, const double2 upper_d2y) {
    const unsigned int num_order = 2;
    const unsigned int poly_order = 2*(num_order+1); // this is actually the order of polynomial plus one
            
    // real part
    double cx[poly_order];
    double AMx[poly_order*(poly_order+1)] = \
                                            {20*pow(upper_x,3), 12*pow(upper_x,2),        6*upper_x,              2,       0, 0, upper_d2y.x,\
                                              5*pow(upper_x,4),  4*pow(upper_x,3), 3*pow(upper_x,2),      2*upper_x,       1, 0, upper_dy.x,\
                                                pow(upper_x,5),    pow(upper_x,4),   pow(upper_x,3), pow(upper_x,2), upper_x, 1, upper_y.x,\
                                                             0,                 0,                0,              2,       0, 0, lower_d2y.x,\
                                                             0,                 0,                0,              0,       1, 0, lower_dy.x,\
                                                             0,                 0,                0,              0,       0, 1, lower_y.x};
    gaussianElimination(AMx, poly_order, cx);

    // imaginary part
    double cy[poly_order];
    double AMy[poly_order*(poly_order+1)] = \
                                            {20*pow(upper_x,3), 12*pow(upper_x,2),        6*upper_x,              2,       0, 0, upper_d2y.y,\
                                              5*pow(upper_x,4),  4*pow(upper_x,3), 3*pow(upper_x,2),      2*upper_x,       1, 0, upper_dy.y,\
                                                pow(upper_x,5),    pow(upper_x,4),   pow(upper_x,3), pow(upper_x,2), upper_x, 1, upper_y.y,\
                                                             0,                 0,                0,              2,       0, 0, lower_d2y.y,\
                                                             0,                 0,                0,              0,       1, 0, lower_dy.y,\
                                                             0,                 0,                0,              0,       0, 1, lower_y.y};
    gaussianElimination(AMy, poly_order, cy);

    double2 target_y;
    target_y.x = 0; target_y.y = 0;
    for (unsigned int i=0; i<poly_order; i++) {
        target_y.x += cx[i]*pow(target_x,double(poly_order-1-i));
        target_y.y += cy[i]*pow(target_x,double(poly_order-1-i));
    }

    return target_y;
}

__device__ double2 call_order3(const double target_x, const double upper_x,
                               const double2 lower_y, const double2 upper_y,
                               const double2 lower_dy, const double2 upper_dy,
                               const double2 lower_d2y, const double2 upper_d2y,
                               const double2 lower_d3y, const double2 upper_d3y) {
    const unsigned int num_order = 3;
    const unsigned int poly_order = 2*(num_order+1); // this is actually the order of polynomial plus one
            
    // real part
    double cx[poly_order];
    double AMx[poly_order*(poly_order+1)] = \
                                            {210*pow(upper_x,4), 120*pow(upper_x,3), 60*pow(upper_x,2),        24*upper_x,                6,              0,       0, 0, upper_d3y.x,\
                                              42*pow(upper_x,5),  30*pow(upper_x,4), 20*pow(upper_x,3), 12*pow(upper_x,2),        6*upper_x,              2,       0, 0, upper_d2y.x,\
                                               7*pow(upper_x,6),   6*pow(upper_x,5),  5*pow(upper_x,4),  4*pow(upper_x,3), 3*pow(upper_x,2),      2*upper_x,       1, 0, upper_dy.x,\
                                                 pow(upper_x,7),     pow(upper_x,6),    pow(upper_x,5),    pow(upper_x,4),   pow(upper_x,3), pow(upper_x,2), upper_x, 1, upper_y.x,\
                                                              0,                  0,                 0,                 0,                6,              0,       0, 0, lower_d3y.x,\
                                                              0,                  0,                 0,                 0,                0,              2,       0, 0, lower_d2y.x,\
                                                              0,                  0,                 0,                 0,                0,              0,       1, 0, lower_dy.x,\
                                                              0,                  0,                 0,                 0,                0,              0,       0, 1, lower_y.x};
    gaussianElimination(AMx, poly_order, cx);

    // imaginary part
    double cy[poly_order];
    double AMy[poly_order*(poly_order+1)] = \
                                            {210*pow(upper_x,4), 120*pow(upper_x,3), 60*pow(upper_x,2),        24*upper_x,                6,              0,       0, 0, upper_d3y.y,\
                                              42*pow(upper_x,5),  30*pow(upper_x,4), 20*pow(upper_x,3), 12*pow(upper_x,2),        6*upper_x,              2,       0, 0, upper_d2y.y,\
                                               7*pow(upper_x,6),   6*pow(upper_x,5),  5*pow(upper_x,4),  4*pow(upper_x,3), 3*pow(upper_x,2),      2*upper_x,       1, 0, upper_dy.y,\
                                                 pow(upper_x,7),     pow(upper_x,6),    pow(upper_x,5),    pow(upper_x,4),   pow(upper_x,3), pow(upper_x,2), upper_x, 1, upper_y.y,\
                                                              0,                  0,                 0,                 0,                6,              0,       0, 0, lower_d3y.y,\
                                                              0,                  0,                 0,                 0,                0,              2,       0, 0, lower_d2y.y,\
                                                              0,                  0,                 0,                 0,                0,              0,       1, 0, lower_dy.y,\
                                                              0,                  0,                 0,                 0,                0,              0,       0, 1, lower_y.y};
    gaussianElimination(AMy, poly_order, cy);

    double2 target_y;
    target_y.x = 0; target_y.y = 0;
    for (unsigned int i=0; i<poly_order; i++) {
        target_y.x += cx[i]*pow(target_x,double(poly_order-1-i));
        target_y.y += cy[i]*pow(target_x,double(poly_order-1-i));
    }

    return target_y;
}

__device__ double2 call_order4(const double target_x, const double upper_x,
                               const double2 lower_y, const double2 upper_y,
                               const double2 lower_dy, const double2 upper_dy,
                               const double2 lower_d2y, const double2 upper_d2y,
                               const double2 lower_d3y, const double2 upper_d3y,
                               const double2 lower_d4y, const double2 upper_d4y) {
    const unsigned int num_order = 4;
    const unsigned int poly_order = 2*(num_order+1); // this is actually the order of polynomial plus one
            
    // real part
    double cx[poly_order];
    double AMx[poly_order*(poly_order+1)] = \
                                            {3024*pow(upper_x,5), 1680*pow(upper_x,4), 840*pow(upper_x,3), 360*pow(upper_x,2),       120*upper_x,                24,                0,              0,       0, 0, upper_d4y.x,\
                                              504*pow(upper_x,6),  336*pow(upper_x,5), 210*pow(upper_x,4), 120*pow(upper_x,3), 60*pow(upper_x,2),        24*upper_x,                6,              0,       0, 0, upper_d3y.x,\
                                               72*pow(upper_x,7),   56*pow(upper_x,6),  42*pow(upper_x,5),  30*pow(upper_x,4), 20*pow(upper_x,3), 12*pow(upper_x,2),        6*upper_x,              2,       0, 0, upper_d2y.x,\
                                                9*pow(upper_x,8),    8*pow(upper_x,7),   7*pow(upper_x,6),   6*pow(upper_x,5),  5*pow(upper_x,4),  4*pow(upper_x,3), 3*pow(upper_x,2),      2*upper_x,       1, 0, upper_dy.x,\
                                                  pow(upper_x,9),      pow(upper_x,8),     pow(upper_x,7),     pow(upper_x,6),    pow(upper_x,5),    pow(upper_x,4),   pow(upper_x,3), pow(upper_x,2), upper_x, 1, upper_y.x,\
                                                               0,                   0,                  0,                  0,                 0,                24,                0,              0,       0, 0, lower_d4y.x,\
                                                               0,                   0,                  0,                  0,                 0,                 0,                6,              0,       0, 0, lower_d3y.x,\
                                                               0,                   0,                  0,                  0,                 0,                 0,                0,              2,       0, 0, lower_d2y.x,\
                                                               0,                   0,                  0,                  0,                 0,                 0,                0,              0,       1, 0, lower_dy.x,\
                                                               0,                   0,                  0,                  0,                 0,                 0,                0,              0,       0, 1, lower_y.x};
    gaussianElimination(AMx, poly_order, cx);

    // imaginary part
    double cy[poly_order];
    double AMy[poly_order*(poly_order+1)] = \
                                            {3024*pow(upper_x,5), 1680*pow(upper_x,4), 840*pow(upper_x,3), 360*pow(upper_x,2),       120*upper_x,                24,                0,              0,       0, 0, upper_d4y.y,\
                                              504*pow(upper_x,6),  336*pow(upper_x,5), 210*pow(upper_x,4), 120*pow(upper_x,3), 60*pow(upper_x,2),        24*upper_x,                6,              0,       0, 0, upper_d3y.y,\
                                               72*pow(upper_x,7),   56*pow(upper_x,6),  42*pow(upper_x,5),  30*pow(upper_x,4), 20*pow(upper_x,3), 12*pow(upper_x,2),        6*upper_x,              2,       0, 0, upper_d2y.y,\
                                                9*pow(upper_x,8),    8*pow(upper_x,7),   7*pow(upper_x,6),   6*pow(upper_x,5),  5*pow(upper_x,4),  4*pow(upper_x,3), 3*pow(upper_x,2),      2*upper_x,       1, 0, upper_dy.y,\
                                                  pow(upper_x,9),      pow(upper_x,8),     pow(upper_x,7),     pow(upper_x,6),    pow(upper_x,5),    pow(upper_x,4),   pow(upper_x,3), pow(upper_x,2), upper_x, 1, upper_y.y,\
                                                               0,                   0,                  0,                  0,                 0,                24,                0,              0,       0, 0, lower_d4y.y,\
                                                               0,                   0,                  0,                  0,                 0,                 0,                6,              0,       0, 0, lower_d3y.y,\
                                                               0,                   0,                  0,                  0,                 0,                 0,                0,              2,       0, 0, lower_d2y.y,\
                                                               0,                   0,                  0,                  0,                 0,                 0,                0,              0,       1, 0, lower_dy.y,\
                                                               0,                   0,                  0,                  0,                 0,                 0,                0,              0,       0, 1, lower_y.y};
    gaussianElimination(AMy, poly_order, cy);

    double2 target_y;
    target_y.x = 0; target_y.y = 0;
    for (unsigned int i=0; i<poly_order; i++) {
        target_y.x += cx[i]*pow(target_x,double(poly_order-1-i));
        target_y.y += cy[i]*pow(target_x,double(poly_order-1-i));
    }

    return target_y;
}

__device__ double2 call_order5(const double target_x, const double upper_x,
                               const double2 lower_y, const double2 upper_y,
                               const double2 lower_dy, const double2 upper_dy,
                               const double2 lower_d2y, const double2 upper_d2y,
                               const double2 lower_d3y, const double2 upper_d3y,
                               const double2 lower_d4y, const double2 upper_d4y,
                               const double2 lower_d5y, const double2 upper_d5y) {
    const unsigned int num_order = 5;
    const unsigned int poly_order = 2*(num_order+1); // this is actually the order of polynomial plus one
            
    // real part
    double cx[poly_order];
    double AMx[poly_order*(poly_order+1)] = \
                                            {55440*pow(upper_x,6), 30240*pow(upper_x,5), 15120*pow(upper_x,4), 6720*pow(upper_x,3), 2520*pow(upper_x,2),        720*upper_x,               120,                 0,                0,              0,       0, 0, upper_d5y.x,\
                                              7920*pow(upper_x,7),  5040*pow(upper_x,6),  3024*pow(upper_x,5), 1680*pow(upper_x,4),  840*pow(upper_x,3), 360*pow(upper_x,2),       120*upper_x,                24,                0,              0,       0, 0, upper_d4y.x,\
                                               990*pow(upper_x,8),   720*pow(upper_x,7),   504*pow(upper_x,6),  336*pow(upper_x,5),  210*pow(upper_x,4), 120*pow(upper_x,3), 60*pow(upper_x,2),        24*upper_x,                6,              0,       0, 0, upper_d3y.x,\
                                               110*pow(upper_x,9),    90*pow(upper_x,8),    72*pow(upper_x,7),   56*pow(upper_x,6),   42*pow(upper_x,5),  30*pow(upper_x,4), 20*pow(upper_x,3), 12*pow(upper_x,2),        6*upper_x,              2,       0, 0, upper_d2y.x,\
                                               11*pow(upper_x,10),    10*pow(upper_x,9),     9*pow(upper_x,8),    8*pow(upper_x,7),    7*pow(upper_x,6),   6*pow(upper_x,5),  5*pow(upper_x,4),  4*pow(upper_x,3), 3*pow(upper_x,2),      2*upper_x,       1, 0, upper_dy.x,\
                                                  pow(upper_x,11),      pow(upper_x,10),       pow(upper_x,9),      pow(upper_x,8),      pow(upper_x,7),     pow(upper_x,6),    pow(upper_x,5),    pow(upper_x,4),   pow(upper_x,3), pow(upper_x,2), upper_x, 1, upper_y.x,\
                                                                0,                    0,                    0,                   0,                   0,                  0,               120,                 0,                0,              0,       0, 0, lower_d5y.x,\
                                                                0,                    0,                    0,                   0,                   0,                  0,                 0,                24,                0,              0,       0, 0, lower_d4y.x,\
                                                                0,                    0,                    0,                   0,                   0,                  0,                 0,                 0,                6,              0,       0, 0, lower_d3y.x,\
                                                                0,                    0,                    0,                   0,                   0,                  0,                 0,                 0,                0,              2,       0, 0, lower_d2y.x,\
                                                                0,                    0,                    0,                   0,                   0,                  0,                 0,                 0,                0,              0,       1, 0, lower_dy.x,\
                                                                0,                    0,                    0,                   0,                   0,                  0,                 0,                 0,                0,              0,       0, 1, lower_y.x};
    gaussianElimination(AMx, poly_order, cx);

    // imaginary part
    double cy[poly_order];
    double AMy[poly_order*(poly_order+1)] = \
                                            {55440*pow(upper_x,6), 30240*pow(upper_x,5), 15120*pow(upper_x,4), 6720*pow(upper_x,3), 2520*pow(upper_x,2),        720*upper_x,               120,                 0,                0,              0,       0, 0, upper_d5y.y,\
                                              7920*pow(upper_x,7),  5040*pow(upper_x,6),  3024*pow(upper_x,5), 1680*pow(upper_x,4),  840*pow(upper_x,3), 360*pow(upper_x,2),       120*upper_x,                24,                0,              0,       0, 0, upper_d4y.y,\
                                               990*pow(upper_x,8),   720*pow(upper_x,7),   504*pow(upper_x,6),  336*pow(upper_x,5),  210*pow(upper_x,4), 120*pow(upper_x,3), 60*pow(upper_x,2),        24*upper_x,                6,              0,       0, 0, upper_d3y.y,\
                                               110*pow(upper_x,9),    90*pow(upper_x,8),    72*pow(upper_x,7),   56*pow(upper_x,6),   42*pow(upper_x,5),  30*pow(upper_x,4), 20*pow(upper_x,3), 12*pow(upper_x,2),        6*upper_x,              2,       0, 0, upper_d2y.y,\
                                               11*pow(upper_x,10),    10*pow(upper_x,9),     9*pow(upper_x,8),    8*pow(upper_x,7),    7*pow(upper_x,6),   6*pow(upper_x,5),  5*pow(upper_x,4),  4*pow(upper_x,3), 3*pow(upper_x,2),      2*upper_x,       1, 0, upper_dy.y,\
                                                  pow(upper_x,11),      pow(upper_x,10),       pow(upper_x,9),      pow(upper_x,8),      pow(upper_x,7),     pow(upper_x,6),    pow(upper_x,5),    pow(upper_x,4),   pow(upper_x,3), pow(upper_x,2), upper_x, 1, upper_y.y,\
                                                                0,                    0,                    0,                   0,                   0,                  0,               120,                 0,                0,              0,       0, 0, lower_d5y.y,\
                                                                0,                    0,                    0,                   0,                   0,                  0,                 0,                24,                0,              0,       0, 0, lower_d4y.y,\
                                                                0,                    0,                    0,                   0,                   0,                  0,                 0,                 0,                6,              0,       0, 0, lower_d3y.y,\
                                                                0,                    0,                    0,                   0,                   0,                  0,                 0,                 0,                0,              2,       0, 0, lower_d2y.y,\
                                                                0,                    0,                    0,                   0,                   0,                  0,                 0,                 0,                0,              0,       1, 0, lower_dy.y,\
                                                                0,                    0,                    0,                   0,                   0,                  0,                 0,                 0,                0,              0,       0, 1, lower_y.y};
    gaussianElimination(AMy, poly_order, cy);

    double2 target_y;
    target_y.x = 0; target_y.y = 0;
    for (unsigned int i=0; i<poly_order; i++) {
        target_y.x += cx[i]*pow(target_x,double(poly_order-1-i));
        target_y.y += cy[i]*pow(target_x,double(poly_order-1-i));
    }

    return target_y;
}

__device__ double2 call_order6(const double target_x, const double upper_x,
                               const double2 lower_y, const double2 upper_y,
                               const double2 lower_dy, const double2 upper_dy,
                               const double2 lower_d2y, const double2 upper_d2y,
                               const double2 lower_d3y, const double2 upper_d3y,
                               const double2 lower_d4y, const double2 upper_d4y,
                               const double2 lower_d5y, const double2 upper_d5y,
                               const double2 lower_d6y, const double2 upper_d6y) {
    const unsigned int num_order = 6;
    const unsigned int poly_order = 2*(num_order+1); // this is actually the order of polynomial plus one
            
    // real part
    double cx[poly_order];
    double AMx[poly_order*(poly_order+1)] = \
                                            {1235520*pow(upper_x,7), 665280*pow(upper_x,6), 332640*pow(upper_x,5), 151200*pow(upper_x,4), 60480*pow(upper_x,3), 20160*pow(upper_x,2),        5040*upper_x,                720,                 0,                 0,                0,              0,       0, 0, upper_d6y.x,\
                                              154440*pow(upper_x,8),  95040*pow(upper_x,7),  55440*pow(upper_x,6),  30240*pow(upper_x,5), 15120*pow(upper_x,4),  6720*pow(upper_x,3), 2520*pow(upper_x,2),        720*upper_x,               120,                 0,                0,              0,       0, 0, upper_d5y.x,\
                                               17160*pow(upper_x,9),  11880*pow(upper_x,8),   7920*pow(upper_x,7),   5040*pow(upper_x,6),  3024*pow(upper_x,5),  1680*pow(upper_x,4),  840*pow(upper_x,3), 360*pow(upper_x,2),       120*upper_x,                24,                0,              0,       0, 0, upper_d4y.x,\
                                               1716*pow(upper_x,10),   1320*pow(upper_x,9),    990*pow(upper_x,8),    720*pow(upper_x,7),   504*pow(upper_x,6),   336*pow(upper_x,5),  210*pow(upper_x,4), 120*pow(upper_x,3), 60*pow(upper_x,2),        24*upper_x,                6,              0,       0, 0, upper_d3y.x,\
                                                156*pow(upper_x,11),   132*pow(upper_x,10),    110*pow(upper_x,9),     90*pow(upper_x,8),    72*pow(upper_x,7),    56*pow(upper_x,6),   42*pow(upper_x,5),  30*pow(upper_x,4), 20*pow(upper_x,3), 12*pow(upper_x,2),        6*upper_x,              2,       0, 0, upper_d2y.x,\
                                                 13*pow(upper_x,12),    12*pow(upper_x,11),    11*pow(upper_x,10),     10*pow(upper_x,9),     9*pow(upper_x,8),     8*pow(upper_x,7),    7*pow(upper_x,6),   6*pow(upper_x,5),  5*pow(upper_x,4),  4*pow(upper_x,3), 3*pow(upper_x,2),      2*upper_x,       1, 0, upper_dy.x,\
                                                    pow(upper_x,13),       pow(upper_x,12),       pow(upper_x,11),       pow(upper_x,10),       pow(upper_x,9),       pow(upper_x,8),      pow(upper_x,7),     pow(upper_x,6),    pow(upper_x,5),    pow(upper_x,4),   pow(upper_x,3), pow(upper_x,2), upper_x, 1, upper_y.x,\
                                                                  0,                     0,                     0,                     0,                    0,                    0,                   0,                720,                 0,                 0,                0,              0,       0, 0, lower_d6y.x,\
                                                                  0,                     0,                     0,                     0,                    0,                    0,                   0,                  0,               120,                 0,                0,              0,       0, 0, lower_d5y.x,\
                                                                  0,                     0,                     0,                     0,                    0,                    0,                   0,                  0,                 0,                24,                0,              0,       0, 0, lower_d4y.x,\
                                                                  0,                     0,                     0,                     0,                    0,                    0,                   0,                  0,                 0,                 0,                6,              0,       0, 0, lower_d3y.x,\
                                                                  0,                     0,                     0,                     0,                    0,                    0,                   0,                  0,                 0,                 0,                0,              2,       0, 0, lower_d2y.x,\
                                                                  0,                     0,                     0,                     0,                    0,                    0,                   0,                  0,                 0,                 0,                0,              0,       1, 0, lower_dy.x,\
                                                                  0,                     0,                     0,                     0,                    0,                    0,                   0,                  0,                 0,                 0,                0,              0,       0, 1, lower_y.x};
    gaussianElimination(AMx, poly_order, cx);

    // imaginary part
    double cy[poly_order];
    double AMy[poly_order*(poly_order+1)] = \
                                            {1235520*pow(upper_x,7), 665280*pow(upper_x,6), 332640*pow(upper_x,5), 151200*pow(upper_x,4), 60480*pow(upper_x,3), 20160*pow(upper_x,2),        5040*upper_x,                720,                 0,                 0,                0,              0,       0, 0, upper_d6y.y,\
                                              154440*pow(upper_x,8),  95040*pow(upper_x,7),  55440*pow(upper_x,6),  30240*pow(upper_x,5), 15120*pow(upper_x,4),  6720*pow(upper_x,3), 2520*pow(upper_x,2),        720*upper_x,               120,                 0,                0,              0,       0, 0, upper_d5y.y,\
                                               17160*pow(upper_x,9),  11880*pow(upper_x,8),   7920*pow(upper_x,7),   5040*pow(upper_x,6),  3024*pow(upper_x,5),  1680*pow(upper_x,4),  840*pow(upper_x,3), 360*pow(upper_x,2),       120*upper_x,                24,                0,              0,       0, 0, upper_d4y.y,\
                                               1716*pow(upper_x,10),   1320*pow(upper_x,9),    990*pow(upper_x,8),    720*pow(upper_x,7),   504*pow(upper_x,6),   336*pow(upper_x,5),  210*pow(upper_x,4), 120*pow(upper_x,3), 60*pow(upper_x,2),        24*upper_x,                6,              0,       0, 0, upper_d3y.y,\
                                                156*pow(upper_x,11),   132*pow(upper_x,10),    110*pow(upper_x,9),     90*pow(upper_x,8),    72*pow(upper_x,7),    56*pow(upper_x,6),   42*pow(upper_x,5),  30*pow(upper_x,4), 20*pow(upper_x,3), 12*pow(upper_x,2),        6*upper_x,              2,       0, 0, upper_d2y.y,\
                                                 13*pow(upper_x,12),    12*pow(upper_x,11),    11*pow(upper_x,10),     10*pow(upper_x,9),     9*pow(upper_x,8),     8*pow(upper_x,7),    7*pow(upper_x,6),   6*pow(upper_x,5),  5*pow(upper_x,4),  4*pow(upper_x,3), 3*pow(upper_x,2),      2*upper_x,       1, 0, upper_dy.y,\
                                                    pow(upper_x,13),       pow(upper_x,12),       pow(upper_x,11),       pow(upper_x,10),       pow(upper_x,9),       pow(upper_x,8),      pow(upper_x,7),     pow(upper_x,6),    pow(upper_x,5),    pow(upper_x,4),   pow(upper_x,3), pow(upper_x,2), upper_x, 1, upper_y.y,\
                                                                  0,                     0,                     0,                     0,                    0,                    0,                   0,                720,                 0,                 0,                0,              0,       0, 0, lower_d6y.y,\
                                                                  0,                     0,                     0,                     0,                    0,                    0,                   0,                  0,               120,                 0,                0,              0,       0, 0, lower_d5y.y,\
                                                                  0,                     0,                     0,                     0,                    0,                    0,                   0,                  0,                 0,                24,                0,              0,       0, 0, lower_d4y.y,\
                                                                  0,                     0,                     0,                     0,                    0,                    0,                   0,                  0,                 0,                 0,                6,              0,       0, 0, lower_d3y.y,\
                                                                  0,                     0,                     0,                     0,                    0,                    0,                   0,                  0,                 0,                 0,                0,              2,       0, 0, lower_d2y.y,\
                                                                  0,                     0,                     0,                     0,                    0,                    0,                   0,                  0,                 0,                 0,                0,              0,       1, 0, lower_dy.y,\
                                                                  0,                     0,                     0,                     0,                    0,                    0,                   0,                  0,                 0,                 0,                0,              0,       0, 1, lower_y.y};
    gaussianElimination(AMy, poly_order, cy);

    double2 target_y;
    target_y.x = 0; target_y.y = 0;
    for (unsigned int i=0; i<poly_order; i++) {
        target_y.x += cx[i]*pow(target_x,double(poly_order-1-i));
        target_y.y += cy[i]*pow(target_x,double(poly_order-1-i));
    }

    return target_y;
}