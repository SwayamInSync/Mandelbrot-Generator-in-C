#define PY_SSIZE_T_CLEAN

#include<Python.h>

#include <stdio.h>
#include <sleef.h>
#include <sleefquad.h>
#include <omp.h>
#include <cmath>
#include <type_traits>
#include <string>
#include <limits>
#include <typeinfo>
#include <quadmath.h>

typedef struct {
    unsigned char r, g, b;
} Color;

template <typename T>
struct Complex
{
    T real;
    T img;
};

template<typename T>
struct MandelbrotResult
{
    T smooth_iter;
    Complex<T> final_z;
};

// template <typename T>
// struct MathOps{};

template <typename T>
struct MathOps
{
    static T add(T a, T b) { return a + b; }
    static T sub(T a, T b) { return a - b; }
    static T mul(T a, T b) { return a * b; }
    static T div(T a, T b) { return a / b; }
    static T sqrt(T a) { return std::sqrt(a); }
    static T log(T a) { return std::log(a); }
    static T abs(T a) { return std::abs(a); }
    static bool greater(T a, T b) { return a > b; }
    static T cast_from_double(double a) { return static_cast<T>(a); }
    static double cast_to_double(T a) { return static_cast<double>(a); }
    static T from_int(int a) { return static_cast<T>(a); }
    static constexpr T pi() { return static_cast<T>(M_PI); }
    static constexpr T ln2() { return static_cast<T>(M_LN2); }

    static const char* dtype_name() {
        return typeid(T).name();
    }
};

template<>
struct MathOps<Sleef_quad>
{
    static Sleef_quad add(Sleef_quad a, Sleef_quad b) { return Sleef_addq1_u05(a, b); }
    static Sleef_quad sub(Sleef_quad a, Sleef_quad b) { return Sleef_subq1_u05(a, b); }
    static Sleef_quad mul(Sleef_quad a, Sleef_quad b) { return Sleef_mulq1_u05(a, b); }
    static Sleef_quad div(Sleef_quad a, Sleef_quad b) { return Sleef_divq1_u05(a, b); }
    static Sleef_quad sqrt(Sleef_quad a) { return Sleef_sqrtq1_u05(a); }
    static Sleef_quad log(Sleef_quad a) { return Sleef_logq1_u10(a); }
    static Sleef_quad abs(Sleef_quad a) { return Sleef_fabsq1(a); }
    static bool greater(Sleef_quad a, Sleef_quad b) { return Sleef_icmpgtq1(a, b); }
    static Sleef_quad cast_from_double(double a) { return Sleef_cast_from_doubleq1(a); }
    static double cast_to_double(Sleef_quad a) { return Sleef_cast_to_doubleq1(a); }
    static Sleef_quad from_int(int a) { return Sleef_cast_from_int64q1(a); }
    static Sleef_quad pi() { return SLEEF_M_PIq; }
    static Sleef_quad ln2() { return SLEEF_M_LN2q; }

    static const char* dtype_name() {
        return "Sleef_quad";
    }
};


 template<>
 struct MathOps<__float128>
 {
     static __float128 add(__float128 a, __float128 b) { return a + b; }
     static __float128 sub(__float128 a, __float128 b) { return a - b; }
     static __float128 mul(__float128 a, __float128 b) { return a * b; }
     static __float128 div(__float128 a, __float128 b) { return a / b; }
     static __float128 sqrt(__float128 a) { return sqrtq(a); }
     static __float128 log(__float128 a) { return logq(a); }
     static __float128 abs(__float128 a) { return fabsq(a); }
     static bool greater(__float128 a, __float128 b) { return a > b; }
     static __float128 cast_from_double(double a) { return (__float128)a; }
     static double cast_to_double(__float128 a) { return (double)a; }
     static __float128 from_int(int a) { return (__float128)a; }
     static __float128 pi() { return M_PIq; }
     static __float128 ln2() { return M_LN2q; }

     static const char* dtype_name() {
         return "__float128";
     }
 };


template<typename T>
T get_abs_value(const Complex<T>& z)
{
    return MathOps<T>::add(MathOps<T>::mul(z.real, z.real), MathOps<T>::mul(z.img, z.img));
}

template<typename T>
Color get_color(T t, T interior_t) {
    Color color;
    const double epsilon = 1e-10;

    // Interior coloring
    if (std::abs(MathOps<T>::cast_to_double(t) - 1.0) < epsilon) {
        unsigned char value = (unsigned char)(255 * MathOps<T>::cast_to_double(interior_t));
        color.r = value;
        color.g = value;
        color.b = value;
        return color;
    }

    // Exterior coloring
    double t_double = MathOps<T>::cast_to_double(t);
    t_double = std::pow(t_double, 0.5);
    t_double = std::fmod(t_double * 20, 1.0);

    if (t_double < 0.16) {
        color.r = 0;
        color.g = (unsigned char)(255 * (t_double / 0.16));
        color.b = (unsigned char)(128 + 127 * (t_double / 0.16));
    }
    else if (t_double < 0.33) {
        color.r = 0;
        color.g = 255;
        color.b = (unsigned char)(255 * (1 - (t_double - 0.16) / 0.17));
    }
    else if (t_double < 0.5) {
        color.r = (unsigned char)(255 * ((t_double - 0.33) / 0.17));
        color.g = 255;
        color.b = 0;
    }
    else if (t_double < 0.66) {
        color.r = 255;
        color.g = (unsigned char)(255 * (1 - (t_double - 0.5) / 0.16));
        color.b = 0;
    }
    else if (t_double < 0.83) {
        color.r = 255;
        color.g = 0;
        color.b = (unsigned char)(255 * ((t_double - 0.66) / 0.17));
    }
    else {
        color.r = (unsigned char)(255 * (1 - (t_double - 0.83) / 0.17));
        color.g = 0;
        color.b = (unsigned char)(128 * ((t_double - 0.83) / 0.17));
    }

    return color;
}

template<typename T>
struct MathOps<Complex<T>>
{
    static Complex<T> add(Complex<T> a, Complex<T> b)
    {
        return {MathOps<T>::add(a.real, b.real), MathOps<T>::add(a.img, b.img)};
    }

    static Complex<T> mul(Complex<T> a, Complex<T> b)
    {
        return {
            MathOps<T>::sub(MathOps<T>::mul(a.real, b.real), MathOps<T>::mul(a.img, b.img)),
            MathOps<T>::add(MathOps<T>::mul(a.real, b.img), MathOps<T>::mul(a.img, b.real))
        };
    }

    static Complex<T> pow(Complex<T> z, int n)
    {
        Complex<T> res = {MathOps<T>::cast_from_double(1.0), MathOps<T>::cast_from_double(0.0)};
        for(int i=0; i < n; i++)
        {
            res = mul(res, z);
        }
        return res;
    }

    static Complex<T> derivative(Complex<T> z, Complex<T> c, int n)
    {
        if(n == 0)
            return {MathOps<T>::cast_from_double(1.0), MathOps<T>::cast_from_double(0.0)};

        Complex<T> prev = derivative(z, c, n-1);
        return mul(prev, {MathOps<T>::cast_from_double(2.0), MathOps<T>::cast_from_double(0.0)});
    }
};

template <typename T>
struct Derivatives
{
    Complex<T> z;
    Complex<T> dz;
    Complex<T> dc;
    Complex<T> dzdz;
};

template<typename T>
Derivatives<T> iterate_and_compute_derivatives(Complex<T> c, int max_iter) {
    Complex<T> z = {MathOps<T>::cast_from_double(0.0), MathOps<T>::cast_from_double(0.0)};
    Complex<T> dz = {MathOps<T>::cast_from_double(1.0), MathOps<T>::cast_from_double(0.0)};
    Complex<T> dc = {MathOps<T>::cast_from_double(0.0), MathOps<T>::cast_from_double(0.0)};
    Complex<T> dzdz = {MathOps<T>::cast_from_double(0.0), MathOps<T>::cast_from_double(0.0)};
    Complex<T> two = {MathOps<T>::cast_from_double(2.0), MathOps<T>::cast_from_double(0.0)};

    for (int i = 0; i < max_iter; i++) {
        dzdz = MathOps<Complex<T>>::add(MathOps<Complex<T>>::mul(two,  MathOps<Complex<T>>::mul(z, dzdz)), MathOps<Complex<T>>::mul(dz, dz));
        dz = MathOps<Complex<T>>::add(MathOps<Complex<T>>::mul(two, MathOps<Complex<T>>::mul(z, dz)), dc);
        z = MathOps<Complex<T>>::add(MathOps<Complex<T>>::mul(z, z), c);
        dc.real = MathOps<T>::cast_from_double(1.0);
        dc.img = MathOps<T>::cast_from_double(0.0);
    }

    return {z, dz, dc, dzdz};
}

template<typename T>
T estimate_interior_distance(Complex<T> c, int max_iter)
{
    Derivatives<T> d = iterate_and_compute_derivatives(c, max_iter);

    T dz_abs_sq = get_abs_value(d.dz);
    T one = MathOps<T>::cast_from_double(1.0);
    T numerator = MathOps<T>::sub(one, dz_abs_sq);

    Complex<T> denominator_term1 =  MathOps<Complex<T>>::mul(d.dc, d.dz);
    Complex<T> denominator_term2 =  MathOps<Complex<T>>::mul(d.dzdz,  MathOps<Complex<T>>::mul(d.z, d.dc));
    Complex<T> denominator =  MathOps<Complex<T>>::add(denominator_term1, denominator_term2);

    T denominator_abs = MathOps<T>::sqrt(get_abs_value(denominator));

    return MathOps<T>::div(numerator, denominator_abs);
}


template <typename T>
T estimate_distance(T cr, T ci, int period)
{
    Complex<T> c = {cr, ci};
    Complex<T> z = c;
    for (int i = 0; i < period; i++)
    {
        z = MathOps<Complex<T>>::add(MathOps<Complex<T>>::pow(z, 2), c);
    }

    Complex<T> dz = MathOps<Complex<T>>::derivative(z, c, period);
    Complex<T> dc = {MathOps<T>::cast_from_double(1.0), MathOps<T>::cast_from_double(0.0)};
    for (int i = 0; i < period; i++)
    {
        dc = MathOps<Complex<T>>::add(MathOps<Complex<T>>::mul(dc, {MathOps<T>::cast_from_double(2.0), MathOps<T>::cast_from_double(0.0)}),
            {MathOps<T>::cast_from_double(1.0), MathOps<T>::cast_from_double(0.0)});
    }

    T abs_dz = MathOps<T>::sqrt(get_abs_value(dz));
    T abs_dc = MathOps<T>::sqrt(get_abs_value(dc));

    return MathOps<T>::div(MathOps<T>::mul(abs_dc, MathOps<T>::log(abs_dz)), abs_dz);
}

// void find_interesting_point(Sleef_quad *center_r, Sleef_quad *center_i, Sleef_quad zoom)
// {
//     Sleef_quad step = Sleef_divq1_u05(Sleef_cast_from_doubleq1(1.0), zoom);
//     Sleef_quad min_distance = Sleef_cast_from_doubleq1(INFINITY);
//     Sleef_quad best_r = *center_r;
//     Sleef_quad best_i = *center_i;

//     for (int i = -10; i <= 10; i++)
//     {
//         for (int j = -10; j <= 10; j++)
//         {
//             Sleef_quad r = Sleef_addq1_u05(*center_r, Sleef_mulq1_u05(Sleef_cast_from_int64q1(i), step));
//             Sleef_quad i = Sleef_addq1_u05(*center_i, Sleef_mulq1_u05(Sleef_cast_from_int64q1(j), step));
//             Sleef_quad distance = estimate_distance(r, i, 100);  // Assuming a max period of 100
//             if (Sleef_icmpltq1(distance, min_distance))
//             {
//                 min_distance = distance;
//                 best_r = r;
//                 best_i = i;
//             }
//         }
//     }

//     *center_r = best_r;
//     *center_i = best_i;
// }

// Sleef_quad interior_distance(Sleef_quad z_real, Sleef_quad z_img) {
//     return Sleef_sqrtq1_u05(get_abs_value(&z_real, &z_img));
// }

template<typename T>
MandelbrotResult<T> mandelbrot(const Complex<T>& c, int MAX_ITER, const T& RADIUS2)
{
    Complex<T> z = {MathOps<T>::cast_from_double(0.0), MathOps<T>::cast_from_double(0.0)};
    T two = MathOps<T>::cast_from_double(2.0);

    MandelbrotResult<T> result;

    for (int i = 0; i < MAX_ITER; i++)
    {
        T z_real2 = MathOps<T>::add(MathOps<T>::sub(MathOps<T>::mul(z.real, z.real), MathOps<T>::mul(z.img, z.img)), c.real);
        T z_img2 = MathOps<T>::add(MathOps<T>::mul(two, MathOps<T>::mul(z.real, z.img)), c.img);

        if (MathOps<T>::greater(get_abs_value((Complex<T>){z_real2, z_img2}), RADIUS2))
        {
            T log_zn = MathOps<T>::log(get_abs_value((Complex<T>){z_real2, z_img2}));
            T nu = MathOps<T>::div(MathOps<T>::log(MathOps<T>::div(log_zn, MathOps<T>::ln2())), MathOps<T>::ln2());

            result.smooth_iter = MathOps<T>::sub(MathOps<T>::add(MathOps<T>::cast_from_double(i), MathOps<T>::cast_from_double(1.0)), nu);
            result.final_z = {z_real2, z_img2};
            return result;
        }

        z.real = z_real2;
        z.img = z_img2;
    }

    result.smooth_iter = MathOps<T>::cast_from_double(MAX_ITER);
    result.final_z = z;
    return result;
}


template<typename T>
PyObject* mandelbrot_set(int width, int height, int max_iter, T center_r, T center_i, T zoom)
{
    unsigned char* img = (unsigned char *)malloc(width * height * 3);
    if (!img)
    {
        PyErr_NoMemory();
        return NULL;
    }

    T RADIUS = MathOps<T>::cast_from_double(2.0);
    T RADIUS2 = MathOps<T>::mul(RADIUS, RADIUS);
    T zoom_q = MathOps<T>::div(MathOps<T>::cast_from_double(1.0), zoom);

    // int total_pixels = 0;
#pragma omp parallel for
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            T cr = MathOps<T>::mul(MathOps<T>::div(MathOps<T>::sub(MathOps<T>::from_int(x), MathOps<T>::from_int(width / 2)), MathOps<T>::from_int(width / 2)), RADIUS);
            T ci = MathOps<T>::mul(MathOps<T>::div(MathOps<T>::sub(MathOps<T>::from_int(y), MathOps<T>::from_int(height / 2)), MathOps<T>::from_int(height / 2)), RADIUS);

            cr = MathOps<T>::add(MathOps<T>::mul(cr, zoom_q), center_r);
            ci = MathOps<T>::add(MathOps<T>::mul(ci, zoom_q), center_i);

            MandelbrotResult<T> result = mandelbrot(Complex<T>{cr, ci}, max_iter, RADIUS2);

            Color color;
            if (MathOps<T>::cast_to_double(result.smooth_iter) == max_iter)
            {
                // Interior point
                Complex<T> c = {cr, ci};
                T distance = estimate_interior_distance(c, max_iter);
                T interior_t = MathOps<T>::sub(distance, MathOps<T>::mul(MathOps<T>::cast_from_double(floor(MathOps<T>::cast_to_double(distance))), MathOps<T>::cast_from_double(1.0)));
                color = get_color(1.0, MathOps<T>::cast_to_double(interior_t));
            }
            else
            {
                T t = MathOps<T>::div(result.smooth_iter, MathOps<T>::cast_from_double(max_iter));
                color = get_color(MathOps<T>::cast_to_double(t), 0.0);
            }

            img[(y * width + x) * 3 + 0] = color.r;
            img[(y * width + x) * 3 + 1] = color.g;
            img[(y * width + x) * 3 + 2] = color.b;
            // total_pixels++;

        }
    }
    // if (total_pixels != width * height) {
    //     PyErr_SetString(PyExc_RuntimeError, "Calculation did not complete for all pixels");
    //     free(img);
    //     return NULL;
    // }

    // int non_zero_pixels = 0;
    // for (int i = 0; i < width * height * 3; i++) {
    //     if (img[i] != 0) {
    //         non_zero_pixels++;
    //     }
    // }
    // printf("Non-zero pixels in img buffer: %d out of %d\n", non_zero_pixels, width * height * 3);

    PyObject* result = PyBytes_FromStringAndSize((char *)img, width * height * 3);

    // if (result)
    // {
    //     Py_ssize_t result_size = PyBytes_Size(result);
    //     printf("Created bytes object size: %zd\n", result_size);
    //     if (result_size == 0)
    //     {
    //         printf("Warning: Created bytes object has zero size\n");
    //     }
    // }


    if (!result) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create bytes object");
        return NULL;
    }
    // printf("Bytes object created successfully\n");
    free(img);
    return result;
}

static PyObject * mandelbrot_set_wrapper(PyObject * self, PyObject * args)
{
    int width, height, max_iter;
    const char * center_r_str, *center_i_str, *zoom_str;
    int use_quad_precision;

    if(!PyArg_ParseTuple(args, "iiisssi", &width, &height, &max_iter, &center_r_str, &center_i_str, &zoom_str, &use_quad_precision))
        return NULL;

    PyObject * result;
    const char * precision_used;

    if(use_quad_precision == 1)
    {
        Sleef_quad center_r = Sleef_strtoq(center_r_str, NULL);
        Sleef_quad center_i = Sleef_strtoq(center_i_str, NULL);
        Sleef_quad zoom = Sleef_strtoq(zoom_str, NULL);

        result = mandelbrot_set<Sleef_quad>(width, height, max_iter, center_r, center_i, zoom);
        precision_used = MathOps<Sleef_quad>::dtype_name();
    }
    else if(use_quad_precision == 2)
    {
        __float128 center_r = strtoflt128(center_r_str, NULL);
        __float128 center_i = strtoflt128(center_i_str, NULL);
        __float128 zoom = strtoflt128(zoom_str, NULL);
        result = mandelbrot_set<__float128>(width, height, max_iter, center_r, center_i, zoom);
        precision_used = MathOps<__float128>::dtype_name();
    }
    else if(use_quad_precision == 3)
    {
        double center_r = std::stod(center_r_str);
        double center_i = std::stod(center_i_str);
        double zoom = std::stod(zoom_str);
        result = mandelbrot_set<double>(width, height, max_iter, center_r, center_i, zoom);
        precision_used = MathOps<double>::dtype_name();
    }
    else
    {
        long double center_r = std::stold(center_r_str);
        long double center_i = std::stold(center_i_str);
        long double zoom = std::stold(zoom_str);
        result = mandelbrot_set<long double>(width, height, max_iter, center_r, center_i, zoom);
        precision_used = MathOps<long double>::dtype_name();
    }

    PyObject* return_tuple = PyTuple_New(2);
    PyTuple_SetItem(return_tuple, 0, result);
    PyTuple_SetItem(return_tuple, 1, PyUnicode_FromString(precision_used));

    return return_tuple;
}

static PyMethodDef MandelbrotMethods[] =
{
    {"mandelbrot_set", mandelbrot_set_wrapper, METH_VARARGS, "Calculate Mandelbrot set"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef mandelbrotmodule =
{
    PyModuleDef_HEAD_INIT,
    "mandelbrot",
    NULL,
    -1,
    MandelbrotMethods
};

PyMODINIT_FUNC PyInit_mandelbrot(void)
{
    return PyModule_Create(&mandelbrotmodule);
}
