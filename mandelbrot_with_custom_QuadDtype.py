import numpy as np
from quaddtype import QuadPrecDType, QuadPrecision
import matplotlib.pyplot as plt


def get_color(t, interior_t):
    epsilon = QuadPrecision("1e-10")

    if abs(t - QuadPrecision(1.0)) < epsilon:
        value = int(255 * interior_t)
        return np.array([value, value, value], dtype=np.uint8)

    t = np.sqrt(t)
    t = np.mod(t * 20, QuadPrecision(1.0))

    if t < 0.16:
        return np.array([0, np.int64(255 * (t / 0.16)), np.int64(128 + 127 * (t / 0.16))], dtype=np.uint8)
    elif t < 0.33:
        return np.array([0, 255, np.int64(255 * (1 - (t - 0.16) / 0.17))], dtype=np.uint8)
    elif t < 0.5:
        return np.array([np.int64(255 * ((t - 0.33) / 0.17)), 255, 0], dtype=np.uint8)
    elif t < 0.66:
        return np.array([255, np.int64(255 * (1 - (t - 0.5) / 0.16)), 0], dtype=np.uint8)
    elif t < 0.83:
        return np.array([255, 0, np.int64(255 * ((t - 0.66) / 0.17))], dtype=np.uint8)
    else:
        return np.array([np.int64(255 * (1 - (t - 0.83) / 0.17)), 0, np.int64(128 * ((t - 0.83) / 0.17))], dtype=np.uint8)


def iterate_and_compute_derivatives(cr, ci, max_iter):
    zr, zi = 0.0, 0.0
    dzr, dzi = 1.0, 0.0
    dcr, dci = 0.0, 0.0
    dzdzr, dzdzi = 0.0, 0.0

    for _ in range(max_iter):
        dzdzr_new = 2 * (zr * dzdzr - zi * dzdzi + dzr * dzr - dzi * dzi)
        dzdzi_new = 2 * (zr * dzdzi + zi * dzdzr + 2 * dzr * dzi)
        dzr_new = 2 * (zr * dzr - zi * dzi) + dcr
        dzi_new = 2 * (zr * dzi + zi * dzr) + dci
        zr_new = zr * zr - zi * zi + cr
        zi_new = 2 * zr * zi + ci

        dzdzr, dzdzi = dzdzr_new, dzdzi_new
        dzr, dzi = dzr_new, dzi_new
        zr, zi = zr_new, zi_new
        dcr, dci = 1.0, 0.0

    return zr, zi, dzr, dzi, dcr, dci, dzdzr, dzdzi


def estimate_interior_distance(cr, ci, max_iter):
    zr, zi, dzr, dzi, dcr, dci, dzdzr, dzdzi = iterate_and_compute_derivatives(
        cr, ci, max_iter)

    dz_abs_sq = dzr * dzr + dzi * dzi
    numerator = 1 - dz_abs_sq

    denominator = np.abs((dcr * dzr + dci * dzi) +
                         (dzdzr * zr + dzdzi * zi) * dcr)

    return numerator / denominator


iters = []
values = []


def mandelbrot(cr, ci, max_iter, radius2):
    zr, zi = QuadPrecision(0.0), QuadPrecision(0.0)
    for i in range(max_iter):
        zr_new = zr * zr - zi * zi + cr
        zi_new = QuadPrecision(2) * zr * zi + ci
        zr, zi = zr_new, zi_new
        if zr * zr + zi * zi > radius2:
            iters.append(i)
            values.append(np.float64(zr * zr + zi * zi))
            log_zn = np.log(zr * zr + zi * zi) / QuadPrecision(2)
            nu = np.log(log_zn / np.log(QuadPrecision(2))) / \
                np.log(QuadPrecision(2))
            return i + QuadPrecision(1) - nu, zr, zi
    return max_iter, zr, zi


def calculate_complex_coordinates(width, height, RADIUS, zoom_q, center_r, center_i):
    cr = np.zeros((height, width), dtype=QuadPrecDType)
    ci = np.zeros((height, width), dtype=QuadPrecDType)

    for y in range(height):
        for x in range(width):
            cr[y, x] = ((x - width / 2) / (width / 2)) * RADIUS
            ci[y, x] = ((y - height / 2) / (height / 2)) * RADIUS
            cr[y, x] = cr[y, x] * zoom_q + center_r
            ci[y, x] = ci[y, x] * zoom_q + center_i

    return cr, ci


def mandelbrot_set(width, height, max_iter, center_r, center_i, zoom):
    radius = QuadPrecision(2.0)
    radius2 = radius * radius
    zoom = 1 / zoom

    # x = np.linspace(np.float64(center_r - radius / zoom),
    #                 np.float64(center_r + radius / zoom),
    #                 width)
    # y = np.linspace(np.float64(center_i - radius / zoom),
    #                 np.float64(center_i + radius / zoom),
    #                 height)
    # cr, ci = np.meshgrid(x, y)
    cr, ci = calculate_complex_coordinates(
        width, height, radius, zoom, center_r, center_i)

    smooth_iter = np.zeros((height, width), dtype=QuadPrecDType)
    final_zr = np.zeros((height, width), dtype=QuadPrecDType)
    final_zi = np.zeros((height, width), dtype=QuadPrecDType)

    for i in range(height):
        for j in range(width):
            smooth_iter[i, j], final_zr[i, j], final_zi[i, j] = mandelbrot(
                cr[i, j], ci[i, j], max_iter, radius2)
        print(f"Row {i}/{height} is completed")

    img = np.zeros((height, width, 3), dtype=np.uint8)

    interior_mask = smooth_iter == QuadPrecision(max_iter)
    interior_cr = cr[interior_mask]
    interior_ci = ci[interior_mask]
    interior_distance = np.array([estimate_interior_distance(
        cr, ci, max_iter) for cr, ci in zip(interior_cr, interior_ci)])
    interior_t = (interior_distance -
                  np.floor(interior_distance)).astype(np.float64)

    exterior_mask = ~interior_mask
    t = smooth_iter[exterior_mask] / QuadPrecision(max_iter)

    interior_colors = np.array(
        [get_color(QuadPrecision(1.0), t) for t in interior_t])
    exterior_colors = np.array([get_color(t, QuadPrecision(0.0)) for t in t])

    img[interior_mask] = interior_colors
    img[exterior_mask] = exterior_colors

    return img


def plot_mandelbrot(width, height, max_iter, center_r, center_i, zoom):
    center_rq = QuadPrecision(center_r)
    center_iq = QuadPrecision(center_i)
    zoom_q = QuadPrecision(zoom)

    img_array = mandelbrot_set(
        width, height, max_iter, center_rq, center_iq, zoom_q)

    plt.figure(figsize=(10, 10))
    plt.imshow(img_array)
    plt.axis('off')
    plt.title(
        f'Mandelbrot Set (zoom: {zoom}, center: {center_r} + {center_i}i, iterations: {max_iter}, dtype: {type(center_rq)})')
    plt.show()


if __name__ == "__main__":
    width = 800
    height = 800
    max_iter = 3000
    center_r = -2.0
    center_i = 0.0
    zoom = 1e15

    plot_mandelbrot(width, height, max_iter, center_r, center_i, zoom)
