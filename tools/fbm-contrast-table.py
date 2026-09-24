#!/usr/bin/env python3
"""Writes the contrast curve of src/include/fbm_noise.h.

One octave of gradient noise swings about seven tenths as far from its middle
as one octave of the value noise it replaced, so a picture drawn with the same
numbers came out paler. The curve takes the values of the one onto the
distribution of the other: it is the quantile map F_old^-1(F_new(v)), measured
over eight million points, laid out as a monotone cubic Hermite table
(Fritsch and Carlson) in the distance from the middle, which is what keeps it
monotone -- a polynomial fitted to the same map turned back on itself near the
ends and left nought to one by a thousandth -- and written out over the value
itself, one cubic a step, which is what keeps it cheap.

Both noises are written here as the C++ writes them: the same hash, the same
sixteen directions, the same quintic and the same gain for the new one, and
randsc's old value noise -- four hashed corners and a smoothstep -- for the
other. The seed of the sampling is fixed, so the table comes out the same
every time.

    python tools/fbm-contrast-table.py

writes src/include/fbm_contrast_table.h, with what it measured at the top.
Needs numpy.
"""
import os
import sys

import numpy as np

KNOTS = 32  # knots over the distance from the middle, nought to a half

np.seterr(over="ignore")
U = np.uint64


def mix(h):
    h ^= h >> U(33)
    h *= U(0xFF51AFD7ED558CCD)
    h ^= h >> U(33)
    h *= U(0xC4CEB9FE1A85EC53)
    h ^= h >> U(33)
    return h


def corner(cx, cy, key):
    return mix((cx.astype(np.int64).view(U) * U(0x9E3779B97F4A7C15)) ^
               (cy.astype(np.int64).view(U) * U(0xC2B2AE3D27D4EB4F)) ^ U(key))


ANGLES = (np.arange(16) + 0.5) * np.pi / 8
DIRS = np.stack([np.cos(ANGLES), np.sin(ANGLES)], 1)
GAIN = np.sqrt(0.5)


def gradient(x, y, key):
    cx, cy = np.floor(x), np.floor(y)
    u, v = x - cx, y - cy
    fu = u * u * u * (u * (u * 6 - 15) + 10)
    fv = v * v * v * (v * (v * 6 - 15) + 10)

    def slope(i, j):
        d = DIRS[(corner(cx + i, cy + j, key) >> U(60)).astype(int)]
        return d[:, 0] * (u - i) + d[:, 1] * (v - j)

    a, b, c, d = slope(0, 0), slope(1, 0), slope(0, 1), slope(1, 1)
    lo = a + (b - a) * fu
    hi = c + (d - c) * fu
    return 0.5 + (lo + (hi - lo) * fv) * GAIN


def value(x, y, key):
    cx, cy = np.floor(x), np.floor(y)
    u, v = x - cx, y - cy
    su, sv = u * u * (3 - 2 * u), v * v * (3 - 2 * v)
    unit = lambda h: (h >> U(11)).astype(np.float64) / 9007199254740992.0
    a = unit(corner(cx, cy, key))
    b = unit(corner(cx + 1, cy, key))
    c = unit(corner(cx, cy + 1, key))
    d = unit(corner(cx + 1, cy + 1, key))
    lo = a + (b - a) * su
    hi = c + (d - c) * su
    return lo + (hi - lo) * sv


def main():
    rng = np.random.default_rng(11)
    old, new = [], []
    for k in range(8):
        p = rng.uniform(-20000, 20000, size=(2, 1000000))
        old.append(np.abs(value(p[0], p[1], 1000 + k) - 0.5))
        new.append(np.abs(gradient(p[0], p[1], 1000 + k) - 0.5))
    old, new = np.concatenate(old), np.concatenate(new)

    # the quantile map, on the distance from the middle, at the knots
    step = 0.5 / KNOTS
    knots = np.arange(KNOTS + 1) * step
    share = np.searchsorted(np.sort(new), knots) / new.size
    vals = np.quantile(old, np.clip(share, 0, 1))
    vals[0], vals[-1] = 0.0, 0.5
    for i in range(1, KNOTS + 1):
        vals[i] = max(vals[i], vals[i - 1] + 1e-9)

    # Fritsch-Carlson: secants, their averages, pulled in where they would
    # overshoot, which is what makes the cubic between two knots monotone
    d = np.diff(vals) / step
    m = np.empty(KNOTS + 1)
    m[1:-1] = (d[:-1] + d[1:]) / 2
    m[0], m[-1] = d[0], d[-1]
    for i in range(KNOTS):
        a, b = m[i] / d[i], m[i + 1] / d[i]
        if a * a + b * b > 9:
            tau = 3 / np.sqrt(a * a + b * b)
            m[i], m[i + 1] = tau * a * d[i], tau * b * d[i]

    def curve(t):
        s = np.minimum(t, 0.5) / step
        i = np.minimum(np.floor(s).astype(int), KNOTS - 1)
        x = s - i
        return ((2 * x ** 3 - 3 * x ** 2 + 1) * vals[i] +
                (x ** 3 - 2 * x ** 2 + x) * step * m[i] +
                (-2 * x ** 3 + 3 * x ** 2) * vals[i + 1] +
                (x ** 3 - x ** 2) * step * m[i + 1])

    mapped = curve(new)
    q = np.linspace(0.001, 0.999, 999)
    rms = lambda t: np.sqrt((t * t).mean())

    # Written out over the value itself, nought to one in 2 * KNOTS steps, as
    # the cubic in the place x within each step: c0 + c1 x + c2 x^2 + c3 x^3.
    # A step of the value is a step of the distance, mirrored below the
    # middle, so each is one piece of the Hermite curve expanded. Looking the
    # value up directly spares the distance, the side and the Hermite basis at
    # run time, which were most of what the curve cost.
    coeffs = []
    for k in range(2 * KNOTS):
        j = k - KNOTS if k >= KNOTS else KNOTS - 1 - k
        y0, y1 = vals[j], vals[j + 1]
        m0, m1 = m[j] * step, m[j + 1] * step
        a = 3 * (y1 - y0) - 2 * m0 - m1
        b = 2 * (y0 - y1) + m0 + m1
        if k >= KNOTS:
            # above the middle: 0.5 + H(x)
            coeffs.append((0.5 + y0, m0, a, b))
        else:
            # below it: 0.5 - H(1 - x), H expanded about the other end
            coeffs.append((0.5 - (y0 + m0 + a + b), m0 + 2 * a + 3 * b,
                           -(a + 3 * b), b))
    coeffs = np.array(coeffs)

    def table(v):
        s = v * 2 * KNOTS
        k = np.minimum(s.astype(int), 2 * KNOTS - 1)
        x = s - k
        c = coeffs[k]
        return c[:, 0] + x * (c[:, 1] + x * (c[:, 2] + x * c[:, 3]))

    fine = table(np.linspace(0, 1, 200001))
    ends = [c[0] + c[1] + c[2] + c[3] for c in coeffs]
    meet = max(abs(ends[k] - coeffs[k + 1][0]) for k in range(2 * KNOTS - 1))
    worst = np.abs(np.quantile(mapped, q) - np.quantile(old, q)).max()
    out = [
        "/* The contrast curve of fbm_noise.h. Written by",
        " * tools/fbm-contrast-table.py, which measures it; change that and",
        " * run it rather than editing this.",
        " *",
        " * Over the value, nought to one in %d steps: in each, the cubic in"
        % (2 * KNOTS),
        " * the place x within the step, lowest power first.",
        " *",
        " * Measured: %.4f from the middle, root mean square, against %.4f for"
        % (rms(mapped), rms(old)),
        " * value noise and %.4f before the curve; the quantiles of the two"
        % rms(new),
        " * agree to %.4f at worst; the curve runs from %.3g to %.17g and"
        % (worst, fine[0], fine[-1]),
        " * never steps back -- its smallest step over two hundred thousand",
        " * is %.2g, forward -- and its pieces meet to %.1g. */"
        % (np.diff(fine).min(), meet),
        "",
        "static const double FBM_NOISE_CURVE[%d][4] = {" % (2 * KNOTS),
        ",\n".join("    {%.17g, %.17g, %.17g, %.17g}" % tuple(c)
                   for c in coeffs) + "};",
    ]
    here = os.path.dirname(os.path.abspath(__file__))
    target = os.path.join(here, "..", "src", "include",
                          "fbm_contrast_table.h")
    with open(target, "w", newline="\n") as f:
        f.write("\n".join(out) + "\n")
    print("\n".join(out[7:13]), file=sys.stderr)
    print("-> " + os.path.normpath(target), file=sys.stderr)


if __name__ == "__main__":
    main()
