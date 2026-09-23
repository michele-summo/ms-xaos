#!/usr/bin/env python3
"""Writes src/sffe/randsctile_tables.h, the tilings randsctile draws.

    python tools/randsctile-tables.py

Every periodic tiling is built here from edge length one, checked -- twenty
thousand points at random, each of which must fall in exactly one tile -- and
scaled so that a tile has unit area on average, as the rest of the randsc
family lays one cell over each unit of area. What is written is the period of
each tiling, the vertices of the tiles of one period, and for each tile its
middle, the radius of the largest circle it holds, and its bounding box.

The tilings that do not repeat are not tables: the multigrid ones are worked
out from a handful of numbers, written here too, and the three made by
substitution are code in sffe_cmplx_gsl.cpp, with only their scale written
here. Needs numpy.
"""
import itertools
import math
import os
import sys

import numpy as np

S2, S3 = math.sqrt(2), math.sqrt(3)
PHI = (1 + math.sqrt(5)) / 2


# --- geometry ------------------------------------------------------------------

def ngon(n, cx, cy, first_deg):
    """a regular n-gon of edge one round (cx, cy), its first corner at that angle"""
    R = 1 / (2 * math.sin(math.pi / n))
    return [(cx + R * math.cos(math.radians(first_deg + 360 * k / n)),
             cy + R * math.sin(math.radians(first_deg + 360 * k / n))) for k in range(n)]


def signed_area(poly):
    return 0.5 * sum(poly[i][0] * poly[(i + 1) % len(poly)][1] -
                     poly[(i + 1) % len(poly)][0] * poly[i][1] for i in range(len(poly)))


def ccw(poly):
    return poly if signed_area(poly) > 0 else poly[::-1]


def centroid(poly):
    A = cx = cy = 0
    for i in range(len(poly)):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % len(poly)]
        c = x0 * y1 - x1 * y0
        A += c
        cx += (x0 + x1) * c
        cy += (y0 + y1) * c
    return (cx / (3 * A), cy / (3 * A))


def inside_convex(poly, p, eps=1e-9):
    for i in range(len(poly)):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % len(poly)]
        if (x1 - x0) * (p[1] - y0) - (y1 - y0) * (p[0] - x0) < -eps:
            return False
    return True


def edge_distance(poly, p):
    best = 1e30
    for i in range(len(poly)):
        a = np.array(poly[i])
        b = np.array(poly[(i + 1) % len(poly)])
        d = b - a
        t = max(0, min(1, np.dot(np.array(p) - a, d) / np.dot(d, d)))
        best = min(best, np.linalg.norm(np.array(p) - (a + t * d)))
    return best


def reach(poly):
    """how far the point of a tile furthest from its edge stands from it: the
    radius of the largest circle the tile holds, found on a grid and then
    refined round the best of it. For a regular polygon it is the apothem; for
    a concave one, a star or an L, it is not the distance from the middle,
    which can sit close to a notch."""
    xs, ys = [p[0] for p in poly], [p[1] for p in poly]
    best, at = -1, None
    step = max(max(xs) - min(xs), max(ys) - min(ys)) / 60
    for _ in range(6):
        cx, cy = (at if at else ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2))
        span = 30 if at is None else 4
        for i in range(-span, span + 1):
            for j in range(-span, span + 1):
                q = (cx + i * step, cy + j * step)
                if not inside_any(poly, q):
                    continue
                d = edge_distance(poly, q)
                if d > best:
                    best, at = d, q
        step /= 4
    return best


def inside_any(poly, p):
    x, y = p
    c = False
    for i in range(len(poly)):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % len(poly)]
        if (y0 > y) != (y1 > y) and x < x0 + (y - y0) * (x1 - x0) / (y1 - y0):
            c = not c
    return c


def translates(faces, a, b, r=2):
    return [[(x + i * a[0] + j * b[0], y + i * a[1] + j * b[1]) for x, y in f]
            for i in range(-r, r + 1) for j in range(-r, r + 1) for f in faces]


def lattice_key(a, b, p):
    inv = np.linalg.inv(np.array([a, b]).T)
    l = inv @ np.array(p)
    return tuple(np.round(l - np.floor(l + 1e-9), 6)), np.floor(l + 1e-9)


def fill_triangles(faces, a, b):
    """the unit equilateral triangles left between the faces given: every three
    corners one apart whose middle no face covers, one of each translation"""
    allf = translates(faces, a, b)
    V = list({(round(x, 7), round(y, 7)): (x, y) for f in allf for x, y in f}.values())
    P = np.array(V)
    unit = np.abs(np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1)) - 1) < 1e-7
    tris, seen = [], set()
    for i in range(len(V)):
        for j, k in itertools.combinations(np.where(unit[i])[0], 2):
            if j < i or k < i or not unit[j, k]:
                continue
            t = ccw([V[i], V[j], V[k]])
            c = centroid(t)
            if any(inside_convex(f, c) for f in allf):
                continue
            key, shift = lattice_key(a, b, c)
            if key in seen:
                continue
            seen.add(key)
            d = shift[0] * np.array(a) + shift[1] * np.array(b)
            tris.append([(x - d[0], y - d[1]) for x, y in t])
    return tris


def dual(faces, a, b):
    """the Laves dual: a tile for every corner, its corners the middles of the
    tiles round that corner, taken in order of angle"""
    allf = translates(faces, a, b)
    mids = [centroid(f) for f in allf]
    corners = {}
    for f in faces:
        for v in f:
            corners.setdefault(lattice_key(a, b, v)[0], v)
    out = []
    for v in corners.values():
        ring = [m for f, m in zip(allf, mids)
                if any(abs(x - v[0]) < 1e-7 and abs(y - v[1]) < 1e-7 for x, y in f)]
        ring.sort(key=lambda m: math.atan2(m[1] - v[1], m[0] - v[0]))
        out.append(ccw(ring))
    return out


def star(n, k, cx, cy, R, first_deg):
    """the star {n/k} with its points on the circle of radius R"""
    rho = R * math.cos(math.pi * k / n) / math.cos(math.pi * (k - 1) / n)
    pts = []
    for i in range(n):
        a = math.radians(first_deg + 360 * i / n)
        b = math.radians(first_deg + 360 * (i + 0.5) / n)
        pts.append((cx + R * math.cos(a), cy + R * math.sin(a)))
        pts.append((cx + rho * math.cos(b), cy + rho * math.sin(b)))
    return pts


def round_middle(poly):
    c = (sum(p[0] for p in poly) / len(poly), sum(p[1] for p in poly) / len(poly))
    return sorted(poly, key=lambda p: math.atan2(p[1] - c[1], p[0] - c[0]))


def starred(tiling, big, k):
    """Every polygon of `big` sides or more made a star {n/k} with its points on
    the polygon's corners. The thin triangle between a notch of the star and the
    polygon's side goes to the tile across that side: two of them across a side
    two stars share make a rhombus, and a small polygon takes one on each side it
    shares with a star, which makes it a star of its own."""
    a, b, faces = tiling
    allf = translates(faces, a, b)

    def key(p):
        return (round(p[0], 6), round(p[1], 6))

    notches = {}
    for f in allf:
        n = len(f)
        if n < big:
            continue
        c = centroid(f)
        R = math.hypot(f[0][0] - c[0], f[0][1] - c[1])
        rho = R * math.cos(math.pi * k / n) / math.cos(math.pi * (k - 1) / n)
        for i in range(n):
            P, Q = f[i], f[(i + 1) % n]
            m = ((P[0] + Q[0]) / 2 - c[0], (P[1] + Q[1]) / 2 - c[1])
            L = math.hypot(*m)
            notches.setdefault(tuple(sorted((key(P), key(Q)))), []).append(
                (c[0] + m[0] / L * rho, c[1] + m[1] / L * rho))
    out, seen = [], set()
    for f in faces:
        n = len(f)
        if n >= big:
            c = centroid(f)
            R = math.hypot(f[0][0] - c[0], f[0][1] - c[1])
            out.append(ccw(star(n, k, c[0], c[1], R,
                                math.degrees(math.atan2(f[0][1] - c[1], f[0][0] - c[0])))))
            for i in range(n):
                P, Q = f[i], f[(i + 1) % n]
                e = tuple(sorted((key(P), key(Q))))
                if len(notches.get(e, [])) != 2:
                    continue
                lk = lattice_key(a, b, ((P[0] + Q[0]) / 2, (P[1] + Q[1]) / 2))[0]
                if lk in seen:
                    continue
                seen.add(lk)
                out.append(ccw(round_middle([P, notches[e][0], Q, notches[e][1]])))
            continue
        pts = []
        for i in range(n):
            P, Q = f[i], f[(i + 1) % n]
            pts.append(P)
            e = tuple(sorted((key(P), key(Q))))
            if e in notches:
                pts.append(notches[e][0])
        out.append(ccw(pts))
    return (a, b, out)


def hexagrams():
    """six-pointed stars on a triangular lattice, point to point, and the regular
    hexagon each gap between three of them turns out to be"""
    r = 1.0
    a, b = (2 * r, 0.0), (r, r * S3)
    s = star(6, 2, 0, 0, r, 0)

    def gap(c):
        pts = []
        for sx, sy in [(0, 0), a, b, (a[0] + b[0], a[1] + b[1])]:
            for q in s:
                P = (sx + q[0], sy + q[1])
                if math.hypot(P[0] - c[0], P[1] - c[1]) < 0.6 * r and \
                        all(math.hypot(P[0] - Q[0], P[1] - Q[1]) > 1e-9 for Q in pts):
                    pts.append(P)
        return ccw(round_middle(pts))

    return (a, b, [ccw(s), gap(((a[0] + b[0]) / 3, (a[1] + b[1]) / 3)),
                   gap((2 * (a[0] + b[0]) / 3, 2 * (a[1] + b[1]) / 3))])


def rows(pattern):
    """rows of squares (S) and of triangles (T), one above the other"""
    faces, y, x0, h = [], 0.0, 0.0, S3 / 2
    for r in pattern:
        if r == "S":
            faces.append([(x0, y), (x0 + 1, y), (x0 + 1, y + 1), (x0, y + 1)])
            y += 1
        else:
            faces.append([(x0, y), (x0 + 1, y), (x0 + 0.5, y + h)])
            faces.append([(x0 + 1, y), (x0 + 1.5, y + h), (x0 + 0.5, y + h)])
            x0 += 0.5
            y += h
    return ((1.0, 0.0), (x0, y), faces)


def hexagons_among_triangles():
    """the triangular tiling with the six triangles round every point of a
    sublattice three apart merged into a hexagon"""
    ea, eb = np.array((1.0, 0.0)), np.array((0.5, S3 / 2))
    A, B = 3 * ea, 3 * eb
    inv = np.linalg.inv(np.array([A, B]).T)
    faces = [ngon(6, 0, 0, 0)]
    for i in range(-8, 9):
        for j in range(-8, 9):
            o = i * ea + j * eb
            for t in ([o, o + ea, o + eb], [o + ea, o + ea + eb, o + eb]):
                c = sum(t) / 3
                l = inv @ c
                if not (0 <= l[0] < 1 and 0 <= l[1] < 1):
                    continue
                near = np.round(l)
                if min(np.linalg.norm(c - ((near[0] + di) * A + (near[1] + dj) * B))
                       for di in (-1, 0, 1) for dj in (-1, 0, 1)) < 0.99:
                    continue
                faces.append([tuple(v) for v in t])
    return (tuple(A), tuple(B), faces)


# --- the periodic tilings, in the order randsctile numbers them ----------------

def build():
    T = []
    add = lambda name, t: T.append((name, t))

    add("squares", ((1, 0), (0, 1), [[(0, 0), (1, 0), (1, 1), (0, 1)]]))
    add("triangles", ((1, 0), (0.5, S3 / 2),
                      [[(0, 0), (1, 0), (0.5, S3 / 2)], [(1, 0), (1.5, S3 / 2), (0.5, S3 / 2)]]))
    add("hexagons", ((S3, 0), (S3 / 2, 1.5), [ngon(6, 0, 0, 30)]))

    p = 1 + S2
    t488 = ((p, 0), (0, p), [ngon(8, 0, 0, 22.5), ngon(4, p / 2, p / 2, 0)])
    add("octagons and squares, 4.8.8", t488)
    a, b = (2, 0), (1, S3)
    t3636 = (a, b, [ngon(6, 0, 0, 0)] + fill_triangles([ngon(6, 0, 0, 0)], a, b))
    add("trihexagonal, 3.6.3.6", t3636)
    L = 1 + S3
    a, b = (L * S3 / 2, L / 2), (0, L)
    hexa = ngon(6, 0, 0, 0)
    sq = []
    for nd in [30, 90, 150]:
        n = (math.cos(math.radians(nd)), math.sin(math.radians(nd)))
        v0, v1 = hexa[(nd // 60) % 6], hexa[(nd // 60 + 1) % 6]
        sq.append(ccw([v0, v1, (v1[0] + n[0], v1[1] + n[1]), (v0[0] + n[0], v0[1] + n[1])]))
    t3464 = (a, b, [hexa] + sq + fill_triangles([hexa] + sq, a, b))
    add("rhombitrihexagonal, 3.4.6.4", t3464)
    L = 2 + S3
    a, b = (L, 0), (L / 2, L * S3 / 2)
    t31212 = (a, b, [ngon(12, 0, 0, 15)] + fill_triangles([ngon(12, 0, 0, 15)], a, b))
    add("truncated hexagonal, 3.12.12", t31212)
    L = 3 + S3
    a, b = (L, 0), (L / 2, L * S3 / 2)
    dod = ngon(12, 0, 0, 15)
    sqs = []
    for nd in [0, 60, 120]:
        n = (math.cos(math.radians(nd)), math.sin(math.radians(nd)))
        k0 = int(((nd - 15) % 360 - 15) / 30) % 12
        v0, v1 = dod[k0], dod[(k0 + 1) % 12]
        sqs.append(ccw([v0, v1, (v1[0] + n[0], v1[1] + n[1]), (v0[0] + n[0], v0[1] + n[1])]))
    c1 = ((a[0] + b[0]) / 3, (a[1] + b[1]) / 3)
    c2 = (2 * (a[0] + b[0]) / 3, 2 * (a[1] + b[1]) / 3)
    t4612 = (a, b, [dod] + sqs + [ngon(6, c1[0], c1[1], 0), ngon(6, c2[0], c2[1], 0)])
    add("truncated trihexagonal, 4.6.12", t4612)
    t33344 = ((1, 0), (0.5, 1 + S3 / 2),
              [[(0, 0), (1, 0), (1, 1), (0, 1)], [(0, 1), (1, 1), (0.5, 1 + S3 / 2)],
               [(1, 1), (1.5, 1 + S3 / 2), (0.5, 1 + S3 / 2)]])
    add("elongated triangular, 3.3.3.4.4", t33344)
    p = math.sqrt(2 + S3)
    a, b = (p, 0), (0, p)
    sqa, sqb = ngon(4, 0, 0, 60), ngon(4, p / 2, p / 2, 30)
    t33434 = (a, b, [sqa, sqb] + fill_triangles([sqa, sqb], a, b))
    add("snub square, 3.3.4.3.4", t33434)
    a = (2.5, S3 / 2)
    b = (a[0] * 0.5 - a[1] * S3 / 2, a[0] * S3 / 2 + a[1] * 0.5)
    t33336 = (a, b, [ngon(6, 0, 0, 0)] + fill_triangles([ngon(6, 0, 0, 0)], a, b))
    add("snub hexagonal, 3.3.3.3.6", t33336)

    for name, t in [("tetrakis square", t488), ("rhombille", t3636),
                    ("deltoidal trihexagonal", t3464), ("triakis triangular", t31212),
                    ("kisrhombille", t4612), ("prismatic pentagonal", t33344),
                    ("Cairo pentagonal", t33434), ("floret pentagonal", t33336)]:
        add(name, (t[0], t[1], dual(t[2], t[0], t[1])))

    add("bricks", ((2, 0), (1, 1), [[(0, 0), (2, 0), (2, 1), (0, 1)]]))
    add("Flemish bond", ((3, 0), (1.5, 1),
                         [[(0, 0), (2, 0), (2, 1), (0, 1)], [(2, 0), (3, 0), (3, 1), (2, 1)]]))
    add("herringbone", ((1, 1), (2, -2),
                        [[(0, 0), (2, 0), (2, 1), (0, 1)], [(0, 1), (1, 1), (1, 3), (0, 3)]]))
    add("basketweave", ((4, 0), (2, 2),
                        [[(0, 0), (2, 0), (2, 1), (0, 1)], [(0, 1), (2, 1), (2, 2), (0, 2)],
                         [(2, 0), (3, 0), (3, 2), (2, 2)], [(3, 0), (4, 0), (4, 2), (3, 2)]]))
    add("Pythagorean", ((2, 1), (-1, 2),
                        [[(0, 0), (2, 0), (2, 2), (0, 2)], [(2, 0), (3, 0), (3, 1), (2, 1)]]))
    add("chevrons", ((2, 0), (0, 1),
                     [[(0, 0), (1, 1), (1, 2), (0, 1)], [(1, 1), (2, 0), (2, 1), (1, 2)]]))
    h = S2 / 2
    add("squares and rhombi", ((1, 0), (h, 1 + h),
                               [[(0, 0), (1, 0), (1, 1), (0, 1)],
                                [(0, 1), (1, 1), (1 + h, 1 + h), (h, 1 + h)]]))
    add("houses", ((1, 0), (0, 2.5),
                   [[(0, 0), (1, 0), (1, 1), (0.5, 1.5), (0, 1)],
                    [(1, 1), (1.5, 1.5), (1.5, 2.5), (0.5, 2.5), (0.5, 1.5)]]))
    add("rows square, square, triangle", rows("SST"))
    add("rows square, triangle, triangle", rows("STT"))
    add("hexagons among triangles", hexagons_among_triangles())
    add("Greek crosses", ((2, 1), (-1, 2),
                          [[(0.5, -0.5), (1.5, -0.5), (1.5, 0.5), (0.5, 0.5), (0.5, 1.5),
                            (-0.5, 1.5), (-0.5, 0.5), (-1.5, 0.5), (-1.5, -0.5), (-0.5, -0.5),
                            (-0.5, -1.5), (0.5, -1.5)]]))
    add("T tetrominoes", ((4, 0), (2, 2),
                          [[(0, 0), (3, 0), (3, 1), (2, 1), (2, 2), (1, 2), (1, 1), (0, 1)],
                           [(3, 0), (4, 0), (4, 1), (5, 1), (5, 2), (2, 2), (2, 1), (3, 1)]]))

    k = S2 - 1
    s8 = [(S2, 0), (1, k), (1, 1), (k, 1), (0, S2), (-k, 1), (-1, 1), (-1, k),
          (-S2, 0), (-1, -k), (-1, -1), (-k, -1), (0, -S2), (k, -1), (1, -1), (1, -k)]
    d = 2 * S2

    def chain(cx, cy, pts):
        return [(cx + x, cy + y) for x, y in pts]

    cross = ([(S2, 0)] + chain(0, 0, [(1, k), (1, 1), (k, 1)]) + [(0, S2)] +
             chain(0, d, [(k, -1), (1, -1), (1, -k)]) + [(S2, d)] +
             chain(d, d, [(-1, -k), (-1, -1), (-k, -1)]) + [(d, S2)] +
             chain(d, 0, [(-k, 1), (-1, 1), (-1, k)]))
    add("eight-pointed stars and crosses", ((d, 0), (0, d), [s8, ccw(cross)]))
    add("six-pointed stars and hexagons", hexagrams())
    add("eight-pointed stars, 4.8.8", starred(t488, 8, 3))
    add("twelve-pointed stars, 3.12.12", starred(t31212, 12, 5))
    add("twelve-pointed stars, 4.6.12", starred(t4612, 12, 4))
    return T


def normalise(t):
    """to unit area a tile, on average"""
    a, b, faces = t
    faces = [ccw([(float(x), float(y)) for x, y in f]) for f in faces]
    cell = abs(a[0] * b[1] - a[1] * b[0])
    k = math.sqrt(len(faces) / cell)
    return ((a[0] * k, a[1] * k), (b[0] * k, b[1] * k),
            [[(x * k, y * k) for x, y in f] for f in faces])


def check(name, t, n=20000):
    a, b, faces = t
    cell = abs(a[0] * b[1] - a[1] * b[0])
    tot = sum(abs(signed_area(f)) for f in faces)
    assert abs(tot - cell) < 1e-9 * cell, f"{name}: tiles {tot} over a period {cell}"
    rg = np.random.default_rng(1)
    s, u = rg.random(n), rg.random(n)
    x, y = s * a[0] + u * b[0], s * a[1] + u * b[1]
    count = np.zeros(n, int)
    for f in translates(faces, a, b, 3):
        c = np.zeros(n, bool)
        for i in range(len(f)):
            x0, y0 = f[i]
            x1, y1 = f[(i + 1) % len(f)]
            if y0 == y1:
                continue
            c ^= ((y0 > y) != (y1 > y)) & (x < x0 + (y - y0) * (x1 - x0) / (y1 - y0))
        count += c
    assert (count == 1).all(), f"{name}: {(count != 1).sum()} points not in exactly one tile"


def candidates(t):
    """the tiles, with their offset in periods, that can hold a point of the
    period [0,1)^2 in lattice coordinates -- by bounding box, so a few too many"""
    a, b, faces = t
    inv = np.linalg.inv(np.array([a, b]).T)
    out = []
    for fi, f in enumerate(faces):
        for i in range(-3, 4):
            for j in range(-3, 4):
                l = np.array([(x + i * a[0] + j * b[0], y + i * a[1] + j * b[1]) for x, y in f]) @ inv.T
                if l[:, 0].max() > 0 and l[:, 0].min() < 1 and l[:, 1].max() > 0 and l[:, 1].min() < 1:
                    out.append((fi, i, j))
    return out


# --- the multigrids ----------------------------------------------------------

GRIDS = [  # families, and the offset of each: generic, and for five summing to one
    (5, [0.1, 0.2, 0.3, 0.15, 0.25]),
    (4, [0.1234, 0.2718, 0.3141, 0.4142]),
    (6, [0.11, 0.23, 0.37, 0.41, 0.53, 0.07]),
    (7, [0.13, 0.29, 0.07, 0.43, 0.31, 0.19, 0.05]),
]


def grid_vectors(n):
    step = 2 * math.pi / n if n % 2 else math.pi / n
    return [(math.cos(j * step), math.sin(j * step)) for j in range(n)]


def grid_mean_area(n):
    """a crossing of families j and k is a rhombus of area |sin| of their angle,
    and they cross |sin| times per unit of grid: the mean is the one sum over the
    other"""
    E = grid_vectors(n)
    s = [abs(E[j][0] * E[k][1] - E[j][1] * E[k][0]) for j in range(n) for k in range(j + 1, n)]
    return sum(x * x for x in s) / sum(s)


def kite_and_dart():
    """a kite and a dart as the kites and darts tiling has them, the legs of the
    kite's halves one: each half with the other mirrored across B C"""
    def whole(A, B, C):
        u = C - B
        t = ((A - B) * u.conjugate()).real / abs(u) ** 2
        F = 2 * (B + t * u) - A
        return [(z.real, z.imag) for z in (A, B, F, C)]
    B = 0j
    A = complex(math.cos(math.radians(-18)), math.sin(math.radians(-18)))
    C = complex(math.cos(math.radians(18)), math.sin(math.radians(18)))
    kite = whole(A, B, C)
    # a half dart as a half kite of legs PHI cuts one off: R Q B
    A2, C2 = A * PHI, C * PHI
    Q = A2 + (B - A2) / PHI
    R = B + (C2 - B) / PHI
    dart = whole(R, Q, B)
    return kite, dart


def kites_mean_area():
    """Robinson's triangles for kites and darts, from a wheel of radius one
    subdivided ten times: the wheel's area over half the number of triangles,
    rescaled to the legs of one the tiles have after as many levels as the wheel
    was scaled up by"""
    tris = []
    for i in range(10):
        U = complex(math.cos((2 * i - 1) * math.pi / 10), math.sin((2 * i - 1) * math.pi / 10))
        V = complex(math.cos((2 * i + 1) * math.pi / 10), math.sin((2 * i + 1) * math.pi / 10))
        if i % 2 == 0:
            U, V = V, U
        tris.append((0, U, 0j, V))
    levels = 10
    for _ in range(levels):
        nxt = []
        for c, A, B, C in tris:
            if c == 0:
                Q = A + (B - A) / PHI
                R = B + (C - B) / PHI
                nxt += [(1, R, Q, B), (0, Q, A, R), (0, C, A, R)]
            else:
                P = C + (A - C) / PHI
                nxt += [(1, B, P, A), (0, P, C, B)]
        tris = nxt
    wheel = 10 * 0.5 * math.sin(math.pi / 5)
    return wheel * PHI ** (2 * levels) / (len(tris) / 2)


# --- writing it ----------------------------------------------------------------

def c_double(x):
    s = repr(float(x))
    return s if ("e" in s or "." in s) else s + ".0"


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    target = os.path.join(here, "..", "src", "sffe", "randsctile_tables.h")
    tilings = []
    for name, t in build():
        t = normalise(t)
        check(name, t)
        tilings.append((name, t))
        print(f"{name:34s} {len(t[2]):2d} tiles a period, covered once", file=sys.stderr)

    verts, faces, cands, rows = [], [], [], []
    for name, (a, b, fs) in tilings:
        first_face, first_cand = len(faces), len(cands)
        for f in fs:
            c = centroid(f)
            xs, ys = [p[0] for p in f], [p[1] for p in f]
            faces.append((len(verts), len(f), c, reach(f),
                          (min(xs), min(ys), max(xs), max(ys))))
            verts += f
        for fi, i, j in candidates((a, b, fs)):
            cands.append((first_face + fi, i, j))
        inv = np.linalg.inv(np.array([a, b]).T)
        rows.append((name, a, b, inv, first_cand, len(cands) - first_cand))

    out = []
    w = out.append
    w("/* The tilings randsctile draws. Written by tools/randsctile-tables.py, which")
    w(" * builds and checks them; change that and run it rather than editing this.")
    w(" *")
    w(" * Every periodic tiling is scaled to unit area a tile on average. For each:")
    w(" * the period, its inverse (a point to lattice coordinates), and the tiles that")
    w(" * can hold a point of one period, each a tile of the table and how many")
    w(" * periods it is moved by. A tile is its vertices counter-clockwise -- it may")
    w(" * be concave -- its middle, the radius of the largest circle it holds, and")
    w(" * its bounding box. */")
    w("")
    w("static const double RANDSCTILE_VERT[][2] = {")
    for x, y in verts:
        w(f"    {{{c_double(x)}, {c_double(y)}}},")
    w("};")
    w("")
    w("static const struct {")
    w("    unsigned short first;")
    w("    unsigned char count;")
    w("    double mx, my, reach, x0, y0, x1, y1;")
    w("} RANDSCTILE_FACE[] = {")
    for v0, n, c, far, bb in faces:
        w(f"    {{{v0}, {n}, {c_double(c[0])}, {c_double(c[1])}, {c_double(far)}, "
          f"{c_double(bb[0])}, {c_double(bb[1])}, {c_double(bb[2])}, {c_double(bb[3])}}},")
    w("};")
    w("")
    w("static const struct {")
    w("    unsigned short face;")
    w("    signed char i, j;")
    w("} RANDSCTILE_CAND[] = {")
    for fi, i, j in cands:
        w(f"    {{{fi}, {i}, {j}}},")
    w("};")
    w("")
    w("static const struct {")
    w("    double ax, ay, bx, by; /* the period */")
    w("    double ia, ib, ja, jb; /* lattice coordinates: ia*x + ib*y, ja*x + jb*y */")
    w("    unsigned short first, count; /* into RANDSCTILE_CAND */")
    w("} RANDSCTILE_PERIODIC[] = {")
    for name, a, b, inv, c0, cn in rows:
        w(f"    /* {name} */")
        w(f"    {{{c_double(a[0])}, {c_double(a[1])}, {c_double(b[0])}, {c_double(b[1])},")
        w(f"     {c_double(inv[0, 0])}, {c_double(inv[0, 1])}, {c_double(inv[1, 0])}, {c_double(inv[1, 1])},")
        w(f"     {c0}, {cn}}},")
    w("};")
    w(f"#define RANDSCTILE_NPERIODIC {len(rows)}")
    w("")
    w("/* The multigrids: n families of lines, the unit vector across each and the")
    w(" * offset of each, and the edge a rhombus has at unit area on average. */")
    w("#define RANDSCTILE_GRID_MAX 7")
    w("static const struct {")
    w("    int n;")
    w("    double e[RANDSCTILE_GRID_MAX][2];")
    w("    double g[RANDSCTILE_GRID_MAX];")
    w("    double edge;")
    w("} RANDSCTILE_GRID[] = {")
    for n, g in GRIDS:
        E = grid_vectors(n)
        es = ", ".join(f"{{{c_double(x)}, {c_double(y)}}}" for x, y in E)
        gs = ", ".join(c_double(x) for x in g)
        w(f"    {{{n}, {{{es}}},")
        w(f"     {{{gs}}}, {c_double(1 / math.sqrt(grid_mean_area(n)))}}},")
    w("};")
    w("")
    w("/* The edge of a Robinson triangle's legs at unit area a tile on average,")
    w(" * kites and darts, and in legs of one the radius of the largest circle a")
    w(" * kite holds and a dart holds. */")
    w(f"#define RANDSCTILE_KITES_EDGE {c_double(1 / math.sqrt(kites_mean_area()))}")
    kite, dart = kite_and_dart()
    for name, q in (("KITE", kite), ("DART", dart)):
        w(f"#define RANDSCTILE_{name}_REACH {c_double(reach(ccw(q)))}")
    with open(target, "w", newline="\n") as f:
        f.write("\n".join(out) + "\n")
    print(f"{len(rows)} periodic tilings, {len(faces)} tiles, {len(verts)} vertices, "
          f"{len(cands)} candidates -> {os.path.normpath(target)}", file=sys.stderr)


if __name__ == "__main__":
    main()
