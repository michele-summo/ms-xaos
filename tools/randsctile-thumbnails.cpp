/* Draws the thumbnails in the Tilings tab of the formula reference.
 *
 *     cmake --build <build directory> --target randsctile-thumbnails
 *
 * writes into src/ui/images/tilings a picture of every tiling randsctile
 * takes, drawn by randsctile itself through the parser, together with the
 * resource file that puts them in the binary and the fingerprints that
 * formula_help_test checks them against. Run it after adding a tiling or
 * changing one: until then that test fails and says which.
 *
 * The pictures are drawn here once rather than each time the tab is opened,
 * since they never change between one opening and the next. What they show
 * is in randsctile-thumbnails.h, which the test reads too.
 */

#include <cmath>
#include <cstdio>
#include <vector>

#include <QDir>
#include <QFile>
#include <QImage>
#include <QTextStream>

#include "config.h"
#include "number_math.h"
#include "sffe.h"
#include "sffe_cmplx_gsl.h"
#include "misc-f.h"

#include "randsctile-thumbnails.h"

const char *qt_gettext(const char * /*context*/, const char *text)
{
    return text;
}

/* A tile is coloured from its value along the ramp the tilings were first
 * drawn with -- ink, vermilion, yellow, paper, blue, navy, sage -- and washed
 * a third of the way towards the paper, so that the outlines over it read. */
static const int RAMP[][3] = {{20, 20, 30},    {200, 60, 40},  {240, 200, 90},
                              {230, 225, 210}, {60, 120, 160}, {25, 40, 70},
                              {150, 170, 120}, {20, 20, 30}};
static const int RAMP_STOPS = sizeof(RAMP) / sizeof(RAMP[0]);
static const double PAPER[3] = {250, 248, 242};
static const double WASH = 0.35;
static const double INK[3] = {40, 40, 48};
static const double FRAME[3] = {150, 150, 150};

/* Outlines are drawn where two neighbouring samples differ, then widened by
 * this many samples: six samples across, a pixel and a half as stored, which
 * is one pixel on a screen at 125%. */
static const int OUTLINE_REACH = 2;
/* The frame round each picture, in stored pixels. */
static const int FRAME_WIDTH = 2;

static void tile_colour(double v, double rgb[3])
{
    double t = v * (RAMP_STOPS - 1);
    int i = (int)std::floor(t);
    if (i < 0)
        i = 0;
    if (i > RAMP_STOPS - 2)
        i = RAMP_STOPS - 2;
    double f = t - i;
    for (int c = 0; c < 3; c++) {
        double ramp = RAMP[i][c] * (1 - f) + RAMP[i + 1][c] * f;
        rgb[c] = ramp * (1 - WASH) + PAPER[c] * WASH;
    }
}

static number_t value_at(sffe *p, double x, double y)
{
    GSL_SET_COMPLEX(&sffe_position, (number_t)x, (number_t)y);
    sffe_iteration = 0;
    cmplx v = sffe_eval(p);
    return GSL_REAL(v);
}

static QImage draw(sffe *p)
{
    const int fine = RANDSCTILE_THUMB_FINE;
    std::vector<double> v((size_t)fine * fine);
    for (int r = 0; r < fine; r++)
        for (int c = 0; c < fine; c++) {
            double x, y;
            randsctile_thumb_point(r, c, &x, &y);
            v[(size_t)r * fine + c] = (double)value_at(p, x, y);
        }

    /* A sample is on an edge when any of its four neighbours is in another
     * tile, which marks both sides and so centres the line on the edge. */
    std::vector<char> edge((size_t)fine * fine, 0);
    for (int r = 0; r < fine; r++)
        for (int c = 0; c < fine; c++) {
            double here = v[(size_t)r * fine + c];
            if ((c > 0 && v[(size_t)r * fine + c - 1] != here) ||
                (c + 1 < fine && v[(size_t)r * fine + c + 1] != here) ||
                (r > 0 && v[(size_t)(r - 1) * fine + c] != here) ||
                (r + 1 < fine && v[(size_t)(r + 1) * fine + c] != here))
                edge[(size_t)r * fine + c] = 1;
        }
    /* Widened by a disc rather than a square, so that a slanting line comes
     * out as wide as an upright one. */
    std::vector<char> ink((size_t)fine * fine, 0);
    for (int r = 0; r < fine; r++)
        for (int c = 0; c < fine; c++) {
            if (!edge[(size_t)r * fine + c])
                continue;
            for (int dr = -OUTLINE_REACH; dr <= OUTLINE_REACH; dr++)
                for (int dc = -OUTLINE_REACH; dc <= OUTLINE_REACH; dc++) {
                    int rr = r + dr, cc = c + dc;
                    if (dr * dr + dc * dc > OUTLINE_REACH * OUTLINE_REACH ||
                        rr < 0 || rr >= fine || cc < 0 || cc >= fine)
                        continue;
                    ink[(size_t)rr * fine + cc] = 1;
                }
        }

    const int side = RANDSCTILE_THUMB_SIDE, n = RANDSCTILE_THUMB_SAMPLES;
    QImage image(side, side, QImage::Format_RGB32);
    for (int py = 0; py < side; py++)
        for (int px = 0; px < side; px++) {
            double sum[3] = {0, 0, 0};
            for (int sy = 0; sy < n; sy++)
                for (int sx = 0; sx < n; sx++) {
                    size_t at = (size_t)(py * n + sy) * fine + (px * n + sx);
                    double rgb[3];
                    if (ink[at])
                        for (int c = 0; c < 3; c++)
                            rgb[c] = INK[c];
                    else
                        tile_colour(v[at], rgb);
                    for (int c = 0; c < 3; c++)
                        sum[c] += rgb[c];
                }
            bool frame = px < FRAME_WIDTH || py < FRAME_WIDTH ||
                         px >= side - FRAME_WIDTH || py >= side - FRAME_WIDTH;
            int out[3];
            for (int c = 0; c < 3; c++)
                out[c] = (int)std::lround(frame ? FRAME[c] : sum[c] / (n * n));
            image.setPixel(px, py, qRgb(out[0], out[1], out[2]));
        }
    return image;
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "usage: %s <src/ui/images/tilings>\n", argv[0]);
        return 2;
    }
    QDir dir(QString::fromLocal8Bit(argv[1]));
    if (!dir.mkpath(".")) {
        fprintf(stderr, "cannot make %s\n", argv[1]);
        return 1;
    }

    /* Every tiling there is: the first number that draws nothing is where
     * they stop, as it is for the test. */
    QString qrc, prints;
    QTextStream q(&qrc), f(&prints);
    q << "<!-- Written by tools/randsctile-thumbnails.cpp: do not edit. -->\n"
      << "<RCC>\n    <qresource prefix=\"/images/tilings\">\n";
    f << "# Written by tools/randsctile-thumbnails.cpp: do not edit.\n"
      << "# Each tiling, then the fingerprint of what randsctile draws over its\n"
      << "# thumbnail -- see randsctile-thumbnails.h.\n";
    int kinds = 0;
    for (int k = 1;; k++) {
        char call[64];
        snprintf(call, sizeof call, RANDSCTILE_THUMB_CALL, k);
        sffe *p = sffe_alloc();
        if (sffe_parse(&p, call)) {
            fprintf(stderr, "cannot parse %s\n", call);
            return 1;
        }
        if (value_at(p, 0.3, 0.7) == 0) {
            sffe_free(&p);
            break;
        }
        QString name = QString("tiling-%1.png").arg(k, 2, 10, QChar('0'));
        if (!draw(p).save(dir.filePath(name), "PNG", 9)) {
            fprintf(stderr, "cannot write %s\n", qPrintable(dir.filePath(name)));
            return 1;
        }
        uint64_t print = randsctile_thumb_fingerprint(
            [p](double x, double y) { return (double)value_at(p, x, y); });
        q << "        <file>" << name << "</file>\n";
        f << k << " " << QString::number((qulonglong)print, 16).rightJustified(16, '0')
          << "\n";
        printf("%s\n", call);
        fflush(stdout);
        sffe_free(&p);
        kinds = k;
    }
    q << "    </qresource>\n</RCC>\n";

    struct {
        const char *name;
        const QString *text;
    } files[] = {{"tilings.qrc", &qrc}, {"fingerprints.txt", &prints}};
    for (auto &out : files) {
        QFile file(dir.filePath(out.name));
        if (!file.open(QIODevice::WriteOnly) ||
            file.write(out.text->toUtf8()) != out.text->toUtf8().size()) {
            fprintf(stderr, "cannot write %s\n", out.name);
            return 1;
        }
    }
    printf("%d tilings drawn into %s\n", kinds, argv[1]);
    return 0;
}
