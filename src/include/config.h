#ifndef CONFIG_H
#define CONFIG_H

// XaoS release.
//
// The two binaries are indistinguishable once they are running -- same title
// bar, same version in the welcome message and the About box -- so the quad
// one says which it is. The name is only what gets displayed; the application
// name Qt is given stays "XaoS" for both, since that is also what decides
// where the settings are stored and the two should go on sharing those.
#define XaoS_VERSION_BASE "1.7"
#ifdef USE_FLOAT128
#define XaoS_NAME "MS XaoS Quad"
#define XaoS_VERSION XaoS_VERSION_BASE "_Quad"
#else
#define XaoS_NAME "MS XaoS"
#define XaoS_VERSION XaoS_VERSION_BASE
#endif

// URLs
//
// Help and Feedback point at this fork, since what they are about -- the
// noise functions, the bailout shapes, the selection zoom, the quad binary --
// is not in the upstream wiki and reporting it there would be a nuisance to
// people who cannot act on it. Fractal types, the web site and the forum are
// upstream's and still describe what they describe.
#define HELP_URL "https://github.com/michele-summo/ms-xaos/blob/master/doc/ms-xaos-guide.md"
#define WEB_URL "http://xaos.sourceforge.net/"
#define DOWNLOAD_URL "https://github.com/michele-summo/ms-xaos/releases"
#define FEEDBACK_URL "https://github.com/michele-summo/ms-xaos/issues"
#define FORUM_URL "https://groups.google.com/d/forum/xaos-users"
#define FRACTALINFO_URL "https://github.com/xaos-project/XaoS/wiki/Fractal-Types#"

// File locations
#ifndef DATAPATH
#define DATAPATH "/usr/share/XaoS"
#endif
#define TUTORIALPATH DATAPATH "/tutorial/"
#define EXAMPLESPATH DATAPATH "/examples/"
#define CATALOGSPATH DATAPATH "/catalogs/"
#define HELPPATH DATAPATH "/help/xaos.hlp"

// Config file name
#ifdef _WIN32
#define CONFIGFILE "XaoS.cfg"
#else
#define CONFIGFILE ".XaoSrc"
#endif

// Optional features
#define USE_PTHREAD

// Numeric type
//
// NUMBER_DIGITS is how many significant decimal digits have to be printed for
// a value to read back as itself. That is ceil(mantissa_bits * log10(2)) + 1,
// the DECIMAL_DIG of the type -- one more than the digits it "carries", which
// is the usual place to go wrong: the code here printed 34 for __float128 and
// 20 for long double, one short in both cases, so a saved position came back
// at a slightly different place than it was saved. Printing more than needed
// only pads with zeros, printing fewer loses the view.
//
// NUMBER_MANTISSA_BITS is what a saved position records about the build that
// wrote it, so that opening it somewhere narrower can say the picture will not
// come back exactly rather than quietly drawing a different one.
#ifdef USE_FLOAT128
typedef __float128 number_t;
#define NUMBER_DIGITS 36 /* 113 bits of mantissa */
#define NUMBER_MANTISSA_BITS 113
#else
#ifdef USE_LONG_DOUBLE
typedef long double number_t;
#define NUMBER_DIGITS 21 /* 64 bits on x87 */
#define NUMBER_MANTISSA_BITS 64
#else
typedef double number_t;
#define NUMBER_DIGITS 17 /* 53 bits */
#define NUMBER_MANTISSA_BITS 53
#endif
#endif

// Supported color depths
#define STRUECOLOR
#define STRUECOLOR16 // required for edge detection and pseudo 3d

// Fractal defaults
#define DEFAULT_MAX_ITER 100
#define DEFAULT_BAILOUT 4
#define MAXSTEP (0.008 * 3)
#define STEP (0.0006 * 3)
#define ROTATIONSPEED 30
#define FRAMERATE 20
#define SPEEDUP 1.05

// Autopilot configuration
#define LOOKSIZE 2 // size explored by autopilot
#define RANGE1 30
#define NGUESSES (RANGE1 * RANGE1 / 2)
#define MAXTIME 10     // maximum zooming time to one direction
#define NGUESSES1 10   // maximum number of guesses using first method
#define NGUESSES2 1000 // maximum number of guesses using second method

// Default user formula
// #define USER_FORMULA "z^log(c)*p"
// #define USER_FORMULA "c^z+im(p)*{0;1}"
/* What the user formula holds before anything is typed into it. The burning
 * ship was there, which is a fractal of its own and a puzzle to meet as a
 * starting point; the Mandelbrot is what a formula is expected to say. */
#define USER_FORMULA "z^2+c"

// Disable optional statistics collection and reporting
//#define STATISTICS
#undef STAT
#ifdef STATISTICS
#define STAT(x) x
#else
#define STAT(x)
#endif

#define NUMBER_BIG ((number_t)INT_MAX)
#endif // CONFIG_H

// How far back a user formula may reach for an earlier value of z: p1 is the
// value it had on the pass before, p9999 the one nine thousand nine hundred
// and ninety-nine passes back. Only the places a formula names are kept, so
// the ceiling costs nothing to raise; it is here to stop a typo like p99999
// from being taken for a variable at all.
#define NUM_P_MAX 9999
