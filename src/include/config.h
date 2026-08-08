#ifndef CONFIG_H
#define CONFIG_H

// XaoS release
#define XaoS_VERSION "4.3.3_MS_HACK5"

// URLs
#define HELP_URL "https://github.com/xaos-project/XaoS/wiki"
#define WEB_URL "http://xaos.sourceforge.net/"
#define DOWNLOAD_URL "https://github.com/xaos-project/XaoS/releases"
#define FEEDBACK_URL "https://github.com/xaos-project/XaoS/issues"
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
#ifdef USE_FLOAT128
typedef __float128 number_t;
#define NUMBER_DIGITS 36 /* 113 bits of mantissa */
#else
#ifdef USE_LONG_DOUBLE
typedef long double number_t;
#define NUMBER_DIGITS 21 /* 64 bits on x87 */
#else
typedef double number_t;
#define NUMBER_DIGITS 17 /* 53 bits */
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
#define USER_FORMULA "(abs(re(z))+i*abs(im(z)))^2+c"

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

// Number of previous z minus 1
#define NUM_P 6
