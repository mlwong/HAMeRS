#ifndef HAMERS_CONFIG_HPP
#define HAMERS_CONFIG_HPP

#include "SAMRAI/SAMRAI_config.h"

/* Enable SIMD */
/* #undef HAMERS_ENABLE_SIMD */

/* Enable assertion checking */
#define HAMERS_DEBUG_CHECK_ASSERTIONS

/* Enable HAMeRS developer assertion checking */
#define HAMERS_DEBUG_CHECK_DEV_ASSERTIONS

/* Define epsilon to prevent divisoin by zero */
#define HAMERS_EPSILON 1.0e-40

#define USE_BOOST

#endif /* HAMERS_CONFIG_HPP */
