#ifndef HAMERS_CONFIG_HPP
#define HAMERS_CONFIG_HPP

#ifdef _OPENMP
/* Enable SIMD */
/* #undef HAMERS_ENABLE_SIMD */
#endif

/* Enable assertion checking */
#define HAMERS_DEBUG_CHECK_ASSERTIONS

/* Enable HAMeRS developer assertion checking */
#define HAMERS_DEBUG_CHECK_DEV_ASSERTIONS

#endif /* HAMERS_CONFIG_HPP */
