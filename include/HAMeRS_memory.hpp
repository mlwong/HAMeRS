#ifndef HAMERS_MEMORY_HPP
#define HAMERS_MEMORY_HPP

#include "HAMeRS_config.hpp"

#ifdef USE_BOOST
#include "boost/enable_shared_from_this.hpp"
#include "boost/make_shared.hpp"
#include "boost/shared_ptr.hpp"
#include "boost/weak_ptr.hpp"
#else /* USE_BOOST */
#include <memory>
#endif /* USE_BOOST */

#ifdef USE_BOOST
#define HAMERS_SHARED_PTR boost::shared_ptr
#define HAMERS_WEAK_PTR boost::weak_ptr
#define HAMERS_MAKE_SHARED boost::make_shared
#define HAMERS_ENABLE_SHARED_FROM_THIS boost::enable_shared_from_this
#define HAMERS_DYNAMIC_POINTER_CAST boost::dynamic_pointer_cast
#define HAMERS_SHARED_PTR_CAST BOOST_CAST
#else /* USE_BOOST */
#define HAMERS_SHARED_PTR std::shared_ptr
#define HAMERS_WEAK_PTR std::weak_ptr
#define HAMERS_MAKE_SHARED std::make_shared
#define HAMERS_ENABLE_SHARED_FROM_THIS std::enable_shared_from_this
#define HAMERS_DYNAMIC_POINTER_CAST std::dynamic_pointer_cast
#define HAMERS_SHARED_PTR_CAST SAMRAI_SHARED_PTR_CAST
#endif /* USE_BOOST */

#endif /* HAMERS_MEMORY_HPP */
