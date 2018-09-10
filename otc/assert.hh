#define BOOST_ENABLE_ASSERT_HANDLER
#include <boost/assert.hpp>
#include "otc/error.h"
#ifdef assert
#undef assert
#endif
#ifdef assert_msg
#undef assert_msg
#endif
#ifdef NDEBUG
#define assert(expr)          (0?(void(expr)):(void(0)))
#define assert_msg(expr,message)  (0?(void(expr)):(void(0)))
#else
#define assert(expr)  BOOST_ASSERT(expr)
#define assert_msg(expr,message)  BOOST_ASSERT_MSG(expr,message)
#endif

#ifndef OTC_ASSERT_H
#define OTC_ASSERT_H
namespace boost
{
    inline void assertion_failed(char const * expr, char const * function, char const * file, long line)
    {
	throw ::otc::OTCError()<<"Assertion ("<<expr<<") failed in '"<<function<<"' at "<<file<<":"<<line;
    }

    inline void assertion_failed_msg(char const * expr, char const * msg, char const * function, char const * file, long line)
    {
	throw ::otc::OTCError()<<"Assertion ("<<expr<<") failed in '"<<function<<"' at "<<file<<":"<<line<<":\n   "<<msg;
    }
}
#endif
