#ifndef HILTY_H
#define HILTY_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/math/distributions/binomial.hpp>

#include <htslib/sam.h>


namespace wallysworld
{

    inline void
    rot(int32_t const n, int32_t& x, int32_t& y, int32_t const rx, int32_t const ry)  {
      if ( ry == 0 ) {
	if ( rx == 1 ) {
	  x = n - 1 - x;
	  y = n - 1 - y;
	}
	int32_t t = x;
	x = y;
	y = t;
      }
    }
  
  
    inline int32_t
    hilbertToPos(int32_t const n, int32_t const x, int32_t const y) {
      int32_t d = 0;
      int32_t inx = x;
      int32_t iny = y;
      for (int32_t s = n / 2; s > 0; s /= 2) {
	int32_t rx = ( inx & s ) > 0;
	int32_t ry = ( iny & s ) > 0;
	d += s * s * ( ( 3 * rx ) ^ ry );
	rot(s, inx, iny, rx, ry);
      }
      return d;
    }

    inline void
    posToHilbert(int32_t const n, int32_t pos, int32_t& x, int32_t& y) {
      int32_t t = pos;
      x = 0;
      y = 0;
      for (int32_t s = 1; s < n; s *= 2) {
	int32_t rx = ( 1 & ( t / 2 ) );
	int32_t ry = ( 1 & ( t ^ rx ) );
	rot(s, x, y, rx, ry);
	x += s * rx;
	y += s * ry;
	t /= 4;
      }
    }
    
}

#endif
