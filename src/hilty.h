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


    template<typename TConfig, typename TCoverage>
    inline void
    drawHilbert(TConfig const& c, Region const& rg, BLContext& img, TCoverage const& covA, TCoverage const& covC, TCoverage const& covG, TCoverage const& covT, TCoverage const& del, std::vector<bool> const& snp) {
      uint32_t maxObsCov = 0;
      for(uint32_t i = 0; i < covA.size(); ++i) {
	uint32_t cumsum = covA[i];
	cumsum += covC[i];
	cumsum += covG[i];
	cumsum += covT[i];
	cumsum += del[i];
	if (cumsum > maxObsCov) maxObsCov = cumsum;
      }

      // Draw relative coverage in grey
      for(int32_t i = 0; i < rg.size; ++i) {
	uint32_t cumsum = covA[i];
	cumsum += covC[i];
	cumsum += covG[i];
	cumsum += covT[i];
	cumsum += del[i];
	int32_t pos = (int32_t) (((double) i / (double) rg.size) * (c.width * c.height));
	double frac = (double) cumsum / (double) maxObsCov;
	int32_t greyval = frac * 255;
	int32_t x = 0;
	int32_t y = 0;
	posToHilbert(c.width, pos, x, y);
	// Just the coverage in grey
	img.fill_rect(BLRectI(x, y, 1, 1), BLRgba32(greyval, greyval, greyval));
      }

      // Overdraw SNPs and deletions
      for(int32_t i = 0; i < rg.size; ++i) {
	if (snp[i]) {
	  // Normalize SNPs to the site coverage
	  uint32_t cumsum = covA[i];
	  cumsum += covC[i];
	  cumsum += covG[i];
	  cumsum += covT[i];
	  cumsum += del[i];
	  double fracA = (double) covA[i] / (double) maxObsCov;
	  double fracC = (double) covC[i] / (double) maxObsCov;
	  double fracG = (double) covG[i] / (double) maxObsCov;
	  double fracT = (double) covT[i] / (double) maxObsCov;
	  double fracDel = (double) del[i] / (double) maxObsCov;
	  int32_t newR = fracA * WALLY_A.r() + fracC * WALLY_C.r() + fracG * WALLY_G.r() + fracT * WALLY_T.r() + fracDel * WALLY_INDEL.r();
	  int32_t newG = fracA * WALLY_A.g() + fracC * WALLY_C.g() + fracG * WALLY_G.g() + fracT * WALLY_T.g() + fracDel * WALLY_INDEL.g();
	  int32_t newB = fracA * WALLY_A.b() + fracC * WALLY_C.b() + fracG * WALLY_G.b() + fracT * WALLY_T.b() + fracDel * WALLY_INDEL.b();
	  int32_t pos = (int32_t) (((double) i / (double) rg.size) * (c.width * c.height));
	  int32_t x = 0;
	  int32_t y = 0;
	  posToHilbert(c.width, pos, x, y);
	  img.fill_circle(x, y, 5, BLRgba32(newR, newG, newB));
	}
      }
    }

}

#endif
