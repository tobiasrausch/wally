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
    drawHilbert(TConfig const& c, Region const& rg, cv::Mat& img, TCoverage const& covA, TCoverage const& covC, TCoverage const& covG, TCoverage const& covT, TCoverage const& del, std::vector<bool> const& snp) {
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
	cv::Rect rect(x, y, 1, 1);
	cv::rectangle(img, rect, cv::Scalar(greyval, greyval, greyval), -1);
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
	  int32_t newR = fracA * WALLY_A.val[2] + fracC * WALLY_C.val[2] + fracG * WALLY_G.val[2] + fracT * WALLY_T.val[2] + fracDel * WALLY_INDEL.val[2];
	  int32_t newG = fracA * WALLY_A.val[1] + fracC * WALLY_C.val[1] + fracG * WALLY_G.val[1] + fracT * WALLY_T.val[1] + fracDel * WALLY_INDEL.val[1];
	  int32_t newB = fracA * WALLY_A.val[0] + fracC * WALLY_C.val[0] + fracG * WALLY_G.val[0] + fracT * WALLY_T.val[0] + fracDel * WALLY_INDEL.val[0];
	  int32_t pos = (int32_t) (((double) i / (double) rg.size) * (c.width * c.height));
	  int32_t x = 0;
	  int32_t y = 0;
	  posToHilbert(c.width, pos, x, y);
	  cv::circle(img, cv::Point(x, y), 5, cv::Scalar(newB, newG, newR), -1);
	}
      }
    }

}

#endif
