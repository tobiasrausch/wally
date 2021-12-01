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
    drawHilbert(TConfig const& c, Region const& rg, cv::Mat& img, TCoverage const& covA, TCoverage const& covC, TCoverage const& covG, TCoverage const& covT, std::vector<bool> const& snp) {
      uint32_t maxObsCov = 0;
      for(uint32_t i = 0; i < covA.size(); ++i) {
	uint32_t cumsum = covA[i];
	cumsum += covC[i];
	cumsum += covG[i];
	cumsum += covT[i];
	if (cumsum > maxObsCov) maxObsCov = cumsum;
      }

      // Draw relative coverage in grey
      for(int32_t i = 0; i < rg.size; ++i) {
	uint32_t cumsum = covA[i];
	cumsum += covC[i];
	cumsum += covG[i];
	cumsum += covT[i];
	int32_t pos = (int32_t) (((double) i / (double) rg.size) * (c.width * c.height));
	double frac = (double) cumsum / (double) maxObsCov;
	int32_t greyval = frac * 255;
	std::cout << pos << ',' << greyval << std::endl;
	int32_t x = 0;
	int32_t y = 0;
	posToHilbert(c.width, pos, x, y);
	cv::Rect rect(x, y, 1, 1);
	cv::rectangle(img, rect, cv::Scalar(greyval, greyval, greyval), -1);
      }

      /*
      // Overdraw SNPs
      px = 0;
      for(int32_t i = 0; i < rg.size; ++i) {
      if (snp[i]) {
      int32_t pxs = px;
      int32_t pxw = c.pxoffset;
	if (pxw < 1) pxw = 1;
	if (c.pxoffset >= WALLY_PX) {
	  pxs = px + 1;
	  pxw = c.pxoffset - 1;
	}
	double frac = (double) covA[i] / (double) maxObsCov;
	int32_t ch = (int32_t) (frac * 2 * c.tlheight);
	int32_t cumCH = ch;
	cv::Rect rectA(pxs, (track-1) * c.tlheight + 2 * c.tlheight - cumCH, pxw, ch);
	cv::rectangle(img, rectA, WALLY_A, -1);
	frac = (double) covC[i] / (double) maxObsCov;
	ch = (int32_t) (frac * 2 * c.tlheight);
	cumCH += ch;
	cv::Rect rectC(pxs, (track-1) * c.tlheight + 2 * c.tlheight - cumCH, pxw, ch);
	cv::rectangle(img, rectC, WALLY_C, -1);
	frac = (double) covG[i] / (double) maxObsCov;
	ch = (int32_t) (frac * 2 * c.tlheight);
	cumCH += ch;
	cv::Rect rectG(pxs, (track-1) * c.tlheight + 2 * c.tlheight - cumCH, pxw, ch);
	cv::rectangle(img, rectG, WALLY_G, -1);
	frac = (double) covT[i] / (double) maxObsCov;
	ch = (int32_t) (frac * 2 * c.tlheight);
	cumCH += ch;
	cv::Rect rectT(pxs, (track-1) * c.tlheight + 2 * c.tlheight - cumCH, pxw, ch);
	cv::rectangle(img, rectT, WALLY_T, -1);
      }
      px += c.pxoffset;
    }

    // Put coverage scale
    if (true) {
      std::string text = "[0-"; 
      text += boost::lexical_cast<std::string>(maxObsCov);
      text += "]";
      double font_scale = 0.4;
      double font_thickness = 1.5;
      int32_t baseline = 0;
      cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_DUPLEX, font_scale, font_thickness, &baseline);
      cv::Rect rect(0,  track * c.tlheight - textSize.height/2, textSize.width, textSize.height);
      cv::rectangle(img, rect, cv::Scalar(0, 255, 255), -1);
      cv::putText(img, text, cv::Point(0, track * c.tlheight + textSize.height/2), cv::FONT_HERSHEY_SIMPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness);
    }
      */
    }

    
}

#endif
