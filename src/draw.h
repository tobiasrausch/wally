#ifndef DRAW_H
#define DRAW_H

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

  

  inline int32_t
    pixelX(int32_t const width, int32_t const rgsz, int32_t const genomicX) {
    return (int32_t) (((double) genomicX / (double) rgsz) * width);
  }

  inline int32_t
    genomicX(int32_t const width, int32_t const rgsz, int32_t const pixelX) {
    return (int32_t) (((double) pixelX / (double) width) * rgsz);
  }

  template<typename TConfig>
  inline void
  drawGenome(TConfig const& c, Region const& rg, cv::Mat& img, int32_t const track) {
    std::string text(boost::lexical_cast<std::string>(rg.end));
    int32_t n = text.length() - 3;
    while (n > 0) {
      text.insert(n, ",");
      n -= 3;
    }
    double font_scale = 0.4;
    double font_thickness = 1.5;
    int32_t baseline = 0;
    cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);

    // Find suitable tick size
    uint32_t modval = 0;
    for(uint32_t i = 10; i < 1000000; i *= 10) {
      if (i * c.pxoffset > textSize.width) {
	modval = i;
	break;
      }
    }
    cv::line(img, cv::Point(0, track * c.tlheight + c.tlheight), cv::Point(c.width, track * c.tlheight + c.tlheight), cv::Scalar(0, 0, 0), 1.8);
    double px = 0;
    for(int32_t i = rg.beg; i < rg.end; ++i) {
      if (i % (modval / 2) == 0) {
	cv::line(img, cv::Point(px, track * c.tlheight), cv::Point(px, track * c.tlheight + c.tlheight), cv::Scalar(0, 0, 0), 1.8);
      }
      if (i % modval == 0) {
	// Font
	text = boost::lexical_cast<std::string>(i);
	n = text.length() - 3;
	while (n > 0) {
	  text.insert(n, ",");
	  n -= 3;
	}
	baseline = 0;
	cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);
	cv::putText(img, text, cv::Point(px - textSize.width/2, (track - 1) * c.tlheight + textSize.height), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness);
      }
      px += c.pxoffset;
    }
  }

  
  template<typename TConfig>
  inline void
  drawRead(TConfig const& c, cv::Mat& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, bool const reverse, bool const tri) {
    cv::Rect rect(x, y, w, h);
    cv::rectangle(img, rect, cv::Scalar(200, 200, 200), -1);
    if (tri) {
      typedef std::vector<cv::Point> TPointVector;
      TPointVector pvec;
      if (reverse) {
	std::vector<cv::Point> pvec{cv::Point(x, y), cv::Point(x, y+h), cv::Point(x-c.pxoffset, y + h/2)};
	cv::polylines(img, pvec, true, cv::Scalar(200, 200, 200), 1);
	cv::fillPoly(img, pvec, cv::Scalar(200, 200, 200));
      } else {
	std::vector<cv::Point> pvec{cv::Point(x+w, y), cv::Point(x+w, y+h), cv::Point(x+w+c.pxoffset, y + h/2)};
	cv::polylines(img, pvec, true, cv::Scalar(200, 200, 200), 1);
	cv::fillPoly(img, pvec, cv::Scalar(200, 200, 200));
      }
    }
  }

  template<typename TConfig>
  inline void
  drawRead(TConfig const& c, Region const& rg, cv::Mat& img, int32_t const track, int32_t const gstart, int32_t const gend, bool const reverse, bool const tri) {
    int32_t px = pixelX(c.width, rg.size, gstart);
    int32_t pxend = pixelX(c.width, rg.size, gend);
    drawRead(c, img, px, track * c.tlheight, pxend - px, c.rdheight, reverse, tri);
  }

  inline void
  drawNuc(cv::Mat& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, char const nuc) {
    // Background rectangle
    cv::Rect rect(x, y, w, h);

    // Font
    std::string text(1, nuc);
    double font_scale = 0.4;
    double font_thickness = 1.5;
    int32_t baseline = 0;
    cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);
    double frac = (double) textSize.width / (double) w;

    // Put nucleotide if there is space
    if (frac < 1) {
      cv::rectangle(img, rect, cv::Scalar(200, 200, 200), -1);
      if ((nuc == 'a') or (nuc == 'A')) {
	cv::putText(img, text, cv::Point(x + w/2 - textSize.width/2, y+h/2+textSize.height/2), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(0, 255, 0), font_thickness);
      }
      else if ((nuc == 'c') or (nuc == 'C')) {
      	cv::putText(img, text, cv::Point(x + w/2 - textSize.width/2, y+h/2+textSize.height/2), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(255, 0, 0), font_thickness);
      }
      else if ((nuc == 'g') or (nuc == 'G')) {
	cv::putText(img, text, cv::Point(x + w/2 - textSize.width/2, y+h/2+textSize.height/2), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(5, 113, 209), font_thickness);
      }
      else if ((nuc == 't') or (nuc == 'T')) {
	cv::putText(img, text, cv::Point(x + w/2 - textSize.width/2, y+h/2+textSize.height/2), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(0, 0, 255), font_thickness);
      }
    } else { 
      if ((nuc == 'a') or (nuc == 'A')) {
	cv::rectangle(img, rect, cv::Scalar(0, 255, 0), -1);
      }
      else if ((nuc == 'c') or (nuc == 'C')) {
	cv::rectangle(img, rect, cv::Scalar(255, 0, 0), -1);
      }
      else if ((nuc == 'g') or (nuc == 'G')) {
	cv::rectangle(img, rect, cv::Scalar(5, 113, 209), -1);
      }
      else if ((nuc == 't') or (nuc == 'T')) {
	cv::rectangle(img, rect, cv::Scalar(0, 0, 255), -1);
      }
    }
  }
  
  template<typename TConfig>
  inline void
  drawNuc(TConfig const& c, Region const& rg, cv::Mat& img, int32_t const track, int32_t const gstart, int32_t const gend, char const nuc) {
    int32_t px = pixelX(c.width, rg.size, gstart);
    int32_t pxend = pixelX(c.width, rg.size, gend);
    drawNuc(img, px, track * c.tlheight, pxend - px, c.rdheight, nuc);
  }

  inline void
  drawDel(cv::Mat& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, int32_t const len) {
    std::string text = boost::lexical_cast<std::string>(len);
    double font_scale = 0.4;
    double font_thickness = 1.5;
    int32_t baseline = 0;
    cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_DUPLEX, font_scale, font_thickness, &baseline);
    cv::line(img, cv::Point(x, y+h/2), cv::Point(x+w, y+h/2), cv::Scalar(0, 0, 0), 1.8);
    double frac = (double) textSize.width / (double) w;
    // Put length if there is space
    if (frac < 0.5) {
      cv::Rect rect(x+w/2-textSize.width/2, y, textSize.width, h);
      cv::rectangle(img, rect, cv::Scalar(255, 255, 255), -1);
      cv::putText(img, text, cv::Point(x+w/2 - textSize.width/2, y + h/2 + textSize.height/2), cv::FONT_HERSHEY_SIMPLEX, font_scale, cv::Scalar(211, 0, 148), font_thickness);
    }
  }

  template<typename TConfig>
  inline void
    drawDel(TConfig const& c, Region const& rg, cv::Mat& img, int32_t const track, int32_t const gstart, int32_t const gend, int32_t const len) {
    int32_t px = pixelX(c.width, rg.size, gstart);
    int32_t pxend = pixelX(c.width, rg.size, gend);
    drawDel(img, px, track * c.tlheight, pxend - px, c.rdheight, len);
  }
  
}

#endif
