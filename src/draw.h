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

  inline void
  drawRead(cv::Mat& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, bool const reverse) {
    cv::Rect rect(x, y, w, h);
    cv::rectangle(img, rect, cv::Scalar(200, 200, 200), -1);
    typedef std::vector<cv::Point> TPointVector;
    TPointVector pvec;
    if (reverse) {
      cv::line(img, cv::Point(x, y), cv::Point(x-2, y + (int)(h/2)), cv::Scalar(150, 150, 150), 1);
      cv::line(img, cv::Point(x, y+h), cv::Point(x-2, y + (int)(h/2)), cv::Scalar(150, 150, 150), 1);
    } else {
      cv::line(img, cv::Point(x+w, y), cv::Point(x+w+2, y + (int)(h/2)), cv::Scalar(150, 150, 150), 1);
      cv::line(img, cv::Point(x+w, y+h), cv::Point(x+w+2, y + (int)(h/2)), cv::Scalar(150, 150, 150), 1);
    }
  }

  template<typename TConfig>
  inline void
    drawRead(TConfig const& c, Region const& rg, cv::Mat& img, int32_t const track, int32_t const gstart, int32_t const gend, bool const reverse) {
    int32_t px = pixelX(c.width, rg.size, gstart);
    int32_t pxend = pixelX(c.width, rg.size, gend);
    drawRead(img, px, track * c.tlheight, pxend - px, c.rdheight, reverse);
  }

  inline void
  drawNuc(cv::Mat& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, char const nuc) {
    cv::Rect rect(x, y, w, h);
    if ((nuc == 'a') or (nuc == 'A')) {
      cv::rectangle(img, rect, cv::Scalar(0, 255, 0), -1);
    }
    else if ((nuc == 'c') or (nuc == 'C')) {
      cv::rectangle(img, rect, cv::Scalar(0, 0, 255), -1);
    }
    else if ((nuc == 'g') or (nuc == 'G')) {
      cv::rectangle(img, rect, cv::Scalar(0, 0, 0), -1);
    }
    else if ((nuc == 't') or (nuc == 'T')) {
      cv::rectangle(img, rect, cv::Scalar(255, 0, 0), -1);
    } else {
      cv::rectangle(img, rect, cv::Scalar(0, 0, 0), -1);
    }
  }
  
  template<typename TConfig>
  inline void
  drawNuc(TConfig const& c, Region const& rg, cv::Mat& img, int32_t const track, int32_t const gstart, int32_t const gend, char const nuc) {
    int32_t px = pixelX(c.width, rg.size, gstart);
    int32_t pxend = pixelX(c.width, rg.size, gend);
    drawNuc(img, px, track * c.tlheight, pxend - px, c.rdheight, nuc);
  }

}

#endif
