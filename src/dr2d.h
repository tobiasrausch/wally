#ifndef DR2D_H
#define DR2D_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/dynamic_bitset.hpp>

#include <htslib/sam.h>

#include "util.h"
#include "bed.h"
#include "draw.h"

namespace wallysworld
{

 

  
  template<typename TConfig>
  inline void
  drawMatch(TConfig const& c, cv::Mat& img, int32_t const x, int32_t const xend, int32_t const y, int32_t const yend, bool const forward) {
    if (forward) {
      cv::line(img, cv::Point(x, y), cv::Point(xend, yend), cv::Scalar(0, 0, 255), 1);
    } else {
      cv::line(img, cv::Point(x, yend), cv::Point(xend, y), cv::Scalar(255, 0, 0), 1);
    }
  }

  template<typename TConfig>
  inline void
  drawMatch(TConfig const& c, Region const& rg1, Region const& rg2, cv::Mat& img, int32_t const gstart, int32_t const gend, int32_t const rstart, int32_t const rend, bool const forward) {
    int32_t px = pixelX(c.width, rg1.size, gstart);
    int32_t pxend = pixelX(c.width, rg1.size, gend);
    int32_t py = pixelX(c.height, rg2.size, rstart);
    int32_t pyend = pixelX(c.height, rg2.size, rend);
    drawMatch(c, img, px, pxend, py, pyend, forward);
  }

  
}

#endif
