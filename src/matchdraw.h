#ifndef MATCHDRAW_H
#define MATCHDRAW_H

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
  drawSplitLine(TConfig const& c, cv::Mat& img, uint32_t const track) {
    // Draw line
    cv::line(img, cv::Point(0, track * c.tlheight), cv::Point(c.width, track * c.tlheight), WALLY_SPLITBORDER, 1);
  }
  
  template<typename TConfig>
  inline void
  drawSplitBorder(TConfig const& c, cv::Mat& img) {
    cv::line(img, cv::Point(0, 0), cv::Point(0, c.height), WALLY_SPLITBORDER, 3);
    cv::line(img, cv::Point(c.width - 1, 0), cv::Point(c.width - 1, c.height), WALLY_SPLITBORDER, 3);
  }

  template<typename TConfig>
  inline uint32_t
  minPixelWidth(TConfig const& c, std::vector<Region> const& rg) {
    uint32_t maxsize = 0;
    for(uint32_t i = 0; i < rg.size(); ++i) {
      std::string textBeg(boost::lexical_cast<std::string>(rg[i].beg));
      insertComma(textBeg);
      std::string textEnd(boost::lexical_cast<std::string>(rg[i].end));
      insertComma(textEnd);
      std::string text = "##" + textBeg + "-" + textEnd + "##";
      double font_scale = c.ftscale;
      double font_thickness = 3 * c.ftscale;
      int32_t baseline = 0;
      cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);
      if (textSize.width > (int32_t) maxsize) maxsize = textSize.width;
    }
    return maxsize * c.splits;
  }

  template<typename TConfig>
  inline void
  drawReadLabel(TConfig const& c, int32_t const track, std::string const& text, cv::Mat& img) {
    double font_scale = c.ftscale;
    double font_thickness = 3 * c.ftscale;
    int32_t baseline = 0;
    cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_DUPLEX, font_scale, font_thickness, &baseline);
    cv::Rect rect(0,  track * c.tlheight, textSize.width, textSize.height);
    cv::rectangle(img, rect, cv::Scalar(0, 255, 255), -1);
    cv::putText(img, text, cv::Point(0, track * c.tlheight + textSize.height), cv::FONT_HERSHEY_SIMPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness, cv::LINE_AA);
  }
  
  template<typename TConfig>
  inline void
  drawCoordinates(TConfig const& c, Region const& rg, std::string const& chr, cv::Mat& img) {
    std::string textBeg(boost::lexical_cast<std::string>(rg.beg));
    insertComma(textBeg);
    std::string textEnd(boost::lexical_cast<std::string>(rg.end));
    insertComma(textEnd);
    std::string text = textBeg + "-" + textEnd;

    double font_scale = c.ftscale;
    double font_thickness = 3 * c.ftscale;
    int32_t baseline = 0;
    cv::Size textSize = cv::getTextSize(chr, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);
    cv::putText(img, chr, cv::Point(c.width / 2 - textSize.width / 2, 0 * c.tlheight + textSize.height), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness, cv::LINE_AA);
    textSize = cv::getTextSize(text, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);
    cv::putText(img, text, cv::Point(c.width / 2 - textSize.width / 2, 1 * c.tlheight + textSize.height), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness, cv::LINE_AA);

    // Draw scale
    int32_t bpscale = 0.5 * c.width * c.bpoffset;
    int32_t mult = 1;
    while (bpscale > 10) {
      bpscale /= 10;
      mult *= 10;
    }
    bpscale *= mult;
    int32_t sep = 5;
    int32_t px = sep;
    int32_t pxend = sep + c.pxoffset * bpscale;
    cv::line(img, cv::Point(px, 2 * c.tlheight + c.tlheight / 2), cv::Point(pxend, 2 * c.tlheight + c.tlheight / 2), cv::Scalar(0, 0, 0), c.lw);
    cv::line(img, cv::Point(px, 2 * c.tlheight + 1), cv::Point(px, 3 * c.tlheight - 1), cv::Scalar(0, 0, 0), c.lw);
    cv::line(img, cv::Point(pxend, 2 * c.tlheight + 1), cv::Point(pxend, 3 * c.tlheight - 1), cv::Scalar(0, 0, 0), c.lw);

    std::string scaletxt(boost::lexical_cast<std::string>(bpscale));
    insertComma(scaletxt);
    scaletxt += "bp";
    textSize = cv::getTextSize(scaletxt, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);
    cv::putText(img, scaletxt, cv::Point(pxend + sep, 2 * c.tlheight + textSize.height), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness, cv::LINE_AA);
  }

  template<typename TConfig>
  inline void
  drawCrossConnect(TConfig const& c, Region const& rg, cv::Mat& img, int32_t const track, Mapping const& prev, Mapping const& mp) {
    bool drawLine = false;
    if ((prev.tid < rg.tid) || ((prev.tid == rg.tid) && (prev.gend < (rg.beg + rg.end) / 2))) {
      if ((mp.tid > rg.tid) || ((mp.tid == rg.tid) && (mp.gstart > (rg.beg + rg.end) / 2))) drawLine = true;
    }
    if ((prev.tid > rg.tid) || ((prev.tid == rg.tid) && (prev.gstart > (rg.beg + rg.end) / 2))) {
      if ((mp.tid < rg.tid) || ((mp.tid == rg.tid) && (mp.gend < (rg.beg + rg.end) / 2))) drawLine = true;
    }
    if (drawLine) {
      int32_t y = track * c.tlheight;
      cv::line(img, cv::Point(0, y), cv::Point(c.width, y), cv::Scalar(0, 0, 0), c.lw);
    }
  }	
  
  template<typename TConfig>
  inline void
  drawBlock(TConfig const& c, Region const& rg, cv::Mat& img, int32_t const track, Mapping const& prev, Mapping const& mp, Mapping const& succ) {

    // Draw rectangle
    int32_t px = pixelX(c.width, rg.size, mp.gstart - rg.beg);
    int32_t pxend = pixelX(c.width, rg.size, mp.gend - rg.beg);
    int32_t x = px;
    int32_t w = pxend - px;
    if (w < 1) w = 1; // Make sure each block is visible
    int32_t y = track * c.tlheight;
    int32_t midpoint = y + c.tlheight/2;

    // Draw connectors
    if (mp.fwd) {
      if (prev.tid != -1) {
	cv::line(img, cv::Point(x - 3, midpoint), cv::Point(x, midpoint), cv::Scalar(0, 0, 0), c.lw);
	cv::line(img, cv::Point(x - 3, y), cv::Point(x - 3, midpoint), cv::Scalar(0, 0, 0), c.lw);
	int32_t prevX = 0;
	if (prev.tid < mp.tid) prevX = 0;
	else if (prev.tid > mp.tid) prevX = c.width;
	else {
	  if (prev.fwd) prevX = pixelX(c.width, rg.size, prev.gend - rg.beg) + 2;
	  else prevX = pixelX(c.width, rg.size, prev.gstart - rg.beg) - 3;
	}
	cv::line(img, cv::Point(x - 3, y), cv::Point(prevX, y), cv::Scalar(0, 0, 0), c.lw);
      }
      if (succ.tid != -1) {
	cv::line(img, cv::Point(x + w, midpoint), cv::Point(x + w + 2, midpoint), cv::Scalar(0, 0, 0), c.lw);
	cv::line(img, cv::Point(x + w + 2, midpoint), cv::Point(x + w + 2, y + c.tlheight), cv::Scalar(0, 0, 0), c.lw);
	int32_t succX = 0;
	if (succ.tid < mp.tid) succX = 0;
	else if (succ.tid > mp.tid) succX = c.width;
	else {
	  if (succ.fwd) succX = pixelX(c.width, rg.size, succ.gstart - rg.beg) - 3;
	  else succX = pixelX(c.width, rg.size, succ.gend - rg.beg) + 2;
	}
	cv::line(img, cv::Point(x + w + 2, y + c.tlheight), cv::Point(succX, y + c.tlheight), cv::Scalar(0, 0, 0), c.lw);
      }
    } else {
      if (succ.tid != -1) {
	cv::line(img, cv::Point(x - 3, midpoint), cv::Point(x, midpoint), cv::Scalar(0, 0, 0), c.lw);
	cv::line(img, cv::Point(x - 3, midpoint), cv::Point(x - 3, y + c.tlheight), cv::Scalar(0, 0, 0), c.lw);
	int32_t succX = 0;
	if (succ.tid < mp.tid) succX = 0;
	else if (succ.tid > mp.tid) succX = c.width;
	else {
	  if (succ.fwd) succX = pixelX(c.width, rg.size, succ.gstart - rg.beg) - 3;
	  else succX = pixelX(c.width, rg.size, succ.gend - rg.beg) + 2;
	}
	cv::line(img, cv::Point(x - 3, y + c.tlheight), cv::Point(succX, y + c.tlheight), cv::Scalar(0, 0, 0), c.lw);
      }
      if (prev.tid != -1) {
	cv::line(img, cv::Point(x + w, midpoint), cv::Point(x + w + 2, midpoint), cv::Scalar(0, 0, 0), c.lw);
	cv::line(img, cv::Point(x + w + 2, y), cv::Point(x + w + 2, midpoint), cv::Scalar(0, 0, 0), c.lw);
	int32_t prevX = 0;
	if (prev.tid < mp.tid) prevX = 0;
	else if (prev.tid > mp.tid) prevX = c.width;
	else {
	  if (prev.fwd) prevX = pixelX(c.width, rg.size, prev.gend - rg.beg) + 2;
	  else prevX = pixelX(c.width, rg.size, prev.gstart - rg.beg) - 3;
	}
	cv::line(img, cv::Point(x + w + 2, y), cv::Point(prevX, y), cv::Scalar(0, 0, 0), c.lw);
      }
    }

    // Draw read block
    cv::Rect rect(x, midpoint - c.rdheight/2, w, c.rdheight);    
    if (mp.fwd) cv::rectangle(img, rect, WALLY_FWDMATCH, -1);
    else cv::rectangle(img, rect, WALLY_REVMATCH, -1);
    
    // Draw read positions
    if (!c.nolabel) {
      std::string rsta = boost::lexical_cast<std::string>(mp.rstart);
      insertComma(rsta);
      std::string	rend = boost::lexical_cast<std::string>(mp.rend);
      insertComma(rend);
      std::string readcoord = rsta + " - " + rend;
      double font_scale = c.ftscale;
      double font_thickness = 3 * c.ftscale;
      int32_t baseline = 0;
      cv::Size textSize = cv::getTextSize(readcoord, cv::FONT_HERSHEY_DUPLEX, font_scale, font_thickness, &baseline);
      if (textSize.width < w) {
	cv::putText(img, readcoord, cv::Point(x + w / 2 - textSize.width / 2, y + c.tlheight/2 + textSize.height/2), cv::FONT_HERSHEY_SIMPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness, cv::LINE_AA);
      }
    }

  }



}

#endif
