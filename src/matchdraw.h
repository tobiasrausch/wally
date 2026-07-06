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
  drawSplitLine(TConfig const& c, BLContext& img, uint32_t const track) {
    // Draw line
    drawLine(img, 0, track * c.tlheight, c.width, track * c.tlheight, WALLY_SPLITBORDER, 1);
  }

  template<typename TConfig>
  inline void
  drawSplitBorder(TConfig const& c, BLContext& img) {
    drawLine(img, 0, 0, 0, c.height, WALLY_SPLITBORDER, 3);
    drawLine(img, c.width - 1, 0, c.width - 1, c.height, WALLY_SPLITBORDER, 3);
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
      TextSize textSize = getTextSize(text, scaleToPx(c.ftscale));
      if (textSize.width > (int32_t) maxsize) maxsize = textSize.width;
    }
    return maxsize * c.splits;
  }

  template<typename TConfig>
  inline void
  drawReadLabel(TConfig const& c, int32_t const track, std::string const& text, BLContext& img) {
    float px = scaleToPx(c.ftscale);
    TextSize textSize = getTextSize(text, px);
    img.fill_rect(BLRectI(0, track * c.tlheight, textSize.width, textSize.height), BLRgba32(255, 255, 0));
    drawText(img, 0, track * c.tlheight + textSize.height, text, BLRgba32(0, 0, 0), px);
  }

  template<typename TConfig>
  inline void
  drawCoordinates(TConfig const& c, Region const& rg, std::string const& chr, BLContext& img) {
    std::string textBeg(boost::lexical_cast<std::string>(rg.beg));
    insertComma(textBeg);
    std::string textEnd(boost::lexical_cast<std::string>(rg.end));
    insertComma(textEnd);
    std::string text = textBeg + "-" + textEnd;

    float px = scaleToPx(c.ftscale);
    TextSize textSize = getTextSize(chr, px);
    drawText(img, c.width / 2 - textSize.width / 2, 0 * c.tlheight + textSize.height, chr, BLRgba32(0, 0, 0), px);
    textSize = getTextSize(text, px);
    drawText(img, c.width / 2 - textSize.width / 2, 1 * c.tlheight + textSize.height, text, BLRgba32(0, 0, 0), px);

    // Draw scale
    int32_t bpscale = 0.5 * c.width * c.bpoffset;
    int32_t mult = 1;
    while (bpscale > 10) {
      bpscale /= 10;
      mult *= 10;
    }
    bpscale *= mult;
    int32_t sep = 5;
    int32_t pxs = sep;
    int32_t pxend = sep + c.pxoffset * bpscale;
    drawLine(img, pxs, 2 * c.tlheight + c.tlheight / 2, pxend, 2 * c.tlheight + c.tlheight / 2, BLRgba32(0, 0, 0), c.lw);
    drawLine(img, pxs, 2 * c.tlheight + 1, pxs, 3 * c.tlheight - 1, BLRgba32(0, 0, 0), c.lw);
    drawLine(img, pxend, 2 * c.tlheight + 1, pxend, 3 * c.tlheight - 1, BLRgba32(0, 0, 0), c.lw);

    std::string scaletxt(boost::lexical_cast<std::string>(bpscale));
    insertComma(scaletxt);
    scaletxt += "bp";
    textSize = getTextSize(scaletxt, px);
    drawText(img, pxend + sep, 2 * c.tlheight + textSize.height, scaletxt, BLRgba32(0, 0, 0), px);
  }

  template<typename TConfig>
  inline void
  drawCrossConnect(TConfig const& c, Region const& rg, BLContext& img, int32_t const track, Mapping const& prev, Mapping const& mp) {
    bool doLine = false;
    if ((prev.tid < rg.tid) || ((prev.tid == rg.tid) && (prev.gend < (rg.beg + rg.end) / 2))) {
      if ((mp.tid > rg.tid) || ((mp.tid == rg.tid) && (mp.gstart > (rg.beg + rg.end) / 2))) doLine = true;
    }
    if ((prev.tid > rg.tid) || ((prev.tid == rg.tid) && (prev.gstart > (rg.beg + rg.end) / 2))) {
      if ((mp.tid < rg.tid) || ((mp.tid == rg.tid) && (mp.gend < (rg.beg + rg.end) / 2))) doLine = true;
    }
    if (doLine) {
      int32_t y = track * c.tlheight;
      drawLine(img, 0, y, c.width, y, BLRgba32(0, 0, 0), c.lw);
    }
  }

  template<typename TConfig>
  inline void
  drawBlock(TConfig const& c, Region const& rg, BLContext& img, int32_t const track, Mapping const& prev, Mapping const& mp, Mapping const& succ) {

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
	drawLine(img, x - 3, midpoint, x, midpoint, BLRgba32(0, 0, 0), c.lw);
	drawLine(img, x - 3, y, x - 3, midpoint, BLRgba32(0, 0, 0), c.lw);
	int32_t prevX = 0;
	if (prev.tid < mp.tid) prevX = 0;
	else if (prev.tid > mp.tid) prevX = c.width;
	else {
	  if (prev.fwd) prevX = pixelX(c.width, rg.size, prev.gend - rg.beg) + 2;
	  else prevX = pixelX(c.width, rg.size, prev.gstart - rg.beg) - 3;
	}
	drawLine(img, x - 3, y, prevX, y, BLRgba32(0, 0, 0), c.lw);
      }
      if (succ.tid != -1) {
	drawLine(img, x + w, midpoint, x + w + 2, midpoint, BLRgba32(0, 0, 0), c.lw);
	drawLine(img, x + w + 2, midpoint, x + w + 2, y + c.tlheight, BLRgba32(0, 0, 0), c.lw);
	int32_t succX = 0;
	if (succ.tid < mp.tid) succX = 0;
	else if (succ.tid > mp.tid) succX = c.width;
	else {
	  if (succ.fwd) succX = pixelX(c.width, rg.size, succ.gstart - rg.beg) - 3;
	  else succX = pixelX(c.width, rg.size, succ.gend - rg.beg) + 2;
	}
	drawLine(img, x + w + 2, y + c.tlheight, succX, y + c.tlheight, BLRgba32(0, 0, 0), c.lw);
      }
    } else {
      if (succ.tid != -1) {
	drawLine(img, x - 3, midpoint, x, midpoint, BLRgba32(0, 0, 0), c.lw);
	drawLine(img, x - 3, midpoint, x - 3, y + c.tlheight, BLRgba32(0, 0, 0), c.lw);
	int32_t succX = 0;
	if (succ.tid < mp.tid) succX = 0;
	else if (succ.tid > mp.tid) succX = c.width;
	else {
	  if (succ.fwd) succX = pixelX(c.width, rg.size, succ.gstart - rg.beg) - 3;
	  else succX = pixelX(c.width, rg.size, succ.gend - rg.beg) + 2;
	}
	drawLine(img, x - 3, y + c.tlheight, succX, y + c.tlheight, BLRgba32(0, 0, 0), c.lw);
      }
      if (prev.tid != -1) {
	drawLine(img, x + w, midpoint, x + w + 2, midpoint, BLRgba32(0, 0, 0), c.lw);
	drawLine(img, x + w + 2, y, x + w + 2, midpoint, BLRgba32(0, 0, 0), c.lw);
	int32_t prevX = 0;
	if (prev.tid < mp.tid) prevX = 0;
	else if (prev.tid > mp.tid) prevX = c.width;
	else {
	  if (prev.fwd) prevX = pixelX(c.width, rg.size, prev.gend - rg.beg) + 2;
	  else prevX = pixelX(c.width, rg.size, prev.gstart - rg.beg) - 3;
	}
	drawLine(img, x + w + 2, y, prevX, y, BLRgba32(0, 0, 0), c.lw);
      }
    }

    // Draw read block
    if (mp.fwd) img.fill_rect(BLRectI(x, midpoint - c.rdheight/2, w, c.rdheight), WALLY_FWDMATCH);
    else img.fill_rect(BLRectI(x, midpoint - c.rdheight/2, w, c.rdheight), WALLY_REVMATCH);

    // Draw read positions
    if (!c.nolabel) {
      std::string rsta = boost::lexical_cast<std::string>(mp.rstart);
      insertComma(rsta);
      std::string	rend = boost::lexical_cast<std::string>(mp.rend);
      insertComma(rend);
      std::string readcoord = rsta + " - " + rend;
      float px = scaleToPx(c.ftscale);
      TextSize textSize = getTextSize(readcoord, px);
      if (textSize.width < w) {
	drawText(img, x + w / 2 - textSize.width / 2, y + c.tlheight/2 + textSize.height/2, readcoord, BLRgba32(0, 0, 0), px);
      }
    }

  }



}

#endif
