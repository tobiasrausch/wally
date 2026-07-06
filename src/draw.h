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
#include <boost/dynamic_bitset.hpp>

#include <htslib/sam.h>

#include "util.h"
#include "bed.h"

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

  inline uint32_t
  findTicks(double const pxoffset, double const textwidth) {
    std::vector<int32_t> bounds{5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 10000000, 25000000, 50000000};
    for(uint32_t i = 0; i < bounds.size(); ++i) {
      if (bounds[i] * pxoffset > textwidth * 1.5) {
	return bounds[i];
	break;
      }
    }
    return 50000000;
  }

  template<typename TConfig>
  inline void
  drawGenome(TConfig const& c, Region const& rg, std::string const& refname, BLContext& img, int32_t const track) {
    std::string text(boost::lexical_cast<std::string>(rg.end));
    int32_t n = text.length() - 3;
    while (n > 0) {
      text.insert(n, ",");
      n -= 3;
    }
    TextSize textSize = getTextSize(text);

    // Find suitable tick size
    uint32_t modval = findTicks(c.pxoffset, textSize.width);
    drawLine(img, 0, track * c.tlheight + c.tlheight, c.width, track * c.tlheight + c.tlheight, BLRgba32(0, 0, 0), 1);
    double px = 0;
    int32_t prevTick = -1;
    int32_t lastTick = -1;
    for(int32_t i = rg.beg; i < rg.end; ++i) {
      if (i % modval == 0) {
	drawLine(img, px - c.pxoffset/2, track * c.tlheight, px - c.pxoffset/2, track * c.tlheight + c.tlheight, BLRgba32(0, 0, 0), 1);
	if (prevTick == -1) prevTick = px;
	else {
	  if (lastTick == -1) {
	    lastTick = px;
	    int32_t midpoint = (prevTick + lastTick) / 2;
	    TextSize textSize = getTextSize(refname);
	    drawText(img, midpoint - textSize.width/2, track * c.tlheight + textSize.height, refname, BLRgba32(0, 0, 0));
	  }
	}
      }
      if (i % modval == 0) {
	// Font
	text = boost::lexical_cast<std::string>(i);
	n = text.length() - 3;
	while (n > 0) {
	  text.insert(n, ",");
	  n -= 3;
	}
	TextSize textSize = getTextSize(text);
	if ((px - c.pxoffset/2 - textSize.width/2 > 0) && (px - c.pxoffset/2 + textSize.width < c.width)) {
	  drawText(img, px - c.pxoffset/2 - textSize.width/2, (track - 1) * c.tlheight + textSize.height, text, BLRgba32(0, 0, 0));
	}
      }
      px += c.pxoffset;
    }
  }

  template<typename TConfig>
  inline void
  drawAnnotation(TConfig const& c, Region const& rg, std::vector<Transcript> const& tr, std::vector<Region> const& anno, BLContext& img, int32_t const track) {

    // Transcripts
    if (!tr.empty()) {
      for(uint32_t i = 0; i < tr.size(); ++i) {
	int32_t px = pixelX(c.width, rg.size, tr[i].rg.beg - rg.beg);
	int32_t pxend = pixelX(c.width, rg.size, tr[i].rg.end - rg.beg);
	drawLine(img, px, track * c.tlheight + c.tlheight/2, pxend, track * c.tlheight + c.tlheight/2, BLRgba32(0, 0, 255), 1);
	int32_t pxi = px + 5;
	if (pxi < 0) pxi = 5;
	if (tr[i].forward) {
	  for(; ((pxi < pxend) && (pxi < (int32_t) c.width)); pxi += 20) {
	    drawLine(img, pxi, track * c.tlheight + 1, pxi+5, track * c.tlheight + c.tlheight/2, BLRgba32(0, 0, 255), 1);
	    drawLine(img, pxi, track * c.tlheight + c.tlheight - 2, pxi+5, track * c.tlheight + c.tlheight/2, BLRgba32(0, 0, 255), 1);
	  }
	} else {
	  for(; ((pxi < pxend) && (pxi < (int32_t) c.width)); pxi += 20) {
	    drawLine(img, pxi, track * c.tlheight + 1, pxi-5, track * c.tlheight + c.tlheight/2, BLRgba32(0, 0, 255), 1);
	    drawLine(img, pxi, track * c.tlheight + c.tlheight - 2, pxi-5, track * c.tlheight + c.tlheight/2, BLRgba32(0, 0, 255), 1);
	  }
	}
      }
    }

    // Annotations (e.g., exons)
    typedef boost::dynamic_bitset<> TBitSet;
    TBitSet blocked(c.width, false);
    if (!anno.empty()) {
      for(uint32_t i = 0; i < anno.size(); ++i) {
	int32_t px = pixelX(c.width, rg.size, anno[i].beg - rg.beg);
	int32_t pxend = pixelX(c.width, rg.size, anno[i].end - rg.beg);
	img.fill_rect(BLRectI(px, track * c.tlheight + 1, pxend - px, c.tlheight - 2), hexToScalar(anno[i].color));

	// Block regions
	if (px < 0) px = 0;
	if (pxend > (int32_t) c.width) pxend = c.width;
	if ((pxend - px > 2) && (!tr.empty())) {
	  for(int32_t k = px; k < pxend; ++k) blocked[k] = true;
	}

	// Add text
	TextSize textSize = getTextSize(anno[i].id);
	int32_t pxmid = (px + pxend) / 2 - textSize.width/2;
	if ((pxmid > px) && (pxmid + textSize.width < pxend)) {
	  drawText(img, pxmid, track * c.tlheight + c.tlheight - 2, anno[i].id, BLRgba32(255, 255, 255));
	}
      }
    }

    // Transcript labels
    if (!tr.empty()) {
      for(uint32_t i = 0; i < tr.size(); ++i) {
	int32_t px = pixelX(c.width, rg.size, tr[i].rg.beg - rg.beg);
	int32_t pxend = pixelX(c.width, rg.size, tr[i].rg.end - rg.beg);
	if (px < 0) px = 0;
	if (pxend > (int32_t) c.width) pxend = c.width;
	int32_t kstart = -1;
	int32_t kend = -1;
	TextSize textSize = getTextSize(tr[i].rg.id);
	for(int32_t k = px; k < pxend; ++k) {
	  if (!blocked[k]) {
	    if (kstart == -1) kstart = k;
	    kend = k;
	    if ((kend - kstart) > 2.5 * textSize.width) {
	      int32_t midpoint = (kstart + kend) / 2;
	      img.fill_rect(BLRectI(midpoint - textSize.width/2, track * c.tlheight, textSize.width, c.tlheight), BLRgba32(255, 255, 255));
	      drawText(img, midpoint - textSize.width/2, track * c.tlheight + c.tlheight - 2, tr[i].rg.id, BLRgba32(0, 0, 255));
	      for(int32_t l = kstart; l < kend; ++l) blocked[l] = true;
	      break;
	    }
	  } else {
	    kstart = -1;
	    kend = -1;
	  }
	}
      }
    }
  }

  template<typename TConfig>
  inline void
  drawReference(TConfig const& c, Region const& rg, BLContext& img, std::string const& ref, int32_t const track) {
    if (c.pxoffset >= WALLY_PX) {
      double px = 0;
      for(int32_t i = 0; i < rg.size; ++i) {
	std::string text(1, ref[i]);
	TextSize textSize = getTextSize(text);
	drawText(img, px + c.pxoffset/2 - textSize.width/2, track * c.tlheight + textSize.height + 2, text, BLRgba32(0, 0, 0));
	px += c.pxoffset;
      }
    }
  }

  template<typename TConfig>
  inline void
  drawBorder(TConfig const& c, BLContext& img) {
    drawLine(img, 0, 0, 0, c.height, WALLY_BORDER, 1);
    drawLine(img, c.width - 1, 0, c.width - 1, c.height, WALLY_BORDER, 1);
  }

  template<typename TConfig, typename TCoverage>
  inline void
  drawCoverage(TConfig const& c, Region const& rg, BLContext& img, TCoverage const& covA, TCoverage const& covC, TCoverage const& covG, TCoverage const& covT, std::vector<bool> const& snp, int32_t const track) {
    uint32_t maxObsCov = 0;
    for(uint32_t i = 0; i < covA.size(); ++i) {
      uint32_t cumsum = covA[i];
      cumsum += covC[i];
      cumsum += covG[i];
      cumsum += covT[i];
      if (cumsum > maxObsCov) maxObsCov = cumsum;
    }
    if (maxObsCov == 0) maxObsCov = 1;

    // Draw horizontal top line
    drawLine(img, 0, (track-1) * c.tlheight, c.width, (track-1) * c.tlheight, WALLY_BORDER, 1);

    // Draw coverage histogram
    double px = 0;
    for(int32_t i = 0; i < rg.size; ++i) {
      if (!snp[i]) {
	uint32_t cumsum = covA[i];
	cumsum += covC[i];
	cumsum += covG[i];
	cumsum += covT[i];
	double frac = (double) cumsum / (double) maxObsCov;
	int32_t ch = (int32_t) (frac * 2 * c.tlheight);
	if (c.pxoffset >= WALLY_PX) {
	  img.fill_rect(BLRectI(px + 1, (track-1) * c.tlheight + 2 * c.tlheight - ch, c.pxoffset - 1, ch), BLRgba32(100, 100, 100));
	} else {
	  img.fill_rect(BLRectI(px, (track-1) * c.tlheight + 2 * c.tlheight - ch, c.pxoffset + 1, ch), BLRgba32(100, 100, 100));
	}
      }
      px += c.pxoffset;
    }

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
	img.fill_rect(BLRectI(pxs, (track-1) * c.tlheight + 2 * c.tlheight - cumCH, pxw, ch), WALLY_A);
	frac = (double) covC[i] / (double) maxObsCov;
	ch = (int32_t) (frac * 2 * c.tlheight);
	cumCH += ch;
	img.fill_rect(BLRectI(pxs, (track-1) * c.tlheight + 2 * c.tlheight - cumCH, pxw, ch), WALLY_C);
	frac = (double) covG[i] / (double) maxObsCov;
	ch = (int32_t) (frac * 2 * c.tlheight);
	cumCH += ch;
	img.fill_rect(BLRectI(pxs, (track-1) * c.tlheight + 2 * c.tlheight - cumCH, pxw, ch), WALLY_G);
	frac = (double) covT[i] / (double) maxObsCov;
	ch = (int32_t) (frac * 2 * c.tlheight);
	cumCH += ch;
	img.fill_rect(BLRectI(pxs, (track-1) * c.tlheight + 2 * c.tlheight - cumCH, pxw, ch), WALLY_T);
      }
      px += c.pxoffset;
    }

    // Put coverage scale
    if (true) {
      std::string text = "[0-";
      text += boost::lexical_cast<std::string>(maxObsCov);
      text += "]";
      TextSize textSize = getTextSize(text);
      img.fill_rect(BLRectI(0, track * c.tlheight - textSize.height/2, textSize.width, textSize.height), BLRgba32(255, 255, 0));
      drawText(img, 0, track * c.tlheight + textSize.height/2, text, BLRgba32(0, 0, 0));
    }
  }

  template<typename TConfig>
  inline void
  drawSampleLabel(TConfig const& c, int32_t const track, std::string const& text, BLContext& img) {
    TextSize textSize = getTextSize(text);
    img.fill_rect(BLRectI(0, track * c.tlheight - textSize.height, textSize.width, textSize.height), BLRgba32(255, 255, 0));
    drawText(img, 0, track * c.tlheight, text, BLRgba32(0, 0, 0));
  }


  template<typename TConfig>
  inline void
  drawRead(TConfig const& c, BLContext& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, bool const reverse, bool const tri, BLRgba32 const& clr) {
    img.fill_rect(BLRectI(x, y, w, h), clr);
    if (tri) {
      if (reverse) {
	BLPoint pvec[3] = {BLPoint(x, y), BLPoint(x, y+h-1), BLPoint(x-c.pxoffset/3, y + h/2)};
	img.fill_polygon(pvec, 3, clr);
      } else {
	BLPoint pvec[3] = {BLPoint(x+w, y), BLPoint(x+w, y+h-1), BLPoint(x+w+c.pxoffset/3, y + h/2)};
	img.fill_polygon(pvec, 3, clr);
      }
    }
  }

  template<typename TConfig>
  inline void
  drawRead(TConfig const& c, Region const& rg, BLContext& img, int32_t const track, int32_t const gstart, int32_t const gend, bool const reverse, bool const tri, BLRgba32 const& clr) {
    if (track == -1) return;
    int32_t px = pixelX(c.width, rg.size, gstart);
    int32_t pxend = pixelX(c.width, rg.size, gend);
    drawRead(c, img, px, track * c.tlheight, pxend - px, c.rdheight, reverse, tri, clr);
  }

  template<typename TConfig>
  inline void
  drawNuc(TConfig const& c, BLContext& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, char const nuc, BLRgba32 const& clr) {
    // Font
    std::string text(1, nuc);
    TextSize textSize = getTextSize(text);

    // Put nucleotide if there is space
    if (c.pxoffset >= WALLY_PX) {
      img.fill_rect(BLRectI(x, y, w, h), clr);
      if (nuc == 'A') {
	drawText(img, x + w/2 - textSize.width/2, y+h/2+textSize.height/2, text, WALLY_A);
      }
      else if (nuc == 'C') {
      	drawText(img, x + w/2 - textSize.width/2, y+h/2+textSize.height/2, text, WALLY_C);
      }
      else if (nuc == 'G') {
	drawText(img, x + w/2 - textSize.width/2, y+h/2+textSize.height/2, text, WALLY_G);
      }
      else if (nuc == 'T') {
	drawText(img, x + w/2 - textSize.width/2, y+h/2+textSize.height/2, text, WALLY_T);
      } else {  // n or N most likely
	drawText(img, x + w/2 - textSize.width/2, y+h/2+textSize.height/2, text, WALLY_N);
      }
    } else {
      int32_t pxw = w;
      if (pxw < 1) pxw = 1; // Make mismatches always visible
      if (nuc == 'A') {
	img.fill_rect(BLRectI(x, y, pxw, h), WALLY_A);
      }
      else if (nuc == 'C') {
	img.fill_rect(BLRectI(x, y, pxw, h), WALLY_C);
      }
      else if (nuc == 'G') {
	img.fill_rect(BLRectI(x, y, pxw, h), WALLY_G);
      }
      else if (nuc == 'T') {
	img.fill_rect(BLRectI(x, y, pxw, h), WALLY_T);
      }
    }
  }

  template<typename TConfig>
  inline void
  drawNuc(TConfig const& c, Region const& rg, BLContext& img, int32_t const track, int32_t const gstart, int32_t const gend, char const nuc, BLRgba32 const& clr) {
    if (track == -1) return;
    int32_t px = pixelX(c.width, rg.size, gstart);
    int32_t pxend = pixelX(c.width, rg.size, gend);
    drawNuc(c, img, px, track * c.tlheight, pxend - px, c.rdheight, nuc, clr);
  }

  template<typename TConfig>
  inline void
  drawModBase(TConfig const& c, Region const& rg, BLContext& img, int32_t const track, int32_t const gstart, int16_t const prob, BLRgba32 const& modClr, BLRgba32 const& unmodClr) {
    if (track == -1) return;
    if (prob < 0) return;
    int32_t px = pixelX(c.width, rg.size, gstart);
    int32_t pxend = pixelX(c.width, rg.size, gstart + 1);
    int32_t w = pxend - px;
    if (w < 1) w = 1;
    double p = (double) prob / 255.0;
    BLRgba32 clr((uint32_t) (unmodClr.r() * (1 - p) + modClr.r() * p), (uint32_t) (unmodClr.g() * (1 - p) + modClr.g() * p), (uint32_t) (unmodClr.b() * (1 - p) + modClr.b() * p));
    img.fill_rect(BLRectI(px, track * c.tlheight, w, c.rdheight), clr);
  }

  inline void
  drawDel(BLContext& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, int32_t const len) {
    int32_t pxw = w;
    if (w < 1) pxw = 1;
    drawLine(img, x, y+h/2, x+pxw, y+h/2, BLRgba32(0, 0, 0), 2);
    if (w > 1) {
      std::string text = boost::lexical_cast<std::string>(len);
      TextSize textSize = getTextSize(text);
      double frac = (double) textSize.width / (double) w;
      // Put length if there is space
      if (frac < 0.5) {
	img.fill_rect(BLRectI(x+w/2-textSize.width/2, y, textSize.width, h), BLRgba32(255, 255, 255));
	drawText(img, x+w/2 - textSize.width/2, y + h/2 + textSize.height/2, text, BLRgba32(148, 0, 211));
      }
    }
  }

  template<typename TConfig>
  inline void
    drawDel(TConfig const& c, Region const& rg, BLContext& img, int32_t const track, int32_t const gstart, int32_t const gend, int32_t const len) {
    if (track == -1) return;
    int32_t px = pixelX(c.width, rg.size, gstart);
    int32_t pxend = pixelX(c.width, rg.size, gend);
    drawDel(img, px, track * c.tlheight, pxend - px, c.rdheight, len);
  }

  inline void
  drawPELine(BLContext& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, BLRgba32 const& clr) {
    drawLine(img, x, y+h/2, x+w, y+h/2, clr, 1);
  }

  template<typename TConfig>
  inline void
  drawPELine(TConfig const& c, Region const& rg, BLContext& img, int32_t const track, int32_t const gstart, int32_t const gend, BLRgba32 const& clr) {
    if ((track != -1) && (gstart < gend)) {
      int32_t px = pixelX(c.width, rg.size, gstart);
      int32_t pxend = pixelX(c.width, rg.size, gend);
      drawPELine(img, px, track * c.tlheight, pxend - px, c.rdheight, clr);
    }
  }

  inline void
  drawRefSkip(BLContext& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h) {
    drawLine(img, x, y+h/2, x+w, y+h/2, BLRgba32(64, 128, 154), 1);
  }

  template<typename TConfig>
  inline void
    drawRefSkip(TConfig const& c, Region const& rg, BLContext& img, int32_t const track, int32_t const gstart, int32_t const gend) {
    if (track == -1) return;
    int32_t px = pixelX(c.width, rg.size, gstart);
    int32_t pxend = pixelX(c.width, rg.size, gend);
    drawRefSkip(img, px, track * c.tlheight, pxend - px, c.rdheight);
  }

  inline void
  drawIns(BLContext& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, int32_t const len) {
    int32_t psw = w/2;
    if (w < 1) psw = -1 + w;
    BLPoint pvec[3] = {BLPoint(x+psw, y), BLPoint(x+psw, y+h), BLPoint(x+w, y + h/2)};
    img.fill_polygon(pvec, 3, BLRgba32(118, 24, 220));
    if (w > 1) {
      std::string text = boost::lexical_cast<std::string>(len);
      TextSize textSize = getTextSize(text);
      // Put length if space takes max. 5 genomic position
      if (textSize.width <= 5 * w) {
	drawText(img, x+w/2 - textSize.width, y + h/2 + textSize.height/2, text, BLRgba32(148, 0, 211));
      }
    }
  }

  template<typename TConfig>
  inline void
    drawIns(TConfig const& c, Region const& rg, BLContext& img, int32_t const track, int32_t const gstart, int32_t const len) {
    if (track == -1) return;
    int32_t px = pixelX(c.width, rg.size, gstart - 1);
    drawIns(img, px, track * c.tlheight, c.pxoffset, c.rdheight, len);
  }

  template<typename TConfig>
  inline void
  drawSC(TConfig const& c, BLContext& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, int32_t const len, bool const leading) {
    if (c.pxoffset < 0.1) return; // Do nothing
    drawLine(img, x + w - 1, y, x + w - 1, y+h, BLRgba32(0, 0, 0), 2);
    if (w < 1) return; // Just draw the soft-clipping line
    std::string text = boost::lexical_cast<std::string>(len);
    TextSize textSize = getTextSize(text);
    // Put length if space takes max. 5 genomic position
    if (textSize.width <= 5 * w) {
      if (leading) {
	drawText(img, x+w - 1, y + h/2 + textSize.height/2, text, BLRgba32(0, 0, 0));
      } else {
	drawText(img, x+w - 1 - textSize.width, y + h/2 + textSize.height/2, text, BLRgba32(0, 0, 0));
      }
    }
  }

  template<typename TConfig>
  inline void
  drawSC(TConfig const& c, Region const& rg, BLContext& img, int32_t const track, int32_t const gstart, int32_t const len, bool const leading) {
    if (track == -1) return;
    if (c.showSoftClip) {
      int32_t px = pixelX(c.width, rg.size, gstart - 1);
      drawSC(c, img, px, track * c.tlheight, c.pxoffset, c.rdheight, len, leading);
    }
  }

}

#endif
