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
    std::vector<int32_t> bounds{5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000};
    for(uint32_t i = 0; i < bounds.size(); ++i) {
      if (bounds[i] * pxoffset > textwidth * 1.5) {
	return bounds[i];
	break;
      }
    }
    return 1000000;
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
    uint32_t modval = findTicks(c.pxoffset, textSize.width);
    cv::line(img, cv::Point(0, track * c.tlheight + c.tlheight), cv::Point(c.width, track * c.tlheight + c.tlheight), cv::Scalar(0, 0, 0), 1.8);
    double px = 0;
    for(int32_t i = rg.beg; i < rg.end; ++i) {
      if (i % modval == 0) {
	cv::line(img, cv::Point(px - c.pxoffset/2, track * c.tlheight), cv::Point(px - c.pxoffset/2, track * c.tlheight + c.tlheight), cv::Scalar(0, 0, 0), 1.8);
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
	if ((px - c.pxoffset/2 - textSize.width/2 > 0) && (px - c.pxoffset/2 + textSize.width < c.width)) {
	  cv::putText(img, text, cv::Point(px - c.pxoffset/2 - textSize.width/2, (track - 1) * c.tlheight + textSize.height), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness);
	}
      }
      px += c.pxoffset;
    }
  }

  template<typename TConfig>
  inline void
  drawAnnotation(TConfig const& c, Region const& rg, std::vector<Transcript> const& tr, std::vector<Region> const& anno, cv::Mat& img, int32_t const track) {

    // Transcripts
    if (!tr.empty()) {
      for(uint32_t i = 0; i < tr.size(); ++i) {
	int32_t px = pixelX(c.width, rg.size, tr[i].rg.beg - rg.beg + 1);
	int32_t pxend = pixelX(c.width, rg.size, tr[i].rg.end - rg.beg + 1);
	cv::line(img, cv::Point(px, track * c.tlheight + c.tlheight/2), cv::Point(pxend, track * c.tlheight + c.tlheight/2), cv::Scalar(255, 0, 0), 1);
	int32_t pxi = px + 5;
	if (pxi < 0) pxi = 5;
	if (tr[i].forward) {
	  for(; ((pxi < pxend) && (pxi < (int32_t) c.width)); pxi += 20) {
	    cv::line(img, cv::Point(pxi, track * c.tlheight + 1), cv::Point(pxi+5, track * c.tlheight + c.tlheight/2), cv::Scalar(255, 0, 0), 1);
	    cv::line(img, cv::Point(pxi, track * c.tlheight + c.tlheight - 2), cv::Point(pxi+5, track * c.tlheight + c.tlheight/2), cv::Scalar(255, 0, 0), 1);
	  }
	} else {
	  for(; ((pxi < pxend) && (pxi < (int32_t) c.width)); pxi += 20) {
	    cv::line(img, cv::Point(pxi, track * c.tlheight + 1), cv::Point(pxi-5, track * c.tlheight + c.tlheight/2), cv::Scalar(255, 0, 0), 1);
	    cv::line(img, cv::Point(pxi, track * c.tlheight + c.tlheight - 2), cv::Point(pxi-5, track * c.tlheight + c.tlheight/2), cv::Scalar(255, 0, 0), 1);
	  }
	}
      }
    }

    // Annotations (e.g., exons)
    double font_scale = 0.4;
    double font_thickness = 1.5;
    int32_t baseline = 0;
    typedef boost::dynamic_bitset<> TBitSet;
    TBitSet blocked(c.width, false);
    if (!anno.empty()) {
      for(uint32_t i = 0; i < anno.size(); ++i) {
	int32_t px = pixelX(c.width, rg.size, anno[i].beg - rg.beg + 1);
	int32_t pxend = pixelX(c.width, rg.size, anno[i].end - rg.beg + 1);
	cv::Rect rect(px, track * c.tlheight + 1, pxend - px, c.tlheight - 2);
	cv::rectangle(img, rect, cv::Scalar(255, 0, 0), -1);

	// Block regions
	if (px < 0) px = 0;
	if (pxend > (int32_t) c.width) pxend = c.width;
	if ((pxend - px > 2) && (!tr.empty())) {
	  for(int32_t k = px; k < pxend; ++k) blocked[k] = true;
	}
	
	// Add text
	cv::Size textSize = cv::getTextSize(anno[i].id, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);
	int32_t pxmid = (px + pxend) / 2 - textSize.width/2;
	if ((pxmid > px) && (pxmid + textSize.width < pxend)) {
	  cv::putText(img, anno[i].id, cv::Point(pxmid, track * c.tlheight + c.tlheight - 2), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(255, 255, 255), font_thickness);
	}
      }
    }

    // Transcript labels
    if (!tr.empty()) {
      for(uint32_t i = 0; i < tr.size(); ++i) {
	int32_t px = pixelX(c.width, rg.size, tr[i].rg.beg - rg.beg + 1);
	int32_t pxend = pixelX(c.width, rg.size, tr[i].rg.end - rg.beg + 1);
	if (px < 0) px = 0;
	if (pxend > (int32_t) c.width) pxend = c.width;
	int32_t kstart = -1;
	int32_t kend = -1;
	cv::Size textSize = cv::getTextSize(tr[i].rg.id, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);
	for(int32_t k = px; k < pxend; ++k) {
	  if (!blocked[k]) {
	    if (kstart == -1) kstart = k;
	    kend = k;
	    if ((kend - kstart) > 2.5 * textSize.width) {
	      int32_t midpoint = (kstart + kend) / 2;
	      cv::Rect rect(midpoint - textSize.width/2, track * c.tlheight, textSize.width, c.tlheight);
	      cv::rectangle(img, rect, cv::Scalar(255, 255, 255), -1);
	      cv::putText(img, tr[i].rg.id, cv::Point(midpoint - textSize.width/2, track * c.tlheight + c.tlheight - 2), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(255, 0, 0), font_thickness);
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
  drawReference(TConfig const& c, Region const& rg, cv::Mat& img, std::string const& ref, int32_t const track) {
    if (c.pxoffset >= WALLY_PX) {
      double px = 0;
      double font_scale = 0.4;
      double font_thickness = 1.5;
      int32_t baseline = 0;
      for(int32_t i = 0; i < rg.size; ++i) {
	std::string text(1, ref[i]);
	cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);
	cv::putText(img, text, cv::Point(px + c.pxoffset/2 - textSize.width/2, track * c.tlheight + textSize.height + 1), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness);
	px += c.pxoffset;
      }
    }
  }

  template<typename TConfig>
  inline void
  drawBorder(TConfig const& c, cv::Mat& img) {
    cv::line(img, cv::Point(0, 0), cv::Point(0, c.height), WALLY_BORDER, 1);
    cv::line(img, cv::Point(c.width - 1, 0), cv::Point(c.width - 1, c.height), WALLY_BORDER, 1);
  }
  
  template<typename TConfig, typename TCoverage>
  inline void
  drawCoverage(TConfig const& c, Region const& rg, cv::Mat& img, TCoverage const& covA, TCoverage const& covC, TCoverage const& covG, TCoverage const& covT, std::vector<bool> const& snp, int32_t const track) {
    uint32_t maxObsCov = 0;
    for(uint32_t i = 0; i < covA.size(); ++i) {
      uint32_t cumsum = covA[i];
      cumsum += covC[i];
      cumsum += covG[i];
      cumsum += covT[i];
      if (cumsum > maxObsCov) maxObsCov = cumsum;
    }

    // Draw horizontal top line
    cv::line(img, cv::Point(0, (track-1) * c.tlheight), cv::Point(c.width, (track-1) * c.tlheight), WALLY_BORDER, 1);
    
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
	  cv::Rect rect(px + 1, (track-1) * c.tlheight + 2 * c.tlheight - ch, c.pxoffset - 1, ch);
	  cv::rectangle(img, rect, cv::Scalar(100, 100, 100), -1);
	} else {
	  cv::Rect rect(px, (track-1) * c.tlheight + 2 * c.tlheight - ch, c.pxoffset + 1, ch);
	  cv::rectangle(img, rect, cv::Scalar(100, 100, 100), -1);
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
	std::vector<cv::Point> pvec{cv::Point(x, y), cv::Point(x, y+h), cv::Point(x-c.pxoffset/3, y + h/2)};
	cv::polylines(img, pvec, true, cv::Scalar(200, 200, 200), 1);
	cv::fillPoly(img, pvec, cv::Scalar(200, 200, 200));
      } else {
	std::vector<cv::Point> pvec{cv::Point(x+w, y), cv::Point(x+w, y+h), cv::Point(x+w+c.pxoffset/3, y + h/2)};
	cv::polylines(img, pvec, true, cv::Scalar(200, 200, 200), 1);
	cv::fillPoly(img, pvec, cv::Scalar(200, 200, 200));
      }
    }
  }

  template<typename TConfig>
  inline void
  drawRead(TConfig const& c, Region const& rg, cv::Mat& img, int32_t const track, int32_t const gstart, int32_t const gend, bool const reverse, bool const tri) {
    if (track == -1) return;
    int32_t px = pixelX(c.width, rg.size, gstart);
    int32_t pxend = pixelX(c.width, rg.size, gend);
    drawRead(c, img, px, track * c.tlheight, pxend - px, c.rdheight, reverse, tri);
  }

  template<typename TConfig>
  inline void
  drawNuc(TConfig const& c, cv::Mat& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, char const nuc) {
    // Font
    std::string text(1, nuc);
    double font_scale = 0.4;
    double font_thickness = 1.5;
    int32_t baseline = 0;
    cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);

    // Put nucleotide if there is space
    if (c.pxoffset >= WALLY_PX) {
      cv::Rect rect(x, y, w, h);
      cv::rectangle(img, rect, cv::Scalar(200, 200, 200), -1);
      if ((nuc == 'a') or (nuc == 'A')) {
	cv::putText(img, text, cv::Point(x + w/2 - textSize.width/2, y+h/2+textSize.height/2), cv::FONT_HERSHEY_DUPLEX, font_scale, WALLY_A, font_thickness);
      }
      else if ((nuc == 'c') or (nuc == 'C')) {
      	cv::putText(img, text, cv::Point(x + w/2 - textSize.width/2, y+h/2+textSize.height/2), cv::FONT_HERSHEY_DUPLEX, font_scale, WALLY_C, font_thickness);
      }
      else if ((nuc == 'g') or (nuc == 'G')) {
	cv::putText(img, text, cv::Point(x + w/2 - textSize.width/2, y+h/2+textSize.height/2), cv::FONT_HERSHEY_DUPLEX, font_scale, WALLY_G, font_thickness);
      }
      else if ((nuc == 't') or (nuc == 'T')) {
	cv::putText(img, text, cv::Point(x + w/2 - textSize.width/2, y+h/2+textSize.height/2), cv::FONT_HERSHEY_DUPLEX, font_scale, WALLY_T, font_thickness);
      }
    } else {
      int32_t pxw = w;
      if (pxw < 1) pxw = 1; // Make mismatches always visible
      cv::Rect rect(x, y, pxw, h);
      if ((nuc == 'a') or (nuc == 'A')) {
	cv::rectangle(img, rect, WALLY_A, -1);
      }
      else if ((nuc == 'c') or (nuc == 'C')) {
	cv::rectangle(img, rect, WALLY_C, -1);
      }
      else if ((nuc == 'g') or (nuc == 'G')) {
	cv::rectangle(img, rect, WALLY_G, -1);
      }
      else if ((nuc == 't') or (nuc == 'T')) {
	cv::rectangle(img, rect, WALLY_T, -1);
      }
    }
  }
  
  template<typename TConfig>
  inline void
  drawNuc(TConfig const& c, Region const& rg, cv::Mat& img, int32_t const track, int32_t const gstart, int32_t const gend, char const nuc) {
    if (track == -1) return;
    int32_t px = pixelX(c.width, rg.size, gstart);
    int32_t pxend = pixelX(c.width, rg.size, gend);
    drawNuc(c, img, px, track * c.tlheight, pxend - px, c.rdheight, nuc);
  }

  inline void
  drawDel(cv::Mat& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, int32_t const len) {
    int32_t pxw = w;
    if (w < 1) pxw = 1;
    cv::line(img, cv::Point(x, y+h/2), cv::Point(x+pxw, y+h/2), cv::Scalar(0, 0, 0), 2);
    if (w > 1) {
      std::string text = boost::lexical_cast<std::string>(len);
      double font_scale = 0.4;
      double font_thickness = 1.5;
      int32_t baseline = 0;
      cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_DUPLEX, font_scale, font_thickness, &baseline);
      double frac = (double) textSize.width / (double) w;
      // Put length if there is space
      if (frac < 0.5) {
	cv::Rect rect(x+w/2-textSize.width/2, y, textSize.width, h);
	cv::rectangle(img, rect, cv::Scalar(255, 255, 255), -1);
	cv::putText(img, text, cv::Point(x+w/2 - textSize.width/2, y + h/2 + textSize.height/2), cv::FONT_HERSHEY_SIMPLEX, font_scale, cv::Scalar(211, 0, 148), font_thickness);
      }
    }
  }

  template<typename TConfig>
  inline void
    drawDel(TConfig const& c, Region const& rg, cv::Mat& img, int32_t const track, int32_t const gstart, int32_t const gend, int32_t const len) {
    if (track == -1) return;
    int32_t px = pixelX(c.width, rg.size, gstart);
    int32_t pxend = pixelX(c.width, rg.size, gend);
    drawDel(img, px, track * c.tlheight, pxend - px, c.rdheight, len);
  }

  inline void
  drawRefSkip(cv::Mat& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h) {
    cv::line(img, cv::Point(x, y+h/2), cv::Point(x+w, y+h/2), cv::Scalar(154, 128, 64), 1);
  }

  template<typename TConfig>
  inline void
    drawRefSkip(TConfig const& c, Region const& rg, cv::Mat& img, int32_t const track, int32_t const gstart, int32_t const gend) {
    if (track == -1) return;
    int32_t px = pixelX(c.width, rg.size, gstart);
    int32_t pxend = pixelX(c.width, rg.size, gend);
    drawRefSkip(img, px, track * c.tlheight, pxend - px, c.rdheight);
  }

  inline void
  drawIns(cv::Mat& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, int32_t const len) {
    int32_t psw = w/2;
    if (w < 1) psw = -1 + w;
    std::vector<cv::Point> pvec{cv::Point(x+psw, y), cv::Point(x+psw, y+h), cv::Point(x+w, y + h/2)};
    cv::polylines(img, pvec, true, cv::Scalar(220, 24, 118), 1);
    cv::fillPoly(img, pvec, cv::Scalar(220, 24, 118));
    if (w > 1) {
      std::string text = boost::lexical_cast<std::string>(len);
      double font_scale = 0.4;
      double font_thickness = 1.5;
      int32_t baseline = 0;
      cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_DUPLEX, font_scale, font_thickness, &baseline);
      // Put length if space takes max. 5 genomic position
      if (textSize.width <= 5 * w) {
	cv::putText(img, text, cv::Point(x+w/2 - textSize.width, y + h/2 + textSize.height/2), cv::FONT_HERSHEY_SIMPLEX, font_scale, cv::Scalar(211, 0, 148), font_thickness);
      }
    }
  }
    
  template<typename TConfig>
  inline void
    drawIns(TConfig const& c, Region const& rg, cv::Mat& img, int32_t const track, int32_t const gstart, int32_t const len) {
    if (track == -1) return;
    int32_t px = pixelX(c.width, rg.size, gstart - 1);
    drawIns(img, px, track * c.tlheight, c.pxoffset, c.rdheight, len);
  }

  inline void
  drawSC(cv::Mat& img, int32_t const x, int32_t const y, int32_t const w, int32_t const h, int32_t const len, bool const leading) {
    if (w < 1) return; // Do nothing
    cv::line(img, cv::Point(x + w - 1, y), cv::Point(x + w - 1, y+h), cv::Scalar(0, 0, 0), 2);
    std::string text = boost::lexical_cast<std::string>(len);
    double font_scale = 0.4;
    double font_thickness = 1.5;
    int32_t baseline = 0;
    cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_DUPLEX, font_scale, font_thickness, &baseline);
    // Put length if space takes max. 5 genomic position
    if (textSize.width <= 5 * w) {
      if (leading) {
	cv::putText(img, text, cv::Point(x+w - 1, y + h/2 + textSize.height/2), cv::FONT_HERSHEY_SIMPLEX, font_scale, cv::Scalar(211, 0, 148), font_thickness);
      } else {
	cv::putText(img, text, cv::Point(x+w - 1 - textSize.width, y + h/2 + textSize.height/2), cv::FONT_HERSHEY_SIMPLEX, font_scale, cv::Scalar(211, 0, 148), font_thickness);
      }
    }
  }

  template<typename TConfig>
  inline void
  drawSC(TConfig const& c, Region const& rg, cv::Mat& img, int32_t const track, int32_t const gstart, int32_t const len, bool const leading) {
    if (track == -1) return;
    if (c.showSoftClip) {
      int32_t px = pixelX(c.width, rg.size, gstart - 1);
      drawSC(img, px, track * c.tlheight, c.pxoffset, c.rdheight, len, leading);
    }
  }
  
}

#endif
