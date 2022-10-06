#ifndef GFADRAW_H
#define GFADRAW_H

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

  struct Segment {
    uint32_t rank;
    uint32_t tid;
    uint32_t pos;
    uint32_t len;

    Segment(uint32_t const rk, uint32_t const t, uint32_t const p, uint32_t const l) : rank(rk), tid(t), pos(p), len(l) {}
  };

  struct Link {
    bool fromrev;
    bool torev;
    uint32_t from;
    uint32_t to;

    Link(bool const fv, bool const tv, uint32_t const fr, uint32_t tos) : fromrev(fv), torev(tv), from(fr), to(tos) {}
  };

  struct Graph {
    std::vector<Segment> segments;
    std::vector<Link> links;
  };


  struct SubGraph {
    std::vector<uint32_t> segments;
    std::vector<uint32_t> links;
  };  


  
  template<typename TConfig>
  inline void
  drawNodes(TConfig const& c, cv::Mat& img, Graph const& g, SubGraph const& gsub, uint32_t const nranks, uint32_t const mnodes) {
    boost::random::mt19937 rng;
    boost::random::uniform_int_distribution<> uf(0,1000);

    uint32_t xoffset = (c.tlwidth - c.nodewidth) / 2;
    uint32_t yoffset = (c.tlheight - c.nodeheight) / 2;

    std::vector<uint32_t> xpos(nranks, 0);
    std::map<uint32_t, uint32_t> xnodemap;
    for(uint32_t i = 0; i < gsub.segments.size(); ++i) {
      uint32_t rk = g.segments[gsub.segments[i]].rank;
      xnodemap.insert(std::make_pair(gsub.segments[i], xpos[rk]));
      
      // Draw node
      cv::Rect rect(xpos[rk] * c.tlwidth + xoffset, rk * c.tlheight + yoffset, c.nodewidth, c.nodeheight);
      cv::rectangle(img, rect, WALLY_FWDMATCH, -1);

      // Draw label
      double font_scale = c.ftscale;
      double font_thickness = 3 * c.ftscale;
      int32_t baseline = 0;

      // Segment coordinate
      std::string posStr(boost::lexical_cast<std::string>(g.segments[gsub.segments[i]].pos));
      insertComma(posStr);
      std::string text(c.chrname[g.segments[gsub.segments[i]].rank][g.segments[gsub.segments[i]].tid]);
      text += ":" + posStr;
      cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_DUPLEX, font_scale, font_thickness, &baseline);
      cv::putText(img, text, cv::Point(xpos[rk] * c.tlwidth + xoffset, rk * c.tlheight + yoffset + textSize.height), cv::FONT_HERSHEY_SIMPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness);
      
      // Segment length
      text = boost::lexical_cast<std::string>(g.segments[gsub.segments[i]].len);
      insertComma(text);
      textSize = cv::getTextSize(text, cv::FONT_HERSHEY_DUPLEX, font_scale, font_thickness, &baseline);
      cv::putText(img, text, cv::Point(xpos[rk] * c.tlwidth + xoffset + c.nodewidth - textSize.width, rk * c.tlheight + yoffset + c.nodeheight - 1), cv::FONT_HERSHEY_SIMPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness);
      
      // Next tile
      ++xpos[rk];
    }

    // Draw links
    for(uint32_t i = 0; i < gsub.links.size(); ++i) {
      uint32_t from = g.links[gsub.links[i]].from;
      bool fromrev = g.links[gsub.links[i]].fromrev;
      uint32_t to = g.links[gsub.links[i]].to;
      bool torev = g.links[gsub.links[i]].torev;

      // Offset lines based on grid position in x- and y-direction
      //float sf = (((float) std::abs((int32_t) xnodemap[from] - (int32_t) xnodemap[to]) / (float) mnodes) + ((float) std::abs((int32_t) g.segments[from].rank - (int32_t) g.segments[to].rank) / (float) nranks)) / 2;
      //double sf = ((double) rand() / (RAND_MAX));
      double sf = (double) uf(rng) / 1000.0;
      std::cerr << sf << std::endl;
      float pry = sf * (yoffset * 0.8) + 1;
      float prx = sf * (xoffset * 0.8) + 1;
      // Offset outgoing edges
      float outg = sf * c.nodeheight / 2;
      
      if (fromrev) {
	// - outgoing edge
	if (torev) {
	  // - incoming edge	      
	  if (g.segments[from].rank == g.segments[to].rank) {
	    // Same rank, - - link
	    if (xnodemap[to] + 1 == xnodemap[from]) {
	      // Consecutive nodes in the same rank
	      cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset, g.segments[from].rank * c.tlheight + c.tlheight / 2), cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth, g.segments[to].rank * c.tlheight + c.tlheight / 2), cv::Scalar(0, 0, 0), c.lw);
	    } else {
	      // Non-consecutive nodes in the same rank
	      pry = pry + yoffset + c.nodeheight;
	      cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Point(xnodemap[from] * c.tlwidth + xoffset, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Scalar(0, 0, 0), c.lw);
	      cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + pry), cv::Point(xnodemap[from] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Scalar(0, 0, 0), c.lw);
	      cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + pry), cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + pry), cv::Scalar(0, 0, 0), c.lw);
	      cv::line(img, cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + pry), cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Scalar(0, 0, 0), c.lw);
	      cv::line(img, cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Scalar(0, 0, 0), c.lw);
	    }
	  } else {
	    // Different rank, - - link
	    if (g.segments[from].rank > g.segments[to].rank) outg = -1 * outg; // Above or below node midline
	    else pry = pry + yoffset + c.nodeheight;
	    cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Point(xnodemap[from] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Point(xnodemap[from] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + pry), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + pry), cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + pry), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + pry), cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Scalar(0, 0, 0), c.lw);	    
	  }
	} else {
	  // + incoming edge
	  if (g.segments[from].rank == g.segments[to].rank) {
	    // Same rank, - + link
	    std::cerr << "-+" << std::endl;
	  } else {
	    // Different rank, - + link
	    if (g.segments[from].rank > g.segments[to].rank) outg = -1 * outg; // Above or below node midline
	    else pry = pry + yoffset + c.nodeheight;
	    cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Point(xnodemap[from] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Point(xnodemap[from] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + pry), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + pry), cv::Point(xnodemap[to] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + pry), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[to] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + pry), cv::Point(xnodemap[to] * c.tlwidth + xoffset - prx, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[to] * c.tlwidth + xoffset - prx, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Point(xnodemap[to] * c.tlwidth + xoffset, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Scalar(0, 0, 0), c.lw);
	  }
	}
      } else {
	// + outgoing edge
	if (torev) {
	  // - incoming edge
	  if (g.segments[from].rank == g.segments[to].rank) {
	    // Same rank, + - link
	    std::cerr << "+-" << std::endl;
	  } else {
	    // Different rank, + - link
	    if (g.segments[from].rank > g.segments[to].rank) outg = -1 * outg; // Above or below node midline
	    else pry = pry + yoffset + c.nodeheight;
	    cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + pry), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + pry), cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + pry), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + pry), cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Point(xnodemap[to] * c.tlwidth + xoffset + c.nodewidth, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Scalar(0, 0, 0), c.lw);
	  }
	} else {
	  // + incoming edge
	  if (g.segments[from].rank == g.segments[to].rank) {
	    // Same rank, + + link
	    if (xnodemap[from] + 1 == xnodemap[to]) {
	      // Consecutive nodes in the same rank
	      cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth, g.segments[from].rank * c.tlheight + c.tlheight / 2), cv::Point(xnodemap[to] * c.tlwidth + xoffset, g.segments[to].rank * c.tlheight + c.tlheight / 2), cv::Scalar(0, 0, 0), c.lw);
	    } else {
	      // Non-consecutive nodes in the same rank
	      cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth, g.segments[from].rank * c.tlheight + c.tlheight / 2 - outg), cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + c.tlheight / 2 - outg), cv::Scalar(0, 0, 0), c.lw);
	      cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + pry), cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + c.tlheight / 2 - outg), cv::Scalar(0, 0, 0), c.lw);
	      cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + pry), cv::Point(xnodemap[to] * c.tlwidth + xoffset - prx, g.segments[to].rank * c.tlheight + pry), cv::Scalar(0, 0, 0), c.lw);
	      cv::line(img, cv::Point(xnodemap[to] * c.tlwidth + xoffset - prx, g.segments[to].rank * c.tlheight + pry), cv::Point(xnodemap[to] * c.tlwidth + xoffset - prx, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Scalar(0, 0, 0), c.lw);
	      cv::line(img, cv::Point(xnodemap[to] * c.tlwidth + xoffset - prx, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Point(xnodemap[to] * c.tlwidth + xoffset, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Scalar(0, 0, 0), c.lw);
	    }
	  } else {
	    // Different rank, + + link
	    if (g.segments[from].rank > g.segments[to].rank) outg = -1 * outg; // Above or below node midline
	    else pry = pry + yoffset + c.nodeheight;
	    cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + c.tlheight / 2 + outg), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + c.tlheight / 2 +  outg), cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + pry), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[from] * c.tlwidth + xoffset + c.nodewidth + prx, g.segments[from].rank * c.tlheight + pry), cv::Point(xnodemap[to] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + pry), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[to] * c.tlwidth + xoffset - prx, g.segments[from].rank * c.tlheight + pry), cv::Point(xnodemap[to] * c.tlwidth + xoffset - prx, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Scalar(0, 0, 0), c.lw);
	    cv::line(img, cv::Point(xnodemap[to] * c.tlwidth + xoffset - prx, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Point(xnodemap[to] * c.tlwidth + xoffset, g.segments[to].rank * c.tlheight + c.tlheight / 2 - outg), cv::Scalar(0, 0, 0), c.lw);
	  }
	}
      }
    }
      
  }



}

#endif
