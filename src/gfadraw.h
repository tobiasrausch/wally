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

    Segment(uint32_t const rk, uint32_t const t, uint32_t const p) : rank(rk), tid(t), pos(p) {}
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
    std::vector<bool> segments;
    std::vector<bool> links;

    SubGraph(Graph const& g) {
      segments.resize(g.segments.size(), false);
      links.resize(g.links.size(), false);
    }
  };  


  
  template<typename TConfig>
  inline void
  drawNodes(TConfig const& c, cv::Mat& img, Graph const& g, SubGraph const& gsub, uint32_t const nranks) {
    uint32_t xoffset = (c.tlwidth - c.nodewidth) / 2;
    uint32_t yoffset = (c.tlheight - c.nodeheight) / 2;

    std::vector<uint32_t> xpos(nranks, 0);
    for(uint32_t i = 0; i < g.segments.size(); ++i) {
      if (gsub.segments[i]) {
	uint32_t rk = g.segments[i].rank;
	if (rk == POS_UNDEF) rk = nranks - 1;

	// Draw node
	cv::Rect rect(xpos[rk] * c.tlwidth + xoffset, rk * c.tlheight + yoffset, c.nodewidth, c.nodeheight);    
	cv::rectangle(img, rect, WALLY_FWDMATCH, -1);

	// Next tile
	++xpos[rk];
      }
    }
  }



}

#endif
