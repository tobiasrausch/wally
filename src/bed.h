#ifndef BED_H
#define BED_H

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
#include <htslib/tbx.h>

#include "util.h"

namespace wallysworld
{

   struct Transcript {
     bool forward;
     Region rg;
     
     Transcript() : forward(true), rg(Region()) {}
     Transcript(Region const& irg, bool const fwd) : forward(fwd), rg(irg) {}
  };


  template<typename TConfig>
  inline bool
  parseAnnotation(bam_hdr_t* hdr, TConfig const& c, Region const& irg, std::vector<Transcript>& tr, std::vector<Region>& itv) {
    if (!c.hasAnnotationFile) return true;
    
    // Load BED annotation file
    htsFile* bedfile = hts_open(c.bedFile.string().c_str(), "r");
    tbx_t* tbx = tbx_index_load(c.bedFile.string().c_str());
    hts_itr_t *itr = tbx_itr_querys(tbx, convertToStr(hdr, irg).c_str());
    kstring_t str = {0,0,0};
    while (tbx_itr_next(bedfile, tbx, itr, &str) >= 0) {
      std::string line = std::string(str.s);
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep("\t");
      Tokenizer tokens(line, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter!=tokens.end()) {
	std::string chrName=*tokIter++;
	if (tokIter != tokens.end()) {
	  std::string start = *tokIter++;
	  if (tokIter != tokens.end()) {
	    std::string end = *tokIter++;
	    std::string str = chrName + ":" + start + "-" + end;
	    if (tokIter != tokens.end()) str += ":" + std::string(*tokIter++);
	    Region tmp;
	    if (!parseRegion(c, hdr, str, tmp)) return false;
	    bool transcript = false;
	    bool fwd = true;
	    if (tokIter != tokens.end()) {
	      std::string column5 = *tokIter++;
	      if (column5 == "transcript") {
		transcript = true;
		if (tokIter != tokens.end()) {
		  std::string strand = *tokIter++;
		  if (strand == "+") fwd = true;
		  else if (strand == "-") fwd = false;
		  else {
		    std::cerr << "Transcript needs forward (+) or reverse (-) in the 6th column!" << std::endl;
		    return false;
		  }
		} else {
		  std::cerr << "Transcript lacks forward (+) or reverse (-) in the 6th column!" << std::endl;
		  return false;
		}
	      } else if ((column5.size()>1) && (column5[0] == '0') && (column5[1] == 'x')) {
		tmp.color = hexColor(column5);
	      }
	    }
	    if (transcript) tr.push_back(Transcript(tmp, fwd));
	    else itv.push_back(tmp);
	  }
	}
      }
    }
    if (itr) hts_itr_destroy(itr);
    tbx_destroy(tbx);
    hts_close(bedfile);
    
    return true;
  }

}

#endif
