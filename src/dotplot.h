#ifndef DOTPLOT_H
#define DOTPLOT_H


#include <iostream>
#include <fstream>

#include <boost/functional/hash.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/timer/timer.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/icl/split_interval_map.hpp>
#include <boost/filesystem.hpp>

#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>

#include <iostream>

#include "version.h"
#include "util.h"
#include "mapping.h"
#include "matchdraw.h"
#include "bed.h"

namespace wallysworld
{

  // Config arguments
  struct ConfigDotplot {
    bool showWindow;
    bool hasReadFile;
    bool hasRegionFile;
    bool storeSequences;
    bool flatten;
    bool flip;
    bool incSelf;
    bool refTop;
    uint32_t minMapQual;
    uint32_t matchlen;
    uint32_t topN;
    uint32_t seqsize;
    uint32_t width;
    uint32_t height;
    uint32_t usedwidth;
    uint32_t usedheight;
    uint32_t tlheight;  // pixel height of a track line
    int32_t format;
    float lw;
    double pxoffset; // 1bp in pixel
    double pyoffset; // 1bp in pixel
    std::string readStr;
    std::string regionStr;
    boost::filesystem::path seqfile;
    boost::filesystem::path readFile;
    boost::filesystem::path genome;
    boost::filesystem::path regionFile;
    boost::filesystem::path file;
  };

  inline char
  complement(char n) {
    switch(n) {
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'G':
      return 'C';
    case 'C':
      return 'G';
    }
    return 'N';
  }

  inline char
  upper(char n) {
    switch(n) {
    case 'a':
      return 'A';
    case 'c':
      return 'C';
    case 'g':
      return 'G';
    case 't':
      return 'T';
    case 'n':
      return 'N';
    }
    return n;
  }

  inline void
  upper(char* nucs) {
    while (*nucs) {
      *nucs = upper(*nucs);
      ++nucs;
    }
  }

  inline void
  revcomplement(char* nucs) {
    char* it = nucs;
    while (*it) {
      *it = complement(*it);
      ++it;
    }
    std::reverse(nucs, nucs + strlen(nucs));
  }
  
  inline uint32_t
  nuc(char const base) {
    return (int(base) & 6) >> 1;
  }

  inline uint64_t
  hashwordShort(std::string const& word) {
    uint32_t j = 0;
    uint64_t h = 0;
    for(int32_t i = word.size() - 1; i>=0; --i, ++j) {
      if (word[i] == 'N') return std::numeric_limits<uint64_t>::max();
      h += (uint64_t) nuc(word[i]) * (uint64_t) std::pow((long double) 4, j);
    }
    return h;
  }

  inline uint64_t
  hashwordLong(std::string const& word) {
    std::size_t seed = hash_string(word.c_str());
    boost::hash<std::string> string_hash;
    boost::hash_combine(seed, string_hash(word));
    //std::cerr << (uint64_t) seed << '\t' << word << std::endl;
    return seed;
  }

  template<typename TConfig, typename THashMap>
  inline void
  hashShort(TConfig const& c, char* seq, int32_t const len, THashMap& hmap, bool forward) {
    typedef typename THashMap::mapped_type TPosVec;
    const uint64_t topmult = 1ULL << (2 * (c.matchlen - 1));
    uint64_t h = 0;
    bool rewind = true;
    for(int32_t k = 0; k < (int32_t) len - (int32_t) c.matchlen + 1; ++k) {
      if (rewind) {
	h = hashwordShort(std::string(seq + k, seq + k + c.matchlen));
	if (h != std::numeric_limits<uint64_t>::max()) rewind = false;
      } else {
	if (seq[k+c.matchlen-1] == 'N') rewind = true;
	else {
	  h -= (uint64_t) nuc(seq[k - 1]) * topmult;
	  h *= 4;
	  h += (uint64_t) nuc(seq[k+c.matchlen-1]);
	}
      }
      if (!rewind) {
	//std::cerr << forward << ',' << h << ',' << k << ',' << std::string(seq + k, seq + k + c.matchlen) << std::endl;
	if (hmap.find(h) == hmap.end()) hmap.insert(std::make_pair(h, TPosVec()));
	if (forward) hmap[h].push_back(k);
	else hmap[h].push_back(len - k - c.matchlen);
      }
    }
  }

  template<typename TConfig, typename THashMap>
  inline void
  hashLong(TConfig const& c, char* seq, int32_t const len, THashMap& hmap, bool forward) {
    typedef typename THashMap::mapped_type TPosVec;
    for(int32_t k = 0; k < (int32_t) len - (int32_t) c.matchlen + 1; ++k) {
      std::string word = std::string(seq + k, seq + k + c.matchlen);
      if (word.find('N') == std::string::npos) {
	uint64_t h = hashwordLong(word);
	if (hmap.find(h) == hmap.end()) hmap.insert(std::make_pair(h, TPosVec()));
	if (forward) hmap[h].push_back(k);
	else hmap[h].push_back(len - k - c.matchlen);
      }
    }
  }

  template<typename TConfig, typename THashMap>
  inline void
  wordMatchShort(TConfig const& c, char* seq, int32_t const xlen, int32_t const ylen, THashMap& fwd, THashMap& rev, BLContext& img) {
    // Find word matches
    const uint64_t topmult = 1ULL << (2 * (c.matchlen - 1));
    uint64_t h = 0;
    bool rewind = true;
    for(int32_t k = 0; k < (int32_t) ylen - (int32_t) c.matchlen + 1; ++k) {
      if (rewind) {
	h = hashwordShort(std::string(seq + k, seq + k + c.matchlen));
	if (h != std::numeric_limits<uint64_t>::max()) rewind = false;
      } else {
	if (seq[k+c.matchlen-1] == 'N') rewind = true;
	else {
	  h -= (uint64_t) nuc(seq[k - 1]) * topmult;
	  h *= 4;
	  h += (uint64_t) nuc(seq[k+c.matchlen-1]);
	}
      }
      if (!rewind) {
	// Forward matches
	if (fwd.find(h) != fwd.end()) {
	  for(uint32_t idx = 0; idx < fwd[h].size(); ++idx) {
	    int32_t px = pixelX(c.usedwidth, xlen, fwd[h][idx]);
	    int32_t pxend = pixelX(c.usedwidth, xlen, fwd[h][idx] + c.matchlen);
	    int32_t py = pixelX(c.usedheight, ylen, k);
	    int32_t pyend = pixelX(c.usedheight, ylen, k + c.matchlen);
	    drawLine(img, px, py, pxend, pyend, BLRgba32(0, 0, 0), c.lw);
	  }
	}
	// Reverse matches
	if (rev.find(h) != rev.end()) {
	  for(uint32_t idx = 0; idx < rev[h].size(); ++idx) {
	    int32_t px = pixelX(c.usedwidth, xlen, rev[h][idx] + c.matchlen);
	    int32_t pxend = pixelX(c.usedwidth, xlen, rev[h][idx]);
	    int32_t py = pixelX(c.usedheight, ylen, k);
	    int32_t pyend = pixelX(c.usedheight, ylen, k + c.matchlen);
	    drawLine(img, px, py, pxend, pyend, BLRgba32(255, 0, 0), c.lw);
	  }
	}
      }
    }
  }
  

  template<typename TConfig, typename THashMap>
  inline void
  wordMatchLong(TConfig const& c, char* seq, int32_t const xlen, int32_t const ylen, THashMap& fwd, THashMap& rev, BLContext& img) {
    // Find word matches
    for(int32_t k = 0; k < (int32_t) ylen - (int32_t) c.matchlen + 1; ++k) {
      std::string word = std::string(seq + k, seq + k + c.matchlen);
      if (word.find('N') == std::string::npos) {
	uint64_t h = hashwordLong(word);
	// Forward matches
	if (fwd.find(h) != fwd.end()) {
	  for(uint32_t idx = 0; idx < fwd[h].size(); ++idx) {
	    int32_t px = pixelX(c.usedwidth, xlen, fwd[h][idx]);
	    int32_t pxend = pixelX(c.usedwidth, xlen, fwd[h][idx] + c.matchlen);
	    int32_t py = pixelX(c.usedheight, ylen, k);
	    int32_t pyend = pixelX(c.usedheight, ylen, k + c.matchlen);
	    drawLine(img, px, py, pxend, pyend, BLRgba32(0, 0, 0), c.lw);
	  }
	}
	// Reverse matches
	if (rev.find(h) != rev.end()) {
	  for(uint32_t idx = 0; idx < rev[h].size(); ++idx) {
	    int32_t px = pixelX(c.usedwidth, xlen, rev[h][idx] + c.matchlen);
	    int32_t pxend = pixelX(c.usedwidth, xlen, rev[h][idx]);
	    int32_t py = pixelX(c.usedheight, ylen, k);
	    int32_t pyend = pixelX(c.usedheight, ylen, k + c.matchlen);
	    drawLine(img, px, py, pxend, pyend, BLRgba32(255, 0, 0), c.lw);
	  }
	}
      }
    }
  }

  template<typename TConfig, typename TReadMapping>
  inline void
  drawXMappings(TConfig const& c, bam_hdr_t* hdr, std::string const& refname, uint32_t const len, std::map<uint32_t, BLRgba32>& cm, TReadMapping& mp, BLContext& img) {
    if (mp.find(refname) != mp.end()) {
      uint32_t runspacer = 4 * c.tlheight;
      for(uint32_t k = 0; k < mp[refname].size(); ++k) {
	int32_t px = pixelX(c.usedwidth, len, mp[refname][k].rstart);
	int32_t pxend = pixelX(c.usedwidth, len, mp[refname][k].rend);
	img.fill_rect(BLRectI(px, c.usedheight + runspacer, pxend - px, c.tlheight), cm[mp[refname][k].tid]);

	if (c.flatten) {
	  drawLine(img, px, c.usedheight + runspacer, px, c.usedheight + runspacer + c.tlheight, BLRgba32(0, 0, 0), 1);
	  drawLine(img, pxend, c.usedheight + runspacer, pxend, c.usedheight + runspacer + c.tlheight, BLRgba32(0, 0, 0), 1);
	} else {
	  std::string text = std::string(hdr->target_name[mp[refname][k].tid]);
	  std::string gstart = boost::lexical_cast<std::string>(mp[refname][k].gstart);
	  insertComma(gstart);
	  std::string gend = boost::lexical_cast<std::string>(mp[refname][k].gend);
	  insertComma(gend);
	  text += ":" + gstart + "-" + gend;
	  if (mp[refname][k].fwd) text += " -->";
	  else text += " <--";
	  TextSize textSize = getTextSize(text);
	  if ((pxend - px) > textSize.width + 10) {
	    drawText(img, (px + pxend) / 2 - textSize.width / 2, c.usedheight + runspacer + textSize.height, text, BLRgba32(0, 0, 0));
	  } else if (px > (int32_t) c.usedwidth / 2) {
	    drawText(img, px - textSize.width - 5, c.usedheight + runspacer + textSize.height, text, BLRgba32(0, 0, 0));
	  } else {
	    drawText(img, pxend + 5, c.usedheight + runspacer + textSize.height, text, BLRgba32(0, 0, 0));
	  }

	  // Next mapping
	  runspacer += c.tlheight;
	}
      }
    }
  }

  template<typename TConfig, typename TReadMapping>
  inline void
  drawYMappings(TConfig const& c, bam_hdr_t* hdr, std::string const& refname, uint32_t const len, std::map<uint32_t, BLRgba32>& cm, TReadMapping& mp, BLContext& img) {
    if (mp.find(refname) != mp.end()) {
      uint32_t runspacer = 4 * c.tlheight;
      for(uint32_t k = 0; k < mp[refname].size(); ++k) {
	int32_t py = pixelX(c.usedheight, len, mp[refname][k].rstart);
	int32_t pyend = pixelX(c.usedheight, len, mp[refname][k].rend);
	img.fill_rect(BLRectI(c.usedwidth + runspacer, py, c.tlheight, pyend - py), cm[mp[refname][k].tid]);

	if (c.flatten) {
	  drawLine(img, c.usedwidth + runspacer, py, c.usedwidth + runspacer + c.tlheight, py, BLRgba32(0, 0, 0), 1);
	  drawLine(img, c.usedwidth + runspacer, pyend, c.usedwidth + runspacer + c.tlheight, pyend, BLRgba32(0, 0, 0), 1);
	} else {
	  std::string text = std::string(hdr->target_name[mp[refname][k].tid]);
	  std::string gstart = boost::lexical_cast<std::string>(mp[refname][k].gstart);
	  insertComma(gstart);
	  std::string gend = boost::lexical_cast<std::string>(mp[refname][k].gend);
	  insertComma(gend);
	  text += ":" + gstart + "-" + gend;
	  if (mp[refname][k].fwd) text += " -->";
	  else text += " <--";
	  TextSize textSize = getTextSize(text);

	  // Rotated text: anchor x is the label column, anchor y is the along-axis start
	  if ((pyend - py) > textSize.width + 10) {
	    drawTextRotated(img, c.usedwidth + runspacer, (py + pyend) / 2 - textSize.width / 2, text, BLRgba32(0, 0, 0));
	  } else if (py > (int32_t) c.usedheight / 2) {
	    drawTextRotated(img, c.usedwidth + runspacer, py - textSize.width - 5, text, BLRgba32(0, 0, 0));
	  } else {
	    drawTextRotated(img, c.usedwidth + runspacer, pyend + 5, text, BLRgba32(0, 0, 0));
	  }

	  // Next mapping
	  runspacer += c.tlheight;
	}
      }
    }
  }


  template<typename TConfig>
  inline void
  drawXScaleDotplot(TConfig const& c, std::string const& refname, uint32_t const len, BLContext& img) {
    std::string text(boost::lexical_cast<std::string>(len));
    insertComma(text);
    TextSize textSize = getTextSize(text);

    // Find suitable tick size
    uint32_t modval = findTicks(c.pxoffset, textSize.width);

    // Scale line
    drawLine(img, 0, c.usedheight, c.usedwidth, c.usedheight, BLRgba32(0, 0, 255), 2);

    // Ticks
    double px = 0;
    for(uint32_t i = 0; i < len; ++i) {
      if (i % modval == 0) {
	drawLine(img, px - c.pxoffset/2, c.usedheight, px - c.pxoffset/2, c.usedheight + 0.75 * c.tlheight, BLRgba32(0, 0, 255), 2);
      }
      if (i % modval == 0) {
	// Font
	text = boost::lexical_cast<std::string>(i);
	insertComma(text);
	TextSize textSize = getTextSize(text);
	if ((px - c.pxoffset/2 - textSize.width/2 > 0) && (px - c.pxoffset/2 + textSize.width < c.usedwidth)) {
	  drawText(img, px - c.pxoffset/2 - textSize.width/2, c.usedheight + c.tlheight + textSize.height, text, BLRgba32(0, 0, 0));
	}
      }
      px += c.pxoffset;
    }

    // Contig name
    if (true) {
      int32_t midpoint = c.usedwidth / 2;
      TextSize textSize = getTextSize(refname);
      drawText(img, midpoint - textSize.width/2, c.usedheight + 2 * c.tlheight + textSize.height, refname, BLRgba32(0, 0, 0));
    }
  }


  // Reference to top
  template<typename TConfig>
  inline void
  drawXScaleDotplotTop(TConfig const& c, std::string const& refname, uint32_t const len, BLContext& img) {
    std::string text(boost::lexical_cast<std::string>(len));
    insertComma(text);
    TextSize textSize = getTextSize(text);
    uint32_t modval = findTicks(c.pxoffset, textSize.width);
    drawLine(img, 0, 0, c.usedwidth, 0, BLRgba32(0, 0, 255), 2);

    double px = 0;
    for(uint32_t i = 0; i < len; ++i) {
      if (i % modval == 0) {
	drawLine(img, px - c.pxoffset/2, -0.55 * c.tlheight, px - c.pxoffset/2, 0, BLRgba32(0, 0, 255), 2);
	std::string txt = boost::lexical_cast<std::string>(i);
	insertComma(txt);
	TextSize ts = getTextSize(txt);
	if ((px - c.pxoffset/2 - ts.width/2 > 0) && (px - c.pxoffset/2 + ts.width < c.usedwidth)) {
	  drawText(img, px - c.pxoffset/2 - ts.width/2, -0.7 * c.tlheight, txt, BLRgba32(0, 0, 0));
	}
      }
      px += c.pxoffset;
    }

    // Contig name
    TextSize ns = getTextSize(refname);
    drawText(img, c.usedwidth / 2 - ns.width / 2, -0.7 * c.tlheight - ns.height - 4, refname, BLRgba32(0, 0, 0));
  }

  template<typename TConfig>
  inline void
  drawYScaleDotplot(TConfig const& c, std::string const& refname, uint32_t const len, BLContext& img) {
    std::string text(boost::lexical_cast<std::string>(len));
    insertComma(text);
    TextSize textSize = getTextSize(text);

    // Find suitable tick size
    uint32_t modval = findTicks(c.pyoffset, textSize.width);

    // Scale line
    drawLine(img, c.usedwidth, 0, c.usedwidth, c.usedheight, BLRgba32(0, 0, 255), 2);

    // Ticks
    double py = 0;
    for(uint32_t i = 0; i < len; ++i) {
      if (i % modval == 0) {
	drawLine(img, c.usedwidth, py - c.pyoffset/2, c.usedwidth + 0.75 * c.tlheight, py - c.pyoffset/2, BLRgba32(0, 0, 255), 2);
      }
      if (i % modval == 0) {
	// Font
	text = boost::lexical_cast<std::string>(i);
	insertComma(text);
	TextSize textSize = getTextSize(text);
	if ((py - c.pyoffset/2 - textSize.width/2 > 0) && (py - c.pyoffset/2 + textSize.width < c.usedheight)) {
	  drawTextRotated(img, c.usedwidth + c.tlheight + textSize.height, py - c.pyoffset/2 - textSize.width/2, text, BLRgba32(0, 0, 0));
	}
      }
      py += c.pyoffset;
    }

    // Contig name
    if (true) {
      TextSize textSize = getTextSize(refname);
      int32_t midpoint = c.usedheight / 2;
      drawTextRotated(img, c.usedwidth + 2 * c.tlheight + textSize.height, midpoint - textSize.width/2, refname, BLRgba32(0, 0, 0));
    }
  }

  
  template<typename TConfigStruct>
  inline int dotplotRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("wally.prof");
#endif
    // Chromosome colors
    BLRgba32 colors[12] = { BLRgba32(51,160,44), BLRgba32(31,120,180), BLRgba32(227,26,28), BLRgba32(255,255,153), BLRgba32(166,206,227), BLRgba32(178,223,138), BLRgba32(251,154,153), BLRgba32(253,191,111), BLRgba32(255,127,0), BLRgba32(202,178,214), BLRgba32(106,61,154), BLRgba32(177,89,40) };

    // Open file handles
    samFile* samfile = NULL;
    bam_hdr_t* hdr = NULL;
    if (c.format == 0) {
      samfile = sam_open(c.file.string().c_str(), "r");
      hts_set_fai_filename(samfile, c.genome.string().c_str());
      hdr = sam_hdr_read(samfile);
    }
    
    // Read mappings
    typedef std::vector<Mapping> TMappings;
    typedef std::map<std::string, TMappings > TReadMappings;
    TReadMappings mp;

    // Make sure the FASTA index is gone
    boost::filesystem::remove(c.seqfile.string());
    boost::filesystem::remove(c.seqfile.string() + ".fai");
    
    // Parse regions
    std::vector<Region> rg; // tid = FASTA index
    std::vector<Region> scanRg; // tid = BAM index
    uint32_t rgcount = 0;
    if ((!c.regionStr.empty()) || (c.hasRegionFile)) {
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Parse regions." << std::endl;
      faidx_t* fai = fai_load(c.genome.string().c_str());
      typedef std::map<std::string, std::pair<int32_t, int32_t> > TChrIdLenMap;
      TChrIdLenMap chrMap;
      for(int32_t refIndex=0; refIndex < faidx_nseq(fai); ++refIndex) {
	std::string tname(faidx_iseq(fai, refIndex));
	chrMap.insert(std::make_pair(tname, std::make_pair(refIndex, faidx_seq_len(fai, tname.c_str()))));
      }
      if (!parseRegions(chrMap, c, rg)) { fai_destroy(fai); return 1; }
      for(uint32_t i = 0; i < rg.size(); ++i) {
	std::string chrn(faidx_iseq(fai, rg[i].tid));
	std::string autoId = chrn + "_" + boost::lexical_cast<std::string>(rg[i].beg + 1) + "_" + boost::lexical_cast<std::string>(rg[i].end);
	if (rg[i].id == autoId) rg[i].id = chrn + ":" + boost::lexical_cast<std::string>(rg[i].beg + 1) + "-" + boost::lexical_cast<std::string>(rg[i].end);
	if ((c.format == 0) && (hdr != NULL)) {
	  int32_t bt = bam_name2id(hdr, chrn.c_str());
	  if (bt >= 0) scanRg.push_back(Region(bt, rg[i].beg, rg[i].end));
	}
      }
      fai_destroy(fai);
    }

    // Parse BAM
    if (c.format == 0) {
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Parse reads." << std::endl;
      typedef std::set<std::string> TReadSet;
      TReadSet reads;
      if (c.topN > 0) {
	if (scanRg.empty()) {
	  std::cerr << "Longest read selection requires a region." << std::endl;
	  return 1;
	}
	selectTopReads(c, scanRg, reads);
      } else if (c.hasReadFile) _parseReads(c, reads);
      else reads.insert(c.readStr);

      // Get read mappings
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Extract read mappings." << std::endl;
      if (!scanRg.empty()) mappingsRegion(c, reads, scanRg, mp);
      else mappings(c, reads, mp);

      // Sort mappings
      for(TReadMappings::iterator it = mp.begin(); it != mp.end(); ++it) std::sort(it->second.begin(), it->second.end());

      // Check size
      if (mp.empty()) {
	std::cerr << "No reads found!" << std::endl;
	std::cerr << "Please check your read names and BAM file!" << std::endl;
	return 1;
      }
    } else if (c.format == 1) {
      // Copy sequence file to append possible reference regions
      boost::filesystem::copy_file(c.file, c.seqfile, boost::filesystem::copy_options::overwrite_existing);
    }

    // Extract region FASTA (append after the read sequences)
    if (!rg.empty()) {
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Extract FASTA for regions." << std::endl;
      faidx_t* fai = fai_load(c.genome.string().c_str());
      std::ofstream sfile;
      sfile.open(c.seqfile.string().c_str(), std::ios_base::app);
      for(int32_t refIndex=0; refIndex < faidx_nseq(fai); ++refIndex) {
	std::string tname(faidx_iseq(fai, refIndex));
	for(uint32_t rgIdx = 0; rgIdx < rg.size(); ++rgIdx) {
	  if (rg[rgIdx].tid == refIndex) {
	    // Load only the region slice
	    int32_t seqlen = -1;
	    char* seq = faidx_fetch_seq(fai, tname.c_str(), rg[rgIdx].beg, rg[rgIdx].end - 1, &seqlen);
	    if (seq != NULL) {
	      sfile << ">" << rg[rgIdx].id << std::endl;
	      sfile << boost::to_upper_copy(std::string(seq)) << std::endl;
	      free(seq);
	      ++rgcount;
	    }
	  }
	}
      }
      sfile.close();
      fai_destroy(fai);
    }
    
    if ((!boost::filesystem::exists(c.seqfile)) || (boost::filesystem::file_size(c.seqfile) == 0)) {
      std::cerr << "Error: No sequences available to plot." << std::endl;
      if (c.format == 0) {
	bam_hdr_destroy(hdr);
	sam_close(samfile);
      }
      return 1;
    }

    // Flip x and y, reverse file
    if (c.flip) {
      faidx_t* fai = fai_load(c.seqfile.c_str());
      int32_t seqpointer = faidx_nseq(fai) - 1;
      std::vector<std::string> seqstore(faidx_nseq(fai));
      std::vector<std::string> seqname(faidx_nseq(fai));
      for(int32_t idx = 0; idx < faidx_nseq(fai); ++idx, --seqpointer) {
	seqname[seqpointer] = std::string(faidx_iseq(fai, idx));
	int32_t seqlen = faidx_seq_len(fai, seqname[seqpointer].c_str());
	int32_t sl = 0;
	char* seq = faidx_fetch_seq(fai, seqname[seqpointer].c_str(), 0, seqlen, &sl);
	seqstore[seqpointer] = std::string(seq);
	free(seq);
      }
      fai_destroy(fai);

      // Clean-up
      boost::filesystem::remove(c.seqfile.string());
      boost::filesystem::remove(c.seqfile.string() + ".fai");

      std::ofstream sfile;
      sfile.open(c.seqfile.string().c_str());
      for(uint32_t idx = 0; idx < seqname.size(); ++idx) {
	sfile << ">" << seqname[idx] << std::endl;
	sfile << seqstore[idx] << std::endl;
      }
      sfile.close();
    }

    // Load sequences from disk for large contigs
    faidx_t* fai = fai_load(c.seqfile.c_str());
    if (fai == NULL) {
      std::cerr << "Error: Unable to index the extracted sequences (" << c.seqfile.string() << ")." << std::endl;
      if (c.format == 0) {
	bam_hdr_destroy(hdr);
	sam_close(samfile);
      }
      return 1;
    }
    int32_t nseq = faidx_nseq(fai);
    bool transpose = (c.refTop && rgcount && (c.format == 0));
    int32_t seqend = nseq;
    if (!c.incSelf) seqend -= 1;
    if (rgcount) seqend = nseq - rgcount;
    int32_t i1beg = 0, i1end = seqend;
    if (transpose) {
      i1beg = nseq - rgcount;
      i1end = nseq;
    }
    int32_t plotCount = 0;
    for(int32_t idx1 = i1beg; idx1 < i1end; ++idx1) {
      std::string seqname1(faidx_iseq(fai, idx1));
      int32_t xlen = faidx_seq_len(fai, seqname1.c_str());
      if (xlen < (int32_t) c.seqsize) continue;
      int32_t sl = 0;
      char* seq1 = faidx_fetch_seq(fai, seqname1.c_str(), 0, xlen, &sl);
      upper(seq1);
      
      // Hash fwd and rev words
      typedef std::vector<uint32_t> TPosVec;
      typedef std::map<uint64_t, TPosVec> TMatchMap;
      TMatchMap fwd;
      if (c.matchlen < 32) hashShort(c, seq1, xlen, fwd, true);
      else hashLong(c, seq1, xlen, fwd, true);
      revcomplement(seq1);
      TMatchMap rev;
      if (c.matchlen < 32) hashShort(c, seq1, xlen, rev, false);
      else hashLong(c, seq1, xlen, rev, false);
      
      // Match 2nd sequence
      int32_t seqoffset = idx1;
      if (!c.incSelf) seqoffset += 1;
      if (rgcount) seqoffset = nseq - rgcount;
      int32_t i2beg = seqoffset, i2end = nseq;
      if (transpose) {
	i2beg = 0;
	i2end = nseq - rgcount;
      }
      for(int32_t idx2 = i2beg; idx2 < i2end; ++idx2) {
	std::string seqname2(faidx_iseq(fai, idx2));
	int32_t ylen = faidx_seq_len(fai, seqname2.c_str());
	if (ylen < (int32_t) c.seqsize) continue;
	sl = 0;
	char* seq2 = faidx_fetch_seq(fai, seqname2.c_str(), 0, ylen, &sl);
	upper(seq2);
	
	// Next plot
	++plotCount;
	std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Plot for " << seqname1 << " and " << seqname2 << std::endl;

	// Estimate image size
	if (transpose) {
	  c.usedwidth = (c.width > 0) ? (int32_t) c.width : 1000;
	  double hh = (double) ylen / (double) xlen * (double) c.usedwidth;
	  if (hh < 1) hh = 1;
	  if (hh > 2.0 * c.usedwidth) hh = 2.0 * c.usedwidth;
	  c.usedheight = (int32_t) hh;
	} else if ((c.width == 0) && (c.height == 0)) {
	  int32_t maxlen = std::max(xlen, ylen);
	  c.usedwidth = std::max(1, (int32_t) ((double) xlen / (double) maxlen * 1024.0));
	  c.usedheight = std::max(1, (int32_t) ((double) ylen / (double) maxlen * 1024.0));
	} else if (c.width == 0) {
	  c.usedwidth = std::max(1, (int32_t) ((double) xlen / (double) ylen * c.height));
	} else if (c.height == 0) {
	  c.usedheight = std::max(1, (int32_t) ((double) ylen / (double) xlen * c.width));
	}

	// Create image
	c.pxoffset = (1.0 / (double) xlen) * (double) (c.usedwidth);
	c.pyoffset = (1.0 / (double) ylen) * (double) (c.usedheight);
	if (transpose) {
	  uint32_t topHeader = 3 * c.tlheight;
	  uint32_t footer = 4 * c.tlheight;
	  uint32_t margin = 4 * c.tlheight;
	  if (c.flatten) margin += c.tlheight;
	  else if (mp.find(seqname2) != mp.end()) margin += (c.tlheight * (mp[seqname2].size() + 1));
	  BLImage img(c.usedwidth + margin, topHeader + c.usedheight + footer, BL_FORMAT_PRGB32);
	  BLContext ctx(img);
	  ctx.fill_all(BLRgba32(255, 255, 255));
	  ctx.translate(0, topHeader);
	  if (c.matchlen < 32) wordMatchShort(c, seq2, xlen, ylen, fwd, rev, ctx);
	  else wordMatchLong(c, seq2, xlen, ylen, fwd, rev, ctx);
	  drawXScaleDotplotTop(c, seqname1, xlen, ctx);
	  drawYScaleDotplot(c, seqname2, ylen, ctx);
	  typedef std::map<uint32_t, BLRgba32> TColorMap;
	  TColorMap cm;
	  uint32_t cidx = 0;
	  if (mp.find(seqname2) != mp.end()) {
	    for(uint32_t k = 0; k < mp[seqname2].size(); ++k) {
	      if (cm.find(mp[seqname2][k].tid) == cm.end()) { cm[mp[seqname2][k].tid] = colors[cidx]; cidx = ((cidx + 1) % 12); }
	    }
	  }
	  drawYMappings(c, hdr, seqname2, ylen, cm, mp, ctx);
	  ctx.end();
	  std::string outfile = seqname2 + "_" + seqname1 + ".png";
	  img.write_to_file(outfile.c_str());
	} else {
	  uint32_t footer = 4 * c.tlheight;
	  if (c.format == 0) {
	    if (c.flatten) footer += c.tlheight;
	    else { if (mp.find(seqname1) != mp.end()) footer += (c.tlheight * (mp[seqname1].size() + 1)); }
	  }
	  uint32_t margin = 4 * c.tlheight;
	  if (c.format == 0) {
	    if (c.flatten) margin += c.tlheight;
	    else { if (mp.find(seqname2) != mp.end()) margin += (c.tlheight * (mp[seqname2].size() + 1)); }
	  }
	  BLImage img(c.usedwidth + margin, c.usedheight + footer, BL_FORMAT_PRGB32);
	  BLContext ctx(img);
	  ctx.fill_all(BLRgba32(255, 255, 255));
	  if (c.matchlen < 32) wordMatchShort(c, seq2, xlen, ylen, fwd, rev, ctx);
	  else wordMatchLong(c, seq2, xlen, ylen, fwd, rev, ctx);
	  drawXScaleDotplot(c, seqname1, xlen, ctx);
	  drawYScaleDotplot(c, seqname2, ylen, ctx);
	  if (c.format == 0) {
	    typedef std::map<uint32_t, BLRgba32> TColorMap;
	    TColorMap cm;
	    uint32_t cidx = 0;
	    if (mp.find(seqname1) != mp.end()) {
	      for(uint32_t k = 0; k < mp[seqname1].size(); ++k) { if (cm.find(mp[seqname1][k].tid) == cm.end()) { cm[mp[seqname1][k].tid] = colors[cidx]; cidx = ((cidx + 1) % 12); } }
	    }
	    if (mp.find(seqname2) != mp.end()) {
	      for(uint32_t k = 0; k < mp[seqname2].size(); ++k) { if (cm.find(mp[seqname2][k].tid) == cm.end()) { cm[mp[seqname2][k].tid] = colors[cidx]; cidx = ((cidx + 1) % 12); } }
	    }
	    drawXMappings(c, hdr, seqname1, xlen, cm, mp, ctx);
	    drawYMappings(c, hdr, seqname2, ylen, cm, mp, ctx);
	  }
	  ctx.end();
	  std::string outfile = seqname1 + "_" + seqname2 + ".png";
	  img.write_to_file(outfile.c_str());
	}

	// Clean-up
	free(seq2);
      }
      free(seq1);
    }
    if (!plotCount) std::cerr << "Warning: No dotplots created. Please provide a second sequence or a region to compare against or use self-alignment." << std::endl;
    fai_destroy(fai);
    if (c.format == 0) {
      bam_hdr_destroy(hdr);
      sam_close(samfile);
    }
    
#ifdef PROFILE
    ProfilerStop();
#endif

    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    return 0;
  }


  int dotplot(int argc, char **argv) {
    ConfigDotplot c;
    c.tlheight = textSize().height + 2;

    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("matchlen,m", boost::program_options::value<uint32_t>(&c.matchlen)->default_value(21), "default match length")
      ("size,s", boost::program_options::value<uint32_t>(&c.seqsize)->default_value(0), "min. sequence size to include")
      ("seqfile,q", boost::program_options::value<boost::filesystem::path>(&c.seqfile), "sequence output file [optional]")      
      ("selfalign,a", "incl. self alignments")
      ("flip,p", "flip x and y-axis")
      ;

    boost::program_options::options_description bammod("BAM mode");
    bammod.add_options()
      ("read,r", boost::program_options::value<std::string>(&c.readStr), "read to display")
      ("rfile,R", boost::program_options::value<boost::filesystem::path>(&c.readFile), "file with reads to display")
      ("map-qual,u", boost::program_options::value<uint32_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("topN,t", boost::program_options::value<uint32_t>(&c.topN)->default_value(0), "plot the N longest reads in a region")
      ("region,e", boost::program_options::value<std::string>(&c.regionStr), "region to display [chrA:35-78]")
      ("reglist,E", boost::program_options::value<boost::filesystem::path>(&c.regionFile), "BED file with regions to display")
      ("flatten,f", "flatten mapping segments")
      ;

    
    boost::program_options::options_description disc("Display options");
    disc.add_options()
      ("linewidth,l", boost::program_options::value<float>(&c.lw)->default_value(1.5), "match line width")
      ("width,x", boost::program_options::value<uint32_t>(&c.width)->default_value(0), "width of the plot [0: best fit]")
      ("height,y", boost::program_options::value<uint32_t>(&c.height)->default_value(0), "height of the plot [0: best fit]")
      ("reftop,T", "transposed layout")
      ;
    
    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.file), "input file")
      ("window,w", "show window")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    // Set the visibility
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(bammod).add(disc).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(bammod).add(disc);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << std::endl;
      std::cout << "Usage: wally " << argv[0] << " [OPTIONS] -g <ref.fa> -R <reads.lst> <sample.bam>" << std::endl;
      std::cout << "       wally " << argv[0] << " [OPTIONS] <sample.fa>" << std::endl;
      std::cout << visible_options << "\n";
      return 0;
    }

    // Show window?
    if (vm.count("window")) c.showWindow = true;
    else c.showWindow = false;

    // Include self alignment
    if (vm.count("selfalign")) c.incSelf = true;
    else c.incSelf = false;

    // Read file
    if (vm.count("rfile")) {
      if (!(boost::filesystem::exists(c.readFile) && boost::filesystem::is_regular_file(c.readFile) && boost::filesystem::file_size(c.readFile))) {
	std::cerr << "File with list of reads is missing: " << c.readFile.string() << std::endl;
	return 1;
      }
      c.hasReadFile = true;
    } else c.hasReadFile = false;

    // Region file
    if (vm.count("reglist")) {
      if (!(boost::filesystem::exists(c.regionFile) && boost::filesystem::is_regular_file(c.regionFile) && boost::filesystem::file_size(c.regionFile))) {
	std::cerr << "BED file with regions is missing: " << c.regionFile.string() << std::endl;
	return 1;
      }
      c.hasRegionFile = true;
    } else c.hasRegionFile = false;

    // Flatten mappings
    if (vm.count("flatten")) c.flatten = true;
    else c.flatten = false;

    // Flip x and y-axis
    if (vm.count("flip")) c.flip = true;
    else c.flip = false;

    // Transposed layout
    if (vm.count("reftop")) c.refTop = true;
    else c.refTop = false;

    // Flip does not work with regions
    if ((c.flip) && ((!c.regionStr.empty()) || (c.hasRegionFile))) {
      std::cerr << "Error: Axis flip cannot be combined with regions." << std::endl;
      return 1;
    }

    // Input format
    c.format = inputType(c.file.string());
    if (c.format == -1) {
      std::cerr << "Unknown input file format!" << std::endl;
    } else if (c.format == 0) {
      // BAM format
      
      // Genome check
      if (vm.count("genome")) {
	if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
	  std::cerr << "Genome file is missing: " << c.genome.string() << std::endl;
	  return 1;
	}
      } else {
	std::cerr << "BAM input requires -g command-line option." << std::endl;
	return 1;
      }

      // Read selection
      bool haveRegion = (vm.count("region") || vm.count("reglist"));
      if (c.topN > 0) {
	if (!haveRegion) {
	  std::cerr << "Read selection requires a region." << std::endl;
	  return 1;
	}
      } else if (vm.count("rfile")) {
	if (!(boost::filesystem::exists(c.readFile) && boost::filesystem::is_regular_file(c.readFile) && boost::filesystem::file_size(c.readFile))) {
	  std::cerr << "File with list of reads is missing: " << c.readFile.string() << std::endl;
	  return 1;
	}
      } else if (!vm.count("read")) {
	if (haveRegion) {
	  c.topN = 10;
	} else {
	  std::cerr << "BAM input requires a read list or a region to auto-select the longest reads." << std::endl;
	  return 1;
	}
      }

      // Store sequences
      c.storeSequences = true;
    } else {
      // Store sequences
      c.storeSequences = false;
    }

    // In case of no automatic estimation
    c.usedwidth = c.width;
    c.usedheight = c.height;

    // Temporary seqfile
    if (!(vm.count("seqfile"))) c.seqfile = boost::filesystem::unique_path().replace_extension(".fa");
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "wally ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    // generate dot plots
    int32_t retVal = dotplotRun(c);

    // Clean-up temporary files
    if (!(vm.count("seqfile"))) {
      boost::filesystem::remove(c.seqfile.string());
      boost::filesystem::remove(c.seqfile.string() + ".fai");
    }

    return retVal;
  }

}

#endif
