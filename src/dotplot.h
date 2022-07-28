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

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

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
    bool storeSequences;
    uint32_t matchlen;
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
    boost::filesystem::path seqfile;
    boost::filesystem::path readFile;
    boost::filesystem::path genome;
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
    case 'N':
      return 'N';
    case 'a':
      return 'T';
    case 't':
      return 'A';
    case 'g':
      return 'C';
    case 'c':
      return 'G';
    case 'n':
      return 'N';
    }   
    return 'N';
  }   

  inline void
  revcomplement(char* nucs) {
    for(uint32_t i = 0; i < strlen(nucs); ++i) nucs[i] = complement(nucs[i]);
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
      if ((word[i] == 'n') || (word[i] == 'N')) return std::numeric_limits<uint64_t>::max();
      h += nuc(word[i]) * std::pow(4, j);
    }
    return h;
  }

  inline uint64_t
  hashwordLong(std::string const& word) {
    boost::hash<std::string> string_hash;
    uint64_t seed = hash_string(word.c_str());
    boost::hash_combine(seed, string_hash(word));
    //std::cerr << seed << '\t' << word << std::endl;
    return seed;
  }

  template<typename TConfig, typename THashMap>
  inline void
  hashShort(TConfig const& c, char* seq, int32_t const len, THashMap& hmap, bool forward) {
    typedef typename THashMap::mapped_type TPosVec;
    uint64_t h = 0;
    bool rewind = true;
    for(int32_t k = 0; k < (int32_t) len - (int32_t) c.matchlen + 1; ++k) {
      if (rewind) {
	h = hashwordShort(std::string(seq + k, seq + k + c.matchlen));
	if (h != std::numeric_limits<uint64_t>::max()) rewind = false;
      } else {
	if ((seq[k+c.matchlen-1] == 'n') || (seq[k+c.matchlen-1] == 'N')) rewind = true;
	else {
	  h -= nuc(seq[k - 1]) * std::pow(4, c.matchlen - 1);
	  h *= 4;
	  h += nuc(seq[k+c.matchlen-1]);
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
      if ((word.find('N') == std::string::npos) && (word.find('n') == std::string::npos)) {
	uint64_t h = hashwordLong(word);
	if (hmap.find(h) == hmap.end()) hmap.insert(std::make_pair(h, TPosVec()));
	if (forward) hmap[h].push_back(k);
	else hmap[h].push_back(len - k - c.matchlen);
      }
    }
  }

  template<typename TConfig, typename THashMap>
  inline void
  wordMatchShort(TConfig const& c, char* seq, int32_t const xlen, int32_t const ylen, THashMap& fwd, THashMap& rev, cv::Mat& img) {
    // Find word matches
    uint64_t h = 0;
    bool rewind = true;
    for(int32_t k = 0; k < (int32_t) ylen - (int32_t) c.matchlen + 1; ++k) {
      if (rewind) {
	h = hashwordShort(std::string(seq + k, seq + k + c.matchlen));
	if (h != std::numeric_limits<uint64_t>::max()) rewind = false;
      } else {
	if ((seq[k+c.matchlen-1] == 'n') || (seq[k+c.matchlen-1] == 'N')) rewind = true;
	else {
	  h -= nuc(seq[k - 1]) * std::pow(4, c.matchlen - 1);
	  h *= 4;
	  h += nuc(seq[k+c.matchlen-1]);
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
	    cv::line(img, cv::Point(px, py), cv::Point(pxend, pyend), cv::Scalar(0, 0, 0), c.lw);
	  }
	}
	// Reverse matches
	if (rev.find(h) != rev.end()) {
	  for(uint32_t idx = 0; idx < rev[h].size(); ++idx) {
	    int32_t px = pixelX(c.usedwidth, xlen, rev[h][idx] + c.matchlen);
	    int32_t pxend = pixelX(c.usedwidth, xlen, rev[h][idx]);
	    int32_t py = pixelX(c.usedheight, ylen, k);
	    int32_t pyend = pixelX(c.usedheight, ylen, k + c.matchlen);
	    cv::line(img, cv::Point(px, py), cv::Point(pxend, pyend), cv::Scalar(0, 0, 255), c.lw);
	  }
	}
      }
    }
  }
  

  template<typename TConfig, typename THashMap>
  inline void
  wordMatchLong(TConfig const& c, char* seq, int32_t const xlen, int32_t const ylen, THashMap& fwd, THashMap& rev, cv::Mat& img) {
    // Find word matches
    for(int32_t k = 0; k < (int32_t) ylen - (int32_t) c.matchlen + 1; ++k) {
      std::string word = std::string(seq + k, seq + k + c.matchlen);
      if ((word.find('N') == std::string::npos) && (word.find('n') == std::string::npos)) {
	uint64_t h = hashwordLong(word);
	// Forward matches
	if (fwd.find(h) != fwd.end()) {
	  for(uint32_t idx = 0; idx < fwd[h].size(); ++idx) {
	    int32_t px = pixelX(c.usedwidth, xlen, fwd[h][idx]);
	    int32_t pxend = pixelX(c.usedwidth, xlen, fwd[h][idx] + c.matchlen);
	    int32_t py = pixelX(c.usedheight, ylen, k);
	    int32_t pyend = pixelX(c.usedheight, ylen, k + c.matchlen);
	    cv::line(img, cv::Point(px, py), cv::Point(pxend, pyend), cv::Scalar(0, 0, 0), c.lw);
	  }
	}
	// Reverse matches
	if (rev.find(h) != rev.end()) {
	  for(uint32_t idx = 0; idx < rev[h].size(); ++idx) {
	    int32_t px = pixelX(c.usedwidth, xlen, rev[h][idx] + c.matchlen);
	    int32_t pxend = pixelX(c.usedwidth, xlen, rev[h][idx]);
	    int32_t py = pixelX(c.usedheight, ylen, k);
	    int32_t pyend = pixelX(c.usedheight, ylen, k + c.matchlen);
	    cv::line(img, cv::Point(px, py), cv::Point(pxend, pyend), cv::Scalar(0, 0, 255), c.lw);
	  }
	}
      }
    }
  }

  template<typename TConfig, typename TReadMapping>
  inline void
  drawXMappings(TConfig const& c, bam_hdr_t* hdr, std::string const& refname, uint32_t const len, std::map<uint32_t, cv::Scalar>& cm, TReadMapping& mp, cv::Mat& img) {
    if (mp.find(refname) != mp.end()) {
      uint32_t runspacer = 4 * c.tlheight;
      for(uint32_t k = 0; k < mp[refname].size(); ++k, runspacer += c.tlheight) {
	int32_t px = pixelX(c.usedwidth, len, mp[refname][k].rstart);
	int32_t pxend = pixelX(c.usedwidth, len, mp[refname][k].rend);
	cv::Rect rect(px, c.usedheight + runspacer, pxend - px, c.tlheight);
	cv::rectangle(img, rect, cm[mp[refname][k].tid], -1);

	std::string text = std::string(hdr->target_name[mp[refname][k].tid]);
	std::string gstart = boost::lexical_cast<std::string>(mp[refname][k].gstart);
	insertComma(gstart);
	std::string gend = boost::lexical_cast<std::string>(mp[refname][k].gend);
	insertComma(gend);
	text += ":" + gstart + "-" + gend;
	if (mp[refname][k].fwd) text += " -->";
	else text += " <--";
	double font_scale = 0.4;
	double font_thickness = 1.5;
	int32_t baseline = 0;
	cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);
	if ((pxend - px) > textSize.width + 10) {
	  cv::putText(img, text, cv::Point((px + pxend) / 2 - textSize.width / 2, c.usedheight + runspacer + textSize.height), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness);
	} else if (px > (int32_t) c.usedwidth / 2) {
	  cv::putText(img, text, cv::Point(px - textSize.width - 5, c.usedheight + runspacer + textSize.height), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness);
	} else {
	  cv::putText(img, text, cv::Point(pxend + 5, c.usedheight + runspacer + textSize.height), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness);
	}
      }
    }
  }

  template<typename TConfig, typename TReadMapping>
  inline void
  drawYMappings(TConfig const& c, bam_hdr_t* hdr, std::string const& refname, uint32_t const len, std::map<uint32_t, cv::Scalar>& cm, TReadMapping& mp, cv::Mat& img) {
    if (mp.find(refname) != mp.end()) {
      uint32_t runspacer = 4 * c.tlheight;
      cv::Mat textimg(img.cols, img.rows, img.type(), cv::Scalar(0, 0, 0, 0)); // rows/cols swapped because of rotation
      for(uint32_t k = 0; k < mp[refname].size(); ++k, runspacer += c.tlheight) {
	int32_t py = pixelX(c.usedheight, len, mp[refname][k].rstart);
	int32_t pyend = pixelX(c.usedheight, len, mp[refname][k].rend);
	cv::Rect rect(c.usedwidth + runspacer, py, c.tlheight, pyend - py);
	cv::rectangle(img, rect, cm[mp[refname][k].tid], -1);

	std::string text = std::string(hdr->target_name[mp[refname][k].tid]);
	std::string gstart = boost::lexical_cast<std::string>(mp[refname][k].gstart);
	insertComma(gstart);
	std::string gend = boost::lexical_cast<std::string>(mp[refname][k].gend);
	insertComma(gend);
	text += ":" + gstart + "-" + gend;
	if (mp[refname][k].fwd) text += " -->";
	else text += " <--";
	double font_scale = 0.4;
	double font_thickness = 1.5;
	int32_t baseline = 0;
	cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);

	if ((pyend - py) > textSize.width + 10) {
	  cv::putText(textimg, text, cv::Point((py + pyend) / 2 - textSize.width / 2, img.cols - (c.usedwidth + runspacer)), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(255, 255, 255), font_thickness);
	} else if (py > (int32_t) c.usedheight / 2) {
	  cv::putText(textimg, text, cv::Point(py - textSize.width - 5, img.cols - (c.usedwidth + runspacer)), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(255, 255, 255), font_thickness);
	} else {
	  cv::putText(textimg, text, cv::Point(pyend + 5, img.cols - (c.usedwidth + runspacer)), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(255, 255, 255), font_thickness);
	}
      }

      // Incorporate text
      cv::rotate(textimg, textimg, cv::ROTATE_90_CLOCKWISE);
      img -= textimg;
    }
  }


  template<typename TConfig>
  inline void
  drawXScaleDotplot(TConfig const& c, std::string const& refname, uint32_t const len, cv::Mat& img) {
    std::string text(boost::lexical_cast<std::string>(len));
    insertComma(text);
    double font_scale = 0.4;
    double font_thickness = 1.5;
    int32_t baseline = 0;
    cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);

    // Find suitable tick size
    uint32_t modval = findTicks(c.pxoffset, textSize.width);

    // Scale line
    cv::line(img, cv::Point(0, c.usedheight), cv::Point(c.usedwidth, c.usedheight), cv::Scalar(255, 0, 0), 2);

    // Ticks
    double px = 0;
    for(uint32_t i = 0; i < len; ++i) {
      if (i % modval == 0) {
	cv::line(img, cv::Point(px - c.pxoffset/2, c.usedheight), cv::Point(px - c.pxoffset/2, c.usedheight + 0.75 * c.tlheight), cv::Scalar(255, 0, 0), 2);
      }
      if (i % modval == 0) {
	// Font
	text = boost::lexical_cast<std::string>(i);
	insertComma(text);
	baseline = 0;
	cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);
	if ((px - c.pxoffset/2 - textSize.width/2 > 0) && (px - c.pxoffset/2 + textSize.width < c.usedwidth)) {
	  cv::putText(img, text, cv::Point(px - c.pxoffset/2 - textSize.width/2, c.usedheight + c.tlheight + textSize.height), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness);
	}
      }
      px += c.pxoffset;
    }

    // Contig name
    if (true) {
      int32_t midpoint = c.usedwidth / 2;
      cv::Size textSize = cv::getTextSize(refname, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);
      cv::putText(img, refname, cv::Point(midpoint - textSize.width/2, c.usedheight + 2 * c.tlheight + textSize.height), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(0, 0, 0), font_thickness);
    }
  }


  template<typename TConfig>
  inline void
  drawYScaleDotplot(TConfig const& c, std::string const& refname, uint32_t const len, cv::Mat& img) {
    std::string text(boost::lexical_cast<std::string>(len));
    insertComma(text);
    double font_scale = 0.4;
    double font_thickness = 1.5;
    int32_t baseline = 0;
    cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);

    // Find suitable tick size
    uint32_t modval = findTicks(c.pyoffset, textSize.width);

    // Scale line
    cv::line(img, cv::Point(c.usedwidth, 0), cv::Point(c.usedwidth, c.usedheight), cv::Scalar(255, 0, 0), 2);

    // Ticks
    cv::Mat textimg(img.cols, img.rows, img.type(), cv::Scalar(0, 0, 0, 0)); // rows/cols swapped because of rotation
    double py = 0;
    for(uint32_t i = 0; i < len; ++i) {
      if (i % modval == 0) {
	cv::line(img, cv::Point(c.usedwidth, py - c.pyoffset/2), cv::Point(c.usedwidth + 0.75 * c.tlheight, py - c.pyoffset/2), cv::Scalar(255, 0, 0), 2);
      }
      if (i % modval == 0) {
	// Font
	text = boost::lexical_cast<std::string>(i);
	insertComma(text);
	baseline = 0;
	cv::Size textSize = cv::getTextSize(text, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);
	if ((py - c.pyoffset/2 - textSize.width/2 > 0) && (py - c.pyoffset/2 + textSize.width < c.usedheight)) {
	  cv::putText(textimg, text, cv::Point(py - c.pyoffset/2 - textSize.width/2, img.cols - (c.usedwidth + c.tlheight + textSize.height)), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(255, 255, 255), font_thickness);
	}
      }
      py += c.pyoffset;
    }

    // Contig name
    if (true) {
      cv::Size textSize = cv::getTextSize(refname, cv::FONT_HERSHEY_SIMPLEX, font_scale, font_thickness, &baseline);
      int32_t midpoint = c.usedheight / 2;
      cv::putText(textimg, refname, cv::Point(midpoint - textSize.width/2, img.cols - (c.usedwidth + 2 * c.tlheight + textSize.height)), cv::FONT_HERSHEY_DUPLEX, font_scale, cv::Scalar(255, 255, 255), font_thickness);
    }

    // Incorporate text
    cv::rotate(textimg, textimg, cv::ROTATE_90_CLOCKWISE);
    img -= textimg;
  }

  
  template<typename TConfigStruct>
  inline int dotplotRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("wally.prof");
#endif
    // Chromosome colors
    cv::Scalar colors[12] = { cv::Scalar(44,160,51), cv::Scalar(180,120,31), cv::Scalar(28,26,227), cv::Scalar(153,255,255), cv::Scalar(227,206,166), cv::Scalar(138,223,178), cv::Scalar(153,154,251), cv::Scalar(111,191,253), cv::Scalar(0,127,255), cv::Scalar(214,178,202), cv::Scalar(154,61,106), cv::Scalar(40,89,177) };

    
    // Read mappings
    typedef std::vector<Mapping> TMappings;
    typedef std::map<std::string, TMappings > TReadMappings;
    TReadMappings mp;

    // Parse BAM
    std::string filename = c.seqfile.string();
    if (c.format == 0) {
      // Make sure the FASTA index is gone
      boost::filesystem::remove(c.seqfile.string());
      boost::filesystem::remove(c.seqfile.string() + ".fai");

      // Parse reads
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Parse reads." << std::endl;
      typedef std::set<std::string> TReadSet;
      TReadSet reads;
      _parseReads(c, reads);
    
      // Get read mappings
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Extract read mappings." << std::endl;
      mappings(c, reads, mp);

      // Sort mappings
      for(TReadMappings::iterator it = mp.begin(); it != mp.end(); ++it) {
	std::sort(it->second.begin(), it->second.end(), SortMappings<Mapping>());
      }
    } else if (c.format == 1) filename = c.file.string();


    // Open file handles
    samFile* samfile = NULL;
    bam_hdr_t* hdr = NULL;
    if (c.format == 0) {
      samfile = sam_open(c.file.string().c_str(), "r");
      hts_set_fai_filename(samfile, c.genome.string().c_str());
      hdr = sam_hdr_read(samfile);
    }
    
    // Load sequences from disk for large contigs
    faidx_t* fai = fai_load(filename.c_str());
    for(int32_t idx1 = 0; idx1 < faidx_nseq(fai) - 1; ++idx1) {
      std::string seqname1(faidx_iseq(fai, idx1));
      int32_t xlen = faidx_seq_len(fai, seqname1.c_str());
      if (xlen < (int32_t) c.seqsize) continue;
      int32_t sl = 0;
      char* seq1 = faidx_fetch_seq(fai, seqname1.c_str(), 0, xlen, &sl);

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
      for(int32_t idx2 = idx1 + 1; idx2 < faidx_nseq(fai); ++idx2) {
	std::string seqname2(faidx_iseq(fai, idx2));
	int32_t ylen = faidx_seq_len(fai, seqname2.c_str());
	if (ylen < (int32_t) c.seqsize) continue;
	sl = 0;
	char* seq2 = faidx_fetch_seq(fai, seqname2.c_str(), 0, ylen, &sl);

	// Next plot
	std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Plot for " << seqname1 << " and " << seqname2 << std::endl;

	// Estimate image size
	if ((c.width == 0) && (c.height == 0)) {
	  int32_t maxlen = std::max(xlen, ylen);
	  c.usedwidth = (int32_t) ((double) xlen / (double) maxlen * 1024.0);
	  c.usedheight = (int32_t) ((double) ylen / (double) maxlen * 1024.0);
	} else if (c.width == 0) {
	  c.usedwidth = (int32_t) ((double) xlen / (double) ylen * c.height);
	} else if (c.height == 0) {
	  c.usedheight = (int32_t) ((double) ylen / (double) xlen * c.width);
	}

	// Create image
	uint32_t footer = 4 * c.tlheight;
	if (c.format == 0) {
	  if (mp.find(seqname1) != mp.end()) footer += (c.tlheight * (mp[seqname1].size() + 1));
	}
	uint32_t margin = 4 * c.tlheight;
	if (c.format == 0) {
	  if (mp.find(seqname2) != mp.end()) margin += (c.tlheight * (mp[seqname2].size() + 1));
	}
	c.pxoffset = (1.0 / (double) xlen) * (double) (c.usedwidth);
	c.pyoffset = (1.0 / (double) ylen) * (double) (c.usedheight);
	cv::Mat img(c.usedheight + footer, c.usedwidth + margin, CV_8UC3, cv::Scalar(255, 255, 255));

	// Compute word matches
	if (c.matchlen < 32) wordMatchShort(c, seq2, xlen, ylen, fwd, rev, img);
	else wordMatchLong(c, seq2, xlen, ylen, fwd, rev, img);

	// Scales
	drawXScaleDotplot(c, seqname1, xlen, img);
	drawYScaleDotplot(c, seqname2, ylen, img);

	// Include mappings
	if (c.format == 0) {
	  // Assign colors to chromosomes
	  typedef std::map<uint32_t, cv::Scalar> TColorMap;
	  TColorMap cm;
	  uint32_t cidx = 0;
	  if (mp.find(seqname1) != mp.end()) {
	    for(uint32_t k = 0; k < mp[seqname1].size(); ++k) {
	      if (cm.find(mp[seqname1][k].tid) == cm.end()) {
		cm[mp[seqname1][k].tid] = colors[cidx];
		cidx = ((cidx + 1) % 12);
	      }
	    }
	  }
	  if (mp.find(seqname2) != mp.end()) {
	    for(uint32_t k = 0; k < mp[seqname2].size(); ++k) {
	      if (cm.find(mp[seqname2][k].tid) == cm.end()) {
		cm[mp[seqname2][k].tid] = colors[cidx];
		cidx = ((cidx + 1) % 12);
	      }
	    }
	  }

	  // Draw mappings
	  drawXMappings(c, hdr, seqname1, xlen, cm, mp, img);
	  drawYMappings(c, hdr, seqname2, ylen, cm, mp, img);
	}

	// Store image (comment this for valgrind, png encoder seems leaky)
	std::string outfile = seqname1;
	outfile += "_";
	outfile += seqname2;
	outfile	+= ".png";
	cv::imwrite(outfile.c_str(), img);
	if (c.showWindow) {
	  cv::imshow(outfile.c_str(), img);
	  cv::waitKey(0);
	}

	// Clean-up
	free(seq2);
      }
      free(seq1);
    }
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
      ("matchlen,m", boost::program_options::value<uint32_t>(&c.matchlen)->default_value(11), "default match length")
      ("size,s", boost::program_options::value<uint32_t>(&c.seqsize)->default_value(0), "min. sequence size to include")
      ("seqfile,q", boost::program_options::value<boost::filesystem::path>(&c.seqfile)->default_value("seq.fa"), "output sequence file in BAM mode")
      ("rfile,R", boost::program_options::value<boost::filesystem::path>(&c.readFile), "file with reads to display")
      ;
    
    boost::program_options::options_description disc("Display options");
    disc.add_options()
      ("linewidth,l", boost::program_options::value<float>(&c.lw)->default_value(1.5), "match line width")
      ("width,x", boost::program_options::value<uint32_t>(&c.width)->default_value(0), "width of the plot [0: best fit]")
      ("height,y", boost::program_options::value<uint32_t>(&c.height)->default_value(0), "height of the plot [0: best fit]")
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
    cmdline_options.add(generic).add(disc).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(disc);
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
	std::cerr << "BAM input requires -g and -R command-line options." << std::endl;
	return 1;
      }

      // Read file check
      if (vm.count("rfile")) {
	if (!(boost::filesystem::exists(c.readFile) && boost::filesystem::is_regular_file(c.readFile) && boost::filesystem::file_size(c.readFile))) {
	  std::cerr << "File with list of reads is missing: " << c.readFile.string() << std::endl;
	  return 1;
	}
      } else {
	std::cerr << "BAM input requires -g and -R command-line options." << std::endl;
	return 1;
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
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "wally ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return dotplotRun(c);
  }

}

#endif
