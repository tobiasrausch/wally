#ifndef MATCHES_H
#define MATCHES_H


#include <iostream>
#include <fstream>

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
#include "matchdraw.h"
#include "bed.h"

namespace wallysworld
{

  // Config arguments
  struct ConfigMatches {
    bool showWindow;
    bool hasReadFile;
    bool separatePlots;
    bool storeSequences;
    uint16_t splits;
    int32_t winsize;
    uint32_t minMatches;
    uint32_t width;
    uint32_t height;
    uint32_t tlheight;  // pixel height of a track line
    uint32_t rdheight;  // pixel height of a single read
    double pxoffset; // 1bp in pixel
    double bpoffset; // 1pixel in bp
    std::string readStr;
    boost::filesystem::path outfile;
    boost::filesystem::path seqfile;
    boost::filesystem::path readFile;
    boost::filesystem::path genome;
    boost::filesystem::path file;
  };

  template<typename TConfig, typename TReads>
  inline void
  _parseReads(TConfig const& c, TReads& reads) {
    typedef typename TReads::value_type TValue;
    std::ifstream readFile(c.readFile.string().c_str(), std::ifstream::in);
    if (readFile.is_open()) {
      while (readFile.good()) {
	std::string line;
	getline(readFile, line);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t");
	Tokenizer tokens(line, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  reads.insert(boost::lexical_cast<TValue>(*tokIter));
	}
      }
      readFile.close();
    }
  }

  template<typename TConfig, typename TReadMappings>
  inline void
  clusterRegions(TConfig const& c, int32_t const nchr, TReadMappings& mp, std::vector<Region>& rg) {
    for(int32_t refIndex = 0; refIndex < nchr; ++refIndex) {

      // Collect intervals for this chromosomes
      typedef boost::icl::interval_set<int32_t> TChrIntervals;
      typedef typename TChrIntervals::interval_type TIVal;
      TChrIntervals cint;
      for(typename TReadMappings::const_iterator it = mp.begin(); it != mp.end(); ++it) {
	for(uint32_t i = 0; i < it->second.size(); ++i) {
	  if (it->second[i].tid == refIndex) {
	    int32_t gstart = it->second[i].gstart - c.winsize;
	    int32_t gend = it->second[i].gend + c.winsize;
	    cint.insert(TIVal::right_open(gstart, gend));
	  }
	}
      }

      // Keep clustered interval with 10% buffer space
      for(typename TChrIntervals::iterator it = cint.begin(); it != cint.end(); ++it) {
	int32_t realstart = it->lower() + c.winsize;
	int32_t realend = it->upper() - c.winsize;
	int32_t offset = (int32_t) (0.1 * (double) (realend - realstart));
	int32_t rgStart = 0;
	if (offset < realstart) rgStart = realstart - offset;
	int32_t rgEnd = realend + offset;
	if (c.minMatches > 1) {
	  uint32_t matchCount = 0;
	  for(typename TReadMappings::const_iterator it = mp.begin(); it != mp.end(); ++it) {
	    for(uint32_t i = 0; i < it->second.size(); ++i) {
	      if (it->second[i].tid == refIndex) {
		if ((it->second[i].gstart >= rgStart) && (it->second[i].gend <= rgEnd)) ++matchCount;
	      }
	    }
	  }
	  if (matchCount >= c.minMatches) rg.push_back(Region(refIndex, rgStart, rgEnd));
	} else {
	  rg.push_back(Region(refIndex, rgStart, rgEnd));
	}
      }
    }

    // Remove mappings if min. matches is used
    if (c.minMatches > 1) {
      for(typename TReadMappings::iterator it = mp.begin(); it != mp.end(); ++it) {
	bool rerun = true;
	while (rerun) {
	  rerun = false;
	  for(uint32_t i = 0; i < it->second.size(); ++i) {
	    bool found = false;
	    for(uint32_t k = 0; k < rg.size(); ++k) {
	      if ((it->second[i].tid == rg[k].tid) && (it->second[i].gstart >= rg[k].beg) && (it->second[i].gend <= rg[k].end)) {
		found = true;
		break;
	      }
	    }
	    if (!found) {
	      it->second.erase(it->second.begin() + i);
	      rerun = true;
	      break;
	    }
	  }
	}
      }
    }
  }
  

  template<typename TConfig>
  inline void
  mappings(TConfig const& c, std::set<std::string> const& reads, std::map<std::string, std::vector<Mapping> >& mp) {
    // Sequence file
    std::ofstream sfile;
    if (c.storeSequences) sfile.open(c.seqfile.string().c_str());
    
    // Open file handles
    samFile* samfile = sam_open(c.file.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.file.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Parse BAM
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Processing... " << hdr->target_name[refIndex] << std::endl;

      // Iterate alignments 
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY)) continue;
	std::string qname = bam_get_qname(rec);
	if (reads.find(qname) != reads.end()) {
      
	  // Get read sequence
	  std::string sequence;
	  sequence.resize(rec->core.l_qseq);
	  uint8_t* seqptr = bam_get_seq(rec);
	  for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	  // Store primary sequence?
	  if (c.storeSequences) {
	    if (!(rec->core.flag & (BAM_FSUPPLEMENTARY))) {
	      sfile << ">" << qname << std::endl;
	      sfile << sequence << std::endl;
	    }
	  }
	  
	  // Parse CIGAR
	  uint32_t* cigar = bam_get_cigar(rec);
	  int32_t gp = rec->core.pos; // Genomic position
	  int32_t gpStart = -1; //Match start
	  int32_t gpEnd = -1; //Match end
	  int32_t sp = 0; // Sequence position
	  int32_t seqStart = -1;  // Match start
	  int32_t seqEnd = -1; // Match end
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	      if (seqStart == -1) {
		seqStart = sp;
		gpStart = gp;
	      }
	      gp += bam_cigar_oplen(cigar[i]);
	      sp += bam_cigar_oplen(cigar[i]);
	      seqEnd = sp;
	      gpEnd = gp;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	      if (seqStart == -1) {
		seqStart = sp;
		gpStart = gp;
	      }
	      sp += bam_cigar_oplen(cigar[i]);
	      seqEnd = sp;
	      gpEnd = gp;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	      if (seqStart == -1) {
		seqStart = sp;
		gpStart = gp;
	      }
	      gp += bam_cigar_oplen(cigar[i]);
	      seqEnd = sp;
	      gpEnd = gp;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	      bool dir = true;
	      if (rec->core.flag & BAM_FREVERSE) {
		dir = false;
		int32_t seqTmp = seqStart;
		seqStart = sp - seqEnd;
		seqEnd = sp - seqTmp;
	      }
	      if (gpStart < gpEnd) {
		if (mp.find(qname) == mp.end()) mp[qname] = std::vector<Mapping>();
		mp[qname].push_back(Mapping(rec->core.tid, gpStart, gpEnd, seqStart, seqEnd, dir, rec->core.qual));
	      }
	      gp += bam_cigar_oplen(cigar[i]);
	      // Reset match
	      gpStart = -1;
	      gpEnd = -1;
	      seqStart = -1;
	      seqEnd = -1;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else {
	      std::cerr << "Warning: Unknown Cigar options!" << std::endl;
	    }
	  }	  
	  bool dir = true;
	  if (rec->core.flag & BAM_FREVERSE) {
	    dir = false;
	    int32_t seqTmp = seqStart;
	    seqStart = sp - seqEnd;
	    seqEnd = sp - seqTmp;
	  }
	  if (gpStart < gpEnd) {
	    if (mp.find(qname) == mp.end()) mp[qname] = std::vector<Mapping>();
	    mp[qname].push_back(Mapping(rec->core.tid, gpStart, gpEnd, seqStart, seqEnd, dir, rec->core.qual));
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    // Close sequence file
    if (c.storeSequences) sfile.close();
  }
  

  template<typename TConfigStruct>
  inline int matchRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("wally.prof");
#endif

    // Parse reads
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Parse reads." << std::endl;
    typedef std::set<std::string> TReadSet;
    TReadSet reads;
    if (c.hasReadFile) _parseReads(c, reads);
    else reads.insert(c.readStr);
    
    // Get read mappings
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Extract read mappings." << std::endl;
    typedef std::vector<Mapping> TMappings;
    typedef std::map<std::string, TMappings > TReadMappings;
    TReadMappings mp;
    mappings(c, reads, mp);

    // Sort mappings
    for(TReadMappings::iterator it = mp.begin(); it != mp.end(); ++it) {
      std::sort(it->second.begin(), it->second.end(), SortMappings<Mapping>());

      // Debug
      //for(uint32_t i = 0; i < it->second.size(); ++i) std::cerr << hdr->target_name[it->second[i].tid] << ':' << it->second[i].gstart << '-' << it->second[i].gend << '\t' << it->second[i].rstart << '-' << it->second[i].rend << '(' << (int) it->second[i].fwd << ')' << '\t' << it->first << std::endl;
    }

    // Number of plots
    uint32_t numPlots = 1;
    uint32_t source_c_width = c.width;
    uint32_t source_c_height = c.height;
    TReadMappings source_mp;
    TReadSet source_reads;
    if (c.separatePlots) {
      numPlots = reads.size();
      source_mp = mp;
      source_reads = reads;
    }

    // Open file handles
    samFile* samfile = sam_open(c.file.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Iterate reads
    typename TReadSet::iterator itread = source_reads.begin();
    for(uint32_t plotk = 0; plotk < numPlots; ++plotk) {
      if (c.separatePlots) {
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Plot for " << *itread << std::endl;
	if (source_mp.find(*itread) != source_mp.end()) {
	  mp.clear();
	  reads.clear();
	  mp.insert(std::make_pair(*itread, source_mp[*itread]));
	  reads.insert(*itread);
	  c.width = source_c_width;
	  c.height = source_c_height;
	}
      }

      // Determine regions
      std::vector<Region> rg;
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Cluster nearby matches." << std::endl;
      clusterRegions(c, hdr->n_targets, mp, rg);
      if (rg.empty()) {
	std::cerr << "Error: No regions to display!" << std::endl;
	return 1;
      }
      c.splits = rg.size();
      // Debug
      //for(uint32_t i = 0; i < rg.size(); ++i) std::cerr << hdr->target_name[rg[i].tid] << ':' << rg[i].beg << '-' << rg[i].end << std::endl;

      // Count number of matches
      int32_t matchCount = 0;
      if (c.separatePlots) matchCount = mp[*itread].size();
      else {
	for(typename TReadMappings::const_iterator it = mp.begin(); it != mp.end(); ++it) matchCount += it->second.size();
      }	
      if (!matchCount) {
	std::cerr << "Error: No mappings found! Are the read names correct?" << std::endl;
	return 1;
      }      

      // Check image height
      int32_t headerTracks = 3;
      if (c.height == 0) c.height = (matchCount + headerTracks + reads.size()) * c.tlheight;
      else if ((matchCount + headerTracks + reads.size()) * c.tlheight > c.height) {
	std::cerr << "Warning: Image height is too small to display all matches!" << std::endl;
	c.height = (matchCount + headerTracks + reads.size()) * c.tlheight;
	std::cerr << "Warning: Adjusting image height to " << c.height << std::endl;
      }


      
      // Check image width
      uint32_t estwidth = minPixelWidth(c, rg);
      if (c.width == 0) c.width = estwidth;
      else if (c.width < estwidth) {
	std::cerr << "Warning: Image width is too small to display all matches!" << std::endl;
	c.width = estwidth;
	std::cerr << "Warning: Adjusting image width to " << c.width << std::endl;
      }
      c.width /= c.splits;

      // Store images
      std::vector<cv::Mat> imageStore;

      // Split region
      for(uint32_t rgIdx = 0; rgIdx < rg.size(); ++rgIdx) {
	// Get pixel width of 1bp
	c.pxoffset = (1.0 / (double) rg[rgIdx].size) * (double) c.width;
	// Get bp of 1 pixel
	c.bpoffset = (1.0 / (double) c.width) * (double) rg[rgIdx].size;
      
	// Generate image
	cv::Mat bg( c.height, c.width, CV_8UC3, cv::Scalar(255, 255, 255));
      
	// Block header tracks
	int32_t maxTracks = c.height / c.tlheight;
	std::vector<int32_t> taken(maxTracks, WALLY_UNBLOCK);
	for(int32_t i = 0; i < headerTracks; ++i) taken[i] = WALLY_BLOCKED;
	
	// Header
	drawCoordinates(c, rg[rgIdx], hdr->target_name[rg[rgIdx].tid], bg);
	
	// Split line
	drawSplitLine(c, bg, 3);
	
	// Draw split border
	drawSplitBorder(c, bg);
      
	// Parse mappings
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Region " << hdr->target_name[rg[rgIdx].tid] << ':' << rg[rgIdx].beg << '-' << rg[rgIdx].end << std::endl;

	// Iterate all reads
	int32_t prevReadOffset = headerTracks;
	for(TReadMappings::iterator it = mp.begin(); it != mp.end(); ++it) {
	  int32_t trackIdx = prevReadOffset;
	  ++trackIdx; // Space for read name
	  bool prevDraw = false;
	  for(uint32_t i = 0; i < it->second.size(); ++i) {
	    if ((it->second[i].tid == rg[rgIdx].tid) && (it->second[i].gstart >= rg[rgIdx].beg) && (it->second[i].gend <= rg[rgIdx].end)) {
	      // Draw match in the current region
	      if ((i > 0) && (i + 1 < it->second.size())) {
		drawBlock(c, rg[rgIdx], bg, trackIdx, it->second[i-1], it->second[i], it->second[i+1]);
	      } else {
		if ((i == 0) && (it->second.size() == 1)) {
		  drawBlock(c, rg[rgIdx], bg, trackIdx, Mapping(), it->second[i], Mapping());
		} else {
		  if (i == 0) drawBlock(c, rg[rgIdx], bg, trackIdx, Mapping(), it->second[i], it->second[i+1]);
		  else drawBlock(c, rg[rgIdx], bg, trackIdx, it->second[i-1], it->second[i], Mapping());
		}
	      }
	      prevDraw = true;
	    } else {
	      // This match is outside but the connection might go through
	      if ((i > 0) && (!prevDraw)) drawCrossConnect(c, rg[rgIdx], bg, trackIdx, it->second[i-1], it->second[i]);
	      prevDraw = false;
	    }	  
	    ++trackIdx;
	  }	
	  // Draw read name
	  if (rgIdx == 0) drawSampleLabel(c, prevReadOffset + 1, it->first, bg);
	  // Draw split line
	  drawSplitLine(c, bg, trackIdx);
	  prevReadOffset = trackIdx;
	}
          
	// Store image (comment this for valgrind, png encoder seems leaky)
	if (c.splits == 1) {
	  std::string outfile = c.outfile.string();
	  if (c.separatePlots) outfile = (*itread) + ".png";
	  cv::imwrite(outfile.c_str(), bg);
	  if (c.showWindow) {
	    cv::imshow(convertToStr(hdr, rg[rgIdx]).c_str(), bg);
	    cv::waitKey(0);
	  }
	} else {
	  imageStore.push_back(bg);
	  // Concatenate
	  if (imageStore.size() == c.splits) {
	    cv::Mat dst;
	    cv::hconcat(imageStore[0], imageStore[1], dst);
	    for(uint32_t i = 2; i < imageStore.size(); ++i) {
	      cv::Mat tdst;
	      cv::hconcat(dst, imageStore[i], tdst);
	      dst = tdst;
	    }
	    std::string outfile = c.outfile.string();
	    if (c.separatePlots) outfile = (*itread) + ".png";
	    cv::imwrite(outfile.c_str(), dst);
	    if (c.showWindow) {
	      cv::imshow(convertToStr(hdr, rg[rgIdx]).c_str(), dst);
	      cv::waitKey(0);
	    }
	    imageStore.clear();
	  }
	}
      }

      // Next read
      if (c.separatePlots) ++itread;
    }

    // Clean-up
    bam_hdr_destroy(hdr);
    sam_close(samfile);

#ifdef PROFILE
    ProfilerStop();
#endif
  
    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    return 0;
  }


  int matches(int argc, char **argv) {
    ConfigMatches c;
    c.tlheight = 15;
    c.rdheight = 12;
    
    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("read,r", boost::program_options::value<std::string>(&c.readStr)->default_value("read_name"), "read to display")
      ("rfile,R", boost::program_options::value<boost::filesystem::path>(&c.readFile), "file with reads to display")
      ("seqfile,q", boost::program_options::value<boost::filesystem::path>(&c.seqfile), "sequence output file [optional]")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.png"), "plot output file")
      ("separate,s", "create one plot for each input read")
      ;
    
    boost::program_options::options_description disc("Display options");
    disc.add_options()
      ("winsize,n", boost::program_options::value<int32_t>(&c.winsize)->default_value(10000), "window size to cluster nearby matches")
      ("matches,m", boost::program_options::value<uint32_t>(&c.minMatches)->default_value(1), "min. number of matches per region")
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
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) { 
      std::cout << std::endl;
      std::cout << "Usage: wally " << argv[0] << " [OPTIONS] -g <ref.fa> <input.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 0;
    }

    // Show window?
    if (vm.count("window")) c.showWindow = true;
    else c.showWindow = false;

    // Separate plots?
    if (vm.count("separate")) c.separatePlots = true;
    else c.separatePlots = false;

    // Store sequences
    if (vm.count("seqfile")) c.storeSequences = true;
    else c.storeSequences = false;

    // Read file
    if (vm.count("rfile")) {
      if (!(boost::filesystem::exists(c.readFile) && boost::filesystem::is_regular_file(c.readFile) && boost::filesystem::file_size(c.readFile))) {
	std::cerr << "File with list of reads is missing: " << c.readFile.string() << std::endl;
	return 1;
      }
      c.hasReadFile = true;
    } else c.hasReadFile = false;

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "wally ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return matchRun(c);
  }

}

#endif
