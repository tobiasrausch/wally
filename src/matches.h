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
    bool hasAnnotationFile;
    bool showSoftClip;
    uint16_t splits;
    int32_t winsize;
    uint32_t width;
    uint32_t height;
    uint32_t tlheight;  // pixel height of a track line
    uint32_t rdheight;  // pixel height of a single read
    double pxoffset; // 1bp in pixel
    double bpoffset; // 1pixel in bp
    boost::filesystem::path bedFile;
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
  clusterRegions(TConfig const& c, int32_t const nchr, TReadMappings const& mp, std::vector<Region>& rg) {
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
	rg.push_back(Region(refIndex, rgStart, rgEnd));
      }
    }
  }
  

  template<typename TConfig>
  inline int32_t
  mappings(TConfig const& c, std::set<std::string> const& reads, std::map<std::string, std::vector<Mapping> >& mp) {
    uint32_t matchCount = 0;
    
    // Open file handles
    samFile* samfile = sam_open(c.file.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.file.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Parse BAM
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing... " << hdr->target_name[refIndex] << std::endl;

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
	      gp += bam_cigar_oplen(cigar[i]);
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
	    ++matchCount;
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

    return matchCount;
  }

  

  template<typename TConfigStruct>
  inline int matchRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("wally.prof");
#endif

    // Parse reads
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Parse reads." << std::endl;
    std::set<std::string> reads;
    _parseReads(c, reads);
    
    // Get read mappings
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Extract read mappings." << std::endl;
    typedef std::vector<Mapping> TMappings;
    typedef std::map<std::string, TMappings > TReadMappings;
    TReadMappings mp;
    int32_t matchCount = mappings(c, reads, mp);
    if (!matchCount) {
      std::cerr << "Error: No mappings found! Are the read names correct?" << std::endl;
      return 1;
    }
    
    // Check image height
    int32_t headerTracks = 4;
    if ((matchCount + headerTracks) * c.tlheight > c.height) {
      std::cerr << "Warning: Image height is too small to display all matches!" << std::endl;
      c.height = (matchCount + headerTracks) * c.tlheight;
      std::cerr << "Warning: Adjusting image height to " << c.height << std::endl;
    }

    // Open file handles
    samFile* samfile = sam_open(c.file.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.file.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Sort mappings
    for(TReadMappings::iterator it = mp.begin(); it != mp.end(); ++it) {
      std::sort(it->second.begin(), it->second.end(), SortMappings<Mapping>());

      // Debug
      //for(uint32_t i = 0; i < it->second.size(); ++i) std::cerr << hdr->target_name[it->second[i].tid] << ':' << it->second[i].gstart << '-' << it->second[i].gend << '\t' << it->second[i].rstart << '-' << it->second[i].rend << '(' << (int) it->second[i].fwd << ')' << '\t' << it->first << std::endl;
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

    // Load genome
    faidx_t* fai = fai_load(c.genome.string().c_str());
    int32_t seqlen;
    std::string oldchr("None");
    char* seq = NULL;

    // Check image width
    uint32_t estwidth = minPixelWidth(c, rg);
    if (c.width < estwidth) {
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
      
      // Lazy loading of genome
      std::string chrName = hdr->target_name[rg[rgIdx].tid];
      if (chrName != oldchr) {
	if (seq != NULL) free(seq);
	seq = faidx_fetch_seq(fai, hdr->target_name[rg[rgIdx].tid], 0, hdr->target_len[rg[rgIdx].tid], &seqlen);
	oldchr = chrName;
      }

      // Header
      drawCoordinates(c, rg[rgIdx], hdr->target_name[rg[rgIdx].tid], bg, 1);

      // Reference
      drawReference(c, rg[rgIdx], bg, boost::to_upper_copy(std::string(seq + rg[rgIdx].beg, seq + rg[rgIdx].end)), 2);

      // Draw split border
      drawSplitBorder(c, bg);
      
      // Parse mappings
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Region " << hdr->target_name[rg[rgIdx].tid] << ':' << rg[rgIdx].beg << '-' << rg[rgIdx].end << std::endl;

      // Iterate all reads
      for(TReadMappings::iterator it = mp.begin(); it != mp.end(); ++it) {
	int32_t trackIdx = headerTracks;
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
	if (rgIdx == 0) drawSampleLabel(c, headerTracks, it->first, bg);
      }
      
      // Store image (comment this for valgrind, png encoder seems leaky)
      if (c.splits == 1) {
	std::string outfile = "test.png";
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
	  std::string outfile = "test.png";
	  cv::imwrite(outfile.c_str(), dst);
	  if (c.showWindow) {
	    cv::imshow(convertToStr(hdr, rg[rgIdx]).c_str(), dst);
	    cv::waitKey(0);
	  }
	  imageStore.clear();
	}
      }
    }

    // Clean-up
    if (seq != NULL) free(seq);
    fai_destroy(fai);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
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
    c.showSoftClip = true;
    
    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("bed,b", boost::program_options::value<boost::filesystem::path>(&c.bedFile), "BED annotation file (optional)")
      ;
    
    boost::program_options::options_description disc("Graphics options");
    disc.add_options()
      ("reads,r", boost::program_options::value<boost::filesystem::path>(&c.readFile), "file with reads to display")
      ("winsize,w", boost::program_options::value<int32_t>(&c.winsize)->default_value(10000), "window size to cluster nearby matches")
      ;
    
    boost::program_options::options_description geno("Display options");
    geno.add_options()
      ("width,x", boost::program_options::value<uint32_t>(&c.width)->default_value(1024), "width of the plot")
      ("height,y", boost::program_options::value<uint32_t>(&c.height)->default_value(1024), "height of the plot")
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
    cmdline_options.add(generic).add(disc).add(geno).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(disc).add(geno);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) { 
      std::cout << std::endl;
      std::cout << "Usage: wally " << argv[0] << " [OPTIONS] -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> ..." << std::endl;
      std::cout << visible_options << "\n";
      return 0;
    }

    // Show window?
    if (vm.count("window")) c.showWindow = true;
    else c.showWindow = false;

    // BED file
    if (vm.count("bed")) {
      if (!(boost::filesystem::exists(c.bedFile) && boost::filesystem::is_regular_file(c.bedFile) && boost::filesystem::file_size(c.bedFile))) {
	std::cerr << "Genome annotation BED file is missing: " << c.bedFile.string() << std::endl;
	return 1;
      }
      c.hasAnnotationFile = true;
    } else c.hasAnnotationFile = false;
    
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
