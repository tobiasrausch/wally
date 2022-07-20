#ifndef DOTPLOT_H
#define DOTPLOT_H


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
  struct ConfigDotplot {
    bool showWindow;
    uint32_t matchlen;
    uint32_t width;
    uint32_t height;
    int32_t format;
    float lw;
    boost::filesystem::path readFile;
    boost::filesystem::path genome;
    boost::filesystem::path file;
  };

  template<typename TConfig>
  inline void
  sequences(TConfig const& c, std::set<std::string> const& reads, std::map<std::string, std::string>& seqmap) {
    
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
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
	std::string qname = bam_get_qname(rec);
	if (reads.find(qname) != reads.end()) {
      
	  // Get read sequence
	  std::string sequence;
	  sequence.resize(rec->core.l_qseq);
	  uint8_t* seqptr = bam_get_seq(rec);
	  for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	  // Store
	  seqmap.insert(std::make_pair(qname, sequence));
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }

  

  template<typename TConfigStruct>
  inline int dotplotRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("wally.prof");
#endif

    typedef std::map<std::string, std::string> TSequences;
    TSequences seqmap;

    if (c.format == 0) {
      // Parse reads
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Parse reads." << std::endl;
      typedef std::set<std::string> TReadSet;
      TReadSet reads;
      _parseReads(c, reads);
    
      // Get read mappings
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Extract reads." << std::endl;
      sequences(c, reads, seqmap);
    } else if (c.format == 1) {
      faidx_t* fai = fai_load(c.file.string().c_str());
      for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
	std::string seqname(faidx_iseq(fai, refIndex));
	uint32_t seqlen = faidx_seq_len(fai, seqname.c_str());
	int32_t sl = 0;
	char* seq = faidx_fetch_seq(fai, seqname.c_str(), 0, seqlen, &sl);
	seqmap.insert(std::make_pair(seqname, std::string(seq)));
	free(seq);
      }
      fai_destroy(fai);
    }
      
    // Find word matches
    for(typename TSequences::const_iterator it = seqmap.begin(); it != seqmap.end(); ++it) {
      int32_t xlen = it->second.size();
      typename TSequences::const_iterator itNext = it;
      ++itNext;
      for(;itNext != seqmap.end(); ++itNext) {
	int32_t ylen = itNext->second.size();
	std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Plot for " << it->first << " and " << itNext->first << std::endl;

	// Generate image, y-axis: itNext, x-axis: it
	int32_t usedwidth = c.width;
	int32_t usedheight = c.height;
	if ((c.width == 0) && (c.height == 0)) {
	  int32_t maxlen = std::max(xlen, ylen);
	  usedwidth = (int32_t) ((double) xlen / (double) maxlen * 1024.0);
	  usedheight = (int32_t) ((double) ylen / (double) maxlen * 1024.0);
	} else if (c.width == 0) {
	  usedwidth = (int32_t) ((double) xlen / (double) ylen * c.height);
	} else if (c.height == 0) {
	  usedheight = (int32_t) ((double) ylen / (double) xlen * c.width);
	}
	cv::Mat img(usedheight, usedwidth, CV_8UC3, cv::Scalar(255, 255, 255));

	// Forward words
	typedef std::vector<uint32_t> TPosVec;
	typedef std::map<std::string, TPosVec> TMatchMap;
	TMatchMap fwd;
	for(int32_t k = 0; k < (int32_t) ylen - (int32_t) c.matchlen; ++k) {
	  std::string word = itNext->second.substr(k, c.matchlen);
	  if (fwd.find(word) == fwd.end()) fwd.insert(std::make_pair(word, TPosVec()));
	  fwd[word].push_back(k);
	}
	// Reverse words
	std::string rc = itNext->second;
	reverseComplement(rc);
	TMatchMap rev;
	for(int32_t k = 0; k < (int32_t) ylen - (int32_t) c.matchlen; ++k) {
	  std::string word = rc.substr(k, c.matchlen);
	  if (rev.find(word) == rev.end()) rev.insert(std::make_pair(word, TPosVec()));
	  rev[word].push_back(ylen - k - c.matchlen);
	}

	// Find word matches
	for(int32_t k = 0; k < (int32_t) xlen - (int32_t) c.matchlen; ++k) {
	  std::string word = it->second.substr(k, c.matchlen);
	  // Forward matches
	  if (fwd.find(word) != fwd.end()) {
	    for(uint32_t idx = 0; idx < fwd[word].size(); ++idx) {
	      int32_t px = pixelX(usedwidth, xlen, k);
	      int32_t pxend = pixelX(usedwidth, xlen, k + c.matchlen);
	      int32_t py = pixelX(usedheight, ylen, fwd[word][idx]);
	      int32_t pyend = pixelX(usedheight, ylen, fwd[word][idx] + c.matchlen);
	      cv::line(img, cv::Point(px, py), cv::Point(pxend, pyend), cv::Scalar(195, 142, 153), c.lw);
	    }
	  }
	  // Reverse matches
	  if (rev.find(word) != rev.end()) {
	    for(uint32_t idx = 0; idx < rev[word].size(); ++idx) {
	      int32_t px = pixelX(usedwidth, xlen, k);
	      int32_t pxend = pixelX(usedwidth, xlen, k + c.matchlen);
	      int32_t py = pixelX(usedheight, ylen, rev[word][idx] + c.matchlen);
	      int32_t pyend = pixelX(usedheight, ylen, rev[word][idx]);
	      cv::line(img, cv::Point(px, py), cv::Point(pxend, pyend), cv::Scalar(64, 163, 241), c.lw);
	    }
	  }
	}

	
	// Store image (comment this for valgrind, png encoder seems leaky)
	std::string outfile = it->first;
	outfile += "_";
	outfile += itNext->first;
	outfile	+= ".png";
	cv::imwrite(outfile.c_str(), img);
	if (c.showWindow) {
	  cv::imshow(outfile.c_str(), img);
	  cv::waitKey(0);
	}
      }
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

    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("rfile,R", boost::program_options::value<boost::filesystem::path>(&c.readFile), "file with reads to display")
      ;
    
    boost::program_options::options_description disc("Display options");
    disc.add_options()
      ("matchlen,m", boost::program_options::value<uint32_t>(&c.matchlen)->default_value(10), "default match length")
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
      if (vm.count("rfile")) {
	if (!(boost::filesystem::exists(c.readFile) && boost::filesystem::is_regular_file(c.readFile) && boost::filesystem::file_size(c.readFile))) {
	  std::cerr << "File with list of reads is missing: " << c.readFile.string() << std::endl;
	  return 1;
	}
      } else {
	std::cerr << "BAM input requires -g and -R command-line options." << std::endl;
	return 1;
      }
    }

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
