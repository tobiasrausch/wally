#ifndef DOTPLOT_H
#define DOTPLOT_H


#include <iostream>
#include <fstream>

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
    uint32_t usedwidth;
    uint32_t usedheight;
    int32_t format;
    float lw;
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
    for(int32_t i = word.size() - 1; i>=0; --i, ++j) h += nuc(word[i]) * std::pow(4, j);
    return h;
  }

  inline uint64_t
  hashwordLong(std::string const& word) {
    boost::hash<std::string> string_hash;
    uint64_t seed = hash_string(word.c_str());
    boost::hash_combine(seed, string_hash(word));
    return seed;
  }

  template<typename TConfig>
  inline void
  sequences(TConfig const& c, std::string const& filename, std::set<std::string> const& reads) {

    std::ofstream sfile(filename.c_str());
    
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

	  // Write to tmp file
	  sfile << ">" << qname << std::endl;
	  sfile << sequence << std::endl;
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    // Close file
    sfile.close();
  }


  template<typename TConfig>
  inline void
  wordMatchShort(TConfig const& c, char* seq1, char* seq2, int32_t const xlen, int32_t const ylen, cv::Mat& img) {
    // Forward words
    typedef std::vector<uint32_t> TPosVec;
    typedef std::map<uint64_t, TPosVec> TMatchMap;
    TMatchMap fwd;	
    uint64_t h = 0;
    for(int32_t k = 0; k < (int32_t) ylen - (int32_t) c.matchlen; ++k) {
      if (k) {
	h -= nuc(seq2[k - 1]) * std::pow(4, c.matchlen - 1);
	h *= 4;
	h += nuc(seq2[k+c.matchlen-1]);
      } else h = hashwordShort(std::string(seq2 + k, seq2 + k + c.matchlen));
      if (fwd.find(h) == fwd.end()) fwd.insert(std::make_pair(h, TPosVec()));
      fwd[h].push_back(k);
    }
    // Reverse words
    revcomplement(seq2);
    TMatchMap rev;
    h = 0;
    for(int32_t k = 0; k < (int32_t) ylen - (int32_t) c.matchlen; ++k) {
      if (k) {
	h -= nuc(seq2[k - 1]) * std::pow(4, c.matchlen - 1);
	h *= 4;
	h += nuc(seq2[k+c.matchlen-1]);
      } else h = hashwordShort(std::string(seq2 + k, seq2 + k + c.matchlen));
      if (rev.find(h) == rev.end()) rev.insert(std::make_pair(h, TPosVec()));
      rev[h].push_back(ylen - k - c.matchlen);
    }
    
    // Find word matches
    h = 0;
    for(int32_t k = 0; k < (int32_t) xlen - (int32_t) c.matchlen; ++k) {
      if (k) {
	h -= nuc(seq1[k - 1]) * std::pow(4, c.matchlen - 1);
	h *= 4;
	h += nuc(seq1[k+c.matchlen-1]);
      } else h = hashwordShort(std::string(seq1 + k, seq1 + k + c.matchlen));
      // Forward matches
      if (fwd.find(h) != fwd.end()) {
	for(uint32_t idx = 0; idx < fwd[h].size(); ++idx) {
	  int32_t px = pixelX(c.usedwidth, xlen, k);
	  int32_t pxend = pixelX(c.usedwidth, xlen, k + c.matchlen);
	  int32_t py = pixelX(c.usedheight, ylen, fwd[h][idx]);
	  int32_t pyend = pixelX(c.usedheight, ylen, fwd[h][idx] + c.matchlen);
	  cv::line(img, cv::Point(px, py), cv::Point(pxend, pyend), cv::Scalar(0, 0, 0), c.lw);
	}
      }
      // Reverse matches
      if (rev.find(h) != rev.end()) {
	for(uint32_t idx = 0; idx < rev[h].size(); ++idx) {
	  int32_t px = pixelX(c.usedwidth, xlen, k);
	  int32_t pxend = pixelX(c.usedwidth, xlen, k + c.matchlen);
	  int32_t py = pixelX(c.usedheight, ylen, rev[h][idx] + c.matchlen);
	  int32_t pyend = pixelX(c.usedheight, ylen, rev[h][idx]);
	  cv::line(img, cv::Point(px, py), cv::Point(pxend, pyend), cv::Scalar(0, 0, 255), c.lw);
	}
      }
    }
  }
  

  template<typename TConfig>
  inline void
  wordMatchLong(TConfig const& c, char* seq1, char* seq2, int32_t const xlen, int32_t const ylen, cv::Mat& img) {
    // Forward words
    typedef std::vector<uint32_t> TPosVec;
    typedef std::map<uint64_t, TPosVec> TMatchMap;
    TMatchMap fwd;
    for(int32_t k = 0; k < (int32_t) ylen - (int32_t) c.matchlen; ++k) {
      uint64_t h = hashwordLong(std::string(seq2 + k, seq2 + k + c.matchlen));
      if (fwd.find(h) == fwd.end()) fwd.insert(std::make_pair(h, TPosVec()));
      fwd[h].push_back(k);
    }
    // Reverse words
    revcomplement(seq2);
    TMatchMap rev;
    for(int32_t k = 0; k < (int32_t) ylen - (int32_t) c.matchlen; ++k) {
      uint64_t h = hashwordLong(std::string(seq2 + k, seq2 + k + c.matchlen));
      if (rev.find(h) == rev.end()) rev.insert(std::make_pair(h, TPosVec()));
      rev[h].push_back(ylen - k - c.matchlen);
    }
    
    // Find word matches
    for(int32_t k = 0; k < (int32_t) xlen - (int32_t) c.matchlen; ++k) {
      uint64_t h = hashwordLong(std::string(seq1 + k, seq1 + k + c.matchlen));
      // Forward matches
      if (fwd.find(h) != fwd.end()) {
	for(uint32_t idx = 0; idx < fwd[h].size(); ++idx) {
	  int32_t px = pixelX(c.usedwidth, xlen, k);
	  int32_t pxend = pixelX(c.usedwidth, xlen, k + c.matchlen);
	  int32_t py = pixelX(c.usedheight, ylen, fwd[h][idx]);
	  int32_t pyend = pixelX(c.usedheight, ylen, fwd[h][idx] + c.matchlen);
	  cv::line(img, cv::Point(px, py), cv::Point(pxend, pyend), cv::Scalar(0, 0, 0), c.lw);
	}
      }
      // Reverse matches
      if (rev.find(h) != rev.end()) {
	for(uint32_t idx = 0; idx < rev[h].size(); ++idx) {
	  int32_t px = pixelX(c.usedwidth, xlen, k);
	  int32_t pxend = pixelX(c.usedwidth, xlen, k + c.matchlen);
	  int32_t py = pixelX(c.usedheight, ylen, rev[h][idx] + c.matchlen);
	  int32_t pyend = pixelX(c.usedheight, ylen, rev[h][idx]);
	  cv::line(img, cv::Point(px, py), cv::Point(pxend, pyend), cv::Scalar(0, 0, 255), c.lw);
	}
      }
    }
  }
  
  template<typename TConfigStruct>
  inline int dotplotRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("wally.prof");
#endif

    typedef std::map<std::string, std::string> TSequences;
    TSequences seqmap;

    // Sequences
    boost::uuids::uuid uuid = boost::uuids::random_generator()();
    std::string filename = "sequences." + boost::lexical_cast<std::string>(uuid) + ".fa";
    if (c.format == 0) {
      // Parse reads
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Parse reads." << std::endl;
      typedef std::set<std::string> TReadSet;
      TReadSet reads;
      _parseReads(c, reads);
    
      // Get read mappings
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Extract reads." << std::endl;
      sequences(c, filename, reads);
    } else if (c.format == 1) filename = c.file.string();

    // Load sequences from disk for large contigs
    faidx_t* fai = fai_load(filename.c_str());
    for(int32_t idx1 = 0; idx1 < faidx_nseq(fai); ++idx1) {
      std::string seqname1(faidx_iseq(fai, idx1));
      int32_t xlen = faidx_seq_len(fai, seqname1.c_str());
      int32_t sl = 0;
      char* seq1 = faidx_fetch_seq(fai, seqname1.c_str(), 0, xlen, &sl);
      for(int32_t idx2 = idx1 + 1; idx2 < faidx_nseq(fai); ++idx2) {
	std::string seqname2(faidx_iseq(fai, idx2));
	int32_t ylen = faidx_seq_len(fai, seqname2.c_str());
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
	cv::Mat img(c.usedheight, c.usedwidth, CV_8UC3, cv::Scalar(255, 255, 255));

	// Compute word matches
	if (c.matchlen < 32) wordMatchShort(c, seq1, seq2, xlen, ylen, img);
	else wordMatchLong(c, seq1, seq2, xlen, ylen, img);
	
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
    
#ifdef PROFILE
    ProfilerStop();
#endif

    // Remove temporary file
    if (c.format == 0) {
      boost::filesystem::remove(filename);
      boost::filesystem::remove(filename + ".fai");
    }
    
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
      ("matchlen,m", boost::program_options::value<uint32_t>(&c.matchlen)->default_value(11), "default match length")
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
