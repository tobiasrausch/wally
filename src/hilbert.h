#ifndef HILBERT_H
#define HILBERT_H


#include <iostream>
#include <fstream>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/timer/timer.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
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
#include "draw.h"
#include "hilty.h"

namespace wallysworld
{

  // Config arguments
  struct ConfigHilbert {
    bool showWindow;
    bool showSupplementary;
    bool hasRegionFile;
    uint32_t minMapQual;
    uint32_t order;
    uint32_t width;
    uint32_t height;
    std::string regionStr;
    boost::filesystem::path genome;
    boost::filesystem::path regionFile;
    boost::filesystem::path file;
  };

  template<typename TConfigStruct>
  inline int hilbertRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("wally.prof");
#endif

    // Debug hilbert
    //for(uint32_t d = 3; d < 10; ++d) {
    //int32_t n = std::pow(2, d);
    //for(uint32_t k = 0; k < n * n; ++k) {
    //int32_t x = 0;
    //	int32_t y = 0;
    //	posToHilbert(n, k, x, y);
    //	std::cerr << d << '\t' << k << '\t' << '(' << x << ',' << y << ')' << std::endl;
    //	int32_t hilt = hilbertToPos(n, x, y);
    //	std::cerr << d << '\t' << hilt << '\t' << '(' << x << ',' << y << ')' << std::endl;
    //}
    //}
    
    // Open file handle
    samFile*  samfile = sam_open(c.file.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.file.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Process regions
    std::vector<Region> rg;
    if (!parseRegions(hdr, c, rg)) return 1;
    for(uint32_t rgIdx = 0; rgIdx < rg.size(); ++rgIdx) {
      // Generate image
      cv::Mat hm( c.width, c.height, CV_8UC3, cv::Scalar(255, 255, 255));

      // Parse BAM files
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Region " << rg[rgIdx].id << std::endl;
	
      // Coverage
      //uint32_t maxCoverage = std::numeric_limits<uint16_t>::max();

      hts_itr_t* iter = sam_itr_queryi(idx, rg[rgIdx].tid, rg[rgIdx].beg, rg[rgIdx].end);	
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY)) continue;
	if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	if ((!c.showSupplementary) && (rec->core.flag & BAM_FSUPPLEMENTARY)) continue;
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);


      // Store image (comment this for valgrind, png encoder seems leaky)	
      std::string outfile = rg[rgIdx].id;
      outfile += ".png";
      cv::imwrite(outfile.c_str(), hm);
      if (c.showWindow) {
	cv::imshow(convertToStr(hdr, rg[rgIdx]).c_str(), hm);
	cv::waitKey(0);
      }
    }

    // Clean-up
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


  int hilbert(int argc, char **argv) {
    ConfigHilbert c;
    
    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ;
    
    boost::program_options::options_description disc("Graphics options");
    disc.add_options()
      ("map-qual,q", boost::program_options::value<uint32_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("region,r", boost::program_options::value<std::string>(&c.regionStr)->default_value("chrA:35-78,chrB:40-80"), "region to display (at least 2)")
      ("rfile,R", boost::program_options::value<boost::filesystem::path>(&c.regionFile), "BED file with regions to display")
      ("supplementary,u", "show supplementary alignments")
      ;
    
    boost::program_options::options_description geno("Display options");
    geno.add_options()
      ("order,d", boost::program_options::value<uint32_t>(&c.order)->default_value(7), "image size is 2^order * 2^order")
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
      std::cout << "Usage: wally " << argv[0] << " [OPTIONS] -g <ref.fa> <sample.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 0;
    }

    // Show window?
    if (vm.count("window")) c.showWindow = true;
    else c.showWindow = false;

    // Supplementary
    if (vm.count("supplementary")) c.showSupplementary = true;
    else c.showSupplementary = false;

    // Order
    if (c.order <= 1) {
      std::cerr << "Order has to be greater than 1!" << std::endl;
      return 1;
    }
    c.width = std::pow(2, c.order);
    c.height = std::pow(2, c.order);
    
    // Region file
    if (vm.count("rfile")) {
      if (!(boost::filesystem::exists(c.regionFile) && boost::filesystem::is_regular_file(c.regionFile) && boost::filesystem::file_size(c.regionFile))) {
	std::cerr << "BED file with regions is missing: " << c.regionFile.string() << std::endl;
	return 1;
      }
      c.hasRegionFile = true;
    } else c.hasRegionFile = false;

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "wally ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return hilbertRun(c);
  }

}

#endif
