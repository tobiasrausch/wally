#ifndef HEATMAP_H
#define HEATMAP_H


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
#include "bed.h"

namespace wallysworld
{

  // Config arguments
  struct ConfigHeatmap {
    bool showWindow;
    bool showSupplementary;
    bool hasRegionFile;
    uint32_t minMapQual;
    uint32_t width;
    uint32_t height;
    double pxoffset; // 1bp in pixel    
    std::string regionStr;
    boost::filesystem::path genome;
    boost::filesystem::path regionFile;
    boost::filesystem::path file;
  };

  // Alignment
  struct Alignment {
    uint32_t alnstart;
    uint32_t alnend;
    bool reverse;
    Alignment(uint32_t const as, uint32_t const ae, bool const rev) : alnstart(as), alnend(ae), reverse(rev) {}
  };
		   

  template<typename TConfigStruct>
  inline int heatmapRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("wally.prof");
#endif

    // Open file handle
    samFile*  samfile = sam_open(c.file.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.file.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Read alignments
    typedef std::vector<Alignment> TAlignments;
    typedef std::map<unsigned, TAlignments> TReadAlign;
    TReadAlign rmap;	
    
    // Process regions
    std::vector<Region> rg;
    if (!parseRegions(hdr, c, rg)) return 1;
    if (rg.size() % 2 != 0) {
      std::cerr << "Number of regions needs to be an even number!" << std::endl;
      return 1;
    }
    for(uint32_t rgIdx = 0; rgIdx < rg.size(); ++rgIdx) {
      // Clear up map
      if (rgIdx % 2 == 0) rmap.clear();
      
      // Get pixel width of 1bp
      c.pxoffset = (1.0 / (double) rg[rgIdx].size) * (double) c.width;

      // Generate image
      cv::Mat hm( c.height, c.width, CV_8UC3, cv::Scalar(255, 255, 255));

      // Parse BAM files
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Region " << rg[rgIdx].id << std::endl;
	
      // Coverage
      uint32_t maxCoverage = std::numeric_limits<uint16_t>::max();

      hts_itr_t* iter = sam_itr_queryi(idx, rg[rgIdx].tid, rg[rgIdx].beg, rg[rgIdx].end);	
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY)) continue;
	if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	if ((!c.showSupplementary) && (rec->core.flag & BAM_FSUPPLEMENTARY)) continue;

	unsigned seed = hash_string(bam_get_qname(rec));
	if (rgIdx % 2 == 0) {
	  bool reverse = rec->core.flag & BAM_FREVERSE;
	  if (rec->core.flag & BAM_FREAD2) reverse = !reverse;  // Flip for read2
	  rmap[seed].push_back(Alignment(rec->core.pos, rec->core.pos + alignmentLength(rec), reverse));
	} else {
	  if (rmap.find(seed) != rmap.end()) {
	    // Common fragment
	    uint32_t alnend = rec->core.pos + alignmentLength(rec);
	    bool reverse = (rec->core.flag & BAM_FREVERSE);
	    if (rec->core.flag & BAM_FREAD2) reverse = !reverse;  // Flip for read2
	    for(uint32_t i = 0; i < rmap[seed].size(); ++i) {
	      std::cout << hdr->target_name[rg[rgIdx-1].tid] << ',' << rmap[seed][i].alnstart << ',' << rmap[seed][i].alnend << '-' << hdr->target_name[rg[rgIdx].tid] << ',' << rec->core.pos << ',' <<  alnend << std::endl;
	    }
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);

      // Store image (comment this for valgrind, png encoder seems leaky)
      /*
      if (c.splits == 1) {
	std::string outfile = rg[rgIdx].id;
	outfile += ".png";
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
	  std::string outfile = rg[rgIdx].id;
	  outfile += ".png";
	  cv::imwrite(outfile.c_str(), dst);
	  if (c.showWindow) {
	    cv::imshow(convertToStr(hdr, rg[rgIdx]).c_str(), dst);
	    cv::waitKey(0);
	  }
	  imageStore.clear();
	}
      }
      */
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


  int heatmap(int argc, char **argv) {
    ConfigHeatmap c;
    
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
    
    return heatmapRun(c);
  }

}

#endif
