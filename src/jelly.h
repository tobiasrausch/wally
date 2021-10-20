#ifndef JELLY_H
#define JELLY_H


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

#include "lodepng.h"
#include <iostream>

#include "version.h"

namespace jellynas
{

  // Config arguments
  struct Config {
    uint32_t minMapQual;
    uint32_t maxCov;
    uint32_t width;
    uint32_t height;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
    std::vector<boost::filesystem::path> files;
  };


  void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
    //Encode the image
    unsigned error = lodepng::encode(filename, image, width, height);

    //if there's an error, display it
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
  }

  
  template<typename TConfigStruct>
  inline int jellyRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("delly.prof");
#endif

    //NOTE: this sample will overwrite the file or test.png without warning!
    const char* filename = "test.png";

    //generate some image
    std::vector<unsigned char> image;
    image.resize(c.width * c.height * 4);
    for(uint32_t y = 0; y < c.height; y++)
      for(uint32_t x = 0; x < c.width; x++) {
	if (x == 100) {
	  image[4 * c.width * y + 4 * x + 0] = 255;
	  image[4 * c.width * y + 4 * x + 1] = 255;
	  image[4 * c.width * y + 4 * x + 2] = 255;
	  image[4 * c.width * y + 4 * x + 3] = 255;
	} else {
	  image[4 * c.width * y + 4 * x + 0] = 0;
	  image[4 * c.width * y + 4 * x + 1] = 0;
	  image[4 * c.width * y + 4 * x + 2] = 0;
	  image[4 * c.width * y + 4 * x + 3] = 255;
	}
      }
   
    encodeOneStep(filename, image, c.width, c.height);

#ifdef PROFILE
    ProfilerStop();
#endif
  
    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    return 0;
  }


  int jelly(int argc, char **argv) {
    Config c;
    
    // Define generic options
    std::string svtype;
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("region.png"), "output file")
      ;
    
    boost::program_options::options_description disc("Read selection options");
    disc.add_options()
      ("map-qual,q", boost::program_options::value<uint32_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("max-cov,m", boost::program_options::value<uint32_t>(&c.maxCov)->default_value(100), "max. coverage to display")
      ;
    
    boost::program_options::options_description geno("Graphics options");
    geno.add_options()
      ("width,x", boost::program_options::value<uint32_t>(&c.width)->default_value(512), "width of the plot")
      ("height,y", boost::program_options::value<uint32_t>(&c.height)->default_value(512), "height of the plot")
      ;

    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
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
      std::cout << "Usage: delly " << argv[0] << " [OPTIONS] -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> ..." << std::endl;
      std::cout << visible_options << "\n";
      return 0;
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "jelly ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return jellyRun(c);
  }

}

#endif
