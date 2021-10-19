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

#include "version.h"

namespace jellynas
{

  // Config arguments
  struct Config {
    bool islr;
    uint16_t minMapQual;
    uint16_t minTraQual;
    uint16_t minGenoQual;
    uint16_t madCutoff;
    uint16_t madNormalCutoff;
    int32_t nchr;
    int32_t minimumFlankSize;
    int32_t indelsize;
    uint32_t graphPruning;
    uint32_t minRefSep;
    uint32_t maxReadSep;
    uint32_t minClip;
    uint32_t maxGenoReadCount;
    uint32_t minCliqueSize;
    float flankQuality;
    bool hasExcludeFile;
    bool hasVcfFile;
    bool isHaplotagged;
    bool hasDumpFile;
    bool svtcmd;
    std::set<int32_t> svtset;
    boost::filesystem::path outfile;
    boost::filesystem::path vcffile;
    boost::filesystem::path genome;
    boost::filesystem::path exclude;
    boost::filesystem::path dumpfile;
    std::vector<boost::filesystem::path> files;
    std::vector<std::string> sampleName;
  };


  template<typename TConfigStruct>
  inline int jellyRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("delly.prof");
#endif

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
    c.isHaplotagged = false;
    c.madNormalCutoff = 5;
    c.islr = false;
    
    // Define generic options
    std::string svtype;
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("svtype,t", boost::program_options::value<std::string>(&svtype)->default_value("ALL"), "SV type to compute [DEL, INS, DUP, INV, BND, ALL]")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("exclude,x", boost::program_options::value<boost::filesystem::path>(&c.exclude), "file with regions to exclude")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("sv.bcf"), "SV BCF output file")
      ;
    
    boost::program_options::options_description disc("Discovery options");
    disc.add_options()
      ("map-qual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(1), "min. paired-end (PE) mapping quality")
      ("qual-tra,r", boost::program_options::value<uint16_t>(&c.minTraQual)->default_value(20), "min. PE quality for translocation")
      ("mad-cutoff,s", boost::program_options::value<uint16_t>(&c.madCutoff)->default_value(9), "insert size cutoff, median+s*MAD (deletions only)")
      ("minclip,c", boost::program_options::value<uint32_t>(&c.minClip)->default_value(25), "min. clipping length")
      ("min-clique-size,z", boost::program_options::value<uint32_t>(&c.minCliqueSize)->default_value(2), "min. PE/SR clique size")
      ("minrefsep,m", boost::program_options::value<uint32_t>(&c.minRefSep)->default_value(25), "min. reference separation")
      ("maxreadsep,n", boost::program_options::value<uint32_t>(&c.maxReadSep)->default_value(40), "max. read separation")
      ;
    
    boost::program_options::options_description geno("Genotyping options");
    geno.add_options()
      ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF/BCF file for genotyping")
      ("geno-qual,u", boost::program_options::value<uint16_t>(&c.minGenoQual)->default_value(5), "min. mapping quality for genotyping")
      ("dump,d", boost::program_options::value<boost::filesystem::path>(&c.dumpfile), "gzipped output file for SV-reads (optional)")
      ;

    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
      ("pruning,j", boost::program_options::value<uint32_t>(&c.graphPruning)->default_value(1000), "PE graph pruning cutoff")
      ("max-geno-count,a", boost::program_options::value<uint32_t>(&c.maxGenoReadCount)->default_value(250), "max. number of reads aligned for SR genotyping")
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
