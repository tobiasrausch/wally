#ifndef GFA_H
#define GFA_H


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
  struct ConfigGfa {
    bool showWindow;
    std::vector<std::string> chrname;
    std::map<std::string, int32_t> nchr;
    boost::filesystem::path gfafile;
    boost::filesystem::path seqfile;
    boost::filesystem::path genome;
  };

  struct Segment {
    uint32_t tid;
    uint32_t pos;
    uint32_t rank;
    std::string name;

    Segment(uint32_t const t, uint32_t const p, uint32_t const rk, std::string const& str) : tid(t), pos(p), rank(rk), name(str) {}
  };

  template<typename TConfig>
  inline bool
  parseGfa(TConfig& c, std::vector<Segment>& segments) {
    // Segment FASTA sequences
    std::ofstream sfile;
    sfile.open(c.seqfile.string().c_str());
    uint32_t id_counter = 0;
    std::ifstream gfaFile(c.gfafile.string().c_str(), std::ifstream::in);
    if (gfaFile.is_open()) {
      while (gfaFile.good()) {
	std::string gline;
	getline(gfaFile, gline);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep("\t");
	Tokenizer tokens(gline, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter != tokens.end()) {
	  // What element
	  if (*tokIter == "#") continue;
	  else if (*tokIter == "S") {
	    // Segment
	    ++tokIter;
	    if (tokIter != tokens.end()) {
	      // Name
	      std::string segname = *tokIter;
	      ++tokIter;
	      if (tokIter != tokens.end()) {
		// Sequence
		std::string chrname;
		uint32_t pos = POS_UNDEF;
		uint32_t rank = POS_UNDEF;
		for(;tokIter != tokens.end(); ++tokIter) {
		  // Optional fields
		  boost::char_separator<char> kvsep(":");
		  Tokenizer tokens(*tokIter, kvsep);
		  Tokenizer::iterator tikv = tokens.begin();
		  if (*tikv == "SN") {
		    ++tikv; ++tikv;
		    chrname = *tikv;
		  } else if (*tikv == "SO") {
		    ++tikv; ++tikv;
		    pos = boost::lexical_cast<uint32_t>(*tikv);
		  } else if (*tikv == "SR") {
		    ++tikv; ++tikv;
		    rank = boost::lexical_cast<uint32_t>(*tikv);
		  }
		}
		uint32_t tid = POS_UNDEF;
		if (c.nchr.find(chrname) != c.nchr.end()) tid = c.nchr[chrname];

		// New segment
		segments.push_back(Segment(tid, pos, rank, segname));
		// Store sequence
		sfile << ">" << id_counter << " " << segname << " " << chrname << ":" << pos << ":" << rank << std::endl;
		sfile << *tokIter << std::endl;
		++id_counter;
	      } else {
		std::cerr << "S segment lacks sequence information!" << std::endl;
		return false;
	      }
	    } else {
	      std::cerr << "S line lacks segment name!" << std::endl;
	      return false;
	    }
	  }
	  else if (*tokIter == "L") {
	    // Link
	  } else {
	    // Todo
	    std::cerr << "Warning: Unknown line " << *tokIter << std::endl;
	    continue;
	  }
	}
      }
      gfaFile.close();
    }
    sfile.close();
    return true;
  }
  
  template<typename TConfigStruct>
  inline int gfaRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("wally.prof");
#endif

    std::vector<Segment> segments;
    parseGfa(c, segments);
    
#ifdef PROFILE
    ProfilerStop();
#endif

    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    return 0;
  }


  int gfa(int argc, char **argv) {
    ConfigGfa c;

    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("fasta,f", boost::program_options::value<boost::filesystem::path>(&c.seqfile)->default_value("seq.fa"), "fasta file for segment sequences")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ;

    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.gfafile), "input file")
      ("window,w", "show window")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    // Set the visibility
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
      std::cout << std::endl;
      std::cout << "Usage: wally " << argv[0] << " [OPTIONS] -g <linear.genome.fa> <input.gfa>" << std::endl;
      std::cout << visible_options << "\n";
      return 0;
    }

    // Show window?
    if (vm.count("window")) c.showWindow = true;
    else c.showWindow = false;

    // Fill baseline genome map
    if (c.nchr.empty()) {
      faidx_t* fai = fai_load(c.genome.string().c_str());
      c.chrname.resize(faidx_nseq(fai));
      for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
	std::string chrName = faidx_iseq(fai, refIndex);
	c.nchr.insert(std::make_pair(chrName, refIndex));
	c.chrname[refIndex] = chrName;
      }
      fai_destroy(fai);
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "wally ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return gfaRun(c);
  }

}

#endif
