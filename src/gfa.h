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
#include "gfadraw.h"
#include "bed.h"

namespace wallysworld
{

  // Config arguments
  struct ConfigGfa {
    typedef std::vector<std::string> TChrNames;
    
    bool showWindow;
    uint32_t width;
    uint32_t height;
    uint32_t tlheight;  // pixel height of a tile
    uint32_t tlwidth;  // pixel width of a tile
    uint32_t nodeheight;  // pixel height of a node
    uint32_t nodewidth;  // pixel width of a node
    float ftscale;  // Font scale
    double pxoffset; // 1bp in pixel
    double bpoffset; // 1pixel in bp
    double lw; // line width
    std::vector<TChrNames> chrname;
    boost::filesystem::path gfafile;
    boost::filesystem::path seqfile;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
  };

  inline uint32_t
  numRanks(Graph const& g, SubGraph const& gsub) {
    std::set<uint32_t> ranks;
    ranks.insert(0); // The reference rank is always present
    for(uint32_t i = 0; i < gsub.segments.size(); ++i) {
      ranks.insert(g.segments[gsub.segments[i]].rank);
    }
    return ranks.size();
  }

  inline void
  subgraphTiles(Graph const& g, SubGraph const& gsub, uint32_t& numranks, uint32_t& mnodes) {
    numranks = numRanks(g, gsub);
    std::vector<uint32_t> ni(numranks, 0);
    for(uint32_t i = 0; i < gsub.segments.size(); ++i) ++ni[g.segments[gsub.segments[i]].rank];
    std::sort(ni.begin(), ni.end(), std::greater<uint32_t>());
    if (!ni.empty()) mnodes = ni[0];
  }
      
  
  template<typename TConfig>
  inline bool
  parseGfa(TConfig& c, Graph& g) {
    typedef std::map<std::string, uint32_t> TChrMap;
    typedef std::vector<TChrMap> TGraphChrMap;
    TGraphChrMap chrmap(1, TChrMap());
    typedef std::map<uint32_t, uint32_t> TRankMap;
    TRankMap rmap;
    rmap.insert(std::make_pair(0, 0));  // Always map rank 0 to rank 0 (reference) even if there are no nodes
    for(uint32_t refIndex = 0; refIndex < c.chrname[0].size(); ++refIndex) chrmap[0].insert(std::make_pair(c.chrname[0][refIndex], refIndex));
    
    // Segment map
    typedef std::map<std::string, uint32_t> TSegmentIdMap;
    TSegmentIdMap smap;
    
    // Segment FASTA sequences
    boost::filesystem::remove(c.seqfile.string());
    boost::filesystem::remove(c.seqfile.string() + ".fai");
    std::ofstream sfile;
    sfile.open(c.seqfile.string().c_str());
    uint64_t seqsize = 0;

    // Parse GFA
    uint32_t id_counter = 0;
    std::ifstream gfaFile(c.gfafile.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(gfaFile);
    std::istream instream(&dataIn);
    std::string gline;
    while(std::getline(instream, gline)) {
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
	      std::string sequence = *tokIter;
	      std::string chrn;
	      uint32_t rank = POS_UNDEF;
	      uint32_t tid = POS_UNDEF;
	      uint32_t pos = POS_UNDEF;
	      for(;tokIter != tokens.end(); ++tokIter) {
		// Optional fields
		boost::char_separator<char> kvsep(":");
		Tokenizer tokens(*tokIter, kvsep);
		Tokenizer::iterator tikv = tokens.begin();
		if (*tikv == "SN") {
		  ++tikv; ++tikv;
		  chrn = *tikv;
		} else if (*tikv == "SO") {
		  ++tikv; ++tikv;
		  pos = boost::lexical_cast<uint32_t>(*tikv);
		} else if (*tikv == "SR") {
		  ++tikv; ++tikv;
		  rank = boost::lexical_cast<uint32_t>(*tikv);
		}
	      }
	      // Remap ranks
	      if (rmap.find(rank) == rmap.end()) {
		// Create new rank level
		uint32_t newrank = rmap.size();
		rmap.insert(std::make_pair(rank, newrank));
	      }
	      rank = rmap[rank];
	      // Set chromosome names to NA and pos to 0 if not present
	      if (chrn.empty()) chrn = "NA";
	      if (pos == POS_UNDEF) pos = 0;
	      // rGFA
	      if ((rank < chrmap.size()) && (chrmap[rank].find(chrn) != chrmap[rank].end())) tid = chrmap[rank][chrn];
	      else {
		if (rank == 0) {
		  std::cerr << "Genome file does not match rank 0 chromosome names!" << std::endl;
		  std::cerr << chrn << " is not present in your genome file!" << std::endl;
		  return false;
		} else {
		  // Insert new chromosome
		  if (rank >= chrmap.size()) {
		    chrmap.resize(rank + 1, TChrMap());
		    c.chrname.resize(rank + 1, ConfigGfa::TChrNames());
		  }
		  tid = chrmap[rank].size();
		  chrmap[rank].insert(std::make_pair(chrn, tid));
		  c.chrname[rank].push_back(chrn);
		}
	      }
	      
	      // New segment
	      g.segments.push_back(Segment(rank, tid, pos, sequence.size()));
	      // Store sequence
	      sfile << ">" << id_counter << " " << segname << " " << chrn << ":" << pos << ":" << rank << std::endl;
	      sfile << sequence << std::endl;
	      seqsize += sequence.size();
	      // Keep segment name <-> id relationship
	      smap.insert(std::make_pair(segname, id_counter));
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
	  ++tokIter;
	  if (tokIter != tokens.end()) {
	    // From
	    if (smap.find(*tokIter) == smap.end()) {
	      std::cerr << "Link with unknown from segment! " << *tokIter << std::endl;
	      return false;
	    }
	    uint32_t fromId = smap[*tokIter];
	    ++tokIter;
	    if (tokIter != tokens.end()) {
	      // FromOrient
	      bool fromrev = false;
	      if (*tokIter == "-") fromrev = true;
	      ++tokIter;
	      if (tokIter != tokens.end()) {
		// To
		if (smap.find(*tokIter) == smap.end()) {
		  std::cerr << "Link with unknown to segment! " << *tokIter << std::endl;
		  return false;
		}
		uint32_t toId = smap[*tokIter];
		++tokIter;
		if (tokIter != tokens.end()) {
		  // ToOrient
		  bool torev = false;
		  if (*tokIter == "-") torev = true;
		  ++tokIter;
		  if (tokIter != tokens.end()) {
		    // Overlap CIGAR
		    if (*tokIter != "0M") {
		      std::cerr << "Currently only 0M links are supported!" << std::endl;
		      return false;
		    }
		    g.links.push_back(Link(fromrev, torev, fromId, toId));
		  }
		}
	      }
	    }
	  }
	} else {
	  // Todo
	  std::cerr << "Warning: Unknown line " << *tokIter << std::endl;
	  continue;
	}
      }
    }
    dataIn.pop();
    dataIn.pop();

    std::cerr << "Parsed: " << g.segments.size() << " segments, " << g.links.size() << " links" << std::endl;
    std::cerr << "Total sequence size: " << seqsize << std::endl;
    
    // Close FASTA file
    sfile.close();

    // Build index
    if (fai_build(c.seqfile.string().c_str())) {
      std::cerr << "Could not build FASTA index!" << std::endl;
      return false;
    }
    
    return true;
  }

  template<typename TConfig>
  inline void
  writeGfa(TConfig& c, Graph& g) {
    std::string filename = "test.out.gfa";
    
    // Output rGFA
    std::ofstream sfile;
    sfile.open(filename.c_str());

    // Output segments
    faidx_t* fai = fai_load(c.seqfile.string().c_str());
    for(uint32_t i = 0; i < g.segments.size(); ++i) {
      std::string seqid = boost::lexical_cast<std::string>(i);
      int32_t seqlen;
      char* seq = faidx_fetch_seq(fai, seqid.c_str(), 0, faidx_seq_len(fai, seqid.c_str()), &seqlen);
      sfile << "S\ts" << (i+1) << "\t" << seq;
      //sfile << "S\t" << (i+1) << "\t" << seq;
      sfile << "\tLN:i:" << g.segments[i].len;
      sfile << "\tSN:Z:" << c.chrname[g.segments[i].rank][g.segments[i].tid];
      sfile << "\tSO:i:" << g.segments[i].pos;
      sfile << "\tSR:i:" << g.segments[i].rank;
      sfile << std::endl;
      free(seq);
    }
    fai_destroy(fai);

    // Output links
    for(uint32_t i = 0; i < g.links.size(); ++i) {
      sfile << "L\ts" << (g.links[i].from+1);
      //sfile << "L\t" << (g.links[i].from+1);
      if (g.links[i].fromrev) sfile << "\t-";
      else sfile << "\t+";
      sfile << "\ts" << (g.links[i].to+1);
      //sfile << "\t" << (g.links[i].to+1);
      if (g.links[i].torev) sfile << "\t-";
      else sfile << "\t+";
      sfile << "\t0M";
      // From and to rank
      sfile << "\tFR:i:" << g.segments[g.links[i].from].rank;
      sfile << "\tTR:i:" << g.segments[g.links[i].to].rank << std::endl;
    }
    sfile.close();
  }
  
  
  template<typename TConfigStruct>
  inline int gfaRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("wally.prof");
#endif

    // Parse pangenome graph
    Graph g;
    if (!parseGfa(c, g)) {
      std::cerr << "GFA parsing failed!" << std::endl;
      return 1;
    }

    // Define subgraph for plotting
    SubGraph gsub;
    for(uint32_t i = 0; i < g.segments.size(); ++i) gsub.segments.push_back(i);
    for(uint32_t i = 0; i < g.links.size(); ++i) gsub.links.push_back(i);


    // Number of x- and y-tiles
    uint32_t numranks = 0;
    uint32_t mnodes = 0;
    subgraphTiles(g, gsub, numranks, mnodes); 
    // Generate image
    if (c.width == 0) {
      c.width = c.tlwidth * mnodes;
      c.height = c.tlheight * numranks;
    }
    cv::Mat img( c.height, c.width, CV_8UC3, cv::Scalar(255, 255, 255));

    // Draw nodes
    drawNodes(c, img, g, gsub, numranks);

    // Store image (comment this for valgrind, png encoder seems leaky)
    cv::imwrite(c.outfile.string().c_str(), img);
    if (c.showWindow) {
      cv::imshow(c.outfile.string().c_str(), img);
      cv::waitKey(0);
    }

    writeGfa(c, g);
    
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
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.png"), "output png file")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("trackheight,t", boost::program_options::value<uint32_t>(&c.tlheight)->default_value(100), "tile height in pixels")
      ("trackwidth,u", boost::program_options::value<uint32_t>(&c.tlwidth)->default_value(300), "tile width in pixels")
      ("ftscale,f", boost::program_options::value<float>(&c.ftscale)->default_value(0.4), "font scale")
      ("width,x", boost::program_options::value<uint32_t>(&c.width)->default_value(0), "width of the plot [0: best fit]")
      ("height,y", boost::program_options::value<uint32_t>(&c.height)->default_value(0), "height of the plot [0: best fit]")
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
      std::cout << "Usage: wally " << argv[0] << " [OPTIONS] -g <linear.genome.fa> <input.gfa.gz>" << std::endl;
      std::cout << visible_options << "\n";
      return 0;
    }

    // Set node parameters
    c.nodeheight = (int) (0.5 * c.tlheight);
    c.nodewidth = (int) (0.8 * c.tlwidth);

    // Set line width
    c.lw = 0.005 * c.tlheight;
    if (c.lw < 1) c.lw = 1;
    
    // Show window?
    if (vm.count("window")) c.showWindow = true;
    else c.showWindow = false;

    // Fill baseline genome map
    if (c.chrname.empty()) {
      faidx_t* fai = fai_load(c.genome.string().c_str());
      c.chrname.resize(1, ConfigGfa::TChrNames()); // Reference chromosomes, rank 0
      c.chrname[0].resize(faidx_nseq(fai));
      for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) c.chrname[0][refIndex] = std::string(faidx_iseq(fai, refIndex));
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
