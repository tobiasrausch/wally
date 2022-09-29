#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>

#define BOOST_DISABLE_ASSERTS

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

#include "version.h"
#include "region.h"
#include "matches.h"
#include "dotplot.h"
#include "heatmap.h"
#include "hilbert.h"
#include "gfa.h"

using namespace wallysworld;

inline void
displayUsage() {
  std::cout << "Usage: wally <command> <arguments>" << std::endl;
  std::cout << std::endl;
  std::cout << "    region         plot genomic region" << std::endl;
  std::cout << "    matches        plot read or contig alignments" << std::endl;
  std::cout << "    dotplot        plot pairwise alignments" << std::endl;
  //std::cout << "    gfa            GFA plotting" << std::endl;
  std::cout << "    hilbert        plot genomic region as hilbert curve" << std::endl; 
  std::cout << std::endl;
  std::cout << std::endl;
}

int main(int argc, char **argv) {
    if (argc < 2) { 
      printTitle("Wally");
      displayUsage();
      return 0;
    }

    if ((std::string(argv[1]) == "version") || (std::string(argv[1]) == "--version") || (std::string(argv[1]) == "--version-only") || (std::string(argv[1]) == "-v")) {
      std::cout << "Wally version: v" << wallyVersionNumber << std::endl;
      std::cout << " using Boost: v" << BOOST_VERSION / 100000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100 << std::endl;
      std::cout << " using HTSlib: v" << hts_version() << std::endl;
      return 0;
    }
    else if ((std::string(argv[1]) == "help") || (std::string(argv[1]) == "--help") || (std::string(argv[1]) == "-h") || (std::string(argv[1]) == "-?")) {
      printTitle("Wally");
      displayUsage();
      return 0;
    }
    else if ((std::string(argv[1]) == "warranty") || (std::string(argv[1]) == "--warranty") || (std::string(argv[1]) == "-w")) {
      displayWarranty();
      return 0;
    }
    else if ((std::string(argv[1]) == "license") || (std::string(argv[1]) == "--license") || (std::string(argv[1]) == "-l")) {
      bsd();
      return 0;
    }
    else if ((std::string(argv[1]) == "region")) {
      return region(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "matches")) {
      return matches(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "dotplot")) {
      return dotplot(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "gfa")) {
      return gfa(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "shared")) {
      return heatmap(argc-1,argv+1);
    }
    else if ((std::string(argv[1]) == "hilbert")) {
      return hilbert(argc-1,argv+1);
    }
    std::cerr << "Unrecognized command " << std::string(argv[1]) << std::endl;
    return 1;
}

