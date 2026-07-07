#include <iostream>
#include <string>
#include <sstream>
#include <emscripten.h>

#include "version.h"
#include "region.h"

using namespace wallysworld;

extern "C" {

  // Region subcommand
  EMSCRIPTEN_KEEPALIVE
  int wally_region(const char* bams, const char* genome, const char* region, int width, int height, int paired, int showSoftClip, int showSupplementary, int showCoverage, int delsize) {
    ConfigRegion c;
    c.showWindow = false;
    c.showSoftClip = (showSoftClip != 0);
    c.showSupplementary = (showSupplementary != 0);
    c.showPairs = (paired != 0);
    c.hasRegionFile = false;
    c.hasAnnotationFile = false;
    c.showCoverage = (showCoverage != 0);
    c.modType = WALLY_MOD_NONE;
    c.delsize = (delsize > 0) ? delsize : 1000;
    c.splits = 1;
    c.minMapQual = 1;
    c.width = (uint32_t) width;
    c.height = (uint32_t) height;
    c.tlheight = 14;
    c.rdheight = 12;
    c.snvcov = 10;
    c.snvvaf = 0.2;
    c.regionStr = std::string(region) + ":wallyplot";
    c.genome = boost::filesystem::path(genome);

    // Split alignment list into individual files
    std::stringstream bamStream(bams);
    std::string bamPath;
    while (std::getline(bamStream, bamPath)) {
      if (!bamPath.empty()) c.files.push_back(boost::filesystem::path(bamPath));
    }
    if (c.files.empty()) {
      std::cerr << "wally error: no alignment files provided" << std::endl;
      return 2;
    }

    try {
      return wallyRun(c);
    } catch (std::exception const& e) {
      std::cerr << "wally error: " << e.what() << std::endl;
      return 2;
    }
  }

}
