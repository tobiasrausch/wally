#include <iostream>
#include <string>
#include <emscripten.h>

#include "version.h"
#include "region.h"

using namespace wallysworld;

extern "C" {

  // Region subcommand
  EMSCRIPTEN_KEEPALIVE
  int wally_region(const char* bam, const char* genome, const char* region, int width, int height, int paired, int showSoftClip, int showSupplementary, int showCoverage) {
    ConfigRegion c;
    c.showWindow = false;
    c.showSoftClip = (showSoftClip != 0);
    c.showSupplementary = (showSupplementary != 0);
    c.showPairs = (paired != 0);
    c.hasRegionFile = false;
    c.hasAnnotationFile = false;
    c.showCoverage = (showCoverage != 0);
    c.modType = WALLY_MOD_NONE;
    c.madCutoff = 9;
    c.splits = 1;
    c.minMapQual = 1;
    c.width = (uint32_t) width;
    c.height = (uint32_t) height;
    c.tlheight = 14;
    c.rdheight = 12;
    c.snvcov = 10;
    c.snvvaf = 0.2;
    // Force a fixed output label so JS knows the filename: "/wallyplot.png"
    c.regionStr = std::string(region) + ":wallyplot";
    c.genome = boost::filesystem::path(genome);
    c.files.push_back(boost::filesystem::path(bam));

    try {
      return wallyRun(c);
    } catch (std::exception const& e) {
      std::cerr << "wally error: " << e.what() << std::endl;
      return 2;
    }
  }

}
