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
  int wally_region(const char* bams, const char* genome, const char* region, int width, int height, int paired, int showSoftClip, int showSupplementary, int showCoverage, int delsize, int mod, int tlheight, int rdheight, int mapQual, const char* bed) {
    ConfigRegion c;
    c.showWindow = false;
    c.autoHeight = (height <= 0);
    c.showSoftClip = (showSoftClip != 0);
    c.showSupplementary = (showSupplementary != 0);
    c.showPairs = (paired != 0);
    c.hasRegionFile = false;
    c.showCoverage = (showCoverage != 0);
    c.modType = ((mod >= WALLY_MOD_NONE) && (mod <= WALLY_MOD_5HMC)) ? mod : WALLY_MOD_NONE;
    c.delsize = (delsize > 0) ? delsize : 1000;
    c.splits = 1;
    c.minMapQual = (mapQual >= 0) ? (uint32_t) mapQual : 1;
    c.width = (uint32_t) width;
    c.height = (uint32_t) height;
    c.tlheight = (tlheight > 0) ? (uint32_t) tlheight : 14;
    c.rdheight = (rdheight > 0) ? (uint32_t) rdheight : 12;
    if (c.rdheight > c.tlheight) c.rdheight = c.tlheight;
    c.snvcov = 10;
    c.snvvaf = 0.2;
    c.regionStr = std::string(region) + ":wallyplot";
    c.genome = boost::filesystem::path(genome);

    // Optional BED annotation
    std::string bedStr(bed ? bed : "");
    if (!bedStr.empty()) {
      c.bedFile = boost::filesystem::path(bedStr);
      c.hasAnnotationFile = true;
    } else {
      c.hasAnnotationFile = false;
    }

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

  // Prefetch
  EMSCRIPTEN_KEEPALIVE
  int wally_slice(const char* bams, const char* genome, const char* region, const char* outbams) {
    std::vector<std::string> inPaths;
    std::vector<std::string> outPaths;
    {
      std::stringstream s(bams);
      std::string line;
      while (std::getline(s, line)) {
	if (!line.empty()) inPaths.push_back(line);
      }
    }
    {
      std::stringstream s(outbams);
      std::string line;
      while (std::getline(s, line)) {
	if (!line.empty()) outPaths.push_back(line);
      }
    }
    if (inPaths.size() != outPaths.size()) return 2;

    for (std::size_t i = 0; i < inPaths.size(); ++i) {
      samFile* in = sam_open(inPaths[i].c_str(), "r");
      if (in == NULL) return 3;
      hts_set_fai_filename(in, genome);
      bam_hdr_t* hdr = sam_hdr_read(in);
      hts_idx_t* idx = sam_index_load(in, inPaths[i].c_str());
      if ((hdr == NULL) || (idx == NULL)) {
	if (idx) hts_idx_destroy(idx);
	if (hdr) bam_hdr_destroy(hdr);
	sam_close(in);
	return 4;
      }
      hts_itr_t* iter = sam_itr_querys(idx, hdr, region);
      if (iter == NULL) {
	hts_idx_destroy(idx);
	bam_hdr_destroy(hdr);
	sam_close(in);
	return 5;
      }

      samFile* out = sam_open(outPaths[i].c_str(), "wb");
      if ((out == NULL) || (sam_hdr_write(out, hdr) < 0)) {
	if (out) sam_close(out);
	hts_itr_destroy(iter);
	hts_idx_destroy(idx);
	bam_hdr_destroy(hdr);
	sam_close(in);
	return 6;
      }
      bam1_t* rec = bam_init1();
      while (sam_itr_next(in, iter, rec) >= 0) {
	if (sam_write1(out, hdr, rec) < 0) break;
      }
      bam_destroy1(rec);
      sam_close(out);
      hts_itr_destroy(iter);
      hts_idx_destroy(idx);
      bam_hdr_destroy(hdr);
      sam_close(in);

      if (sam_index_build(outPaths[i].c_str(), 0) < 0) return 7;
    }
    return 0;
  }

}
