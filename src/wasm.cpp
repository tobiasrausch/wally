#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include <unistd.h>
#include <emscripten.h>

#include "version.h"
#include "region.h"
#include "dotplot.h"

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

  // Dotplot
  EMSCRIPTEN_KEEPALIVE
  int wally_dotplot(const char* bams, const char* genome, const char* region, int numReads, int matchlen, double linewidth, int width, int flatten, int mapQual, int refTop) {
    // Use the first alignment file
    std::string bamPath;
    {
      std::stringstream s(bams);
      std::string line;
      while (std::getline(s, line)) { if (!line.empty()) { bamPath = line; break; } }
    }
    if (bamPath.empty()) {
      std::cerr << "wally error: no alignment file provided" << std::endl;
      return 2;
    }

    ConfigDotplot c;
    c.showWindow = false;
    c.hasReadFile = false;
    c.hasRegionFile = false;
    c.storeSequences = true;
    c.flatten = (flatten != 0);
    c.flip = false;
    c.incSelf = false;
    c.refTop = (refTop != 0);
    c.minMapQual = (mapQual >= 0) ? (uint32_t) mapQual : 1;
    c.matchlen = (matchlen > 0) ? (uint32_t) matchlen : 31;
    c.topN = (numReads > 0) ? (uint32_t) numReads : 10;
    c.seqsize = 0;
    c.width = (width > 0) ? (uint32_t) width : 0;
    c.height = 0;
    c.usedwidth = c.width;
    c.usedheight = c.height;
    c.tlheight = textSize().height + 2;
    c.lw = (linewidth > 0) ? (float) linewidth : 1.5f;
    c.format = 0;
    c.readStr = "";
    c.regionStr = std::string(region);
    c.readFile = boost::filesystem::path("/dpreads.txt");
    c.seqfile = boost::filesystem::path("/dpseq.fa");
    c.genome = boost::filesystem::path(genome);
    c.file = boost::filesystem::path(bamPath);

    // dotplot directory
    boost::filesystem::remove_all("/dot");
    boost::filesystem::create_directories("/dot");
    if (chdir("/dot") != 0) { std::cerr << "wally error: chdir /dot failed" << std::endl; return 7; }

    int32_t rc = 2;
    try {
      rc = dotplotRun(c);
    } catch (std::exception const& e) {
      std::cerr << "wally error: " << e.what() << std::endl;
      rc = 2;
    }
    if (chdir("/") != 0) { }
    return rc;
  }

}
