#ifndef REGION_H
#define REGION_H


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
  struct ConfigRegion {
    bool showWindow;
    bool showSoftClip;
    bool showSupplementary;
    bool showPairs;
    bool hasRegionFile;
    bool hasAnnotationFile;
    uint16_t madCutoff;
    uint16_t splits;
    uint32_t minMapQual;
    uint32_t width;
    uint32_t height;
    uint32_t tlheight;  // pixel height of a track line
    uint32_t rdheight;  // pixel height of a single read
    uint32_t snvcov;
    float snvvaf;
    double pxoffset; // 1bp in pixel    
    std::string regionStr;
    boost::filesystem::path bedFile;
    boost::filesystem::path genome;
    boost::filesystem::path regionFile;
    std::vector<boost::filesystem::path> files;
  };



  inline bool
  readMM(bam1_t* rec) {
    uint8_t *mm = bam_aux_get(rec, "MM");
    if ((!mm) || (mm[0] != 'Z')) {
      std::cerr << "Read lacks MM tag or wrong data type (not 'Z')" << std::endl;
      return false;
    }
    uint8_t *ml = bam_aux_get(rec, "ML");
    if ((!ml) || (ml[0] != 'B') || (ml[1] != 'C')) {
      std::cerr << "Read lacks ML tag or wrong data type (not 'B')" << std::endl;
      return false;
    }
    std::string mmval;
    char* cpMM = (char *)mm+1;
    while(*cpMM) mmval += *cpMM++;

    std::cerr << mmval << std::endl;
    
    return true;
  }
  
  template<typename TConfigStruct>
  inline int wallyRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("wally.prof");
#endif

    // Get library parameters
    typedef std::vector<LibraryInfo> TLibInfo;
    TLibInfo sampleLib(c.files.size());
    getLibraryParams(c, sampleLib);
    //for(uint32_t i = 0; i < c.files.size(); ++i) std::cerr << sampleLib[i].rs << ',' << sampleLib[i].median << ',' << sampleLib[i].mad << std::endl;
    
    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);

    // Load genome
    faidx_t* fai = fai_load(c.genome.string().c_str());
    int32_t seqlen;
    std::string oldchr("None");
    char* seq = NULL;

    // Adjust image width by number of split images?
    c.width /= c.splits;
    if (c.width < 10) {
      std::cerr << "Image width is too small or too many splits!" << std::endl;
      return 1;
    }
    std::vector<cv::Mat> imageStore;
    
    // Split region
    std::vector<Region> rg;
    if (!parseRegions(hdr, c, rg)) return 1;
    for(uint32_t rgIdx = 0; rgIdx < rg.size(); ++rgIdx) {
      // Parse annotation
      std::vector<Transcript> tr;
      std::vector<Region> anno;
      if (!parseAnnotation(hdr, c, rg[rgIdx], tr, anno)) return 1;

      // Debug code
      //for(uint32_t i = 0; i < tr.size(); ++i) std::cerr << hdr->target_name[tr[i].rg.tid] << ':' << tr[i].rg.beg << '-' << tr[i].rg.end << '\t' << tr[i].rg.id << "\t" << tr[i].forward << std::endl;
      //for(uint32_t i = 0; i < anno.size(); ++i) std::cerr << hdr->target_name[anno[i].tid] << ':' << anno[i].beg << '-' << anno[i].end << '\t' << anno[i].id << std::endl;
      
      // Get pixel width of 1bp
      c.pxoffset = (1.0 / (double) rg[rgIdx].size) * (double) c.width;

      // Generate image
      cv::Mat bg( c.height, c.width, CV_8UC3, cv::Scalar(255, 255, 255));

      // Tracks
      int32_t headerTracks = 4;
      int32_t maxTracks = c.height / c.tlheight;
      if (maxTracks < headerTracks + 10 * (int32_t) c.files.size()) {
	std::cerr << "Image height is too small!" << std::endl;
	return 1;
      }

      // Block header tracks
      std::vector<int32_t> taken(maxTracks, WALLY_UNBLOCK);
      for(int32_t i = 0; i < headerTracks; ++i) taken[i] = WALLY_BLOCKED;
      
      // Split by number of BAMs
      int32_t bamTrackSize = (int32_t) ((maxTracks - headerTracks) / c.files.size());
          
      // Lazy loading of genome
      std::string chrName = hdr->target_name[rg[rgIdx].tid];
      if (chrName != oldchr) {
	if (seq != NULL) free(seq);
	seq = faidx_fetch_seq(fai, hdr->target_name[rg[rgIdx].tid], 0, hdr->target_len[rg[rgIdx].tid], &seqlen);
	oldchr = chrName;
	boost::to_upper(seq);
      }

      // Header
      drawGenome(c, rg[rgIdx], hdr->target_name[rg[rgIdx].tid], bg, 1);

      // Reference
      drawReference(c, rg[rgIdx], bg, boost::to_upper_copy(std::string(seq + rg[rgIdx].beg, seq + rg[rgIdx].end)), 2);

      // Genes/Annotation
      drawAnnotation(c, rg[rgIdx], tr, anno, bg, 3);
      
      // Read offset
      int32_t genomicReadOffset  = (2.0 / (double) c.width) * (double) rg[rgIdx].size;
      
      // Iterate files
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Parse BAM files
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Region " << rg[rgIdx].id << "; File " << c.files[file_c].stem().string() << std::endl;
	
	// Reserve tracks for other BAM files
	int32_t lowerBound = file_c * bamTrackSize + headerTracks;
	int32_t upperBound = (file_c + 1) * bamTrackSize + headerTracks;
	for(int32_t i = headerTracks; i < maxTracks; ++i) {
	  if (i < lowerBound) taken[i] = WALLY_BLOCKED;
	  else if (i >= upperBound) taken[i] = WALLY_BLOCKED;
	  else taken[i] = WALLY_UNBLOCK;
	}
	// Reserve coverage track
	for(int32_t i = lowerBound; i < lowerBound + 2; ++i) taken[i] = WALLY_BLOCKED;
	
	// Coverage
	uint32_t maxCoverage = std::numeric_limits<uint16_t>::max();
	std::vector<uint16_t> covA(rg[rgIdx].size, 0);
	std::vector<uint16_t> covC(rg[rgIdx].size, 0);
	std::vector<uint16_t> covG(rg[rgIdx].size, 0);
	std::vector<uint16_t> covT(rg[rgIdx].size, 0);
	
	// Read alignments
	int32_t lastAlignedPos = 0;
	std::set<std::size_t> lastAlignedPosReads;
	std::map<std::size_t, int32_t> assignedTrack;
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], rg[rgIdx].tid, rg[rgIdx].beg, rg[rgIdx].end);	
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY)) continue;
	  if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	  if ((!c.showSupplementary) && (rec->core.flag & BAM_FSUPPLEMENTARY)) continue;
	  
	  // Load sequence
	  std::string sequence;
	  sequence.resize(rec->core.l_qseq);
	  uint8_t* seqptr = bam_get_seq(rec);
	  for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	  // Modified bases
	  //if (!readMM(rec)) return 1;
	  
	  // Find out layout
	  cv::Scalar readCol = WALLY_READ1;
	  if (c.showPairs) {
	    if (rec->core.flag & BAM_FREAD2) readCol = WALLY_READ2;
	    uint8_t pl = layout(rec);
	    uint32_t minSep = 50;
	    if (pl == 0) {
	      if (std::abs((int) rec->core.pos - (int) rec->core.mpos) >= minSep) readCol = cv::Scalar(175, 175, 63);
	    } else if (pl == 1) {
	      if (std::abs((int) rec->core.pos - (int) rec->core.mpos) >= minSep) readCol = cv::Scalar(213, 100, 78);
	    } else if (pl == 2) {
	      if (std::abs(rec->core.isize) > sampleLib[file_c].median + c.madCutoff * sampleLib[file_c].mad) readCol = cv::Scalar(63, 63, 213);
	    } else if (pl == 3) {
	      if (std::abs((int) rec->core.pos - (int) rec->core.mpos) >= minSep) readCol = cv::Scalar(63, 175, 63);
	    } else if (pl == DELLY_SVT_TRANS + 0) {
	      readCol = cv::Scalar(153, 154, 251);
	    } else if (pl == DELLY_SVT_TRANS + 1) {
	      readCol = cv::Scalar(111, 191, 253);
	    } else if (pl == DELLY_SVT_TRANS + 2) {
	      readCol = cv::Scalar(0, 127, 255);
	    } else if (pl == DELLY_SVT_TRANS + 3) {
	      readCol = cv::Scalar(214, 178, 202);
	    }
	  }

	  // Search empty track
	  int32_t trackIdx = -1;
	  int32_t alnend = rec->core.pos + alignmentLength(rec);
	  if ((c.showPairs) && !(rec->core.flag & BAM_FSUPPLEMENTARY)) {
	    unsigned seed = hash_string(bam_get_qname(rec));
	    
	    // Clean-up the read store for identical alignment positions
	    if (rec->core.pos > lastAlignedPos) {
	      lastAlignedPosReads.clear();
	      lastAlignedPos = rec->core.pos;
	    }
	    if (_firstPairObs(rec, seed, lastAlignedPosReads)) {
	      // First observed read
	      trackIdx = firstEmptyTrack(taken, rec->core.pos);
	      if (trackIdx != -1) {
		if ((rec->core.tid == rec->core.mtid) && (rec->core.mpos < rg[rgIdx].end)) {
		  // Block until mate
		  taken[trackIdx] = std::max((int32_t) rec->core.mpos + genomicReadOffset, alnend + genomicReadOffset);
		  assignedTrack[seed] = trackIdx;
		  // draw the paired-end connector line
		  drawPELine(c, rg[rgIdx], bg, trackIdx, (alnend - rg[rgIdx].beg), (rec->core.mpos - rg[rgIdx].beg), readCol);
		} else {
		  // Treat like single-end
		  taken[trackIdx] = alnend + genomicReadOffset;
		}
	      }
	    } else {
	      // Second observed read
	      if (assignedTrack.find(seed) == assignedTrack.end()) {
		// Mate outside window
		trackIdx = firstEmptyTrack(taken, rec->core.pos);  
		if (trackIdx != -1) taken[trackIdx] = alnend + genomicReadOffset;
	      } else {
		// Mate inside window
		trackIdx = assignedTrack[seed];
		assignedTrack.erase(seed);
		taken[trackIdx] = std::max(alnend + genomicReadOffset, taken[trackIdx]);
	      }
	    }
	  } else {
	    trackIdx = firstEmptyTrack(taken, rec->core.pos);
	    if (trackIdx != -1) taken[trackIdx] = alnend + genomicReadOffset;
	  }

	  // Parse CIGAR
	  uint32_t rp = rec->core.pos; // reference pointer
	  uint32_t sp = 0; // sequence pointer
	  bool firstBox = true;
	  uint32_t leadingSC = 0;
	  // Parse the CIGAR
	  uint32_t* cigar = bam_get_cigar(rec);
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	      bool drawTriangle = false;
	      if (firstBox) {
		if (rec->core.flag & BAM_FREVERSE) drawTriangle = true;
		firstBox = false;
	      }
	      if ((rp + bam_cigar_oplen(cigar[i]) == (uint32_t) alnend) && (!(rec->core.flag & BAM_FREVERSE))) drawTriangle = true;
	      drawRead(c, rg[rgIdx], bg, trackIdx, (rp - rg[rgIdx].beg), (rp + bam_cigar_oplen(cigar[i]) - rg[rgIdx].beg), (rec->core.flag & BAM_FREVERSE), drawTriangle, readCol);
	      if (leadingSC > 0) {
		drawSC(c, rg[rgIdx], bg, trackIdx, (rp - rg[rgIdx].beg), leadingSC, true);
		leadingSC = 0;
	      }
	      // match or mismatch
	      for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
		// Increase coverage
		int32_t rpadj = (int32_t) rp - (int32_t) rg[rgIdx].beg;
		if ((rpadj >= 0) && (rpadj < (int32_t) rg[rgIdx].size)) {
		  if (sequence[sp] == 'A') {
		    if (covA[rpadj] < maxCoverage) ++covA[rpadj];
		  } else if (sequence[sp] == 'C') {
		    if (covC[rpadj] < maxCoverage) ++covC[rpadj];
		  } else if (sequence[sp] == 'G') {
		    if (covG[rpadj] < maxCoverage) ++covG[rpadj];
		  } else if (sequence[sp] == 'T') {
		    if (covT[rpadj] < maxCoverage) ++covT[rpadj];
		  }
		}
		// Draw nucleotide for mismatches
		if (rec->core.l_qseq) {
		  if (sequence[sp] != seq[rp]) {
		    drawNuc(c, rg[rgIdx], bg, trackIdx, (rp - rg[rgIdx].beg), (rp + 1 - rg[rgIdx].beg), sequence[sp], readCol);
		  }
		}
		++sp;
		++rp;
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	      drawDel(c, rg[rgIdx], bg, trackIdx, (rp - rg[rgIdx].beg), (rp + bam_cigar_oplen(cigar[i]) - rg[rgIdx].beg), bam_cigar_oplen(cigar[i]));
	      rp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	      drawIns(c, rg[rgIdx], bg, trackIdx, (rp - rg[rgIdx].beg), bam_cigar_oplen(cigar[i]));
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	      if (sp == 0) leadingSC = bam_cigar_oplen(cigar[i]);
	      else drawSC(c, rg[rgIdx], bg, trackIdx, (rp - rg[rgIdx].beg), bam_cigar_oplen(cigar[i]), false);
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if(bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	      if (sp == 0) leadingSC = bam_cigar_oplen(cigar[i]);
	      else drawSC(c, rg[rgIdx], bg, trackIdx, (rp - rg[rgIdx].beg), bam_cigar_oplen(cigar[i]), false);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	      drawRefSkip(c, rg[rgIdx], bg, trackIdx, (rp - rg[rgIdx].beg), (rp + bam_cigar_oplen(cigar[i]) - rg[rgIdx].beg));
	      rp += bam_cigar_oplen(cigar[i]);
	    } else {
	      std::cerr << "Unknown Cigar options" << std::endl;
	      return 1;
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);

	// Fill-in coverage track
	std::vector<bool> snp(rg[rgIdx].size, false);
	for(int32_t rp = rg[rgIdx].beg; rp < rg[rgIdx].end; ++rp) {
	  int32_t rpadj = (int32_t) rp - (int32_t) rg[rgIdx].beg;
	  uint32_t cumsum = covA[rpadj];
	  cumsum += covC[rpadj];
	  cumsum += covG[rpadj];
	  cumsum += covT[rpadj];
	  if (cumsum >= c.snvcov) {
	    if (seq[rp] == 'A') {
	      if (((double) covA[rpadj] / (double) cumsum) < (1 - c.snvvaf)) {
		snp[rpadj] = true;
	      }
	    } else if (seq[rp] == 'C') {
	      if (((double) covC[rpadj] / (double) cumsum) < (1 - c.snvvaf)) {
		snp[rpadj] = true;
	      }
	    } else if (seq[rp] == 'G') {
	      if (((double) covG[rpadj] / (double) cumsum) < (1 - c.snvvaf)) {
		snp[rpadj] = true;
	      }
	    } else if (seq[rp] == 'T') {
	      if (((double) covT[rpadj] / (double) cumsum) < (1 - c.snvvaf)) {
		snp[rpadj] = true;
	      }
	    }
	  }
	}
	drawCoverage(c, rg[rgIdx], bg, covA, covC, covG, covT, snp, lowerBound + 1);
	drawSampleLabel(c, lowerBound + 3, c.files[file_c].stem().string(), bg);
      }
      drawBorder(c, bg);

      
      // Store image (comment this for valgrind, png encoder seems leaky)
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
    }

    // Clean-up
    if (seq != NULL) free(seq);
    fai_destroy(fai);
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }

#ifdef PROFILE
    ProfilerStop();
#endif
  
    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    return 0;
  }


  int region(int argc, char **argv) {
    ConfigRegion c;
    c.tlheight = 14;
    c.rdheight = 12;
    
    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("bed,b", boost::program_options::value<boost::filesystem::path>(&c.bedFile), "BED annotation file (optional)")
      ;
    
    boost::program_options::options_description disc("Graphics options");
    disc.add_options()
      ("map-qual,q", boost::program_options::value<uint32_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("mad-cutoff,m", boost::program_options::value<uint16_t>(&c.madCutoff)->default_value(9), "insert size cutoff, median+m*MAD")
      ("snv-vaf,v", boost::program_options::value<float>(&c.snvvaf)->default_value(0.2), "min. SNV VAF")
      ("snv-cov,t", boost::program_options::value<uint32_t>(&c.snvcov)->default_value(10), "min. SNV coverage")
      ("region,r", boost::program_options::value<std::string>(&c.regionStr)->default_value("chrA:35-78"), "region to display")
      ("rfile,R", boost::program_options::value<boost::filesystem::path>(&c.regionFile), "BED file with regions to display")
      ("paired,p", "paired-end view")
      ("supplementary,u", "show supplementary alignments")
      ("clip,c", "show soft- and hard-clips")
      ;
    
    boost::program_options::options_description geno("Display options");
    geno.add_options()
      ("split,s", boost::program_options::value<uint16_t>(&c.splits)->default_value(1), "number of horizontal images")
      ("width,x", boost::program_options::value<uint32_t>(&c.width)->default_value(1024), "width of the plot")
      ("height,y", boost::program_options::value<uint32_t>(&c.height)->default_value(1024), "height of the plot")
      ;

    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
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
      std::cout << "Usage: wally " << argv[0] << " [OPTIONS] -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> ..." << std::endl;
      std::cout << visible_options << "\n";
      return 0;
    }

    // Show window?
    if (vm.count("window")) c.showWindow = true;
    else c.showWindow = false;

    // Soft-clips
    if (vm.count("clip")) c.showSoftClip = true;
    else c.showSoftClip = false;

    // Supplementary
    if (vm.count("supplementary")) c.showSupplementary = true;
    else c.showSupplementary = false;

    // Paired-end view
    if (vm.count("paired")) c.showPairs = true;
    else c.showPairs = false;

    // Check splits
    if (c.splits < 1) {
      std::cerr << "The number of horizontal images needs to be >=1." << std::endl;
      return 1;
    }
      
    // Region file
    if (vm.count("rfile")) {
      if (!(boost::filesystem::exists(c.regionFile) && boost::filesystem::is_regular_file(c.regionFile) && boost::filesystem::file_size(c.regionFile))) {
	std::cerr << "BED file with regions is missing: " << c.regionFile.string() << std::endl;
	return 1;
      }
      c.hasRegionFile = true;
    } else c.hasRegionFile = false;

        // BED file
    if (vm.count("bed")) {
      if (!(boost::filesystem::exists(c.bedFile) && boost::filesystem::is_regular_file(c.bedFile) && boost::filesystem::file_size(c.bedFile))) {
	std::cerr << "Genome annotation BED file is missing: " << c.bedFile.string() << std::endl;
	return 1;
      }
      c.hasAnnotationFile = true;
    } else c.hasAnnotationFile = false;
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "wally ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return wallyRun(c);
  }

}

#endif
