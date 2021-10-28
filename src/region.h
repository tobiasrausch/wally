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

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include <iostream>

#include "version.h"
#include "util.h"
#include "draw.h"

namespace wallysworld
{

  // Config arguments
  struct Config {
    bool showWindow;
    uint32_t minMapQual;
    uint32_t width;
    uint32_t height;
    uint32_t tlheight;  // pixel height of a track line
    uint32_t rdheight;  // pixel height of a single read
    uint32_t snvcov;
    float snvvaf;
    double pxoffset; // 1bp in pixel    
    std::string regionStr;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
    boost::filesystem::path regionFile;
    std::vector<boost::filesystem::path> files;
  };

  template<typename TConfigStruct>
  inline int wallyRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("wally.prof");
#endif

    cv::Mat bg( c.height, c.width, CV_8UC3, cv::Scalar(255, 255, 255));

    // Tracks
    int32_t maxTracks = c.height / c.tlheight;
    std::vector<int32_t> taken(maxTracks, -2000000);
    if (maxTracks < 10) {
      std::cerr << "Image height is too small!" << std::endl;
      return 1;
    }
    // First 3 tracks: coverage
    taken[0] = 1073741824;
    taken[1] = 1073741824;
    taken[2] = 1073741824;
    taken[3] = 1073741824;
    taken[4] = 1073741824;
    taken[5] = 1073741824;
    
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

    // Parse BAM files
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Parsing BAMs" << std::endl;

    // Split region
    Region rg;
    if (!parseRegion(hdr, c.regionStr, rg)) return 1;
    c.pxoffset = (1.0 / (double) rg.size) * (double) c.width;

    // Load genome
    faidx_t* fai = fai_load(c.genome.string().c_str());
    int32_t seqlen;
    char* seq = faidx_fetch_seq(fai, hdr->target_name[rg.tid], 0, hdr->target_len[rg.tid], &seqlen);
    
    // Header
    drawGenome(c, rg, bg, 1);

    // Reference
    drawReference(c, rg, bg, boost::to_upper_copy(std::string(seq + rg.beg, seq + rg.end)), 2);
    
    // Coverage
    uint32_t maxCoverage = std::numeric_limits<uint16_t>::max();
    std::vector<uint16_t> covA(rg.size, 0);
    std::vector<uint16_t> covC(rg.size, 0);
    std::vector<uint16_t> covG(rg.size, 0);
    std::vector<uint16_t> covT(rg.size, 0);

    // Read offset
    int32_t genomicReadOffset = 0.05 * rg.size;
    if (genomicReadOffset > 5) genomicReadOffset = 5;
    else if (genomicReadOffset < 1) genomicReadOffset = 1;
    
    // Iterate files
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      // Read alignments
      hts_itr_t* iter = sam_itr_queryi(idx[file_c], rg.tid, rg.beg, rg.end);	
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
	if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;

	// Load sequence
	std::string sequence;
	sequence.resize(rec->core.l_qseq);
	uint8_t* seqptr = bam_get_seq(rec);
	for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	// Search empty track
	int32_t trackIdx = -1;
	for(uint32_t i = 0; i < taken.size(); ++i) {
	  if (taken[i] < rec->core.pos) {
	    trackIdx = i;
	    taken[i] = rec->core.pos + alignmentLength(rec) + genomicReadOffset;
	    break;
	  }
	}
	if (trackIdx == -1) continue;
	
	// Parse CIGAR
	uint32_t alnend = rec->core.pos + alignmentLength(rec);
	uint32_t rp = rec->core.pos; // reference pointer
	uint32_t sp = 0; // sequence pointer
	bool firstBox = true;
	// Parse the CIGAR
	uint32_t* cigar = bam_get_cigar(rec);
	for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	  if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	    bool drawTriangle = false;
	    if (firstBox) {
	      if (rec->core.flag & BAM_FREVERSE) drawTriangle = true;
	      firstBox = false;
	    }
	    if ((rp + bam_cigar_oplen(cigar[i]) == alnend) && (!(rec->core.flag & BAM_FREVERSE))) drawTriangle = true;
	    drawRead(c, rg, bg, trackIdx, (rp - rg.beg), (rp + bam_cigar_oplen(cigar[i]) - rg.beg), (rec->core.flag & BAM_FREVERSE), drawTriangle);
	    // match or mismatch
	    for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
	      // Increase coverage
	      int32_t rpadj = (int32_t) rp - (int32_t) rg.beg;
	      if ((rpadj >= 0) && (rpadj < (int32_t) rg.size)) {
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
		  drawNuc(c, rg, bg, trackIdx, (rp - rg.beg), (rp + 1 - rg.beg), sequence[sp]);
		}
	      }
	      ++sp;
	      ++rp;
	    }
	  } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	    drawDel(c, rg, bg, trackIdx, (rp - rg.beg), (rp + bam_cigar_oplen(cigar[i]) - rg.beg), bam_cigar_oplen(cigar[i]));
	    rp += bam_cigar_oplen(cigar[i]);
	  } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	    // ToDo
	    sp += bam_cigar_oplen(cigar[i]);
	  } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	    sp += bam_cigar_oplen(cigar[i]);
	  } else if(bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	  } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
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
      std::vector<bool> snp(rg.size, false);
      for(int32_t rp = rg.beg; rp < rg.end; ++rp) {
	int32_t rpadj = (int32_t) rp - (int32_t) rg.beg;
	uint32_t cumsum = covA[rpadj];
	cumsum += covC[rpadj];
	cumsum += covG[rpadj];
	cumsum += covT[rpadj];
	if (cumsum >= c.snvcov) {
	  if ((seq[rp] == 'a') || (seq[rp] == 'A')) {
	    if (((double) covA[rpadj] / (double) cumsum) < (1 - c.snvvaf)) {
	      snp[rpadj] = true;
	    }
	  } else if ((seq[rp] == 'c') || (seq[rp] == 'C')) {
	    if (((double) covC[rpadj] / (double) cumsum) < (1 - c.snvvaf)) {
	      snp[rpadj] = true;
	    }
	  } else if ((seq[rp] == 'g') || (seq[rp] == 'G')) {
	    if (((double) covG[rpadj] / (double) cumsum) < (1 - c.snvvaf)) {
	      snp[rpadj] = true;
	    }
	  } else if ((seq[rp] == 't') || (seq[rp] == 'T')) {
	    if (((double) covT[rpadj] / (double) cumsum) < (1 - c.snvvaf)) {
	      snp[rpadj] = true;
	    }
	  }
	}
      }
      drawCoverage(c, rg, bg, covA, covC, covG, covT, snp, 4);
    }

    // Store image
    cv::imwrite(c.outfile.string().c_str(), bg);
    if (c.showWindow) {
      std::string str = hdr->target_name[rg.tid];
      str += ":" + boost::lexical_cast<std::string>(rg.beg + 1);
      str += "-" + boost::lexical_cast<std::string>(rg.end + 1);
      cv::imshow(str.c_str(), bg);
      cv::waitKey(0);
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
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    return 0;
  }


  int region(int argc, char **argv) {
    Config c;
    c.tlheight = 14;
    c.rdheight = 12;
    
    // Define generic options
    std::string svtype;
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("region.jpg"), "output file")
      ;
    
    boost::program_options::options_description disc("Graphics options");
    disc.add_options()
      ("map-qual,q", boost::program_options::value<uint32_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("snv-vaf,s", boost::program_options::value<float>(&c.snvvaf)->default_value(0.2), "min. SNV VAF")
      ("snv-cov,t", boost::program_options::value<uint32_t>(&c.snvcov)->default_value(10), "min. SNV coverage")
      ("region,r", boost::program_options::value<std::string>(&c.regionStr)->default_value("chrA:35-78"), "region to display")
      ("rfile,R", boost::program_options::value<boost::filesystem::path>(&c.regionFile), "BED file with regions to display")
      ;
    
    boost::program_options::options_description geno("Display options");
    geno.add_options()
      ("width,x", boost::program_options::value<uint32_t>(&c.width)->default_value(1024), "width of the plot")
      ("height,y", boost::program_options::value<uint32_t>(&c.height)->default_value(1024), "height of the plot")
      ("window,w", "show window")
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
      std::cout << "Usage: wally " << argv[0] << " [OPTIONS] -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> ..." << std::endl;
      std::cout << visible_options << "\n";
      return 0;
    }

    // Show window?
    if (vm.count("window")) c.showWindow = true;
    else c.showWindow = false;
    
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
