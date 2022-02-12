#ifndef MATCHES_H
#define MATCHES_H


#include <iostream>
#include <fstream>

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
#include "draw.h"
#include "bed.h"

namespace wallysworld
{

  // Config arguments
  struct ConfigMatches {
    bool showWindow;
    bool hasAnnotationFile;
    bool showSoftClip;
    uint16_t splits;
    int32_t winsize;
    uint32_t width;
    uint32_t height;
    uint32_t tlheight;  // pixel height of a track line
    uint32_t rdheight;  // pixel height of a single read
    double pxoffset; // 1bp in pixel    
    boost::filesystem::path bedFile;
    boost::filesystem::path readFile;
    boost::filesystem::path genome;
    boost::filesystem::path file;
  };


  struct Mapping {
    int32_t tid;
    int32_t gstart;
    int32_t gend;
    int32_t rstart;
    int32_t rend;
    bool fwd;
    uint16_t qual;

    Mapping(int32_t const t, int32_t const gs, int32_t const ge, int32_t const rs, int32_t const re, bool const val, uint16_t const qval) : tid(t), gstart(gs), gend(ge), rstart(rs), rend(re), fwd(val), qual(qval) {}
  };

  template<typename TMapping>
  struct SortMappings : public std::binary_function<TMapping, TMapping, bool>
  {
    inline bool operator()(TMapping const& mp1, TMapping const& mp2) {
      return ((mp1.rstart < mp2.rstart) || ((mp1.rstart == mp2.rstart) && (mp1.rend < mp2.rend)));
    }
  };

  

  template<typename TConfig, typename TReads>
  inline void
  _parseReads(TConfig const& c, TReads& reads) {
    typedef typename TReads::value_type TValue;
    std::ifstream readFile(c.readFile.string().c_str(), std::ifstream::in);
    if (readFile.is_open()) {
      while (readFile.good()) {
	std::string line;
	getline(readFile, line);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t");
	Tokenizer tokens(line, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  reads.insert(boost::lexical_cast<TValue>(*tokIter));
	}
      }
      readFile.close();
    }
  }

  template<typename TConfig, typename TReadMappings>
  inline void
  clusterRegions(TConfig const& c, int32_t const nchr, TReadMappings const& mp, std::vector<Region>& rg) {
    for(int32_t refIndex = 0; refIndex < nchr; ++refIndex) {

      // Collect intervals for this chromosomes
      typedef boost::icl::interval_set<uint32_t> TChrIntervals;
      typedef typename TChrIntervals::interval_type TIVal;
      TChrIntervals cint;
      for(typename TReadMappings::const_iterator it = mp.begin(); it != mp.end(); ++it) {
	for(uint32_t i = 0; i < it->second.size(); ++i) {
	  if (it->second[i].tid == refIndex) {
	    int32_t gstart = 0;
	    if (c.winsize < it->second[i].gstart) gstart = it->second[i].gstart - c.winsize;
	    int32_t gend = it->second[i].gend + c.winsize;
	    cint.insert(TIVal::right_open(gstart, gend));
	  }
	}
      }

      // Keep clustered interval
      for(typename TChrIntervals::iterator it = cint.begin(); it != cint.end(); ++it) {
	rg.push_back(Region(refIndex, it->lower(), it->upper()));
      }
    }
  }
  

  template<typename TConfig>
  inline void
  mappings(TConfig const& c, std::set<std::string> const& reads, std::map<std::string, std::vector<Mapping> >& mp) {
    // Open file handles
    samFile* samfile = sam_open(c.file.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.file.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Parse BAM
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing... " << hdr->target_name[refIndex] << std::endl;

      // Iterate alignments 
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY)) continue;
	std::string qname = bam_get_qname(rec);
	if (reads.find(qname) != reads.end()) {
      
	  // Get read sequence
	  std::string sequence;
	  sequence.resize(rec->core.l_qseq);
	  uint8_t* seqptr = bam_get_seq(rec);
	  for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	
	  // Parse CIGAR
	  uint32_t* cigar = bam_get_cigar(rec);
	  int32_t gp = rec->core.pos; // Genomic position
	  int32_t gpStart = -1; //Match start
	  int32_t gpEnd = -1; //Match end
	  int32_t sp = 0; // Sequence position
	  int32_t seqStart = -1;  // Match start
	  int32_t seqEnd = -1; // Match end
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	      if (seqStart == -1) {
		seqStart = sp;
		gpStart = gp;
	      }
	      gp += bam_cigar_oplen(cigar[i]);
	      sp += bam_cigar_oplen(cigar[i]);
	      seqEnd = sp;
	      gpEnd = gp;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	      if (seqStart == -1) {
		seqStart = sp;
		gpStart = gp;
	      }
	      sp += bam_cigar_oplen(cigar[i]);
	      seqEnd = sp;
	      gpEnd = gp;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	      if (seqStart == -1) {
		seqStart = sp;
		gpStart = gp;
	      }
	      gp += bam_cigar_oplen(cigar[i]);
	      seqEnd = sp;
	      gpEnd = gp;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	      gp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else {
	      std::cerr << "Warning: Unknown Cigar options!" << std::endl;
	    }
	  }	  
	  bool dir = true;
	  if (rec->core.flag & BAM_FREVERSE) {
	    dir = false;
	    int32_t seqTmp = seqStart;
	    seqStart = sp - seqEnd;
	    seqEnd = sp - seqTmp;
	  }
	  if (gpStart < gpEnd) {
	    if (mp.find(qname) == mp.end()) mp[qname] = std::vector<Mapping>();
	    mp[qname].push_back(Mapping(rec->core.tid, gpStart, gpEnd, seqStart, seqEnd, dir, rec->core.qual));
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }

  

  template<typename TConfigStruct>
  inline int matchRun(TConfigStruct& c) {
#ifdef PROFILE
    ProfilerStart("wally.prof");
#endif

    // Parse reads
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Parse reads." << std::endl;
    std::set<std::string> reads;
    _parseReads(c, reads);
    
    // Get read mappings
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Extract read mappings." << std::endl;
    typedef std::vector<Mapping> TMappings;
    typedef std::map<std::string, TMappings > TReadMappings;
    TReadMappings mp;
    mappings(c, reads, mp);
    if (mp.empty()) {
      std::cerr << "Error: No mappings found! Are the read names correct?" << std::endl;
      return 1;
    }

    // Open file handles
    samFile* samfile = sam_open(c.file.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.file.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Sort mappings
    for(TReadMappings::iterator it = mp.begin(); it != mp.end(); ++it) {
      std::sort(it->second.begin(), it->second.end(), SortMappings<Mapping>());

      // Debug
      //for(uint32_t i = 0; i < it->second.size(); ++i) std::cerr << hdr->target_name[it->second[i].tid] << ':' << it->second[i].gstart << '-' << it->second[i].gend << '\t' << it->second[i].rstart << '-' << it->second[i].rend << '(' << (int) it->second[i].fwd << ')' << '\t' << it->first << std::endl;
    }

    // Determine regions
    std::vector<Region> rg;
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Cluster nearby matches." << std::endl;
    clusterRegions(c, hdr->n_targets, mp, rg);
    if (rg.empty()) {
      std::cerr << "Error: No regions to display!" << std::endl;
      return 1;
    }
    c.splits = rg.size();
    // Debug
    //for(uint32_t i = 0; i < rg.size(); ++i) std::cerr << hdr->target_name[rg[i].tid] << ':' << rg[i].beg << '-' << rg[i].end << std::endl;

    // Load genome
    faidx_t* fai = fai_load(c.genome.string().c_str());
    int32_t seqlen;
    std::string oldchr("None");
    char* seq = NULL;

    // Check image width
    c.width /= c.splits;
    if (c.width < 10) {
      std::cerr << "Image width is too small to display all matches, please increase!" << std::endl;
      return 1;
    }
    std::vector<cv::Mat> imageStore;

    // Split region
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
      //if (maxTracks < headerTracks + 10 * (int32_t) c.files.size()) {
      if (maxTracks < headerTracks + 10) {
	std::cerr << "Image height is too small!" << std::endl;
	return 1;
      }

      // Block header tracks
      std::vector<int32_t> taken(maxTracks, WALLY_UNBLOCK);
      for(int32_t i = 0; i < headerTracks; ++i) taken[i] = WALLY_BLOCKED;
      
      // Split by number of BAMs
      //int32_t bamTrackSize = (int32_t) ((maxTracks - headerTracks) / c.files.size());
      int32_t file_c = 0;
      int32_t bamTrackSize = (int32_t) (maxTracks - headerTracks);
          
      // Lazy loading of genome
      std::string chrName = hdr->target_name[rg[rgIdx].tid];
      if (chrName != oldchr) {
	if (seq != NULL) free(seq);
	seq = faidx_fetch_seq(fai, hdr->target_name[rg[rgIdx].tid], 0, hdr->target_len[rg[rgIdx].tid], &seqlen);
	oldchr = chrName;
      }

      // Header
      drawGenome(c, rg[rgIdx], bg, 1);

      // Reference
      drawReference(c, rg[rgIdx], bg, boost::to_upper_copy(std::string(seq + rg[rgIdx].beg, seq + rg[rgIdx].end)), 2);

      // Genes/Annotation
      drawAnnotation(c, rg[rgIdx], tr, anno, bg, 3);
      
      // Read offset
      int32_t genomicReadOffset  = (2.0 / (double) c.width) * (double) rg[rgIdx].size;
      
      // Parse BAM file
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Region " << rg[rgIdx].id << "; File " << c.file.stem().string() << std::endl;
	
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
      hts_itr_t* iter = sam_itr_queryi(idx, rg[rgIdx].tid, rg[rgIdx].beg, rg[rgIdx].end);	
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY)) continue;
	
	// Load sequence
	std::string sequence;
	sequence.resize(rec->core.l_qseq);
	uint8_t* seqptr = bam_get_seq(rec);
	for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	
	// Find out layout
	cv::Scalar readCol = WALLY_READ1;
	
	// Search empty track
	int32_t trackIdx = -1;
	int32_t alnend = rec->core.pos + alignmentLength(rec);
	trackIdx = firstEmptyTrack(taken, rec->core.pos);
	if (trackIdx != -1) taken[trackIdx] = alnend + genomicReadOffset;
	
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
      drawBorder(c, bg);

      
      // Store image (comment this for valgrind, png encoder seems leaky)
      if (c.splits == 1) {
	std::string outfile = "test.png";
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
	  std::string outfile = "test.png";
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
    hts_idx_destroy(idx);
    sam_close(samfile);

#ifdef PROFILE
    ProfilerStop();
#endif
  
    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    return 0;
  }


  int matches(int argc, char **argv) {
    ConfigMatches c;
    c.tlheight = 14;
    c.rdheight = 12;
    c.showSoftClip = true;
    
    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("bed,b", boost::program_options::value<boost::filesystem::path>(&c.bedFile), "BED annotation file (optional)")
      ;
    
    boost::program_options::options_description disc("Graphics options");
    disc.add_options()
      ("reads,r", boost::program_options::value<boost::filesystem::path>(&c.readFile), "file with reads to display")
      ("winsize,w", boost::program_options::value<int32_t>(&c.winsize)->default_value(10000), "window size to cluster nearby matches")
      ;
    
    boost::program_options::options_description geno("Display options");
    geno.add_options()
      ("width,x", boost::program_options::value<uint32_t>(&c.width)->default_value(1024), "width of the plot")
      ("height,y", boost::program_options::value<uint32_t>(&c.height)->default_value(1024), "height of the plot")
      ;

    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.file), "input file")
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
    
    return matchRun(c);
  }

}

#endif
