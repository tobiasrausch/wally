#ifndef MAPPING_H
#define MAPPING_H


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
#include "matchdraw.h"
#include "bed.h"

namespace wallysworld
{

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

  template<typename TConfig>
  inline void
  mappings(TConfig const& c, std::set<std::string> const& reads, std::map<std::string, std::vector<Mapping> >& mp) {
    // Sequence file
    std::ofstream sfile;
    if (c.storeSequences) sfile.open(c.seqfile.string().c_str());
    
    // Open file handles
    samFile* samfile = sam_open(c.file.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.file.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Parse BAM
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Processing... " << hdr->target_name[refIndex] << std::endl;

      // Iterate alignments 
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY)) continue;
	std::string qname = bam_get_qname(rec);
	if (reads.find(qname) != reads.end()) {
	  // Size check
	  if (sequenceLength(rec) < c.seqsize) continue;

	  // Get read sequence
	  if (c.storeSequences) {
	    if (!(rec->core.flag & (BAM_FSUPPLEMENTARY))) {
	      std::string sequence;
	      sequence.resize(rec->core.l_qseq);
	      uint8_t* seqptr = bam_get_seq(rec);
	      for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	      // Reverse?
	      if (rec->core.flag & BAM_FREVERSE) reverseComplement(sequence);
	      sfile << ">" << qname << std::endl;
	      sfile << sequence << std::endl;
	    }
	  }
	  
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
	      gp += bam_cigar_oplen(cigar[i]);
	      // Reset match
	      gpStart = -1;
	      gpEnd = -1;
	      seqStart = -1;
	      seqEnd = -1;
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

    // Close sequence file
    if (c.storeSequences) sfile.close();
  }
  
}

#endif
