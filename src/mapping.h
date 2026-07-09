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
  mappings(TConfig const& c, std::set<std::string> const& reads, std::map<std::string, std::vector<Mapping> >& mp, int32_t const indelSplit = 0) {
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
	  int32_t readlen = sequenceLength(rec); // Full read length
	  int32_t seqStart = -1;  // Match start
	  int32_t seqEnd = -1; // Match end
	  
	  // Emit the current match
	  auto flushSegment = [&]() {
	    if (gpStart < gpEnd) {
	      bool dir = true;
	      int32_t rs = seqStart;
	      int32_t re = seqEnd;
	      if (rec->core.flag & BAM_FREVERSE) {
		dir = false;
		rs = readlen - seqEnd;
		re = readlen - seqStart;
	      }
	      if (mp.find(qname) == mp.end()) mp[qname] = std::vector<Mapping>();
	      mp[qname].push_back(Mapping(rec->core.tid, gpStart, gpEnd, rs, re, dir, rec->core.qual));
	    }
	    gpStart = -1;
	    gpEnd = -1;
	    seqStart = -1;
	    seqEnd = -1;
	  };
	  
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
	      int32_t oplen = bam_cigar_oplen(cigar[i]);
	      if ((indelSplit > 0) && (oplen > indelSplit)) {
		// Large insertion
		flushSegment();
		sp += oplen;
	      } else {
		if (seqStart == -1) {
		  seqStart = sp;
		  gpStart = gp;
		}
		sp += oplen;
		seqEnd = sp;
		gpEnd = gp;
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	      int32_t oplen = bam_cigar_oplen(cigar[i]);
	      if ((indelSplit > 0) && (oplen > indelSplit)) {
		// Large deletion
		flushSegment();
		gp += oplen;
	      } else {
		if (seqStart == -1) {
		  seqStart = sp;
		  gpStart = gp;
		}
		gp += oplen;
		seqEnd = sp;
		gpEnd = gp;
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	      // Reference skip
	      flushSegment();
	      gp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else {
	      std::cerr << "Warning: Unknown Cigar options!" << std::endl;
	    }
	  }
	  flushSegment();
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

  // Aligned segments of primary record
  inline void
  emitPrimarySegments(bam1_t* rec, std::string const& qname, int32_t const readlen, std::map<std::string, std::vector<Mapping> >& mp) {
    uint32_t* cigar = bam_get_cigar(rec);
    int32_t gp = rec->core.pos;
    int32_t gpStart = -1;
    int32_t gpEnd = -1;
    int32_t sp = 0;
    int32_t seqStart = -1;
    int32_t seqEnd = -1;
    bool rev = (rec->core.flag & BAM_FREVERSE);
    auto flush = [&]() {
      if (gpStart < gpEnd) {
	bool dir = true;
	int32_t rs = seqStart;
	int32_t re = seqEnd;
	if (rev) {
	  dir = false;
	  rs = readlen - seqEnd;
	  re = readlen - seqStart;
	}
	if (mp.find(qname) == mp.end()) mp[qname] = std::vector<Mapping>();
	mp[qname].push_back(Mapping(rec->core.tid, gpStart, gpEnd, rs, re, dir, rec->core.qual));
      }
      gpStart = gpEnd = seqStart = seqEnd = -1;
    };
    for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
      int32_t op = bam_cigar_op(cigar[i]);
      int32_t oplen = bam_cigar_oplen(cigar[i]);
      if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) {
	if (seqStart == -1) {
	  seqStart = sp;
	  gpStart = gp;
	}
	gp += oplen;
	sp += oplen;
	seqEnd = sp;
	gpEnd = gp;
      } else if (op == BAM_CINS) {
	if (seqStart == -1) {
	  seqStart = sp;
	  gpStart = gp;
	}
	sp += oplen;
	seqEnd = sp;
	gpEnd = gp;
      } else if (op == BAM_CDEL) {
	if (seqStart == -1) {
	  seqStart = sp;
	  gpStart = gp;
	}
	gp += oplen;
	seqEnd = sp;
	gpEnd = gp;
      } else if (op == BAM_CSOFT_CLIP) {
	sp += oplen;
      }
      else if (op == BAM_CREF_SKIP) {
	flush();
	gp += oplen;
      }
      else if (op == BAM_CHARD_CLIP) {
	sp += oplen;
      }
    }
    flush();
  }

  // Parse SA tag
  inline void
  parseSATag(bam1_t* rec, bam_hdr_t* hdr, int32_t const readlen, std::string const& qname, std::map<std::string, std::vector<Mapping> >& mp) {
    uint8_t* satag = bam_aux_get(rec, "SA");
    if (satag == NULL) return;
    char* saz = bam_aux2Z(satag);
    if (saz == NULL) return;
    typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
    boost::char_separator<char> semi(";");
    std::string sa(saz);
    Tokenizer entries(sa, semi);
    for(Tokenizer::iterator eit = entries.begin(); eit != entries.end(); ++eit) {
      std::string entry = *eit;
      if (entry.empty()) continue;
      // fields: rname,pos,strand,CIGAR,mapQ,NM
      boost::char_separator<char> comma(",");
      Tokenizer fields(entry, comma);
      std::vector<std::string> f;
      for(Tokenizer::iterator fit = fields.begin(); fit != fields.end(); ++fit) f.push_back(*fit);
      if (f.size() < 5) continue;
      int32_t tid = bam_name2id(hdr, f[0].c_str());
      if (tid < 0) continue;
      int32_t gpStart = 0;
      try {
	gpStart = boost::lexical_cast<int32_t>(f[1]) - 1;
      } catch (...) {
	continue;
      }
      if (gpStart < 0) gpStart = 0;
      char strand = f[2].empty() ? '+' : f[2][0];
      uint16_t mapq = 0;
      try {
	mapq = (uint16_t) boost::lexical_cast<int32_t>(f[4]);
      } catch (...) {
	mapq = 0;
      }
      // Walk the CIGAR string
      int32_t refspan = 0, qspan = 0, lead = 0;
      bool aligned = false;
      uint32_t num = 0;
      bool haveNum = false;
      const std::string& cig = f[3];
      for(std::size_t k = 0; k < cig.size(); ++k) {
	char ch = cig[k];
	if ((ch >= '0') && (ch <= '9')) {
	  num = num * 10 + (uint32_t) (ch - '0');
	  haveNum = true;
	} else {
	  if (!haveNum) continue;
	  int32_t len = (int32_t) num; num = 0; haveNum = false;
	  switch (ch) {
	  case 'M': case '=': case 'X': refspan += len; qspan += len; aligned = true; break;
	  case 'I': qspan += len; aligned = true; break;
	  case 'D': case 'N': refspan += len; aligned = true; break;
	  case 'S': case 'H': if (!aligned) lead += len; break;
	  default: break;
	  }
	}
      }
      int32_t gpEnd = gpStart + refspan;
      int32_t seqStart = lead;
      int32_t seqEnd = lead + qspan;
      bool dir; int32_t rs, re;
      if (strand == '-') {
	dir = false;
	rs = readlen - seqEnd;
	re = readlen - seqStart;
      } else {
	dir = true;
	rs = seqStart;
	re = seqEnd;
      }
      if (rs < 0) rs = 0;
      if ((gpStart < gpEnd) && (rs < re)) {
	if (mp.find(qname) == mp.end()) mp[qname] = std::vector<Mapping>();
	mp[qname].push_back(Mapping(tid, gpStart, gpEnd, rs, re, dir, mapq));
      }
    }
  }

  // Region-limited read mappings + SA tag
  template<typename TConfig>
  inline void
  mappingsRegion(TConfig const& c, std::set<std::string> const& reads, std::vector<Region> const& scanRg, std::map<std::string, std::vector<Mapping> >& mp) {
    std::ofstream sfile;
    if (c.storeSequences) sfile.open(c.seqfile.string().c_str());

    samFile* samfile = sam_open(c.file.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.file.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    std::set<std::string> processed;  // a read may overlap several query regions
    for(uint32_t rgIdx = 0; rgIdx < scanRg.size(); ++rgIdx) {
      hts_itr_t* iter = sam_itr_queryi(idx, scanRg[rgIdx].tid, scanRg[rgIdx].beg, scanRg[rgIdx].end);
      if (iter == NULL) continue;
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	// Primary records only; supplementary segments come from the SA tag
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
	std::string qname = bam_get_qname(rec);
	if (reads.find(qname) == reads.end()) continue;
	if (processed.find(qname) != processed.end()) continue;
	if (sequenceLength(rec) < c.seqsize) continue;
	processed.insert(qname);
	int32_t readlen = sequenceLength(rec);

	// Store read sequence (original orientation)
	if (c.storeSequences) {
	  std::string sequence;
	  sequence.resize(rec->core.l_qseq);
	  uint8_t* seqptr = bam_get_seq(rec);
	  for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	  if (rec->core.flag & BAM_FREVERSE) reverseComplement(sequence);
	  sfile << ">" << qname << std::endl;
	  sfile << sequence << std::endl;
	}
	emitPrimarySegments(rec, qname, readlen, mp);
	parseSATag(rec, hdr, readlen, qname, mp);
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
    if (c.storeSequences) sfile.close();
  }

  // Select the topN longest reads
  template<typename TConfig>
  inline void
  selectTopReads(TConfig const& c, std::vector<Region> const& scanRg, std::set<std::string>& reads) {
    samFile* samfile = sam_open(c.file.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.file.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    std::map<std::string, int32_t> len;
    for(uint32_t rgIdx = 0; rgIdx < scanRg.size(); ++rgIdx) {
      hts_itr_t* iter = sam_itr_queryi(idx, scanRg[rgIdx].tid, scanRg[rgIdx].beg, scanRg[rgIdx].end);
      if (iter == NULL) continue;
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
	if ((rec->core.qual < c.minMapQual) || (rec->core.tid < 0)) continue;
	std::string qname = bam_get_qname(rec);
	int32_t sl = (int32_t) rec->core.l_qseq;
	std::map<std::string, int32_t>::iterator it = len.find(qname);
	if (it == len.end()) len.insert(std::make_pair(qname, sl));
	else if (sl > it->second) it->second = sl;
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    std::vector<std::pair<int32_t, std::string> > rl;
    for(std::map<std::string, int32_t>::const_iterator it = len.begin(); it != len.end(); ++it) rl.push_back(std::make_pair(it->second, it->first));
    std::sort(rl.begin(), rl.end(), [](std::pair<int32_t, std::string> const& a, std::pair<int32_t, std::string> const& b) { return a.first > b.first; });

    std::ofstream ofile;
    bool writeFile = !c.readFile.empty();
    if (writeFile) ofile.open(c.readFile.string().c_str());
    uint32_t cnt = 0;
    for(std::size_t i = 0; (i < rl.size()) && (cnt < c.topN); ++i, ++cnt) {
      reads.insert(rl[i].second);
      if (writeFile) ofile << rl[i].second << std::endl;
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Read " << (cnt + 1) << ": " << rl[i].second << " (" << rl[i].first << " bp)" << std::endl;
    }
    if (writeFile) ofile.close();
  }

}

#endif
