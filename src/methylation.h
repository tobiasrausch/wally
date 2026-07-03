#ifndef METHYLATION_H
#define METHYLATION_H

#include <boost/algorithm/string.hpp>

#include <htslib/sam.h>

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cctype>

#include "util.h"

namespace wallysworld
{

  #ifndef WALLY_MOD_NONE
  #define WALLY_MOD_NONE 0
  #endif
  #ifndef WALLY_MOD_5MC
  #define WALLY_MOD_5MC 1
  #endif
  #ifndef WALLY_MOD_5HMC
  #define WALLY_MOD_5HMC 2
  #endif

  struct ModHit {
    int32_t pos;
    char code;
    uint8_t prob;
    bool rev;
    char base;

    ModHit(int32_t const p, char const c, uint8_t const r, bool const s, char const b) : pos(p), code(c), prob(r), rev(s), base(b) {}
  };

  inline char
  complementBase(char b) {
    switch (std::toupper(static_cast<unsigned char>(b))) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    default: return b;
    }
  }

  inline char
  modCode(int32_t const modType) {
    if (modType == WALLY_MOD_5MC) return 'm';
    if (modType == WALLY_MOD_5HMC) return 'h';
    return 0;
  }

  inline bool
  buildModProb(bam1_t* rec, char const targetCode, std::vector<int16_t>& methProb) {
    int32_t l = rec->core.l_qseq;
    methProb.assign(l, -1);
    if (l == 0) return false;

    uint8_t* mm_aux = bam_aux_get(rec, "MM");
    if ((!mm_aux) || (*mm_aux != 'Z')) return false;
    bool readRev = (rec->core.flag & BAM_FREVERSE);

    // Read sequence
    const uint8_t* seqptr = bam_get_seq(rec);
    std::string sequence(l, 'N');
    for (int32_t i = 0; i < l; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

    // Fwd sequence
    std::string fwdseq = sequence;
    if (readRev) reverseComplement(fwdseq);
    std::unordered_map<char, std::vector<int32_t> > basepos;
    for (int32_t i = 0; i < l; ++i) basepos[std::toupper(static_cast<unsigned char>(fwdseq[i]))].push_back(i);

    // Parse MM tag
    const char* mmstr = reinterpret_cast<const char*>(mm_aux + 1);
    std::string mmString(mmstr);
    std::vector<ModHit> modhits;
    std::vector<std::string> tokens;
    boost::split(tokens, mmString, boost::is_any_of(";"));
    for (const auto& tok : tokens) {
      if (tok.empty()) continue;
      std::size_t idx = 0;
      char base = tok[idx++];
      if (idx >= tok.size()) continue;
      char strand = tok[idx++];
      bool revMod = (strand == '-');
      std::string mod_codes;
      while ((idx < tok.size()) && (tok[idx] != ',')) {
	char ch = tok[idx++];
	if ((ch == '?') || (ch == '.')) continue;
	else if (std::isalpha(static_cast<unsigned char>(ch))) mod_codes.push_back(ch);
      }
      if ((idx < tok.size()) && (tok[idx] == ',')) {
	std::string pos_str = tok.substr(idx + 1);
	if (!pos_str.empty()) {
	  std::vector<std::string> pos_tokens;
	  boost::split(pos_tokens, pos_str, boost::is_any_of(","));
	  int32_t current = -1;
	  for (const auto& pt : pos_tokens) {
	    if (pt.empty()) continue;
	    int32_t delta = 0;
	    try {
	      delta = std::stoi(pt);
	    } catch (...) {
	      break;
	    }
	    current += delta + 1;
	    for (char mc : mod_codes) modhits.emplace_back(current, mc, static_cast<uint8_t>(0), revMod, base);
	  }
	}
      }
    }

    // Parse ML tag
    uint8_t* ml_aux = bam_aux_get(rec, "ML");
    if ((ml_aux) && (*ml_aux == 'B')) {
      char subtype = *reinterpret_cast<char*>(ml_aux + 1);
      if (subtype == 'C') {
	const uint8_t* p = ml_aux + 2;
	int32_t n = (int32_t)(p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24));
	const uint8_t* data = p + 4;
	int32_t assign = std::min<int32_t>(n, (int32_t)modhits.size());
	for (int32_t i = 0; i < assign; ++i) modhits[i].prob = data[i];
      }
    }

    // Map mods onto fwd positions
    for (const auto& mh : modhits) {
      if (mh.code != targetCode) continue;
      char ub = std::toupper(static_cast<unsigned char>(mh.base));
      char target_base = mh.rev ? complementBase(ub) : ub;
      auto it = basepos.find(target_base);
      if (it == basepos.end()) continue;
      const auto& occs = it->second;
      if ((mh.pos < 0) || (static_cast<std::size_t>(mh.pos) >= occs.size())) continue;
      int32_t fwdPos = occs[mh.pos];
      methProb[fwdPos] = (int16_t) mh.prob;
    }

    return true;
  }

}

#endif
