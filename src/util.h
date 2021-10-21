#ifndef UTIL_H
#define UTIL_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/math/distributions/binomial.hpp>

#include <htslib/sam.h>


namespace wallysworld
{

  inline unsigned hash_string(const char *s) {
    unsigned h = 37;
    while (*s) {
      h = (h * 54059) ^ (s[0] * 76963);
      s++;
    }
    return h;
  }

  inline std::size_t hash_read(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    boost::hash_combine(seed, (rec->core.flag & BAM_FREAD2));
    return seed;
  }
  
  inline std::size_t hash_pair(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    return seed;
  }

  inline std::size_t hash_pair_mate(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    return seed;
  }

  inline bool is_gff3(boost::filesystem::path const& f) {
    std::ifstream in(f.string().c_str());
    if (!in) return false;
    in.close();

    std::ifstream file(f.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    std::getline(instream, gline);
    bool gff = false;
    if ((gline.size()>=5) && (gline.substr(0,5) == "##gff")) gff = true;
    file.close();
    return gff;
  }
    
  
  inline bool is_gz(boost::filesystem::path const& f) {
    std::ifstream bfile(f.string().c_str(), std::ios_base::binary | std::ios::ate);
    bfile.seekg(0, std::ios::beg);
    char byte1;
    bfile.read(&byte1, 1);
    char byte2;
    bfile.read(&byte2, 1);
    bfile.close();
    if ((byte1 == '\x1F') && (byte2 == '\x8B')) return true;
    else return false;
  }
    

  // F+ 0
  // F- 1
  // R+ 2
  // R- 3
  inline uint8_t layout(bam1_t const* rec) {
    if (rec->core.flag & BAM_FREAD1) {
      if (!(rec->core.flag & BAM_FREVERSE)) {
	if (!(rec->core.flag & BAM_FMREVERSE)) return (rec->core.pos < rec->core.mpos) ? 0 : 1;
	else return (rec->core.pos < rec->core.mpos) ? 2 : 3;
      } else {
	if (!(rec->core.flag & BAM_FMREVERSE)) return (rec->core.pos > rec->core.mpos) ? 2 : 3;
	else return (rec->core.pos > rec->core.mpos) ? 0 : 1;
      }
    } else {
      if (!(rec->core.flag & BAM_FREVERSE)) {
	if (!(rec->core.flag & BAM_FMREVERSE)) return (rec->core.pos < rec->core.mpos) ? 1 : 0;
	else return (rec->core.pos < rec->core.mpos) ? 2 : 3;
      } else {
	if (!(rec->core.flag & BAM_FMREVERSE)) return (rec->core.pos > rec->core.mpos) ? 2 : 3;
	else return (rec->core.pos > rec->core.mpos) ? 1 : 0;
      }
    }
  }
  
  inline uint32_t alignmentLength(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t alen = 0;
    for (uint32_t i = 0; i < rec->core.n_cigar; ++i)
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF) || (bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) alen += bam_cigar_oplen(cigar[i]);
    return alen;
  }

  inline uint32_t
  lastAlignedPosition(bam1_t const* rec) {
    return rec->core.pos + alignmentLength(rec);
  }


  inline bool
  getSMTag(std::string const& header, std::string const& fileName, std::string& sampleName) {
    std::set<std::string> smIdentifiers;
    std::string delimiters("\n");
    typedef std::vector<std::string> TStrParts;
    TStrParts lines;
    boost::split(lines, header, boost::is_any_of(delimiters));
    TStrParts::const_iterator itH = lines.begin();
    TStrParts::const_iterator itHEnd = lines.end();
    bool rgPresent = false;
    for(;itH!=itHEnd; ++itH) {
      if (itH->find("@RG")==0) {
	std::string delim("\t ");
	TStrParts keyval;
	boost::split(keyval, *itH, boost::is_any_of(delim));
	TStrParts::const_iterator itKV = keyval.begin();
	TStrParts::const_iterator itKVEnd = keyval.end();
	for(;itKV != itKVEnd; ++itKV) {
	  size_t sp = itKV->find(":");
	  if (sp != std::string::npos) {
	    std::string field = itKV->substr(0, sp);
	    if (field == "SM") {
	      rgPresent = true;
	      std::string rgSM = itKV->substr(sp+1);
	      smIdentifiers.insert(rgSM);
	    }
	  }
	}
      }
    }
    if (!rgPresent) {
      sampleName = fileName;
      return true;
    } else if (smIdentifiers.size() == 1) {
      sampleName = *(smIdentifiers.begin());
      return true;
    } else {
      sampleName = "";
      return false;
    }
  }

}

#endif
