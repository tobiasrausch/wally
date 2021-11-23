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

  #ifndef WALLY_PX
  #define WALLY_PX 8
  #endif
  
  #ifndef WALLY_A
  #define WALLY_A cv::Scalar(0,150,0)
  #endif

  #ifndef WALLY_C
  #define WALLY_C cv::Scalar(255,0,0)
  #endif

  #ifndef WALLY_G
  #define WALLY_G cv::Scalar(5,113,209)
  #endif

  #ifndef WALLY_T
  #define WALLY_T cv::Scalar(0,0,255)
  #endif

  #ifndef WALLY_BORDER
  #define WALLY_BORDER cv::Scalar(120,120,120)
  #endif

  #ifndef WALLY_BLOCKED
  #define WALLY_BLOCKED 1073741824
  #endif

  #ifndef WALLY_UNBLOCK
  #define WALLY_UNBLOCK -1073741824
  #endif
  
  struct Region {
    int32_t tid;
    int32_t beg;
    int32_t end;
    int32_t size;
    std::string id;
    
    Region() : tid(0), beg(0), end(0), size(0), id("") {}
  };


  struct LibraryInfo {
    int32_t rs;
    int32_t median;
    int32_t mad;

    LibraryInfo() : rs(0), median(0), mad(0) {}
  };

  template<typename TConfig, typename TSampleLibrary>
  inline void
  getLibraryParams(TConfig const& c, TSampleLibrary& sampleLib) {
    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<bam_hdr_t*> TSamHeader;
    TSamFile samfile(c.files.size());
    TSamHeader hdr(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      hdr[file_c] = sam_hdr_read(samfile[file_c]);
    }

    // Iterate all samples
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
      uint32_t maxAlignmentsScreened=10000;
      uint32_t minAlignmentsScreened=1000;
      uint32_t alignmentCount=0;
      typedef std::vector<uint32_t> TSizeVector;
      TSizeVector vecISize;
      TSizeVector readSize;

      bam1_t* rec = bam_init1();
      while (sam_read1(samfile[file_c], hdr[file_c], rec) >= 0) {
	if (!(rec->core.flag & BAM_FREAD2)) {
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	  if (alignmentCount > maxAlignmentsScreened) break;
	  ++alignmentCount;
	  readSize.push_back(rec->core.l_qseq);
	  // Paired-end?
	  if ((rec->core.flag & BAM_FPAIRED) && !(rec->core.flag & BAM_FMUNMAP) && (rec->core.tid==rec->core.mtid)) vecISize.push_back(abs(rec->core.isize));
	}
      }
      bam_destroy1(rec);

      // Get library parameters
      if (alignmentCount > minAlignmentsScreened) {
	std::sort(readSize.begin(), readSize.end());
	sampleLib[file_c].rs = readSize[readSize.size() / 2];

	if (vecISize.size() > (minAlignmentsScreened / 2)) {
	  std::sort(vecISize.begin(), vecISize.end());
	  int32_t median = vecISize[vecISize.size() / 2];
	  std::vector<uint32_t> absDev;
	  for(uint32_t i = 0; i < vecISize.size(); ++i) absDev.push_back(std::abs((int32_t) vecISize[i] - median));
	  std::sort(absDev.begin(), absDev.end());
	  int32_t mad = absDev[absDev.size() / 2];
	  sampleLib[file_c].median = median;
	  sampleLib[file_c].mad = mad;
	}
      } else {
	// Some defaults
	sampleLib[file_c].rs = 100;
	sampleLib[file_c].median = 200;
	sampleLib[file_c].mad = 30;
      }
    }

    // Clean-up
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      bam_hdr_destroy(hdr[file_c]);
      sam_close(samfile[file_c]);
    }
  }
  
  inline std::string
  convertToStr(bam_hdr_t* hdr, Region const& irg) {
    std::string str = hdr->target_name[irg.tid];
    str += ":" + boost::lexical_cast<std::string>(irg.beg + 1);
    str += "-" + boost::lexical_cast<std::string>(irg.end + 1);
    return str;
  }

  inline bool
  parseRegion(bam_hdr_t* hdr, std::string const& regionStr, Region& rg) {
    std::size_t pos = regionStr.find(":");
    if (pos == std::string::npos) {
      std::cerr << "Invalid region " << regionStr << std::endl;
      std::cerr << "No chromosome separator found ':'" << std::endl;
      return false;
    }
    std::string chrName(regionStr.substr(0, pos));
    std::string tmp = regionStr.substr(pos+1);
    pos = tmp.find("-");
    if (pos == std::string::npos) {
      std::cerr << "Invalid region " << regionStr << std::endl;
      std::cerr << "No position separator found '-'" << std::endl;
      return false;
    }
    rg.beg = boost::lexical_cast<int32_t>(tmp.substr(0, pos));
    tmp = tmp.substr(pos+1);
    pos = tmp.find(":");
    if (pos == std::string::npos) {
      rg.end = boost::lexical_cast<int32_t>(tmp);
      rg.id = chrName + "_" + boost::lexical_cast<std::string>(rg.beg) + "_" + boost::lexical_cast<std::string>(rg.end);
    } else {
      rg.end = boost::lexical_cast<int32_t>(tmp.substr(0, pos));
      rg.id = tmp.substr(pos+1);
    }
    rg.tid = bam_name2id(hdr, chrName.c_str());
    if (rg.tid < 0) {
      std::cerr << "Invalid region " << regionStr << std::endl;
      std::cerr << "Chromosome not found in BAM file." << std::endl;
      return false;
    }
    if (rg.beg >= rg.end) {
      std::cerr << "Invalid region " << regionStr << std::endl;
      std::cerr << "Region begin has to be smaller than region end." << std::endl;
      return false;
    }
    if (rg.end - rg.beg > 100000)  {
      std::cerr << "Invalid region " << regionStr << std::endl;
      std::cerr << "Region is larger than 100kbp." << std::endl;
      return false;
    }
    // Regions are 1-based, offset
    if (rg.beg > 0) {
      --rg.beg;
      --rg.end;
    } else {
      std::cerr << "Invalid region " << regionStr << std::endl;
      std::cerr << "Regions are 1-based." << std::endl;
      return false;
    }
    rg.size = rg.end - rg.beg;
    return true;
  }

  template<typename TConfig>
  inline bool
  parseRegions(bam_hdr_t* hdr, TConfig const& c, std::vector<Region>& rg) {
    if (c.hasRegionFile) {
      std::ifstream regionFile(c.regionFile.string().c_str(), std::ifstream::in);
      if (regionFile.is_open()) {
	while (regionFile.good()) {
	  std::string gline;
	  getline(regionFile, gline);
	  typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	  boost::char_separator<char> sep(" \t,;");
	  Tokenizer tokens(gline, sep);
	  Tokenizer::iterator tokIter = tokens.begin();
	  if (tokIter!=tokens.end()) {
	    std::string chrName=*tokIter++;
	    if (tokIter != tokens.end()) {
	      std::string start = *tokIter++;
	      if (tokIter != tokens.end()) {
		std::string end = *tokIter++;
		std::string str = chrName + ":" + start + "-" + end;
		if (tokIter != tokens.end()) str += ":" + std::string(*tokIter++);
		Region tmp;
		if (!parseRegion(hdr, str, tmp)) return false;
		rg.push_back(tmp);
	      }
	    }
	  }
	}
	regionFile.close();
      }
    } else {
      // Split multiple command-line regions
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep(",");
      Tokenizer tokens(c.regionStr, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      for(;tokIter != tokens.end(); ++tokIter) {
	std::string regStr = *tokIter;
	Region tmp;
	if (!parseRegion(hdr, regStr, tmp)) return false;
	rg.push_back(tmp);
      }
    }
    return true;
  }

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
    
  inline bool
  _translocation(bam1_t const* rec) {
    return (rec->core.tid != rec->core.mtid);
  }
  
  // F+ 0
  // F- 1
  // R+ 2
  // R- 3
  inline uint8_t
  layout(bam1_t const* rec) {
    if ((!_translocation(rec)) && (rec->core.flag & BAM_FPAIRED)) {
      if (!(rec->core.flag & BAM_FREVERSE)) {
	if (!(rec->core.flag & BAM_FMREVERSE)) return (rec->core.pos < rec->core.mpos) ? 0 : 1;
	else return (rec->core.pos < rec->core.mpos) ? 2 : 3;
      } else {
	if (!(rec->core.flag & BAM_FMREVERSE)) return (rec->core.pos > rec->core.mpos) ? 2 : 3;
	else return (rec->core.pos > rec->core.mpos) ? 0 : 1;
      }
    } else {
      if (_translocation(rec)) return 255;
      else return 255;
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
