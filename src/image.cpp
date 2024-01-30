#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <opencv2/core.hpp>
#include <opencv2/core/bindings_utils.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <emscripten.h>
#include <emscripten/bind.h>

#include <boost/algorithm/string.hpp>

#include <htslib/faidx.h>

namespace ems = emscripten;

int main() {
  std::string genome = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa";
  faidx_t* fai = fai_load(genome.c_str());
  int32_t seqlen;
  char* seq = faidx_fetch_seq(fai, "chr1", 0, 5, &seqlen);
  std::string test = boost::to_upper_copy(seq + 0, seq + seqlen);
  if (seq != NULL) free(seq);
  fai_destroy(fai);
  std::cout << test << std::endl;
  return 0;
}
