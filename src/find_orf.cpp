#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector internal_find_orf(const std::string seq, const IntegerVector start, const std::vector< std::string > stop_codon )
{
  IntegerVector collected(0);
  int last_end = 0;
  const int n_start = start.size();
  const int n_stop_codon = stop_codon.size();
  const int seq_len = seq.size();
  for(int j = 0; j < n_start; ++j){
    int codon_start = start[j];
    if(codon_start < last_end){
      continue;
    }
    bool no_stop = TRUE;
    while((codon_start + 3 <= seq_len) & no_stop){
      int codon_end = codon_start + 2;
      std::string cur_codon = seq.substr(codon_start, 3);
      //Rcout<<cur_codon<<"\t";
      for(int z = 0; z < n_stop_codon; ++z){
        if(cur_codon == stop_codon[z]){
          collected.push_back(start[j]);
          collected.push_back(codon_end);
          last_end = codon_end;
          no_stop = FALSE;
          break;
        }
      }
      codon_start = codon_start + 3;
    }
    if(no_stop){
      int codon_end = codon_start - 1;
      collected.push_back(start[j]);
      collected.push_back(codon_end);
      last_end = codon_end;
    }
  }
  return collected;
}

// [[Rcpp::export]]
List cxx_find_orf(std::vector< std::string > seqs, List starts, const std::vector< std::string > stop_codon){
  List out = List(seqs.size());
  for(int i = 0; i < seqs.size(); ++i){
    std::string cur_seq = seqs[i];
    IntegerVector cur_start = starts[i];
    out[i]=internal_find_orf(cur_seq, cur_start, stop_codon);
  }
  return out;
}
