/*
  example modified from http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html

  compile and run with:
  g++ vcf.cpp -I../htslib/ ../htslib/libhts.a -lz -D__WITH_MAIN__

./a.out my.bcf
 */

#include <stdio.h>
#include <vector>
#include <htslib/vcf.h>
#include <cmath>
#include <limits>
#include <string>

void usage() {
        puts(
        "NAME\n"
        "    03_vcf - High quality calls for single sample\n"
        "SYNOPSIS\n"
        "    03_vcf vcf_file sample\n"
        "DESCRIPTION\n"
        "    Given a <vcf file>, extract all calls for <sample> and filter for\n"
        "    high quality, homozygous SNPs. This will omit any positions\n"
        "    that are homozygous ref (0/0) or heterozygous.  The exact filter\n"
        "    used is\n"
        "        FI == 1 & GQ > 20 & GT != '0/0'\n"
        "        [NOTE: FI == 1 implies homozygous call]\n"
        "        \n"
        "    The returned format is\n"
        "        chrom pos[0-based]  REF  ALT GQ|DP\n"
        "    and can be used for Marei's personalizer.py\n"
                );
}

// https://github.com/samtools/htslib/blob/57fe419344cb03e2ea46315443abd242489c32f2/vcf.c#L53
// uint32_t bcf_float_missing    = 0x7F800001;
// uint32_t bcf_float_vector_end = 0x7F800002;

const int PHREDMAX=256;
float pl2ln[PHREDMAX];

// double pl2ln[256];
float pl2ln_f(int32_t & val){
  if(val>=PHREDMAX){
    return log(pow(10.0,-0.1*val));
  } else {
    return pl2ln[val];
  }
    
}

template <class T>
bool same_val_vcf(T a, T b) {
  return std::fabs(a - b) < std::numeric_limits<T>::epsilon();  
}


// https://en.cppreference.com/w/c/numeric/math/isnan
bool is_nan_vcf(double x) { return x != x; }

//from angsd
double emFrequency(double *loglike,int numInds, int iter,double start,char *keep,int keepInd){

  if(keepInd == 0)
    return 0.0;
  
     
  float W0;
  float W1;
  float W2;
  // fprintf(stderr,"start=%f\n",start);
  float p=(float)start;
  float temp_p=(float)start;
  double accu=0.00001;
  double accu2=0;
  float sum;

  
  int it=0;
  
  for(it=0;it<iter;it++){
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=exp(loglike[i*3+0])*pow(1-p,2);
      W1=exp(loglike[i*3+1])*2*p*(1-p);
      W2=exp(loglike[i*3+2])*(pow(p,2));
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      if(std::isnan(sum))
	fprintf(stderr,"PRE[%d]:gls:(%f,%f,%f) W(%f,%f,%f) sum=%f\n",i,loglike[i*3],loglike[i*3+1],loglike[i*3+2],W0,W1,W2,sum);
    }
    
    p=sum/keepInd;
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
    temp_p=p;
  }
  
  if(std::isnan(p)){
    fprintf(stderr,"[%s] caught nan will not exit\n",__FUNCTION__);
    fprintf(stderr,"logLike (3*nInd). nInd=%d\n",numInds);
    //print_array(stderr,loglike,3*numInds);
    fprintf(stderr,"keepList (nInd)\n");
    //print_array(stderr,keep,numInds);
    fprintf(stderr,"used logLike (3*length(keep))=%d\n",keepInd);
    
    for(int ii=0;1&&ii<numInds;ii++){
      if(keep!=NULL && keep[ii]==1){
	//	fprintf(stderr,"1\t");
	for(int gg=0;gg<3;gg++)
	  fprintf(stderr,"%f\t",loglike[ii*3+gg]);
	fprintf(stderr,"\n");
      }
    }
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=exp(loglike[i*3+0])*pow(1-p,2);
      W1=exp(loglike[i*3+1])*2*p*(1-p);
      W2=exp(loglike[i*3+2])*(pow(p,2));
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      fprintf(stderr,"[%s.%s():%d] p=%f W %f\t%f\t%f sum=%f loglike: %f\n",__FILE__,__FUNCTION__,__LINE__,p,W0,W1,W2,sum,exp(loglike[i*3+2])*pow(1-p,2));
      break;
    }
    p=-999;
    assert(p!=999);
    return p;
  }

  return(p);
}


size_t getgls(char*fname,std::vector<double *> &mygl, std::vector<double> &freqs,int minind,double minfreq, std::string &vcf_format_field, std::string &vcf_allele_field,std::vector<char *> &posinfo){
  for(int i=0;i<PHREDMAX;i++){    
    pl2ln[i] = log(pow(10.0,-0.1*i));
  }
  //   http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html
  // counters
  int n    = 0;  // total number of records in file
  int nsnp = 0;  // number of SNP records in file
  int nseq = 0;  // number of sequences
  int nsamples = 0;

  // pl data for each call
  int npl_arr = 0;
  int npl     = 0;
  int32_t *pl = NULL;

   // gt data for each call
  int32_t ngt_arr = 0;
  int32_t ngt     = 0;
  int32_t *gt     = NULL;

  
  // af1/af data for each call
  int naf_arr = 0;
  int naf     = 0;
  float *af     = NULL;

  
  htsFile * inf = hts_open(fname, "r");
  if (inf == NULL) {
    fprintf(stderr,"bcf_open is null:");
    exit(0);
  }
  
  // read header
  bcf_hdr_t *hdr = bcf_hdr_read(inf);
  nsamples = bcf_hdr_nsamples(hdr);
  fprintf(stderr, "\t-> File %s contains %i samples\n", fname, nsamples);
  const char **seqnames = NULL;
  seqnames = bcf_hdr_seqnames(hdr, &nseq);
#if __WITH_MAIN__

  if (seqnames == NULL) {
    fprintf(stderr," error1\n");
    exit(0);
  }
  //        fprintf(stderr, "Sequence names:\n");
  for (int i = 0;0&& i < nseq; i++) {
    // bcf_hdr_id2name is another way to get the name of a sequence
    fprintf(stderr, "  [%2i] %s (bcf_hdr_id2name -> %s)\n", i, seqnames[i],
	    bcf_hdr_id2name(hdr, i));
  }
#endif
  char *chr;
  // struc for storing each record
  bcf1_t *rec = bcf_init();
  if (rec == NULL) {
    fprintf(stderr,"\t-> problem making bcf1_t\n");
    exit(0);
  }
  
  while (bcf_read(inf, hdr, rec) == 0) {
    n++;
    if (!bcf_is_snp(rec))
      continue;

    nsnp++;

    if(rec->n_allele>=3||rec->n_allele==1)//last case shouldt happen
      continue;

    float ln_gl[3*nsamples];    

    if(vcf_format_field == "PL") {
      npl = bcf_get_format_int32(hdr, rec, "PL", &pl, &npl_arr);
      if(npl<0){
        // return codes: https://github.com/samtools/htslib/blob/bcf9bff178f81c9c1cf3a052aeb6cbe32fe5fdcc/htslib/vcf.h#L667
        // no PL tag is available
        fprintf(stderr, "BAD SITE %s:%d. return code:%d while fetching PL tag\n", bcf_seqname(hdr,rec), rec->pos, npl);
        continue;
      }
      // https://github.com/samtools/bcftools/blob/e9c08eb38d1dcb2b2d95a8241933daa1dd3204e5/plugins/tag2tag.c#L151
      
      for (int i=0; i<npl; i++){
        if ( pl[i]==bcf_int32_missing ){
          bcf_float_set_missing(ln_gl[i]);
        } else if ( pl[i]==bcf_int32_vector_end ){
          bcf_float_set_vector_end(ln_gl[i]);
        } else{
          ln_gl[i] = pl2ln_f(pl[i]);
        }
	//         fprintf(stderr, "%d %f\n", pl[i], ln_gl[i]);
      }
    } else if(vcf_format_field == "GT"){
       int ngts = bcf_get_genotypes(hdr, rec, &gt, &ngt_arr);
       if ( ngts<0 ){
         fprintf(stderr, "BAD SITE %s:%d. return code:%d while fetching GT tag\n", bcf_seqname(hdr,rec), rec->pos, npl);
         continue;
       }
       for(int ns=0; ns<nsamples;ns++){
         int32_t *curr_ptr = gt + ns*2;
         float *ln_gl_ptr = ln_gl + ns*3;
         if ( bcf_gt_is_missing(curr_ptr[0]) ||
              bcf_gt_is_missing(curr_ptr[1]) ){ // obs genotype is missing
           // missing
           ln_gl_ptr[0] = 0;
           ln_gl_ptr[1] = 0;
           ln_gl_ptr[2] = 0;
           
         } else if (bcf_gt_allele(curr_ptr[0])!=bcf_gt_allele(curr_ptr[1])){ // this is obs genotype
           // het
	   ln_gl_ptr[0] = -INFINITY;
	   ln_gl_ptr[1] = 0;
	   ln_gl_ptr[2] = -INFINITY;
	 } else if(bcf_gt_allele(curr_ptr[0])==1){ // this is obs genotype
           // hom alt
	   ln_gl_ptr[0] = -INFINITY;
	   ln_gl_ptr[1] = -INFINITY;
	   ln_gl_ptr[2] = 0;
	 } else{ // this is obs genotype
           // hom ref
	   ln_gl_ptr[0] = 0;
	   ln_gl_ptr[1] = -INFINITY;
	   ln_gl_ptr[2] = -INFINITY;
	 }
       }
    } else {
      fprintf(stderr, "\t\t-> BIG TROUBLE. Can only take one of two tags, GT or PL\n");
    }
    
    int keepInd=0;
    char keep[nsamples];
    double *tmp = new double[3*nsamples];    
    for(int ns=0;ns<nsamples;ns++){
      float *ary= ln_gl+ns*3;
      if ((is_nan_vcf(ary[0]) || is_nan_vcf(ary[1]) || is_nan_vcf(ary[2])) ||(same_val_vcf(ary[0], ary[1]) && same_val_vcf(ary[0], ary[2]))){
        keep[ns]=0;
      }else{
	keep[ns]=1;
	keepInd++;
      }
      tmp[ns*3] = ary[0];
      tmp[ns*3+1] = ary[1];
      tmp[ns*3+2] = ary[2];
      // fprintf(stderr, "TMP: %d %d: %f %f %f\n", ns+1, rec->pos+1, tmp[ns*3], tmp[ns*3+1], tmp[ns*3+2]);
    }
    //    fprintf(stderr,"keepind:%d\n",keepInd);

    
    naf = bcf_get_info_float(hdr, rec, vcf_allele_field.c_str(), &af, &naf_arr);
    // fprintf(stderr,"rec->pos:%d npl:%d ngl:%d naf:%d rec->n_allele:%d\n",rec->pos,npl,ngl,naf,rec->n_allele);    
    //if multiple alt alleles then n_allele>3. We only care about diallelic ref/alt alleless
    //		if(rec->n_allele==4) fprintf(stdout,"\n%s\n",rec->d.allele[2]);
    //ok this is a bit messed up. apparantly sometime the allele is <*> sometimes not.
    // just use the first two alleles now and discard the rest of the alleles.

    double freq;
    if(naf==1){
      freq = af[0];
    }else{
      freq = emFrequency(tmp,nsamples,50,0.05,keep,keepInd);
    }
    //should matter, program should never run on such low freqs, this is just for validation between formats
    if(freq>0.999)
      freq=1;
    if(freq<0.001)
      freq=0;
    //fprintf(stderr,"freq:%f minfreq:%f keepInd:%d minind:%d\n",freq,minfreq,keepInd,minind);
    //filtering
    if(keepInd>minind&&freq>=minfreq && freq<= (1-minfreq)) {

#ifdef __WITH_MAIN__
      fprintf(stdout,"%s\t%i\t%s\t%s\tqual:%f n_info:%d n_allele:%d n_fmt:%d n_sample:%d n_samplws_with_data:%d freq:%f",
	      seqnames[rec->rid],
	      rec->pos+1,
	      rec->d.allele[0],
	      rec->d.allele[1],
	      rec->qual,
	      rec->n_info,
	      rec->n_allele,
	      rec->n_fmt,
	      rec->n_sample,
	      keepInd,
	      freq
	      );
      for(int i=0;i<3*nsamples;i++)
	fprintf(stdout," %f",ln_gl[i]);
      fprintf(stdout,"\n");
#endif
      mygl.push_back(tmp);
      freqs.push_back(freq);
      //populate debug names
#if 1
      char ptmp[1024];
      snprintf(ptmp,1024,"%s_%d", seqnames[rec->rid],rec->pos+1);
      posinfo.push_back(strdup(ptmp));
      //	fprintf(stderr,"pushhh\n");
#endif
    } else {
      delete [] tmp;
    }
    // fprintf(stderr,"rec->pos:%d npl:%d naf:%d rec->n_allele:%d af[0]:%f\n",rec->pos,npl,naf,rec->n_allele,freq);
    // exit(0);
  }
  fprintf(stderr, "\t-> [vcf.cpp] Read %i records %i of which were SNPs number of sites with data:%lu\n", n, nsnp,mygl.size());
  free(pl);
  free(gt);
  bcf_hdr_destroy(hdr);
  bcf_close(inf);
  bcf_destroy(rec);

  return nsamples;
 
}

#ifdef __WITH_MAIN__

int main(int argc, char **argv) {
  if (argc != 2) {
    usage();
    return 1;
  }
  std::vector<double *> gls;
  std::vector<double> freqs;
  std::vector<char *> posinfo;
  std::string pl=std::string("PL");
  std::string fr=std::string("AFngsrelate");
  int nsites = getgls(argv[1],gls,freqs,2,0.04,pl,fr,posinfo);
  return 0;
}

#endif
