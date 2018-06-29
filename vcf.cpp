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

double pl2ln[256];

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


int getgls(char*fname,std::vector<double *> &mygl, std::vector<double> &freqs,int minind,double minfreq){
  fprintf(stderr,"\t-> [getgls] fname:%s minind:%d minfreq:%f\n",fname,minind,minfreq);

  for(int i=0;i<256;i++){
    pl2ln[i] = log(pow(10.0,-0.1*i));
    //    fprintf(stderr,"%d) %f %f\n",i,exp(pl2ln[i]),pl2ln[i]);
  }
 
  // counters
  int n    = 0;  // total number of records in file
  int nsnp = 0;  // number of SNP records in file
  int nseq = 0;  // number of sequences
  int nsamples = 0;

  // pl data for each call
  int npl_arr = 0;
  int npl     = 0;
  int *pl     = NULL;
  
  // gl data for each call
  int ngl_arr = 0;
  int ngl     = 0;
  float *gl     = NULL;

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
#if __WITH_MAIN__
   seqnames = bcf_hdr_seqnames(hdr, &nseq);
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
    
#if 0
    ngl = bcf_get_format_float(hdr, rec, "GL", &gl, &ngl_arr);
#endif
    npl = bcf_get_format_int32(hdr, rec, "PL", &pl, &npl_arr);
    
    naf = bcf_get_info_float(hdr, rec, "AF1", &af, &naf_arr);
    //    fprintf(stderr,"rec->pos:%d npl:%d ngl:%d naf:%d rec->n_allele:%d af[0]:%f\n",rec->pos,npl,ngl,naf,rec->n_allele,af[0]);

    //if multiple alt alleles then n_allele>3. We only care about diallelic ref/alt alleless
    //		if(rec->n_allele==4) fprintf(stdout,"\n%s\n",rec->d.allele[2]);
    //ok this is a bit messed up. apparantly sometime the allele is <*> sometimes not.
    // just use the first two alleles now and discard the rest of the alleles.
    if(rec->n_allele>3||rec->n_allele==1)//last case shouldt happen
      continue;

    int offs = rec->n_allele==2?3:6;
      
    //lets only deal with pl for now.
    double *tmp = new double[3*nsamples];
    int keepInd=0;
    char keep[nsamples];
    //  memset(keep,'0',nsamples);
    for(int n=0;n<nsamples;n++){
      for(int nn=0;nn<3;nn++){
	//fprintf(stderr,"%d\n",pl[n*offs+nn]);
	tmp[n*3+nn] = pl2ln[pl[n*offs+nn]];
	// 
      }
      double *ary= tmp+n*3;
      if(ary[0]==ary[1]&&ary[0]==ary[2])
	keep[n]=0;
      else{
	keep[n]=1;
	keepInd++;
	}
      
    }
    double freq;
    if(naf==1)
      freq = af[0];
    else
      freq= emFrequency(tmp,nsamples,50,0.05,keep,keepInd);
    // fprintf(stdout,"%f %f\n",af[0],freq);;
    //    fprintf(stdout,"%d %f\n",keepInd,freq);
    //exit(0);
    //filtering
    if(keepInd>minind&&freq>=minfreq) {
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
	fprintf(stdout," %f",tmp[i]);
      fprintf(stdout,"\n");
#endif
      mygl.push_back(tmp);
      freqs.push_back(freq);
      //	fprintf(stderr,"pushhh\n");
    }else
      delete [] tmp;
  }
  fprintf(stderr, "\t-> [vcf.cpp] Read %i records %i of which were SNPs number of sites with data:%lu\n", n, nsnp,mygl.size());
   free(pl);
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
  int nsites = getgls(argv[1],gls,freqs,2,0.04);
  return 0;
}

#endif
