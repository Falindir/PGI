#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXLLINE 100000
#define MAXLWORD 1000
#define MAXMOTIFLG 10
#define MAXNBALLELE 10

#define true 1
#define false 0

/* data structure  */


struct microsat{
  int ctg_index;
  int start, end;
  int replen, avlen;
  char motif[MAXMOTIFLG];
  int motiflen;
  int nb_alleles;
  int nb_genotyped_ind;
  double heterozygosity;
  int** genotype;	    //ind,gene_copy
  int** coverage;	    //ind,gene_copy
};

struct contig{
  char* name;
  int index;
  int nb_msat;
  int* msat_id;
};

struct dataset{
  int nb_indiv;
  char** indiv_names;
  int nb_contigs;
  struct contig* ctg;
  int nb_msat; 
  struct microsat* msat;
};




/* check_pt */
void check_pt(void* pt, char mess[]){

  if(!pt) {printf("problem: %s\n", mess); exit(0);}

}



/* usage */
void usage(){

  printf("\n");
  printf("usage: amp infile outfile min_cov_tot min_cov_var max_cov_tot\n");
  printf("\n");
  exit(0);
}


/* read_data */

struct dataset* read_data(FILE* infile){

  struct dataset* data;
  char line[MAXLLINE+1], word[MAXLWORD+1], *tok, *ret, *ctg_name;
  int i, j, cur_ctg, cur_msat, new_ctg, g0, g1, c0, c1, prev_msat;

  data=(struct dataset*)calloc(1, sizeof(struct dataset));

  //parse header, read indiv nb and names
  ret=fgets(line, MAXLLINE, infile);
  check_pt(ret, "cannot read from input file");
  strtok(line, "\t \n"); strtok(NULL, "\t \n"); strtok(NULL, "\t \n");
  data->nb_indiv=0;
  while(1){
    tok=strtok(NULL, "\t \n");
    if(!tok) break;
    data->nb_indiv++;
  }
  if(data->nb_indiv==0) check_pt(NULL, "pb header line");
  data->indiv_names=(char**)calloc(data->nb_indiv, sizeof(char*));
  check_pt(data->indiv_names, "not enough memory");
  for(i=0;i<data->nb_indiv;i++){
    data->indiv_names[i]=(char*)calloc(MAXLWORD+1, sizeof(char));
    check_pt(data->indiv_names[i], "not enough memory");
  }
  rewind(infile);
  fgets(line, MAXLLINE, infile);
  strtok(line, "\t \n"); strtok(NULL, "\t \n"); strtok(NULL, "\t \n");
  for(i=0;i<data->nb_indiv;i++)
    sprintf(data->indiv_names[i], "%s", strtok(NULL, "\t \n"));  

  //first pass: count microsats
  rewind(infile);
  fgets(line, MAXLLINE+1, infile);
  data->nb_msat=0;
  while(1){
    ret=fgets(line, MAXLLINE+1, infile);
    if(!ret) break;
    if(line[0]=='\n') continue;
    data->nb_msat++;
  }

  //allocate
  data->msat=(struct microsat*)calloc(data->nb_msat, sizeof(struct microsat));
  check_pt(data->msat, "not enough memory");
  data->ctg=(struct contig*)calloc(data->nb_msat, sizeof(struct contig));
  check_pt(data->ctg, "not enough memory");

  //second pass: read all
  rewind(infile);
  fgets(line, MAXLLINE+1, infile);
  cur_ctg=-1;
  cur_msat=-1;
  while(1){
    //read line
    ret=fgets(line, MAXLLINE+1, infile);
    if(!ret) break;
    if(line[0]=='\n') continue;
    cur_msat++;
    //word 1 = contig name
    ctg_name=strtok(line, "_ \t");
    new_ctg=0;
    if(cur_ctg==-1) new_ctg=1;
    else if(strcmp(ctg_name, data->ctg[cur_ctg].name)!=0) new_ctg=1;
    if(new_ctg){
      cur_ctg++;
      data->ctg[cur_ctg].name=(char*)calloc(MAXLWORD+1, sizeof(char));
      sprintf(data->ctg[cur_ctg].name, "%s", ctg_name);
      data->ctg[cur_ctg].index=cur_ctg;
      if(cur_ctg==1) prev_msat=0; 
      if(cur_ctg>1) prev_msat=data->ctg[cur_ctg-2].nb_msat;
      if(cur_ctg>=1){
        data->ctg[cur_ctg-1].nb_msat=cur_msat-prev_msat;
        data->ctg[cur_ctg-1].msat_id=(int*)calloc(data->ctg[cur_ctg-1].nb_msat, sizeof(int));
        check_pt(data->ctg[cur_ctg-1].msat_id, "not enough memory");
        for(j=0;j<data->ctg[cur_ctg-1].nb_msat;j++)
          data->ctg[cur_ctg-1].msat_id[j]=j+prev_msat;
      }
    } 
    
    data->msat[cur_msat].ctg_index=cur_ctg;
    //word 2 and 3 = boundaries
    sscanf(strtok(NULL, "_ \t"), "%d", &(data->msat[cur_msat].start));
    sscanf(strtok(NULL, "_ \t"), "%d", &(data->msat[cur_msat].end));
    //word 4 = replen
    sscanf(strtok(NULL, "_ \t"), "%d", &(data->msat[cur_msat].replen));
    //word 5 = motif
    sprintf(data->msat[cur_msat].motif, "%s", strtok(NULL, "_ \t"));
    data->msat[cur_msat].motiflen=(int)strlen(data->msat[cur_msat].motif);
    //genotypes
    data->msat[cur_msat].genotype=(int**)calloc(data->nb_indiv, sizeof(int*));
    check_pt(data->msat[cur_msat].genotype, "not enough memory");
    data->msat[cur_msat].coverage=(int**)calloc(data->nb_indiv, sizeof(int*));
    check_pt(data->msat[cur_msat].coverage, "not enough memory");
    for(j=0;j<data->nb_indiv;j++){
      sprintf(word, "%s", strtok(NULL, "_ \t"));
      data->msat[cur_msat].genotype[j]=(int*)calloc(2, sizeof(int));
      check_pt(data->msat[cur_msat].genotype[j], "not enough memory");
      data->msat[cur_msat].coverage[j]=(int*)calloc(2, sizeof(int));
      check_pt(data->msat[cur_msat].coverage[j], "not enough memory");
      if(word[0]=='.') {g0=g1=c0=c1=-1;}
      else if(word[4]=='|') {
        sscanf(word, "%d", &g0); sscanf(word+2, "%d", &c0);
        sscanf(word+4, "%d", &g1); sscanf(word+6, "%d", &c1);
      }
      else {sscanf(word, "%d", &g0); g1=g0; sscanf(word+2, "%d", &c0); c1=c0/2; c0=c0-c1;}
	// Bogi: outputs allele 1, and allele2 as well as coverage for both alleles for the current locus
      data->msat[cur_msat].genotype[j][0]=g0;
      data->msat[cur_msat].genotype[j][1]=g1;
      data->msat[cur_msat].coverage[j][0]=c0;
      data->msat[cur_msat].coverage[j][1]=c1;
    }
  }
  data->nb_contigs=cur_ctg+1;

  //add microsat list to last contig (not done in main loop above)
  if(cur_ctg==1) prev_msat=0; 
  if(cur_ctg>1) prev_msat=cur_msat-data->ctg[cur_ctg-2].nb_msat;
  if(cur_ctg>=1){
    data->ctg[cur_ctg-1].nb_msat=cur_msat-prev_msat;
    data->ctg[cur_ctg-1].msat_id=(int*)calloc(data->ctg[cur_ctg-1].nb_msat, sizeof(int));
    check_pt(data->ctg[cur_ctg-1].msat_id, "not enough memory");
    for(j=0;j<data->ctg[cur_ctg-1].nb_msat;j++)
      data->ctg[cur_ctg-1].msat_id[j]=j+prev_msat;
  }

  return data;
}



/* which_in_list */

which_in_list(int k, int *list, int n){

  int i;

  for(i=0;i<n;i++)
    if(k==list[i]) return i;

  return -1;
}



/* calculate_msat_attributes */

void calculate_msat_attributes(struct dataset* data, int min_cov_ind, int min_cov_var, int max_cov_ind){

  int i, j, k, ge, al, br;
  int* allele_list, *allele_freq, *genotyped;
  double nbgen;
  
  allele_list=(int*)calloc(MAXNBALLELE, sizeof(int));
  check_pt(allele_list, "not enough memory");
  allele_freq=(int*)calloc(MAXNBALLELE, sizeof(int));
  check_pt(allele_freq, "not enough memory");
  genotyped=(int*)calloc(data->nb_indiv, sizeof(int));
  check_pt(genotyped, "not enough memory");

  for(i=0;i<data->nb_msat;i++){

    br=0;
    for(j=0;j<data->nb_indiv;j++) genotyped[j]=0;

    //nb_genotyped_ind
    data->msat[i].nb_genotyped_ind=0;
    for(j=0;j<data->nb_indiv;j++){
      if(data->msat[i].genotype[j][0]==-1) continue;
	// Bogi: Checks if coverage for a locus is bad
	// 1st line: cov(allele1)+cov(allele2)<min_lim?
	// 2nd line: cov(allele1)<var_lim?
	// 3rd line: cov(allele2)<var_lim?
	// 4th line: cov(allele1)+cov(allele2)>max_lim?
	// Why does he call it variance if he doesn't calculate any variance?
	// Continue here prevents loop from going through other statements if current one isn't satisfied
	// Returns to the begining of the loop
      if(data->msat[i].coverage[j][0]+data->msat[i].coverage[j][1]<min_cov_ind) continue;
      if(data->msat[i].coverage[j][0]<min_cov_var) continue; 
      if(data->msat[i].coverage[j][1]<min_cov_var) continue;
      if(data->msat[i].coverage[j][0]+data->msat[i].coverage[j][1]>max_cov_ind) continue;
      data->msat[i].nb_genotyped_ind++;
      genotyped[j]=1;
    }

    //nb_alleles and allele frequencies
    for(j=0;j<MAXNBALLELE;j++) allele_list[j]=0;
    for(j=0;j<MAXNBALLELE;j++) allele_freq[j]=0;
    data->msat[i].nb_alleles=0;
    for(j=0;j<data->nb_indiv;j++){
      for(k=0;k<2;k++){
        ge=data->msat[i].genotype[j][k];
        if(genotyped[j]){
          al=which_in_list(ge, allele_list, data->msat[i].nb_alleles);
          if(al>=0) allele_freq[al]++;
          else{
            allele_list[data->msat[i].nb_alleles]=ge;
            allele_freq[data->msat[i].nb_alleles]=1;
            data->msat[i].nb_alleles++;
            if(data->msat[i].nb_alleles==MAXNBALLELE){
              al=MAXNBALLELE;
              printf("more than %d alleles in microsat %d (skipped)\n", al, i+1);
              br=1; break;
            }
          }
        }
      }
      if(br) break;
    }

    //average repeat length
    

    //heterozygosity
    if(br) {data->msat[i].heterozygosity=-1.; getchar(); continue;}
    nbgen=2.*data->msat[i].nb_genotyped_ind;
    data->msat[i].heterozygosity=1.;
    for(j=0;j<data->msat[i].nb_alleles;j++)
	// Bogi: why is he multiplying allele_freq[j] with itself?
	// Probably subtracts homozygosity
	// Starts at H=1 and then decreases H with ???
	// nbgen is number of individuals with a given genotype
      data->msat[i].heterozygosity-=(allele_freq[j]/nbgen)*(allele_freq[j]/nbgen);
  }

}



/* data_summary */

void data_summary(struct dataset* data){

  int i;

  printf("%d individuals found\n", data->nb_indiv);
  if(data->nb_indiv<20){
    for(i=0;i<data->nb_indiv;i++)
      printf("%s  ", data->indiv_names[i]);
    printf("\n");
  }

  printf("%d microsatellites foud in %d distinct contigs\n", data->nb_msat, data->nb_contigs);

}



/* data_summary2 */

void data_summary2(struct dataset* data){

  int i, j, *nbgen;

  nbgen=calloc(data->nb_indiv, sizeof(int));

  for(i=0;i<data->nb_msat;i++){
    for(j=1;j<data->nb_indiv-1;j++){
      if(data->msat[i].nb_genotyped_ind>j) nbgen[j]++;
    }
  }

  for(j=1;j<data->nb_indiv-1;j++)
    printf("%d microsat loci are genotyped in >%d individuals\n", nbgen[j], j);

}



/* write_results */

void write_results(struct dataset* data, FILE* outfile){
// Why is he defining ctg_ind when it seems that it doesn't change between iterations?
  int i, ctg_ind;
  struct microsat* msat;

  fprintf(outfile, "msat,contig,start,end,motif,motiflen,reflen,nb_genotyped,nb_allele,heterozygosity\n");

  for(i=0;i<data->nb_msat;i++){
    msat=data->msat+i;
    ctg_ind = msat->ctg_index;
    fprintf(outfile, "%d,", i+1);
    fprintf(outfile, "%s,", data->ctg[ctg_ind].name);
    fprintf(outfile, "%d,", msat->start);
    fprintf(outfile, "%d,", msat->end);
    fprintf(outfile, "%s,", msat->motif);
    fprintf(outfile, "%d,", msat->motiflen);
    fprintf(outfile, "%d,", msat->replen);
    fprintf(outfile, "%d,", msat->nb_genotyped_ind);
    fprintf(outfile, "%d,", msat->nb_alleles);
    fprintf(outfile, "%f\n", msat->heterozygosity);
  }

}


/* main */

main(int argc, char** argv){

  FILE *infile, *outfile;
  struct dataset* data;
  int ret, i, min_cov_ind, min_cov_var, max_cov_ind;

  if(argc!=6) usage();


  //read from file
  infile=fopen(argv[1], "r");
  check_pt(infile, "cannot read file");
  printf("reading data from file: %s...\n", argv[1]);
  data=read_data(infile);
  data_summary(data);

  //read arguments
  ret=sscanf(argv[3], "%d", &min_cov_ind);
  if(ret!=1) usage();
  ret=sscanf(argv[4], "%d", &min_cov_var);
  if(ret!=1) usage();
  ret=sscanf(argv[5], "%d", &max_cov_ind);
  if(ret!=1) usage();

  //calculate
  printf("calculating...\n");
  calculate_msat_attributes(data, min_cov_ind, min_cov_var, max_cov_ind);
  data_summary2(data);

  //write
    outfile=fopen(argv[2], "w");
    printf("writing to file: %s\n", argv[2]);
    write_results(data, outfile);

  printf("done\n");

}



