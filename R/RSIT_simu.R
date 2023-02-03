#' Simulation for RSIT test
#'
#' This function use pre-defined parameters to make the simulation data for the RSIT test (including type I and power test)
#'
#' @param n.sample Numeric, sample size, number of individuals
#' @param geno matrix, genotype matrix

#' @param n_casual Vector, number of causal SNPs for common SNPs and rare SNPs
#' @param main_effect Vector, main effect for common SNPs and rare SNPs
#' @param main_effect Vector, gene-environment interaction effect for common SNPs and rare SNPs
#' @param direction Vector, main effect and gene-environment interaction effect direction
#' @param GEindependent Logical variable, indicating whether gene and environment variable are independent
#' @param power Logical variable, indicating whether include disease genes in the generated genes
#' @param maineffect Logical variable, indicating whether include main effect
#' 
#' 
#' 
#' @return A list object is returned to be used as object for RSIT test
#' @export

simu.binary.phe<-function(n.sample, geno, n_casual=c(2,8),  main_effect=c(0.1,0.2), GE_effect=c(0.3,0.88), 
direction=c("same","same"), GEindependent=T, power=T, maineffect=T,SNPtype="rare_common",main_GE=T){
    geno=as.matrix(geno)
    ss.raw=dim(geno)[1]
    maf=colMeans(geno)/2
    if(SNPtype=="rare"){
    geno=geno[,which(maf<0.05)]
    } else if(SNPtype=="common"){
    geno=geno[,which(maf>=0.05)]
    } else{
    geno=geno}
    maf=colMeans(geno)/2
    x1.raw = rnorm(ss.raw, mean=62.4, sd=11.5)											
    x2.raw = rbinom(ss.raw, prob=0.5, size=1)	
    
    if(min(maf)<0.05){
        G.rare.main=geno[,sample(which(maf<0.05),n_casual[2],replace=F)]
        G.rare.GE=geno[,sample(which(maf<0.05),n_casual[2],replace=F)]
    }else {
        G.rare.main=NULL
        G.rare.GE=NULL
    }

    if(max(maf)>=0.05){
        G.common.main=geno[,sample(which(maf>0.05),n_casual[1],replace=F)]
        G.common.GE=geno[,sample(which(maf>0.05),n_casual[1],replace=F)]
    }else {
        G.common.main=NULL
        G.common.GE=NULL
    }
    SNP=cbind(G.common.main,G.rare.main)
    SNP.GE=cbind(G.common.GE,G.rare.GE)
    if(GEindependent){													
        e1.raw = rbinom(ss.raw, prob=0.5, size=1) 
    }else{
        logitE=-(rowSums(SNP))*0.5
        e1.raw = rbinom(ss.raw, prob=1/(1+exp(-logitE)), size=1) 
    }
    alpha1 = runif(n_casual[1], min=main_effect[1]-0.02, max=main_effect[1]+0.02)
    alpha2 = runif(n_casual[2], min=main_effect[2]-0.02, max=main_effect[2]+0.02)
    beta1 = runif(n_casual[1], min=GE_effect[1]-0.02, max=GE_effect[1]+0.02)
    beta2 = runif(n_casual[2], min=GE_effect[2]-0.02, max=GE_effect[2]+0.02)
    y_common = y_rare = y_common.main = y_rare.main = rep(0,ss.raw)
    if(direction[1]=="same"){
        main_sign=rep(1,sum(n_casual))
    }else if(direction[1]=="oppo"){
        main_sign=rep(c(1,-1),sum(n_casual))
    } else{
        stop("Error: direction shoud be \"same\" or \"oppo\" ")
    }

    if(direction[2]=="same"){
        GE_sign=rep(1,sum(n_casual))
    }else if(direction[2]=="oppo"){
        GE_sign=rep(c(1,-1),sum(n_casual))
    } else{
        stop("Error: direction shoud be \"same\" or \"oppo\" ")
    }
    if(max(maf)>=0.05){
        for (i in 1:n_casual[1]){
            if(main_GE){
                y_common = y_common + beta1[i]*G.common.main[,i]*e1.raw*GE_sign[i]
            } else{y_common = y_common + beta1[i]*G.common.GE[,i]*e1.raw*GE_sign[i]}
            y_common.main = y_common.main + alpha1[i]*G.common.main[,i]*main_sign[i]
        }
    }
    if(min(maf)<0.05){
        for (j in 1:n_casual[2]){
            if(main_GE){
            y_rare = y_rare + beta2[j]*G.rare.main[,j]*e1.raw*GE_sign[j]
            } else{y_rare = y_rare + beta2[j]*G.rare.GE[,j]*e1.raw*GE_sign[j]}
            y_rare.main = y_rare.main + alpha2[j]*G.rare.main[,j]*main_sign[j+n_casual[1]]
        }
    }

    if(power){
        if(maineffect){
            linear = -6.2 + 0.05*x1.raw+0.5*x2.raw+0.64*e1.raw + y_rare +  y_rare.main + y_common + y_common.main
                        variance_prop=c(var(linear),var(y_rare.main + y_common.main),var(y_rare + y_common),dim(geno)[2])
        }else{
            linear = -6.2 + 0.05*x1.raw+0.5*x2.raw+0.64*e1.raw + y_rare +  y_common
                        variance_prop=c(var(linear),0,var(y_rare + y_common),dim(geno)[2])
        }
    }else{
        if(maineffect){
            linear = -6.2 + 0.05*x1.raw+0.5*x2.raw+0.64*e1.raw + y_common.main+  y_rare.main
                        variance_prop=c(var(linear),var(y_rare.main + y_common.main),dim(geno)[2])
        }else{
            linear = -6.2 + 0.05*x1.raw+0.5*x2.raw+0.64*e1.raw
                        variance_prop=c(var(linear),0,dim(geno)[2])
        }
    }

    pi = exp(linear)/(1+exp(linear))
    y.raw = rbinom(ss.raw, prob=pi, size=1)
    pre = mean(pi)
    cat(pre,"\n")
    y.case.index = which(y.raw == 1)
    y.control.index = which(y.raw == 0)
    subcase=sample(y.case.index,n.sample/2)
    subcontrol=sample(y.control.index,n.sample/2)

    y = y.raw[c(subcase,subcontrol)]
    x1 = x1.raw[c(subcase,subcontrol)]
    x2 = x2.raw[c(subcase,subcontrol)]
    e1 = e1.raw[c(subcase,subcontrol)]
    G =  geno[c(subcase,subcontrol),]

    return(list(y=y, x1=x1,x2=x2,e1=e1,G=G,pre=pre,variance_prop=variance_prop));
}

#' Simulate the genotype data
#'
#' This function use pre-defined parameters to simulate the genotype data
#'
#' @param n.sample Numeric, sample size, number of individuals
#' @param mat.count Numeric, number of SNP regions
#' @param common_snp Logical variable, indicating whether  the region include common SNPs
#' @param independent Logical variable, indicating whether SNPs in the same region are independent
#
#' 
#' 
#' 
#' @return A list object is returned to be used as genotype for RSIT test
#' @export

simu_snp_mat<-function(  n.sample, mat.count,common_snp=F,independent=F)
{
  file.snp.hap1 <- system.file("extdata", "skat-test-1.hap.gz", package="RSIT");
  file.snp.pos2 <- system.file("extdata", "skat-test-1.pos", package="RSIT");
  
  snp.hap <- utils::read.table(file.snp.hap1, header=F);
  snp.pos <- utils::read.table(file.snp.pos2, header=T);
  
  snp.maxpos <- max(snp.pos$CHROM_POS);
  snp.minpos <- min(snp.pos$CHROM_POS);
  
  rare.cutoff <- 0.05;
  genelength <- 20*1000
  
  mat.ret <- list();
  n.mat <- 1;
  avoid_range <- c();
  while(n.mat<=mat.count)
  {
    snp.start <- as.integer(runif(1, min=snp.minpos, max=snp.maxpos- genelength));
    snp.ends <- snp.start +  genelength;
    if(!independent){
      p.sel<-which( snp.pos$CHROM_POS>snp.start & snp.pos$CHROM_POS<snp.ends );
    }else{
      p.sel=sample(1:(dim(snp.hap)[2]-2),300,replace=F)
    }
    
    snp.mat1 <- snp.hap[sample(nrow(snp.hap),n.sample,replace=T), p.sel + 2 ];
    snp.mat2 <- snp.hap[sample(nrow(snp.hap),n.sample,replace=T), p.sel + 2 ];
    snp.mat <- snp.mat1 + snp.mat2 -2;
    
    maf <- colMeans(snp.mat)/2;
    m.same <- which( maf==1 | maf==0 );
    if (length(m.same)>0)
      snp.mat <- snp.mat[, -m.same, drop=F ];
    
    #make sure of 2: Minor Allels, 0:Major Allels
    maf <- colMeans(snp.mat)/2;
    m.minor <- which( maf>0.5 );
    if (length(m.minor)>0)
      for(i in 1:length(m.minor))
        snp.mat[,m.minor[i]] <- 2 - snp.mat[,m.minor[i]];
    
    maf <- colMeans(snp.mat)/2;
    m.smallmaf <- which( maf < 0.001);
    if (length(m.smallmaf)>0)
      snp.mat <- snp.mat[, -m.smallmaf,drop=F ];
    
    if (dim(snp.mat)[2]==1)
      next;
    
    maf <- colMeans(snp.mat)/2;
    n.rare <- length(which(maf<rare.cutoff));
    if (n.rare<1)
      next;
    
    m.commonmaf <- which( maf > rare.cutoff)
    if( length(m.commonmaf)>0&(!common_snp))
      snp.mat <- snp.mat[, -m.commonmaf,drop=F ];
    
    if(independent){ 
      if(common_snp){
        id=c(sample(which(maf>rare.cutoff),10,replace=F),sample(which(maf<rare.cutoff),40,replace=F))   
      }else{
        maf <- colMeans(snp.mat)/2
        id=sample(which(maf<rare.cutoff),40,replace=F)
      }
      snp.mat=snp.mat[,id]
    }
    rownames(snp.mat) <- paste("ID", 1:NROW(snp.mat), sep="");
    mat.ret[[n.mat]] <- snp.mat;
    
    n.mat <- n.mat + 1;
  }
  
  return(mat.ret);
}

