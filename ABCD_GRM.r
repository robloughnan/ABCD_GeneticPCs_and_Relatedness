library(SNPRelate)
library(GWASTools)
library(GENESIS)
library(ggplot2)
library(qqman)

# This script runs population and relatedness inference using GENESIS (https://github.com/UW-GAC/GENESIS) 
# on ABCD data. The approach is considered best practice for performing GWAS
# in samples of mixed ancestry with a high degree of relatives and is used by both the TOPMED consortium and PAGE Study.
# The code essentially lifts from the two tutorials in the package documentation 
# (http://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html and
# http://bioconductor.org/packages/devel/bioc/vignettes/GENESIS/inst/doc/assoc_test.html)
# The main contributions this framework has over other approaches are:
# 1) PC-AIR: computing genetic PCs in sample of individuals with relatedness and ancestry admixture.
# 2) PC-Relate: compute genetic relateness in sample of mixed ancestries/population structure.

# Genetics (path to PLINK files)
fn_prefix='/path/to/non_imputed/ABCD_20220428.updated.nodups.curated'

###################### Create GDS file if one does not exist #################
gds.fn=paste0(fn_prefix, '.gds')
# If gds file does not exist create it
if (!file.exists(gds.fn)){
    print('Generating and saving GDS')
    snpgdsBED2GDS(bed.fn = paste0(fn_prefix, '.bed'),
    bim.fn = paste0(fn_prefix, '.bim'),
    fam.fn = paste0(fn_prefix, '.fam'),
    out.gdsfn = gds.fn)
    print('Completed GDS generation')
}

###################### Prunned SNPs and Kinship Matrix #################
pruned_file = paste0(fn_prefix, '_prunned.snps')
kinship_file = paste0(fn_prefix, '_kinship.RData')

gds <- snpgdsOpen(gds.fn)
## Prunned set of SNPs
if (!file.exists(pruned_file)){
    print('Generating prunned set list')
    snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, 
                      ld.threshold=sqrt(0.1), verbose=FALSE)
    pruned <- unlist(snpset, use.names=FALSE)
    print(paste0(length(pruned), ' pruned snps identified'))
    write(pruned, pruned_file)
    print(paste0('Saved to ', pruned_file))
}else{
    pruned = read.csv(pruned_file, header=F)$V1
}
## Kinship matrix
if (!file.exists(kinship_file)){
    print('Generating kinship matrix')
    ibd.robust = snpgdsIBDKING(gds)
    save(ibd.robust, file=kinship_file)
    print(paste0('Saved kinship matrix to ', kinship_file))
}else{
    # Load File
    load(kinship_file)
}
KINGmat = as(ibd.robust$kinship, 'dsyMatrix')
row.names(KINGmat) = ibd.robust$sample.id
colnames(KINGmat) = ibd.robust$sample.id
snpgdsClose(gds)

###################### PCA-AiR #################
# read in GDS data
ABCD_geno <- GdsGenotypeReader(filename = gds.fn)
# create a GenotypeData class object
ABCD_genoData <- GenotypeData(ABCD_geno)
ABCD_genoData

print('Generating pcair')
mypcair <- pcair(ABCD_genoData, kinobj = KINGmat, divobj = KINGmat, 
                 snp.include = pruned)

pcair_file = paste0(fn_prefix, '_pcair.RData')
save(mypcair, file=pcair_file)
print(paste0('Saved pcair ', pcair_file))

pcs = mypcair$vectors
colnames(pcs) = paste0('C', seq(dim(pcs)[2]))
pcair_file = paste0(fn_prefix, '_pcair.tsv')
write.table(pcs, pcair_file, quote=FALSE, sep='\t', col.names=NA)

# Output from PC-AiR
summary(mypcair) 

###################### PCA-Relate #################
pcrelate_file = paste0(fn_prefix, '_pcrelate.RData')
if (!file.exists(pcrelate_file)){
    print('Generating PC-Relate')
    ABCD_genoData <- GenotypeBlockIterator(ABCD_genoData, snpInclude=pruned)
    mypcrelate <- pcrelate(ABCD_genoData, pcs = mypcair$vectors[,1:2], 
                        training.set = mypcair$unrels)
    save(mypcrelate, file=pcrelate_file)
    print(paste0('Saved PC-Relate ', pcrelate_file))
}else{
    load(pcrelate_file)
}

###################### Create GRM #################
print('Generating GRM')
grm = pcrelateToMatrix(mypcrelate)
grm_mat = format(as.matrix(grm), digits=8)
grm_out = paste0(fn_prefix, '_GRM.tsv')
write.table(grm_mat, grm_out, quote=FALSE, sep='\t', col.names=NA)
print(paste0('Saved GRM to ', grm_out))

################ Compute PI_hat and Save GRM Long Form ####
mypcrelate$kinBtwn$k1 = 1 - mypcrelate$kinBtwn$k2 - mypcrelate$kinBtwn$k0
# Calculate PI_HAT as https://www.cog-genomics.org/plink/1.9/ibd
mypcrelate$kinBtwn$PI_HAT = mypcrelate$kinBtwn$k2 + 0.5*mypcrelate$kinBtwn$k1
# Clip values between 0 and 1
mypcrelate$kinBtwn$PI_HAT = pmin(pmax(mypcrelate$kinBtwn$PI_HAT, 0), 1)
mypcrelate$kinBtwn$GeneticRel = 2*mypcrelate$kinBtwn$kin
grm_long_out = paste0(fn_prefix, '_PCRelate_long.tsv')
cols_to_save = c('ID1','ID2','k0','k1','k2','nsnp','PI_HAT','GeneticRel')
write.table(mypcrelate$kinBtwn[, cols_to_save], grm_long_out, quote=FALSE, sep='\t', row.names=FALSE)
print(paste0('Saved long GRM to ', grm_long_out))
