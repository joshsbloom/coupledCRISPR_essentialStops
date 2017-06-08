# FUNCTIONS -----------------------------------------------------------------------------------------------------------

# some global variables
x <- org.Sc.sgdALIAS
mapped_probes <- mappedkeys(x)

gene2alias <- sapply(as.list(x[mapped_probes]), paste, collapse=' ')
t.var=c(0,24,48,72,96)


# setup for GO enrichment analysis ------------------------------------------------------------------------
myGene2GO.full.table=read.delim('/data/eQTL/reference/go_annotations_sgd.txt',  sep='\t', header=F, stringsAsFactors=F)

# data frame of gene then GO
myGene2GO=cbind( sapply(strsplit(myGene2GO.full.table[,11], '\\|'), function(x) x[1]), myGene2GO.full.table[,5])
myGene2GO=na.omit(myGene2GO)

SYS2ORF=unique(cbind(sapply(strsplit(myGene2GO.full.table[,11], '\\|'), function(x) x[1]), myGene2GO.full.table[,3]))
SYS2ORF.key=list()
SYS2ORF.key[SYS2ORF[,1]]=SYS2ORF[,2]
gene2GOList = lapply(unique(myGene2GO[,1]), function(x){myGene2GO[myGene2GO[,1] == x, 2]})
names(gene2GOList) = unique(myGene2GO[,1])

GO2geneList = lapply(unique(myGene2GO[,2]), function(x){myGene2GO[myGene2GO[,2] == x, 1]})
names(GO2geneList) = unique(myGene2GO[,2])


doBIC=function(l, k, n) { -2*l+k*log(n) }

countIndels=function(seqin, mat) { 
    am=apply(mat,2,function(x) x=='-')
    if(is.null(nrow(am))) { counts=sum(am) } else {  counts=rowSums(am) }
    if(is.null(nrow(am))) { pos=which(am) } else {    pos=apply(am,1, function(x) which(x)) }
    return(list(counts=counts, pos=pos))
}
countMM=function(seqin, mat, pos.offset=0) { 
    am=( t(t(mat)==seqin) + 0 + (((mat)!='-') + 0)==1) 
    counts=rowSums(am)
    pos=apply(am,1, function(x) which(x))
    #, collapse=':'))
    return(list(counts=counts, pos=pos))
}

# gos is table split by oligo
fit.oligo.model =function(gos, t.var, tot.counts,cl.cnt=2) {
    registerDoParallel(cl.cnt)
    oligo.fits=foreach( gn=names(gos)  ) %dopar% {
        print(gn)
        ag1=gos[[gn]] #gt.WT.pe[g1,]
        Xg1=data.matrix(ag1[,9:13])
        Xg1[Xg1[,2]==0,2]=1
        b=as.factor(ag1$barcode)
        y1=as.vector(Xg1)
        t1=rep(t.var, each=nrow(Xg1))
        b1= rep(b, ncol(Xg1))
        ofl=rep(tot.counts, each=nrow(Xg1))
        expt=rep(as.factor(ag1$expt), ncol(Xg1))

            if(nrow(Xg1)>1 & length(unique(b))>1 ){
                # random slope and intercept model, test if fixed effect of slope is signficant
                tryCatch({
                         g1=glmer(y1~offset(log(ofl))+t1+(1+t1|b1)+(1|b1), 
                                  family=poisson(link=log),control=glmerControl(optimizer="Nelder_Mead"))
                            }, error=function(e) {print(gn)}


                         )
                g0=glmer(y1~offset(log(ofl))+(1+t1|b1)+(1|b1) , family=poisson(link=log))
                mtest=anova(g0,g1,test='Chisq')
                pval=mtest[,ncol(mtest)][2]
                reffs=ranef(g0)[[1]]
                feffs=fixef(g1)
            } else {   
                g1=glm(y1~offset(log(ofl))+t1 , family=poisson(link=log))
                g0=glm(y1~offset(log(ofl)) , family=poisson(link=log))
                reffs=NA
                feffs=coef(g1)
                pval=NA
            }  
            #cInt=confint(g1, parm='t1', .95, quiet=T)
            cInt=suppressMessages( confint(g1, parm='t1', .95, method='Wald', quiet=T))
            # we want random effects from gN1
            # return(list(g0=g0, g1=g1, p.value=pval, cInt=cInt))
            gc()
            return(list(bc.cnt=nrow(Xg1), reffs=reffs, feffs=feffs, p.value=pval, cInt=cInt))
    }  
   names(oligo.fits)=names(gos) #[1:100]
   return(oligo.fits) 
}


# Visualize poisson regressions 
#gos=gBs
#tot.counts=colSums(gBsA[,9:13])
#gn=names(gos)[27]
#plot(t.var, Xg1[1,], type='p', ylim=c(0,300), ylab='counts', xlab='time')
#points(t.var, predict(glm(Xg1[1,]~offset(log(tot.counts))+t.var, family=poisson(link=log)), type='response'), col=1, type='l')
#for(n in 4:nrow(Xg1)){
#points(t.var, Xg1[n,], type='p', col=n)
#points(t.var, predict(glm(Xg1[n,]~offset(log(tot.counts))+t.var, family=poisson(link=log)), type='response'), col=n, type='l')
#readline()
#}



fit.oligo.model.combined =function(gos, t.var, tot.counts,cl.cnt=2) {
    registerDoParallel(cl.cnt)
  
    #removed reffs because of errors
    oligo.fits=foreach(gn=names(gos)) %dopar% {
        print(gn)
        #gn='chrI:101014:101016:4320:GCATACTTGGAGTTGGCAAG:+:0:-:TGA:130:132:3:44:384:90:92:YAL025C:nonsense:CCA:TGA:P:*:307:263'

        ag1=gos[[gn]] #gt.WT.pe[g1,]
        Xg1=data.matrix(ag1[,9:13])
        Xg1[Xg1[,2]==0,2]=1
        b   =  as.factor(paste0(ag1$barcode, ag1$expt)) #as.factor(ag1$barcode)
        y1  =  as.vector(Xg1)
        t1  =  rep(t.var, each=nrow(Xg1))
        b1  =  rep(b, ncol(Xg1))
        ofl =  rep(tot.counts, each=nrow(Xg1))
        expt = rep(as.factor(ag1$expt), ncol(Xg1))
        tryCatch({

        if(nrow(Xg1)>1 & length(unique(b))>1 ) {

               fixed.coeffs=t(apply(Xg1,1, function(y) coef(glm(y~offset(log(tot.counts))+t.var, family=poisson(link=log)))))
               
               g3= glmer(y1~offset(log(ofl))+t1+(1+t1|b1)+(1|b1), 
                                   family=poisson(link=log),control=glmerControl(optimizer="Nelder_Mead"))
               g0= glmer(y1~offset(log(ofl))+(1+t1|b1)+(1|b1) , family=poisson(link=log))

               # supress this model for now
               #gN= glmer(y1~offset(log(ofl))+(1|b1)+(0+t1|b1) , family=poisson(link=log))
                
               #cor.test(as.vector(Xg1), as.vector(exp(predict(gN))))
             
              if(length(unique(expt))>1) {
                    g2= glmer(y1~offset(log(ofl))+t1+expt+(1+t1|b1)+(1|b1), 
                                  family=poisson(link=log),control=glmerControl(optimizer="Nelder_Mead"))
                    g1= glmer(y1~offset(log(ofl))+expt+(1+t1|b1)+(1|b1), 
                                     family=poisson(link=log),control=glmerControl(optimizer="Nelder_Mead"))
                 mtest.time=anova(g2,g1,test='Chisq')
                 mtest.expt=anova(g2,g3,test='Chisq')
                
                 pval.time=mtest.time[,ncol(mtest.time)][2]
                 pval.expt=mtest.expt[,ncol(mtest.expt)][2]
                
                 #reffs=ranef(g1)[[1]]
                 feffs=fixef(g2)
                 
                 cInt=suppressMessages(confint(g2, parm='t1', .95, method='Wald', quiet=T))
                    
              } else {
                    mtest.time=anova(g3,g0,test='Chisq')
                    pval.time=mtest.time[,ncol(mtest.time)][2]
                    pval.expt=NA
                    #reffs=ranef(g0)[[1]]
                    feffs=fixef(g3)
                    cInt=suppressMessages(confint(g3, parm='t1', .95, method='Wald', quiet=T))
              }
         
       } 
       else {   
                g1=glm(y1~offset(log(ofl))+t1 , family=poisson(link=log))
                fixed.coeffs=coef(g1)
                g0=glm(y1~offset(log(ofl)) , family=poisson(link=log))
                #reffs=NA
                feffs=coef(g1)
                pval.time=NA
                pval.expt=NA
               cInt=suppressMessages(confint(g1, parm='t1', .95, method='Wald', quiet=T))
        }  

         }, error=function(e) {print('HELP ERROR')})

         gc()
         return(list(bc.cnt=nrow(Xg1), feffs=feffs, p.value.time=pval.time, p.value.expt=pval.expt, cInt=cInt, fixed.coeffs=fixed.coeffs))
    }  
   names(oligo.fits)=names(gos) #[1:100]
   return(oligo.fits) 
}
        

# filter on errors in repair template and  minimum number of read counts at t=0  (set to 20, was 35)
addFilters=function(giant.table, cutoff_t0=20) {
    # exact length
    perfect.filter=  giant.table$start.align==1 & giant.table$end.align==101 & giant.table$n.indel==0 & giant.table$n.mismatch==0 & (giant.table$WT_0> cutoff_t0 | giant.table$NMDm_0> cutoff_t0)

    clipping.filter= giant.table$start.align<5 & giant.table$end.align>96 & giant.table$n.indel==0 & giant.table$n.mismatch==0 &    (giant.table$WT_0> cutoff_t0 | giant.table$NMDm_0> cutoff_t0)

    d1.filter= giant.table$start.align<5 & giant.table$end.align>96 & giant.table$pam.indel.cnt==0 & giant.table$pam.mm.cnt==0  & giant.table$u1.indel.cnt==0 &
    giant.table$u1.mm.cnt==0 & giant.table$u.indel.cnt==0 & giant.table$u.mm.cnt==0 & giant.table$d.indel.cnt==0 & giant.table$d.mm.cnt==0 & giant.table$d1.indel.cnt<3 &
    giant.table$d1.mm.cnt<3 & (giant.table$WT_0> cutoff_t0 | giant.table$NMDm_0> cutoff_t0)

    d_d1.filter= giant.table$start.align<5 & giant.table$end.align>96 & giant.table$pam.indel.cnt==0 & giant.table$pam.mm.cnt==0  & giant.table$u1.indel.cnt==0 &
    giant.table$u1.mm.cnt==0 & giant.table$u.indel.cnt==0 & giant.table$u.mm.cnt==0 & giant.table$d.indel.cnt<3 & giant.table$d.mm.cnt<3 & giant.table$d1.indel.cnt<3 &
    giant.table$d1.mm.cnt<3 & (giant.table$WT_0> cutoff_t0 | giant.table$NMDm_0> cutoff_t0)

    u1_d_d1.filter= giant.table$start.align<5 & giant.table$end.align>96 & giant.table$pam.indel.cnt==0 & giant.table$pam.mm.cnt==0  & giant.table$u1.indel.cnt<3 &
    giant.table$u1.mm.cnt<3 & giant.table$u.indel.cnt==0 & giant.table$u.mm.cnt==0 & giant.table$d.indel.cnt<3 & giant.table$d.mm.cnt<3 & giant.table$d1.indel.cnt<3 &
    giant.table$d1.mm.cnt<3 & (giant.table$WT_0> cutoff_t0 | giant.table$NMDm_0> cutoff_t0)


    giant.table$perfect.filter= perfect.filter
    giant.table$clipping.filter= clipping.filter
    giant.table$d1.filter= d1.filter
    giant.table$d_d1.filter= d_d1.filter
    giant.table$u1_d_d1.filter=u1_d_d1.filter

    return(giant.table)
}

# BUILD OLIGO DATA STRUCTURE ---------------------------------------------------------------------------------------------------------
    # Build annotation structures 
    #nPAMs$CDS_length=as.vector(unlist(protein.lengths[nPAMs$GENEID]))/3
    #nPAMs$dist_from_CDS_end=nPAMs$CDS_length-unlist(nPAMs$PROTEINLOC)
    
    # Load various annotation files ----------------------------------------------------------------------------------------
    load('/media/jbloom/d1/coupled_CRISPR/Reference/eStops_oligos.RData')
    # 'oligos' contains information about each oligo
    
    #information in name about position is now incorrect
   # how oligo names were constructed for reference contigs
    oligo.name=apply(oligos,1,function(x) paste(x[c(1:3, 5:23, 25,26)], collapse=':'))
    oligo.name=gsub(' ' , '', oligo.name) 
    oligos$oligo=oligo.name #rownames(oligos)=oligo.name

    # fix discrepant lengths---------------------------------------------------
    aa=read.fasta('/media/jbloom/d1/coupled_CRISPR/Reference/orf_trans.fasta')
    aa.lengths=sapply(aa,length)
    #plengths=sapply(protein.lengths, function(x)x/3)
    #plengths=plengths[na.omit(match(names(aa.lengths), names(plengths)))]
    
    aa.lengths=aa.lengths[na.omit(match(unique(oligos$GENEID), names(aa.lengths)))]
    g.lengths=oligos$CDS_length[match(unique(oligos$GENEID), oligos$GENEID)]
    names(g.lengths)=unique(oligos$GENEID)
    g.lengths=g.lengths[match(names(aa.lengths), names(g.lengths))]

    discrep.lengths=aa.lengths[g.lengths!=aa.lengths]
    ogm=oligos$GENEID[which(oligos$GENEID %in% names(discrep.lengths))]
    oligos$CDS_length[which(oligos$GENEID %in% names(discrep.lengths))]=discrep.lengths[ogm]
    oligos$dist_from_CDS_end[which(oligos$GENEID %in% names(discrep.lengths))]=   oligos$CDS_length[which(oligos$GENEID %in% names(discrep.lengths))]-oligos$PROTEINLOC[which(oligos$GENEID %in% names(discrep.lengths))]
    #--------------------------------------------------------------------------

    rm(aa, aa.lengths, g.lengths, discrep.lengths, ogm)

    # flag oligos as coming from dubious ORFs-----------------------------------------------------------
    load('/media/jbloom/d1/coupled_CRISPR/Reference/dubious_essential_genes.Rdata')
    dubious.e.oligos=which(oligos$GENEID %in% dubious.essential)

    oligos$dubious=FALSE
    oligos$dubious[dubious.e.oligos]=TRUE
    #---------------------------------------------------------------------------------------------------

    # flag oligos that are actually causing synonymous changes -----------------------------------------
    synolis=read.delim(file='/media/jbloom/d1/coupled_CRISPR/Reference/ESS_synonymous_controls.csv', header=T, sep='\t')
    flagsyn=sort(as.vector(unlist(sapply(as.character(synolis$guide), function(x)  grep(x, oligos$guide)))))
    oligos$flagsyn=FALSE
    oligos$flagsyn[flagsyn]=TRUE
    
    # additional genes flagged by Meru as non-essential
    load(file='/media/jbloom/d1/coupled_CRISPR/Reference/additional_genes_to_drop_not_essential.RData') 
    genes.to.drop=as.vector(unlist(genes.to.drop))
    oligos$drop=oligos$GENEID %in% genes.to.drop
    oligos$drop=oligos$drop | oligos$dubious  # modified 02/27/17 | oligos$flagsyn 

     mpam=sapply(oligos$stopPAM_oligos, function(x) s2c(x)[51])
     gtRC=rep('', length(mpam))
     gtRC[as.vector(which(mpam!='G'))]=as.character(reverseComplement(DNAStringSet(oligos$stopPAM_oligos[as.vector(which(mpam!='G'))])))
     gtRC[as.vector(which(mpam=='G'))]=(oligos$stopPAM_oligos[as.vector(which(mpam=='G'))])
     oligos$repair.coding.strand=gtRC
     rm(mpam)
     rm(gtRC)

    oligos$unique.Index=1:nrow(oligos)
    ################################################################################################################################################
    # an experiment to hand check crispr induction for 10 transformed colonies 
    #hand.verified=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/hand_verified.txt', header=F, sep=':', stringsAsFactors=F)
    #checkme=lapply( hand.verified[,1], function(x) grep(x, oligos$stopPAM_oligos))
    #checkme=do.call('c', checkme)
    ################################################################################################################################################
    #save(oligos, file=paste0(paste(out.base.dir, 'processed/RData/', sep=''), 'oligoAnnotations.RData'))
    # load oligo annotations 
    #load(paste0(paste(out.base.dir, 'processed/RData/', sep=''), 'oligoAnnotations.RData'))

    data(fcabundance)
    fcabundance$yORF=as.character(fcabundance$yORF)
    data(gfp)
    gfp$yORF=as.character(gfp$yORF)
    half_lives=read.delim('/media/jbloom/d1/coupled_CRISPR/Experiments/082316/ref/pnas_0605420103_SuppDataSet.txt', header=T ,sep='\t', stringsAsFactors=F)[,1:4]
    names(half_lives)[1]='yORF'

    abundance=merge(gfp, fcabundance, by='yORF', all.x=T)
    abundance=merge(abundance, half_lives, by='yORF', all.x=T)

    # some additional information
    evolvability=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/Evolvability.csv', header=T, sep='\t', stringsAsFactors=F)
    #evolvability$Systematic.name 
    #evolvability$X
    comp_human=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/Marcotte_complementsHuman.csv', header=T, sep='\t', stringsAsFactors=F)
    #comp_human$ScENSP
    #comp_human$Final.CompStatus

    #haploinsufficiency
    haploinf=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/haploinsufficiency.txt',  header=T ,sep='\t', stringsAsFactors=F)
    haploinf$orf=gsub(' ','',haploinf$orf)
    # note this introducing NAs
    haploinf$HET_AV=as.numeric(haploinf$HET_AV)
    haploinf$HET_AV=as.numeric(haploinf$HET_AV)

    # additional information about conservation and protein domains
    all.coding.files=list.files('/media/jbloom/d1/coupled_CRISPR/Reference/coding/', pattern='aa.mfa', full.names=T)
    #domains=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/domains.tab', header=F, sep='\t', stringsAsFactors=F)
    domains=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/pfam.tab', header=F, sep='\t', stringsAsFactors=F)
    domains2=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/pfam.tab2', header=F, sep='\t', stringsAsFactors=F)

    #domains=domains[domains[,4]=='Pfam',]
    dgsplit=split(domains2, domains2[,1]) 

    #retain aa
    #registerDoParallel(70)

    conservation=foreach(gene=unique(oligos$GENEID)) %dopar% {
        print(gene)
        coding.file=all.coding.files[grep(paste0('_', gene, '_'), all.coding.files)][1]
        if(file.exists(coding.file)){
        r=bio3d::read.fasta(coding.file)
        r.conserv=conserv(r)
        no.del=which(r$ali[1,]!='-')
        return(r.conserv[no.del]) } else {return(NULL) }
    }
    names(conservation)=unique(oligos$GENEID)

    # for example/ks
    #plot(oligos[match(unique(oligos$GENEID), oligos$GENEID),]$dubious, sapply(conservation, function(x) sum(x>.9999)/length(x)) )
    frac.aa.perfect.conserved=sapply(conservation, function(x) sum(x>.9999)/length(x))
    names(frac.aa.perfect.conserved)=unique(oligos$GENEID)
    mean.aa.perfect.conserved=sapply(conservation, function(x) sum(x)/length(x))
    names(mean.aa.perfect.conserved)=unique(oligos$GENEID)

    # integrate conservation signal from gene end
    end.conservation=rep(NA, nrow(oligos))
    for(i in 1:nrow(oligos) ) {
        cvec=conservation[[oligos[i,'GENEID']]]
        if(is.null(cvec)) { next 
        } else {
        cvec.end=cvec[oligos[i,'PROTEINLOC']:(length(cvec)-1)]
        #end.conservation[i]=sum(cvec.end)/length(cvec.end)
        end.conservation[i]=sum(cvec.end>.9999)/length(cvec.end)
        }
    }

    dn_ds=read.xls('/media/jbloom/d1/coupled_CRISPR/Reference/ranked_dnds.xlsx', skip=1, sheet=3)[,c(1:8)]
    dn_ds$Gene=as.character(dn_ds$Gene)

    load('/media/jbloom/d1/coupled_CRISPR/Reference/overlapping_genes.RData')
    str(overlapping_features)

    # annotations about viability
    viable_table=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/viable_annotations.txt' ,header=T, sep='\t', stringsAsFactors=F, comment.char='!')

    # low complexity regions 
    # the easy part -------------------------------------------
    lc.fasta=read.fasta('/media/jbloom/d1/coupled_CRISPR/Reference/uniprot-sCerevisiaeProteome_160404.out.fasta',forceDNAtolower = FALSE)
    lc.intervals=lapply(lc.fasta, function(x) {
                       x.ind = seq(1:length(x))
                       x.lc  = grepl('[a-z]', x)
                       xrle=rle(x.lc)
                       xrle$values=as.factor(1:length(xrle$values))
                       boundary.set=split(x.ind[x.lc], inverse.rle(xrle)[x.lc])
                       t(sapply(boundary.set[sapply(boundary.set, length)>0], range))
                    })
    #-----------------------------------------------------------

    # keep original name to double check we haven't screwed up name matching
    # rename list
    uniprot.IDs_from_lc.list=sapply(strsplit(names(lc.intervals), '\\|'), function(x)x[2])

    #straight from here
    #http://bioconductor.org/packages/release/data/annotation/manuals/org.Sc.sgd.db/man/org.Sc.sgd.db.pdf
    x = org.Sc.sgdUNIPROT
    # Get the Systematic ORF IDs that are mapped to a Uniprot ID
    mapped_genes <- mappedkeys(x)
    # Convert to a list
    xx = as.list(x[mapped_genes])

    uniprot2sys=data.frame(uniprot.key=as.vector(unlist(xx)), 
               yname.key=rep(names(xx), sapply(xx, length)), 
               stringsAsFactors=F)

    ygene.matches=match(uniprot.IDs_from_lc.list, uniprot2sys[,1])

    for(n in 1:length(lc.intervals)) {
        attr(lc.intervals[[n]],'fasta.name' )=names(lc.intervals)[n]
        attr(lc.intervals[[n]],'uniprot' )=uniprot.IDs_from_lc.list[n]
        attr(lc.intervals[[n]],'systematic.name' )=uniprot2sys[ygene.matches[n],2] 
    }

    has.lowcomp=sapply(lc.intervals, ncol)>0
    sapply(lc.intervals, nrow)
    # now build data.frame
    lc.start.stop=do.call('rbind', lc.intervals[has.lowcomp])
    gene=as.character(rep(sapply(lc.intervals[has.lowcomp], function(x) attr(x,'systematic.name')), sapply(lc.intervals[has.lowcomp],nrow)))
    low.complexity=data.frame(gene=gene, start=lc.start.stop[,1], stop=lc.start.stop[,2], stringsAsFactors=F)
    low.complexity=na.omit(low.complexity)
    attr(low.complexity, 'na.action')=NULL
    lcsplit=split(low.complexity, low.complexity[,1])
    #save(low.complexity,file='/media/jbloom/d1/coupled_CRISPR/Reference/low.complexity.RData')

    # integrate information on low complexity regions and on domain structures
    osg=oligos
    osg=split(osg, osg$GENEID)
    osg=lapply(osg, function(x) {
               x[order(x$dist_from_CDS_end, decreasing=F),]
                   })
    for(g in names(osg)) {
         da=osg[[g]]$CDS_length-osg[[g]]$dist_from_CDS_end #as.numeric(names(so[[g]])) #which(so[[g]]==1)))
        
        # count of regions of low complexity downstream of oligo
        if(!is.null(lcsplit[[g]])) {
            dg=Intervals(lcsplit[[g]][,c(2,3)])
            low.complexity.table=sapply(dg[,2], function(x) da<x)
            low.complexity.downstream.cnt=rowSums(low.complexity.table)
       } else {low.complexity.downstream.cnt=rep(0, length(da)) }

        #domains
        if(!is.null(dgsplit[[g]])) {
            dg=Intervals(dgsplit[[g]][,c(6,7)])
            udc=sapply(dg[,2], function(x) da<x)
            colnames(udc)=dgsplit[[g]][,4]
            # true or value, are there dowmains downstream
            domain.downstream=rowSums(udc)>0
            for(n in 1:ncol(udc)) {
                        udc[udc[,n]==T,n]=colnames(udc)[n]
                        udc[udc[,n]==F,n]='' #colnames(udc)[n]
            }
            #  list the domains
            domains.downstream=apply(udc, 1, paste, collapse=' ')
        } else {
             domain.downstream=rep('FALSE', length(da))
             domains.downstream=rep('', length(da))
        }
          # overwrite osg
        osg[[g]]=data.frame(osg[[g]],     low.complexity.downstream.cnt=low.complexity.downstream.cnt,
                          domain.downstream=domain.downstream,
                         domains.downstream=domains.downstream, stringsAsFactors=F)
    }

    oligos=rbindlist(osg)
    oligos=data.frame(oligos, stringsAsFactors=F)
    oligos$domain.downstream=as.logical(oligos$domain.downstream)
    oligos$domains.downstream=gsub('^\\s*', '' , oligos$domains.downstream)
    oligos$domains.downstream=gsub('\\s*$', '' , oligos$domains.downstream)
    oligos$domains.downstream=gsub('\\s\\s', ' ' , oligos$domains.downstream)
    oligos=oligos[order(oligos$unique.Index),]


    # merge in abundance information
    oligos=data.frame(oligos, abundance[match(oligos$GENEID, abundance$yORF),], stringsAsFactors=F)

    # conservation info (gene level)
    oligos$mean.aa.perfect.conserved=mean.aa.perfect.conserved[oligos$GENEID]
    oligos$frac.aa.perfect.conserved=frac.aa.perfect.conserved[oligos$GENEID]
    # assuming oligos are sorted by unique.Index (oligo level)
    oligos$end.conservation = end.conservation

    #dn_ds (H0_w) 
    oligos=merge(oligos, dn_ds, by.x='GENEID', by.y='Gene', all.x=T, sort=F)
    oligos=oligos[order(oligos$unique.Index),]

    # haploinsufficiency
    oligos = merge(oligos, haploinf, by.x='GENEID', by.y='orf', all.x=T, sort=F)
    oligos=oligos[order(oligos$unique.Index),]

    # viability
    oligos$viable_annotation=oligos$GENEID %in% viable_table$Gene.Systematic.Name

    #evolvability
    oligos$evolvability=evolvability$X[match(oligos$GENEID, evolvability$Systematic.name)]

    # human complements null
    oligos$human_complements_null=comp_human$Final.CompStatus[match(oligos$GENEID, comp_human$ScENSP)]

    # calculate the number of matches for expected guide RNA vs rest of genome 
    # guide GC content 
    dsog=DNAStringSet(oligos$guide)
    lf=letterFrequency(dsog, letters='ACGT', OR=0)
    oligos$guide.GCcontent=(lf[,2]+lf[,3])/rowSums(lf)
    rm(dsog)
    rm(lf)

    oligos$PAMvariantCNT=apply(cbind(t(sapply(oligos$REFCODON, s2c)),t(sapply(oligos$VARCODON, s2c))), 1,
          function(x) sum(x[1:3]!=x[4:6]))

    og=DNAStringSet(paste0(oligos$guide, 'NGG'))
    guideRC=og
    soi=oligos$guideStrand=='-'
    guideRC[soi]=reverseComplement(guideRC[soi])

    source('/media/jbloom/d1/coupled_CRISPR/code/BLAST.R')
    # run once 
    #bl = makeblast(db="/media/jbloom/d1/coupled_CRISPR/Reference/BLAST/sacCer3.fasta")
    #blast_all_guides=predictBLAST(bl, guideRC, 
    #                          BLAST_args='-evalue=100 -word_size=8 -gapopen=5 -gapextend=2 -reward=2 -penalty=-3 -num_threads=70',
    #                          custom_format='sseq qseq'
    #                          )
    #save(blast_all_guides, file ='/media/jbloom/d1/coupled_CRISPR/Reference/BLAST/BLAST_oligos.RData')
    load('/media/jbloom/d1/coupled_CRISPR/Reference/BLAST/BLAST_oligos.RData')

    gi=paste0('Query_', 1:length(guideRC))
    blast_all_guides$ID=match(as.character(blast_all_guides$QueryID), gi)
    blast_all_guides$strand=oligos$guideStrand[blast_all_guides$ID]
    #8  + NGG
    b2=blast_all_guides[(blast_all_guides$strand=='+' & blast_all_guides$Q.end==23 & grepl('GG$', blast_all_guides$sseq)) |
                                      (blast_all_guides$strand=='-' & blast_all_guides$Q.start==1 & grepl('^CC', blast_all_guides$sseq)),] 
    #now filter on exact match for 8 bases upstream of NGG
    b2=b2[ (b2$strand=='-' &  substr(b2$qseq, 4,11)== substr(b2$sseq, 4,11)) |
                         (b2$strand=='+' &  substr(b2$qseq, 13,20)== substr(b2$sseq, 13,20)),]
    b2=b2[b2$Alignment.Length>18,]
    b2=b2[(23-b2$Alignment.Length+b2$Mismatches)<6,]

    its=split(b2, b2$ID)
    itg=sapply(its, nrow)
    itgn=(as.numeric(names(its)))
    # note that 
    oligos$cnt.offtarget=as.vector(itg[match(oligos$unique.Index, itgn)])
    oligos$cnt.offtarget[is.na(oligos$cnt.offtarget)]=50
    oligos$drop=oligos$drop | oligos$cnt.offtarget==50

    #exactly duplicated
    edup=sapply(its, function(x) 
          if(nrow(x)>1) {
           sum(x$Mismatches==1 & x$Gap.Openings==0 & x$Alignment.Length==23)
          } else {1 } )
    edup=ifelse(edup==1,0,1)
    oligos$dup.gRNA=as.vector(edup[match(oligos$unique.Index, itgn)])
    oligos$dup.gRNA[is.na(oligos$dup.gRNA)]=0

    #http://www.flyrnai.org/evaluateCrispr/
    gE=read.xls('/media/jbloom/d1/coupled_CRISPR/Reference/gRNA_efficiency.xls')
    oligos$U6.Terminator=gE$U6.Terminator
    oligos$Score=gE$Score
#--------------------------------------------------------------------------------------------end building oligos data structure










# HMM model
doHMM =function(big.mm, oligo.stats) {
      
    # get priors from big.mm
    #atog.by.gene.D=split(big.mm[big.mm$GENEID %in% oligo.stats$GENEID[oligo.stats$dubious],],
    #                   big.mm[big.mm$GENEID %in% oligo.stats$GENEID[oligo.stats$dubious],'GENEID'])

    #atog.by.gene.D=lapply(atog.by.gene.D, function(x) {
    #                        x=x[order(x$dist_from_CDS_end),];
    #                        x$tt=as.numeric(as.factor(x$dist_from_CDS_end)); 
    #                        return(x);
    #            })
    #table(as.numeric(unlist(sapply(atog.by.gene.D, function(x) x$DA))))

    # Emission probabilities if actual state is alive (built from barcodes for 'dubious genes')
    #   1    2 
    #3789 2887 =   
    #3789/(3789+2887) probability of emitting alive if alive p=0.5676  1-p=.4324
   

    atog.by.gene=split(big.mm[big.mm$GENEID %in% oligo.stats$GENEID[!oligo.stats$drop],],
                       big.mm[big.mm$GENEID %in% oligo.stats$GENEID[!oligo.stats$drop],'GENEID'])

    atog.by.gene=lapply(atog.by.gene, function(x) {
                            x=x[order(x$dist_from_CDS_end),];
                            x$tt=as.numeric(as.factor(x$dist_from_CDS_end)); 
                            return(x); })

    # Emission probabilities if actual state is dead
    # table(as.numeric(unlist(sapply(atog.by.gene, function(x) x$DA[x$tt==max(x$tt)]))))
    #   1    2 
    # 1646 6823
    # probability emit alive if dead = p=.1944   probability emit dead if dead= p=.8056


    initprobs2=c(.5,.5)
    emissionProbs2=rbind(c(.5676, .4324),
                         c(.1944, .8056))
    transitionProbs2=rbind(c(.5,.5),
                           c(0,  1))

    # 3 - state HMM with error state 
    initprobs3=c(.49,.49,.02)
    transitionProbs3=rbind(c(.5,.5,0),
                           c(0 ,.98, .02),
                           c(0,1,0))
    emissionProbs3=rbind(c(.5676, .4324),
                         c(.1944, .8056),
                         c(.5676, .4324))

    

    hmm.out=list()
    for(gene in names(atog.by.gene)){
        print(gene)
       abg=atog.by.gene[[gene]]
       abg$DA=as.numeric(as.vector(unlist(sapply(split(abg$DA, abg$tt), sort))))
             
       observation=split(abg$DA, abg$tt)

       ihmm2=initHMM.mobs(c('A','D'), c('1', '2'), startProbs=initprobs2, transProbs=transitionProbs2, emissionProbs=emissionProbs2)
       mout2=viterbi.mobs(ihmm2, observation)
       hmm2=t(posterior.mobs(ihmm2, observation))
       pout2=hmm2 #hmm2/rowSums(hmm2)
       colnames(pout2)=paste0('posterior', colnames(pout2))
       ihmm3=initHMM.mobs(c('A','D', 'E'), c('1', '2'), startProbs=initprobs3, transProbs=transitionProbs3, emissionProbs=emissionProbs3)
       mout3=viterbi.mobs(ihmm3, observation)
       hmm3=t(posterior.mobs(ihmm3, observation))
       pout3=hmm3 #/rowSums(hmm3)
       colnames(pout3)=paste0('posterior', colnames(pout3),3)
       oname=abg$oligo[match(unique(abg$tt), abg$tt)]
       dff=data.frame(name=oname, hmm2=mout2,pout2,hmm3=mout3,pout3)
       rownames(dff)=oname
       hmm.out[[gene]]=dff
     }
    hmm.out.df=do.call('rbind', hmm.out)
    return(hmm.out.df)
}



# HMM plots 
makeHMMplots=function(outdir, oligo.stats, big.mm, conservation,dgsplit) {

    #updated plotting
    X11.options(type='cairo')
    t.var=c(0,24,48,72,96)
    #pdf(file='~/Desktop/plots.pdf', width=10,height=12)
    dir.create(outdir)

    # split oligo stats by gene
    osg=split(oligo.stats, oligo.stats$GENEID)
    # split raw data by gene 
    bigmmg=split(big.mm, big.mm$GENEID)

    for(i in 1:length(osg)) {
        gene=names(osg)[i]
        og=osg[[gene]]
        bg=bigmmg[[gene]]
        if(is.null(bg) | is.null(og)) {next;} 
        dub.stat=ifelse(og$dubious[1], 'DUBIOUS', 'ESSENTIAL')
        pdf(file=paste0(outdir, gene, ':', dub.stat, '.pdf'), width=12, height=13)

        conserv.vector=conservation[[gene]]

        cds.length=og$CDS_length[1]
        
        bg.s=split(bg, bg$expt)
        bg.s=lapply(bg.s, function(x) split(x, x$dist_from_CDS_end))
        gWgL=bg.s$WT
        gNgL=bg.s$NMD
        
        
        wt.length=length(gWgL)
        nmd.length=length(gNgL)

        wt.pos=as.numeric(names(gWgL))
        nmd.pos=as.numeric(names(gNgL))

        noligos=unique(sort(c(nmd.pos, wt.pos)))
        total.oligos=length(noligos)

        m=cbind(c(2:(total.oligos+1)),c((total.oligos+2):(total.oligos+total.oligos+1)))
        m=rbind(c(1,1), c(1,1), m)
        layout(m)
        par(mar=c(3,1,1,1), oma=c(2,1,3,1))
        if( is.null(conserv.vector) ) {  # length(coding.file)==0 | is.na(coding.file)){
            plot(0,0, xlim=c(1,cds.length), ylim=c(0,1), type='n') 
        } else {  
               plot(conserv.vector, type='l', xlim=c(1,cds.length), ylim=c(0,1), xaxt='n', ylab='fraction aa conserved', col='#00000066', lwd=2)
        }
        # color ticks by HMM output
        hmm.calls=og$hmm2
        hmm.calls3=og$hmm3

        unique.xaxs.locations=unique(c(wt.pos, nmd.pos))

       # text(cds.length-wt.pos, .5, wt.pos, srt=90, col='red', cex=1.5)
       # text(cds.length-nmd.pos, .5, nmd.pos, srt=90, col='red', cex=1.5)
        axis(1, at=c(cds.length-unique.xaxs.locations), labels = unique.xaxs.locations, cex.axis=1) 
         #at=c(1, cds.length))

        if(!is.null(dgsplit[[gene]])) {
            rect(dgsplit[[gene]][,6], 0, dgsplit[[gene]][,7], 1, col='#0000ff33')
            text(dgsplit[[gene]][,6]+20,.3, dgsplit[[gene]][,5], cex=1.2)
        }
       if(!is.null(hmm.calls)){
            if(length(which(hmm.calls=='A'))>0 ) { }
                yp=cds.length-og$dist_from_CDS_end[which(hmm.calls=='A')]
                points(yp, rep(.1, length(yp)), type='l', col='green', pch=20, cex=2, lwd=4, lty=1)
                points(yp, rep(.1, length(yp)), type='p', col='green', pch=20, cex=2, lwd=4, lty=1)

            if(length(which(hmm.calls=='D'))>0){
                yp=cds.length-og$dist_from_CDS_end[which(hmm.calls=='D')]
                points(yp, rep(.1, length(yp)), type='l', col='red', pch=20, cex=2, lwd=4, lty=1)
                points(yp, rep(.1, length(yp)), type='p', col='red', pch=20, cex=2, lwd=4, lty=1)

            }
        }

       #points(og$PROTEINLOC, ifelse(-log10(og$ALL.p.value)>20, 20, -log10(og$ALL.p.value))/20, cex=3, pch='*' )
       prattempt=seq(0,20,1)
       axis(4, at=prattempt/20, labels=prattempt)

        # extract conservation track get info for scer coords 
        #par(mfrow=c(length(gWgL[[i]]), 1), 
        #par(mar=c(2.5,4,1,0), oma=c(1,1,3,1))
        for(jj in 1:total.oligos) {
            if( noligos[[jj]] %in% wt.pos) {
                j=match(noligos[jj], wt.pos) 
                plot(t.var, gWgL[[j]][1,41:45], type='n', xlab='time', ylab='RPKM', ylim=c(0, max( gWgL[[j]][,41:46])), main=names(gWgL)[j] )
                lim <- par("usr")
                if(!is.na(match(noligos[jj], og$dist_from_CDS_end)) ){
                   mind=match(noligos[jj],og$dist_from_CDS_end )
                    if(!is.na(hmm.calls[[mind]])) {
                     if(hmm.calls[[mind]]=='D' ) {
                         rect(lim[1]-1, lim[3]-1, lim[2]+1, lim[4]+1,  col = ifelse(og$posteriorD[mind]>.9, "#ff000033",  "#ff000011"))
                      } 
                      if(hmm.calls[[mind]]=='A') {
                          rect(lim[1]-1, lim[3]-1, lim[2]+1, lim[4]+1, col = ifelse(og$posteriorA[mind]>.7, "#00ff0033",  "#00ff0011")) #"#00ff0011")
                     }
                    }
                }
                 if(!is.na(match(noligos[jj], og$dist_from_CDS_end)) ){
                     mind=match(noligos[jj], og$dist_from_CDS_end )
                    if(!is.na(hmm.calls3[[mind]])) {
                     
                     if(hmm.calls3[[mind]]=='E') {
                     legend('topright', 'E', text.col='brown', cex=1.5, bty='n')
                    }
                    }
                 }
               # if( sum(gWgL[[j]]$flagsyn)>0) { legend('topright', 'flagsyn', text.col='red') }
                for(k in 1:nrow(gWgL[[j]]) ){ points(t.var, gWgL[[j]][k,41:45], type='b', xlab='time', ylab='RPKM', col=gWgL[[j]][k,'matched.barcode']+1 )          
                                                   points(72+rnorm(1,2),   gWgL[[j]][k,46], col='blue') 
                }
            } else {plot(0,0, type='n', axes=FALSE) }
         }

        for(jj in 1:total.oligos) {
            if( noligos[[jj]] %in% nmd.pos) {
                j=match(noligos[jj], nmd.pos) 
                plot(t.var, gNgL[[j]][1,41:45], type='n', xlab='time', ylab='RPKM', ylim=c(0, max( gNgL[[j]][,41:46])), main=names(gNgL)[j] )
                lim <- par("usr")
                 if(!is.na(match(noligos[jj], og$dist_from_CDS_end)) ){
                   mind=match(noligos[jj],og$dist_from_CDS_end  )
                    if(!is.na(hmm.calls[[mind]])) {

                     if(hmm.calls[[mind]]=='D') {
                         rect(lim[1]-1, lim[3]-1, lim[2]+1, lim[4]+1,  col =ifelse(og$posteriorD[mind]>.9, "#ff000033",  "#ff000011") )# "#ff000011")
                      } 
                      if(hmm.calls[[mind]]=='A') {
                          rect(lim[1]-1, lim[3]-1, lim[2]+1, lim[4]+1, col =  ifelse(og$posteriorA[mind]>.7, "#00ff0033",  "#00ff0011") ) #"#00ff0011")
                     }
                    }
                }
                 if(!is.na(match(noligos[jj], og$dist_from_CDS_end)) ){
                     mind=match(noligos[jj], og$dist_from_CDS_end )
                    if(!is.na(hmm.calls3[[mind]])) {
                    
                     if(hmm.calls3[[mind]]=='E') {
                     legend('topright', 'E', text.col='brown',cex=1.5, bty='n')
                         }
                    }
                 }

                #if( sum(gNgL[[j]]$flagsyn)>0) { legend('topright', 'flagsyn', text.col='red') }
                for(k in 1:nrow(gNgL[[j]]) ){ points(t.var, gNgL[[j]][k,41:45], type='b', xlab='time', ylab='RPKM', col=gNgL[[j]][k,'matched.barcode']+1)
                                                   points(72+rnorm(1,2),   gNgL[[j]][k,46], col='blue')
                 }
            } else {plot(0,0, type='n', axes=FALSE) }

        }


        if(!is.na(match(gene, names(gene2alias))) ) { #null(gene2alias[[gene]])) {
            if(!is.na(SYS2ORF[match(gene,SYS2ORF[,1]),2]) ) { 
            title(paste(gene, ',', dub.stat, ',',SYS2ORF[match(gene,SYS2ORF[,1]),2],
                        '\n', gene2alias[[gene]]), outer=T)
            
            } else {
            title(paste(gene, ',', dub.stat, '\n', gene2alias[[gene]]), outer=T)
            }
        } else {
            title(paste(gene, ',', dub.stat), outer=T) 
        }


        dev.off()

        }
}

# takes a named vector of gene scores
# the other needed variables should be in the global namespace
doGO=function(testGenes) {
    geneSel.fx=function(allScore) { return(allScore<1e9) }
    GO.out=list()
    for(thisOntology in c('BP', 'MF', 'CC') ) {
            print(thisOntology)
            #thisOntology='BP'
            #thisOntology='MF'
            #thisOntology='CC'
            #GOData = new("topGOdata", ontology=thisOntology, allGenes = testGenes, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=5)
            GOData = new("topGOdata", ontology=thisOntology, allGenes = testGenes, geneSel=geneSel.fx, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=5)
            GOresult.dec = runTest(GOData, algorithm="classic", statistic="ks", scoreOrder='decreasing')
            GOresult.inc = runTest(GOData, algorithm="classic", statistic="ks", scoreOrder='increasing')

            gt=GenTable(GOData, GOresult.dec, numChar=140, topNodes = length(score(GOresult.dec)))
            gt2=GenTable(GOData, GOresult.inc, numChar=140, topNodes = length(score(GOresult.inc)))
            gt3=merge(gt,gt2,by='GO.ID')
            min.p=ifelse(as.numeric(gt3$result1.x)<as.numeric(gt3$result1.y), as.numeric(gt3$result1.x), as.numeric(gt3$result1.y))
            effect.direction=ifelse(as.numeric(gt3$result1.x)<as.numeric(gt3$result1.y), 'persisting', 'dead')
            
            gt3$min.p=as.numeric(min.p)
            gt3$effect.direction=effect.direction

            gt3=gt3[order(gt3$min.p, decreasing=F),]
            gt3=gt3[,c(1,2,3,12,13)]
            names(gt3)=c('GO.ID', 'Term', 'Annotated', 'p.value', 'class')
            gene.sets=genesInTerm(GOData, (gt3$GO.ID))
            genes.in.set.sorted=sapply(gene.sets, function(x)     paste(as.character(na.omit(x[match(names(testGenes), x)])), collapse=' ')           )
            gt3$genes.in.set.sorted=genes.in.set.sorted
            GO.out[[thisOntology]]=gt3
            print(head(gt,20))
            print(head(gt2,20))
    }
    return(GO.out)
}

doGO.threshold=function(geneschosen, all.tested) {
    GO.out=list()
    for(thisOntology in c('BP', 'MF', 'CC') ) {
            print(thisOntology)
            #thisOntology='BP'
            #thisOntology='MF'
            #thisOntology='CC'
            #GOData = new("topGOdata", ontology=thisOntology, allGenes = testGenes, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=5)
          
            testGenes= factor(0+(all.tested %in% geneschosen))
            names(testGenes)=all.tested

            GOData = new("topGOdata", ontology=thisOntology, allGenes = testGenes,  annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=5)
            GOresult = runTest(GOData, algorithm="classic", statistic="fisher") #, scoreOrder='decreasing')
            #GOresult.inc = runTest(GOData, algorithm="classic", statistic="ks", scoreOrder='increasing')

            gt=GenTable(GOData, GOresult, numChar=140, topNodes = length(score(GOresult)))
            gt$result1=as.numeric(gt$result1)
            gt=gt[order(gt$result1, decreasing=F),]
            names(gt)[6]='p.value'

            gene.sets=genesInTerm(GOData, (gt$GO.ID))
            genes.in.set=sapply(gene.sets, function(x) paste(x[x %in% geneschosen], collapse= ' '))

            gt$genes.in.set=genes.in.set
            GO.out[[thisOntology]]=gt
            print(head(gt,20))
    }
    return(GO.out)
}

# run GLMER and augment oligo.stats with 
doGLMER.oligo.gene=function(formula.input, data.in, oligo.stats, lab, doGene=TRUE) {
    #optimizer="Nelder_Mead",
    binarized.mmR=glmer(formula.input, family=binomial(link="logit"), data=data.in,
         control=glmerControl(optCtrl = list(maxfun = 1e6)), verbose=T)               
    binarized.ranefR=ranef(binarized.mmR, cond=TRUE)
    dA=cbind(binarized.ranefR$oligo[,1],as.vector(attr(binarized.ranefR$oligo, 'postVar')) )
    rownames(dA)=rownames(binarized.ranefR$oligo)
    colnames(dA)=c('binarized.oligo.blup', 'binarized.oligo.blup.postVar' )
    if(lab=='') { } else{  colnames(dA)=paste(colnames(dA), lab, sep='.')    }

    if(doGene) {
        gr=cbind(binarized.ranefR$GENEID[,1], as.vector(attr(binarized.ranefR$GENEID, 'postVar')))
        colnames(gr)=c( 'binarized.gene.blup', 'binarized.gene.blup.postVar')
        if(lab=='') { } else{   colnames(gr)=paste(colnames(gr), lab, sep='.')      }
        rownames(gr)=rownames(binarized.ranefR$GENEID)
            
        # rewrite without using merge ()
        os= data.frame(oligo.stats, dA[match(oligo.stats$oligo, rownames(dA)),], stringsAsFactors=F)
        os= data.frame(os, gr[match(oligo.stats$GENEID, rownames(gr)),], stringsAsFactors=F)
        return(os)
    } else{
        os= data.frame(oligo.stats, dA[match(oligo.stats$oligo, rownames(dA)),], stringsAsFactors=F)
        return(os)
    }

  }

doLMER.oligo.gene=function(formula.input, data.in, oligo.stats, lab='', doGene=TRUE) {
    #optimizer="Nelder_Mead",
    #slope.mmR=glmer(formula.input, family=binomial(link="logit"), data=data.in,
    #     control=glmerControl(optCtrl = list(maxfun = 1e6)), verbose=T)               
    slope.mmR=lmer(formula.input, data=data.in, verbose=T)
    slope.ranefR=ranef(slope.mmR, cond=TRUE)
    dA=cbind(slope.ranefR$oligo[,1],as.vector(attr(slope.ranefR$oligo, 'postVar')) )
    rownames(dA)=rownames(slope.ranefR$oligo)
    colnames(dA)=c('slope.oligo.blup', 'slope.oligo.blup.postVar' )
    if(lab=='') { } else{
    colnames(dA)=paste(colnames(dA), lab, sep='.')
    }

    if(doGene) {
        gr=cbind(slope.ranefR$GENEID[,1], as.vector(attr(slope.ranefR$GENEID, 'postVar')))
        colnames(gr)=c( 'slope.gene.blup', 'slope.gene.blup.postVar')
        if(lab=='') { } else{  colnames(gr)=paste(colnames(gr), lab, sep='.')    }
        rownames(gr)=rownames(slope.ranefR$GENEID)
            
        # rewrite without using merge ()
        os= data.frame(oligo.stats, dA[match(oligo.stats$oligo, rownames(dA)),], stringsAsFactors=F)
        os= data.frame(os, gr[match(oligo.stats$GENEID, rownames(gr)),], stringsAsFactors=F)
        return(os)
    } else{
        os= data.frame(oligo.stats, dA[match(oligo.stats$oligo, rownames(dA)),], stringsAsFactors=F)
        return(os)
    }

  }


resampler <- function(data) {
    n <- nrow(data)
    resample.rows <- sample(1:n,size=n,replace=TRUE)
    return(data[resample.rows,])
}

spline.estimator <- function(data,m=300) {
    fit <- smooth.spline(x=data[,1],y=data[,2],cv=TRUE)
    eval.grid <- seq(from=min(data[,1]),to=max(data[,1]),length.out=m)
    return(predict(fit,x=eval.grid)$y) # We only want the predicted values
}

spline.cis <- function(data,B,alpha=0.05,m=300) {
    spline.main <- spline.estimator(data,m=m)
    spline.boots <- replicate(B,spline.estimator(resampler(data),m=m))
    cis.lower <- 2*spline.main - apply(spline.boots,1,quantile,probs=1-alpha/2)
    cis.upper <- 2*spline.main - apply(spline.boots,1,quantile,probs=alpha/2)
    return(list(main.curve=spline.main,lower.ci=cis.lower,upper.ci=cis.upper,
    x=seq(from=min(data[,1]),to=max(data[,1]),length.out=m)))
}




#as.vector(edup[match(oligos$unique.Index, itgn)])
#y=oligos$binarized.oligo.blup[oligo.stats$unique.Index %in%itgn]
#itg.b=ifelse(itg==1,0,1)
#t.test(y~itg.b)
#cor.test(y, itgc)
#plot(y, itgc)
#summary(lm(y~itg.b+itgc))
#cor.test(itg.b, itgc)
#summary(lm(y~itg.b+itgc))
#plot(oligo.stats$binarized.oligo.blup[oligo.stats$unique.Index %in%itgn], itg)
#todo refcodon to varcodon count

#t.test(oligo.stats$binarized.oligo.blup~oligo.stats$dup.gRNA)
#
#binot=ifelse(oligos$cnt.offtarget==1,0,1)
#with(oligo.stats, {
#    summary(lm(scale(binarized.oligo.blup)~-1+guide.GCcontent+binot+(PAMvariantCNT)))})
#o2=sapply(split(oligo.stats$binarized.oligo.blup, oligo.stats$PAMvariantCNT), na.omit)
#   Mitg.b+itgc))

#buildSacCer3_genome_dictionaryFR=function(sacCer3, unique.chrs) {
#    sc3.dna=list()
#    for(chr in  unique.chrs ) { 
#        sc3.dna[[chr]]=DNAString(sacCer3[[chr]])
#    }
#    sc3.set=DNAStringSet(sc3.dna)
#}
#
#unique.chrs=paste0('chr', as.roman(1:16))
#sc3.set=buildSacCer3_genome_dictionaryFR(sacCer3, unique.chrs)
#
#build.4base.pdict=function(guide.seq) { 
#    guide.seeds.N.pdict=list()
#    for(base in c('A', 'C', 'G', 'T') ) {
#        print(base)
#        gtable.seq=sapply(as.character(guide.seq), s2c)[1:16,]
#        gtable.seq[14,]=base
#        all.guide.seeds=as.vector(apply(gtable.seq,2, c2s))
#        seed.set=DNAStringSet(unlist(all.guide.seeds))
#        guide.seeds.N.pdict[[base]]=PDict(seed.set, tb.start=15,tb.width=2)
#    }
#    return(guide.seeds.N.pdict)
#}
#
#guideForward=guideRC[!soi]
#gFp=build.4base.pdict(subseq(guideForward, 8,23)[1:10]) #guideForward)
#test=vcountPDict(gFp[[4]], sc3.set, with.indels=F, max.mismatch=2) 
#test=matchPDict(gFp[[4]], sc3.set[[1]], with.indels=F, max.mismatch=2) 
#
##with.indels=F, max.mismatch=2) # with.indels=F)
#test=matchPDict(gFp[[4]][1], sc3.set[[1]], with.indels=F, max.mismatch=2) # with.indels=F)
#
#test=vcountPDict(subseq(guideForward, 8,23)[1:10], sc3.set, min.mismatch=1, max.mismatch=2, with.indels=F) # with.indels=F)
#test=vmatchPDict(subseq(guideForward, 8,23)[1:10], sc3.set, min.mismatch=1, max.mismatch=2, with.indels=F) # with.indels=F)
#
#
#
#
#seed.set.dict=PDict(guideRC[1], tb.start=1, tb.end=13) #,tb.start=13, tb.end=13)
#vmatchPattern(seed.set.dict, sc3.set) # with.indels=F)
#perfect.matches=vcountPDict(seed.set.dict,  sc3.set, min.mismatch=0, max.mismatch=0, with.indels=F)
#
#mat = nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = FALSE, type='DNA')
##989
#test=pairwiseAlignment(pattern=guideRC[1],    
#                       subject=sc3.set[[16]], 
#                       type='local', 
#                       substitutionMatrix=mat, gapOpening=5, gapExtension=2)
#mpas=mclapply(tas, function(x) {
#                     pairwiseAlignment(pattern=DNAStringSet(x), 
#                                        subject=subject.seq,type='local', substitutionMatrix=mat, gapOpening=5, gapExtension=2) }, mc.cores=threads)




# set up for GO enrichment analysis -----------------------------------------------------------------------------------------------
#gcoord.key= build.gcoord.key('/data/eQTL/reference/sacCer3.fasta')


    

#
#        ag1=gBsA #gos[[gn]] #gt.WT.pe[g1,]
#
#        #ag1=do.call('rbind', gos)
#        # for WT
#        Xg1=data.matrix(ag1[,9:13])
#        Xg1[Xg1[,2]==0,2]=1
#        b   =  as.factor(paste0(ag1$barcode, ag1$expt)) #as.factor(ag1$barcode)
#        y1  =  as.vector(Xg1)
#        t1  =  rep(t.var, each=nrow(Xg1))
#        b1  =  rep(b, ncol(Xg1))
#        tot.count2=colSums(Xg1)
#        ofl =  rep(tot.count2, each=nrow(Xg1))
#        expt = rep(as.factor(ag1$expt), ncol(Xg1))
#        olig.fact = rep(as.factor(rownames(ag1)), ncol(Xg1))
#        gene.fact= rep(as.factor(ag1$GENEID), ncol(Xg1))
#
#
#        g3= glmer(y1~offset(log(ofl))+(t1|gene.fact)+(t1|olig.fact)+(t1|b1)+(1|b1), 
#                                   family=poisson(link=log) ,nAGQ=0,
#                                   control=glmerControl(optimizer = "nloptwrap") )
#
#        #,control=glmerControl(optimizer="Nelder_Mead"))
#
#    
#        
#
#
## print(cor.test(y1, exp(predict(g1))))
## print(anova(g0, gN1, refit=FALSE))
# py1=matrix(exp(predict(gN)), nrow(Xg1), ncol(Xg1))
# py1=matrix(exp(pff.p), nrow(Xg1), ncol(Xg1))
#
#
# plot(t1, y1, typ='n', ylim=c(0,300))
#    for(i in 1:nrow(Xg1) ) {
#        points(t.var, Xg1[i,], type='b')
#        points(t.var, py1[i,], type='b', col='blue')
#        readline()
#    }
#

#ioi=c(989,990,1553,3110,3111,3112,3113,5181,5422,5423,6383,6386,6387,6393,6818,7101,7102,7139,8799,9237,9238,9239,9240,9241,9242,9244,9245,9936)
#fisher.test(rbind(c(5,28),c(981,10971)) )


#sample data
#data<-data.frame(x=rnorm(100), y=rnorm(100))

#run and plot
#sp.cis <- spline.cis(data, B=1000,alpha=0.05)
#lines(x=sp.cis$x,y=sp.cis$main.curve)
#lines(x=sp.cis$x,y=sp.cis$lower.ci, lty=2)
#lines(x=sp.cis$x,y=sp.cis$upper.ci, lty=2)




# using a cutoff of cut, calculate sensitivity, specificity, and classification rate
#mod=glm(slope.binarized~cnt+p_intercept+expt+binot+U6.Terminator+Score, family=binomial(link="logit"), data=big.mm)
#mod=lm(slope~cnt+p_intercept+expt+binot+U6.Terminator+Score, data=big.mm)
#
#y=big.mm$slope.binarized
#perf = function(cut, mod, y)
#{
#   yhat = (mod$fit>cut)
#   w = which(y==1)
#   sensitivity = mean( yhat[w] == 1 ) 
#   specificity = mean( yhat[-w] == 0 ) 
#   c.rate = mean( y==yhat ) 
#   d = cbind(sensitivity,specificity)-c(1,1)
#   d = sqrt( d[1]^2 + d[2]^2 ) 
#   out = t(as.matrix(c(sensitivity, specificity, c.rate,d)))
#   colnames(out) = c("sensitivity", "specificity", "c.rate", "distance")
#   return(out)
#}
#
#s = seq(.01,.99,length=1000)
#OUT = matrix(0,1000,4)
#for(i in 1:1000) OUT[i,]=perf(s[i],mod,y)
#plot(s,OUT[,1],xlab="Cutoff",ylab="Value",cex.lab=1.5,cex.axis=1.5,ylim=c(0,1),type="l",lwd=2,axes=FALSE,col=2)
#axis(1,seq(0,1,length=5),seq(0,1,length=5),cex.lab=1.5)
#axis(2,seq(0,1,length=5),seq(0,1,length=5),cex.lab=1.5)
#lines(s,OUT[,2],col="darkgreen",lwd=2)
#lines(s,OUT[,3],col=4,lwd=2)
#lines(s,OUT[,4],col="darkred",lwd=2)
#box()
#legend(0,.25,col=c(2,"darkgreen",4,"darkred"),lwd=c(2,2,2,2),c("Sensitivity","Specificity","Classification Rate","Distance"))








## previous HMM code
# # HMM output 
#    # get prior probs for oligo closest to end --------------------------------------------------------------------------------------
#    #ogs=split(oligos, oligos$GENEID)
#    #og.closest=sapply(ogs, function(x) x$unique.Index[which.min(x$dist_from_CDS_end)])
#    #oblup.ind=rownames(red.effs$oligo)
#    #match(oblup.ind, rownames(oligos))
#    #R> sum(red.effs$oligo[na.omit(match(og.closest, match(oblup.ind, rownames(oligos)))),]>0)
#    #[1] 570
#    #R> length(red.effs$oligo[na.omit(match(og.closest, match(oblup.ind, rownames(oligos)))),]>0)
#    #[1] 909
#    ##570/909
#     #ih2=baumWelch.mobs(ihmm2, observation, maxIterations=5)
#       #posterior.mobs(ih2$hmm, observation)
#       #viterbi.mobs(ih2$hmm, observation)
#       #hmm4=t(posterior.mobs(ih2$hmm, observation))
#       #hmm4/rowSums(hmm4)
#     # prior probabilities of alive and dead at first position closest to end of gene
#    #tda1=sapply(sapply(atog.by.gene, function(x) x$DA[x$tt==1]) , function(x) table(x))
#    #tda2=sapply(sapply(atog.by.gene, function(x) x$DA[x$tt==2]) , function(x) table(x))
#    #tda3=sapply(sapply(atog.by.gene, function(x) x$DA[x$tt==3]) , function(x) table(x))
#    #tda4=sapply(sapply(atog.by.gene, function(x) x$DA[x$tt==4]) , function(x) table(x))
#    #tda5=sapply(sapply(atog.by.gene, function(x) x$DA[x$tt==5]) , function(x) table(x))
#    #tda6=sapply(sapply(atog.by.gene, function(x) x$DA[x$tt==6]) , function(x) table(x))
#    #tda7=sapply(sapply(atog.by.gene, function(x) x$DA[x$tt==7]) , function(x) table(x))
#    #tda8=sapply(sapply(atog.by.gene, function(x) x$DA[x$tt==8]) , function(x) table(x))
#    #tda9=sapply(sapply(atog.by.gene, function(x) x$DA[x$tt==9]) , function(x) table(x))
#    #tda10=sapply(sapply(atog.by.gene, function(x) x$DA[x$tt==10]) , function(x) table(x))
#
#    #table(apply(tda1, 2, which.max))
#    #  1   2 
#    #445 589
#    #p(alive)=.4304  p(dead)=.5696
#      
#    #p(alive|d0) =.4304*.5676+.5696*.4324 
#    #p(dead|d0) = .4304*.1944+.5696*.8056 
#    #table(apply(tda2, 2, which.max))
#    #sum(apply(tda2, 2, which.max)==2 &  apply(tda1, 2, which.max)==1)/sum(apply(tda1,2,which.max)==1)
#    #sum(apply(tda3, 2, which.max)==2 &  apply(tda2, 2, which.max)==1)/sum(apply(tda2, 2, which.max)==1)
#    #sum(apply(tda4, 2, which.max)==2 &  apply(tda3, 2, which.max)==1)/sum(apply(tda3, 2, which.max)==1)
#    #sum(apply(tda5, 2, which.max)==2 &  apply(tda4, 2, which.max)==1)/sum(apply(tda4, 2, which.max)==1)
#    #sum(apply(tda6, 2, which.max)==2 &  apply(tda5, 2, which.max)==1)/sum(apply(tda5, 2, which.max)==1)
#    #sum(apply(tda7, 2, which.max)==2 &  apply(tda6, 2, which.max)==1)/sum(apply(tda6, 2, which.max)==1)
#    #sum(apply(tda8, 2, which.max)==2 &  apply(tda7, 2, which.max)==1)/sum(apply(tda7, 2, which.max)==1)
#    #sum(apply(tda9, 2, which.max)==2 &  apply(tda8, 2, which.max)==1)/sum(apply(tda8, 2, which.max)==1)
#    #o=list()
#    #o[[1]]=rep(c('1','2'), apply(tda1,1,sum))
#    #o[[2]]=rep(c('1','2'), apply(tda2,1,sum))
#    #o[[3]]=rep(c('1','2'), apply(tda3,1,sum))
#    #o[[4]]=rep(c('1','2'), apply(tda4,1,sum))
#    #o[[5]]=rep(c('1','2'), apply(tda5,1,sum))
#    #o[[6]]=rep(c('1','2'), apply(tda6,1,sum))
#    #o[[7]]=rep(c('1','2'), apply(tda7,1,sum))
#    #o[[8]]=rep(c('1','2'), apply(tda8,1,sum))
#    #o[[9]]=rep(c('1','2'), apply(tda9,1,sum))
#    #o[[10]]=rep(c('1','2'), apply(tda10,1,sum))
#    #initprobs2=c(.63,.37)
#    #transitionProbs2=rbind(c(.1,.9),
#    #                       c(0,  1))
#    # top row was 60% and 40%
#    # from g.alive reset to 74% and 26%
#    #emissionProbs2=rbind(c(.74, .26), 
#    #                     c(.26, .74))
#
#    # rows are states, columns are observations
#    # two observable states and 3 hidden 
#    #initprobs3=c(.62, .35, .02)
#    #emissionProbs3=rbind(  c(.74,   .26),
#    #                       c(.26,  .74 ),
#    #                       c(.74,   .26 ))
#                                #A     D     #E
#    #transitionProbs3=rbind( c(.1,   .9,      0),
#    #                        c( 0,   .98,   .02),
#    #                        c(0,     1,     0))
#    #initprobs2=c(.6105, .3985)
#    #initprobs2.new=initprobs2
#    #mout2=doHMM.2state(abg,transitionProbs2,emissionProbs2, initprobs2)
##   viterbi.msm(mout2)
##


