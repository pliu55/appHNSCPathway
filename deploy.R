#!/usr/bin/env Rscript

#
#  pliu 20190124
#
#  prepare input data for app
#

library(data.table)
library(stringr)
library(RColorBrewer)
library(scales)
library(rsconnect)

main <- function() {
    ntwdir = '/ua/pliu/comp/ahlquist/network/'

    prm = list(
        set = 'training',

        greys = brewer.pal(6, 'Greys'),
        red   = brewer.pal(6, 'Reds' )[4],
        blue  = brewer.pal(6, 'Blues')[4],

        fnode = paste0(ntwdir, '03_pathway/node.tsv'),
        fedge = paste0(ntwdir, '03_pathway/edge.tsv'),
        fmut  = paste0(ntwdir, '04_nMut.tsv'),
        fprtb = paste0(ntwdir, '17_nPrtb/n_prtb_199.tsv'),
        fmed  = paste0(ntwdir, '18p_cmpHpvPrtb/median_199.tsv'),
        fpval = paste0(ntwdir, '19_IplOmicP/logrank_p.tsv'),

        fout_med  = './median.tsv',
        fout_node = './node.tsv',
        fout_edge = './edge.tsv'
    )

    nodedt = fread(prm$fnode, header=T)
    edgedt = fread(prm$fedge, header=T)
    mutdt = fread(prm$fmut, header=T, select=c('set', 'hpv', 'pid', 
                  'entity', 'n_mut'))[ set == prm$set ]
    prtbdt = fread(prm$fprtb, header=T)
    pvaldt = getPval(prm$fpval, prm$pval_cutoff)

    mrgdt = merge(pvaldt, prtbdt, by=c('hpv', 'pid', 'entity'), all=T)
    mutdt[, set := NULL]
    mrgmrgdt = merge(mrgdt, mutdt, by=c('hpv', 'pid', 'entity'), all.x=T)

    out_edgedt = selDefEdge(edgedt, mrgmrgdt, prm)
    out_nodedt = selDefNode(nodedt, mrgmrgdt, prm)

    write.table(out_nodedt, prm$fout_node, quote=F, sep="\t", row.names=F)
    write.table(out_edgedt, prm$fout_edge, quote=F, sep="\t", row.names=F)
    file.copy(prm$fmed, prm$fout_med)

    cat('file written:', prm$fout_node, "\n")
    cat('file written:', prm$fout_edge, "\n")
    cat('file written:', prm$fout_med , "\n")


    deployApp( appDir  = './', 
               appName = 'HNSC_pathway', 
               account = 'uwbiostat')
}


selDefNode <- function(in_nodedt, mrgmrgdt, prm) {
    pth_nodedt = in_nodedt[ pid %in% mrgmrgdt$pid ]

    ## entities with the same name, but different type
    ## e.g. PAK1, complex/protein in pid 36885
    uni_nodedt = pth_nodedt[, list(type = paste(type, collapse='/')), 
                            by=list(pid, pth, source, entity)]

    nodedt = merge(mrgmrgdt, uni_nodedt, by=c('pid', 'entity'), all=T)
    setnames(nodedt, 'entity', 'id')

    nodedt[, `:=`( mut_frc  = n_mut/n_patient,
                   prtb_frc = n_prtb/n_patient)]

    nodedt[, `:=`(
        s_omic = ifelse(has_omic == FALSE, '<li>no mutation/CNA/RNA data</li>',
                        paste0('<li>CNA: p=', round(p_cna, digits=4), '</li>',
                               '<li>RNA: p=', round(p_rna, digits=4), '</li>',
                               '</ul>',
                               'mutation: ', n_mut, '/', n_patient, '=',
                               percent(mut_frc)) ),
        shape = ifelse(prtb_frc > 0.5, 'star',
                       ifelse(prtb_frc > 0, 'triangle', 'dot')),
        size  = ifelse(prtb_frc > 0.5, 35, ifelse(prtb_frc > 0, 20, 15)),
        color = ifelse(prtb_frc > 0,
                       ifelse((! is.na(p_ipl)) & (p_ipl < 0.05),
                              ifelse(((! is.na(p_cna)) & (p_cna < 0.05) ) |
                                     ((! is.na(p_rna)) & (p_rna < 0.05) ),
                                     prm$red, prm$blue), prm$greys[2]), 
                       prm$greys[2])
    )]

    nodedt[, `:=`(
        label = ifelse( (color == prm$red) & (prtb_frc > 0.5), id, NA),
        font.size = 45,
        title = paste0(id, ' (', type, ')<br>',
            'perturbed patients: ', n_prtb, '/', n_patient, '=',
            percent(prtb_frc), '<br>',
            'correlation with overall survival:<br>',
            '<ul>', '<li>pathway level: p=', round(p_ipl, digits=3), '</li>',
            s_omic )
    )]

    out_nodedt = nodedt[, .(id, hpv, pth, shape, size, color, label, font.size,
                            title)] 

    return(out_nodedt)
}


selDefEdge <- function(in_edgedt, mrgmrgdt, prm) {
    edgedt = in_edgedt[ pid %in% mrgmrgdt$pid ]

    edgedt[, `:=`(arrows = ifelse(title == 'component>', 'middle;to', 'to'),
                  width  = 2,
                  color  = prm$greys[2],
                  dashes = ifelse(str_sub(title, -1, -1) == '|', T, F) )]

    out_edgedt = edgedt[, .(pth, from, to, title, arrows, width, color, dashes)]
    return(out_edgedt)
}


getPval <- function(fin, pval_cutoff) {
    ## ignore intermediate node (*__$i) created for max_in_degree=5
    indt = fread(fin, header=T)[ ! grepl('__', entity) ]

    ipldt  = indt[ feat == 'ipl' ]
    omicdt = indt[ feat != 'ipl' ]

    setnames(ipldt, c('entity', 'n_repressed', 'n_activated', 'n_normal', 'p'), 
             c('pid_entity', 'n_repressed_ipl', 'n_activated_ipl', 
               'n_normal_ipl', 'p_ipl'))
    ipldt[, `:=`( pid    = as.integer(gsub('(\\d{5})_(.*)', '\\1', pid_entity)),
                  entity = gsub('(\\d{5})_(.*)', '\\2', pid_entity) )]
    ipldt[, c('feat', 'pid_entity') := NULL]

    wide_omicdt = dcast(omicdt, hpv + entity ~ feat, value.var=c( 
                        'n_repressed', 'n_activated', 'n_normal', 'p'))
    wide_omicdt[, has_omic := TRUE]
    mrgdt = merge(ipldt, wide_omicdt, by=c('hpv', 'entity'), all.x=T)
    mrgdt[, has_omic := ifelse(is.na(has_omic), FALSE, has_omic)]

    return(mrgdt)
}


system.time( main() )
