#================================================================
# basic functions
#================================================================
libs = c('dplyr','DBI','RSQLite','dbplyr', 'stringr', 'shiny', 'shinydashboard', 'networkD3', 'igraph', 'ggplot2', 'ggsci', 'DT', 'plotly', 'forestplot', 'data.table', 'htmlTable')
lapply(libs, require, character.only = TRUE) 
options(stringsAsFactors=F, width = 250)

## pattern of number
# e.g., grepl(pattern_number, '1e-5')
pattern_number = "^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$"

## determine if a vec is string number, omit na
is_numeric <- function(x) {
  !any(is.na(suppressWarnings(as.numeric(na.omit(x))))) & is.character(x)
}

## extract element with pattern from vector
# e.g., get('a', c('AS', 'Ba', 'C')) ; get('a', c('AS', 'Ba', 'C'), exact=F) 
get = function(key, vector, exact=F, v=F){
    if (exact==T){out =vector[grep(key, vector, invert=v)]}
    else {out = vector[grep(toupper(key), toupper(vector), invert=v)]}
    return(out)
}
## replace element in df
# e.g., map(data.frame(t=c('a', 'b')), c('a', 'b'), c('xxx', 'yyy'))
map = function(df, keys, values){
  for (i in 1:length(keys)){df[df==keys[i]] = values[i]}
  return(df)
}


#================================================================
# var
#================================================================
file_mrdb = 'data/db_mr.csv'
file_mrdb = 'data/collect/mr_collect.csv'
file_ao = 'data/ao.rdata'
file_pairs = 'data/pairs.rdata'

# mr data
methods = c('ivw', 'egger', 'wmedian', 'wmode') # combine Wald ratio with IVW, when nsnp=1
names = c('exposure', 'outcome', 'nsnp', paste0(rep(methods, each=3), c('_b', '_se', '_p')))

dbs = c('UKB', 'FIN', 'SWE', 'COVID-19 HGI', 'IEU')

# shiny
method_names = c('IVW', 'MR-Egger', 'Weighted median', 'Weighted mode') # in shiny


#================================================================
# process mr db
#================================================================
## combinePvalVector, from netbid2
# e.g., combinePvalVector(c(0.1,1e-3,1e-5))
# 
combinePvalVector = function(pvals, method = 'Stouffer', signed = TRUE, twosided = TRUE) {
  # remove NA pvalues
  pvals <- pvals[!is.na(pvals) & !is.null(pvals)]
  pvals[which(abs(pvals)<=0)] <- .Machine$double.xmin
  if (sum(is.na(pvals)) >= 1) {
    stat <- NA
    pval <- NA
  } else{
    if (twosided & (sum(pvals > 1 | pvals < -1) >= 1))
      stop('pvalues must between 0 and 1!\n')
    if (!twosided & (sum(pvals > 0.5 | pvals < -0.5) >= 1))
      stop('One-sided pvalues must between 0 and 0.5!\n')
    if (!signed) {
      pvals <- abs(pvals)
    }
    signs <- sign(pvals)
    signs[signs == 0] <- 1
    if (grepl('Fisher', method, ignore.case = TRUE)) {
      if (twosided & signed) {
        neg.pvals <- pos.pvals <- abs(pvals) / 2
        pos.pvals[signs < 0] <- 1 - pos.pvals[signs < 0]
        neg.pvals[signs > 0] <- 1 - neg.pvals[signs > 0]
      } else{
        neg.pvals <- pos.pvals <- abs(pvals)
      }
      pvals <-
        c(1, -1) * c(
          pchisq(
            -2 * sum(log(as.numeric(pos.pvals))),
            df = 2 * base::length(pvals),
            lower.tail = FALSE
          ) / 2,
          pchisq(
            -2 * sum(log(as.numeric(neg.pvals))),
            df = 2 * base::length(pvals),
            lower.tail = FALSE
          ) / 2
        )
      pval <- base::min(abs(pvals))[1]
      # if two pvals are equal, pick up the first one
      stat <-
        sign(pvals[abs(pvals) == pval])[1] * qnorm(pval, lower.tail = F)[1]
      pval <- 2 * pval
    }
    else if (grepl('Stou', method, ignore.case = TRUE)) {
      if (twosided) {
        zs <- signs * qnorm(abs(pvals) / 2, lower.tail = FALSE)
        stat <- sum(zs) / sqrt(base::length(zs))
        pval <- 2 * pnorm(abs(stat), lower.tail = FALSE)
      }
      else{
        zs <- signs * qnorm(abs(pvals), lower.tail = FALSE)
        stat <- sum(zs) / sqrt(base::length(zs))
        pval <- pnorm(abs(stat), lower.tail = FALSE)
      }
    }
    else{
      stop('Only Fisher or Stouffer method is supported !!!\n')
    }
  }
  return(c(`Z-statistics` = stat, `P.Value` = pval))
}


#================================================================
# mr functions
#================================================================
# extract link with db, nsnp, and pthres of 4 methods #exposure_db%in%input$db & outcome_db%in%input$db
select_links = function(input, db_mr){
  links = db_mr%>%dplyr::filter(nsnp>=input$thres_nsnp &
    ivw_p<input$pthres_ivw & egger_p<input$pthres_egger & wmedian_p<input$pthres_wmedian & wmode_p<input$pthres_wmode)
  links = as.data.frame(links)
  # print(sprintf('select_links done! N links: %s', nrow(links)))
  return(links)
}

# get expo and outcome traits
extract_traits = function(links){
  traits = unique(unlist(links[, c("exposure_name", "outcome_name")], use.names = F))
  traits = traits[order(!grepl("^FIN", traits), traits)] # here to control the order
  # print(sprintf('extract_traits done! N trait: %s', length(traits)))
  return(traits)
}

# get link start from exposure, with given step
select_graph = function(links, trait=input$trait, depth=input$net_depth){
  # print(sprintf('select_graph with %s', trait))
  traits = c(trait)
  graph = data.frame()
  for (i in 1:depth){
    add = links%>%dplyr::filter(exposure_name%in%traits|outcome_name%in%traits)
    graph = rbind(graph, add)
    traits = c(traits, extract_traits(graph)) # update traits to extract new network
  }
  graph = graph%>%distinct()%>%mutate_if(is_numeric,as.numeric)
  # print(sprintf('select_graph done! N edge: %s', nrow(links)))
  # print(graph[,c('exposure_name', 'outcome_name')])
  return(graph)
}

# graph to igraph
make_igraph = function(graph){
  vertices = extract_traits(graph)
  igraph = graph_from_data_frame(graph%>%dplyr::select(exposure_name, outcome_name), directed = TRUE, vertices = vertices)
  return(igraph)
}

# igraph to d3graph
make_d3graph = function(igraph, links=links, trait=input$trait){
  d3graph = igraph_to_networkD3(igraph)
  # add group
  groups = sapply(unlist(d3graph$nodes), function(x){strsplit(x, ':')[[1]][1]}); names(groups) = NULL
  d3graph$nodes = d3graph$nodes%>%mutate(group=groups)%>%mutate(group=ifelse(name==trait, 'The focus trait', group))
  d3graph$links$value=1
  return(d3graph)
}

# harmo dat
get_mrdat = function(graph, exposure=input$exposure, outcome=input$outcome){
  temp = graph %>% dplyr::filter(exposure_name==!!exposure & outcome_name==!!outcome)
  file = sprintf('data/collect/harmo_dat/%s_TO_%s.rdata', temp$exposure, temp$outcome)
  load(file)
  dat = mr_res[['dat']]
  return(dat)
}

# get mr_res for plot
get_mr_res = function(graph, exposure, outcome){
  mr_res = list()
  dat = get_mrdat(graph, exposure, outcome)
  line = graph%>%dplyr::filter(exposure_name==!!exposure&outcome_name==!!outcome)
  coef = as.data.frame(matrix(unlist(line[paste0(rep(methods, each=3), c('_b', '_se', '_p'))]), ncol=3, byrow=T))
  coef = cbind(method_names, line$nsnp, coef, c(0, line$egger_int, 0, 0))
  names(coef) = c('method', 'snp', 'b', 'se', 'p', 'a') # a is intercept
  het = line[c('ivw_q', 'ivw_q_p', 'egger_q', 'egger_q_p')]
  pleio = line[c('egger_int', 'egger_int_se', 'egger_int_p')]
  dat = dat%>%dplyr::select(-id.exposure, -id.outcome, -exposure, -outcome, -mr_keep)%>%
    rename(A1=effect_allele.exposure, A2=other_allele.exposure, FRQ=eaf.exposure)
  mr_res[['dat']] = dat
  mr_res[['coef']] = coef
  mr_res[['supp']] = c(het, pleio)
  return(mr_res)
}  



mr_heatmap = function(graph){
    bs = c(); ps = c()
    for (method in methods){
        for (outcome in graph$outcome_name){
        b = graph%>%dplyr::filter(outcome_name==!!outcome)%>%pull(paste0(method, '_b'))
        p = graph%>%dplyr::filter(outcome_name==!!outcome)%>%pull(paste0(method, '_p'))
        bs = c(bs, method, outcome, b)
        ps = c(ps, method, outcome, p)
        }
    }
    mat_b = data.frame(matrix(bs, ncol=3, byrow=T))%>%mutate_if(is_numeric,as.numeric)
    mat_p = data.frame(matrix(ps, ncol=3, byrow=T))%>%mutate_if(is_numeric,as.numeric)
    mat_b = map(mat_b, methods, method_names)
    mat_p = map(mat_p, methods, method_names)
    mat_b$method = factor(method_names, levels=method_names)
    mat_p$method = factor(method_names, levels=method_names)
    plot = ggplot(data = mat_b, aes(x=X1, y=X2, fill=X3)) + 
        xlab('Methods') + 
        labs(fill = "Causal effect (β)") +
        geom_tile(width = 0.1)+
        theme_bw()+
        theme(axis.text = element_text(size = 11),
          axis.title = element_text(size = 15),
          axis.title.y = element_blank()) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
    return(plot)
}

mr_forest = function(graph, exposure, outcome){
    method_names = c('IVW', 'MR-Egger', 'Weighted median', 'Weighted mode')
    df = graph%>%dplyr::filter(exposure_name==!!exposure & outcome_name == !!outcome)%>%mutate_if(is_numeric,as.numeric)
    nsnp = df$nsnp 
    b = df[,get('_b', names(df))]
    se = df[,get('_se', names(df))][1:4]
    p = df[,get('_p', names(df))][1:4]
    p = sapply(p, function(x){sprintf("%.2e", x)})
    
    tabletext=list(list(expression(bold('Method')), method_names[1], method_names[2], method_names[3], method_names[4]), 
        list(expression(bold(italic(N)['SNP'])), nsnp, nsnp, nsnp, nsnp), list(expression(bold(italic(P)['MR'])), p[1], p[2], p[3], p[4]))
    
    plot = forestplot(tabletext,
        mean  = cbind(unlist(c(NA, b))), 
        lower = cbind(unlist(c(NA, b-1.96*se))), 
        upper = cbind(unlist(c(NA, b+1.96*se))),
        is.summary=c(FALSE, rep(FALSE, 4)), vertices=T, 
        new_page=T, boxsize=0.15, line.margin=0.15, 
        xlog=F, xlab=expression("Causal Effect" * (italic(β))), 
        col=fpColors(box=c("blue"), line="darkblue", summary="royalblue"),
        txt_gp = fpTxtGp (label=gpar(cex=1.25),ticks=gpar(cex=1.1),xlab=gpar(cex=1.2),title=gpar(cex=1.2)))
    return(plot)
}

mr_scatter = function (res, dat){
    # TwoSampleMR package reverse some effect size, is problematic
    # index = dat$beta.exposure < 0
    # dat$beta.exposure[index] = dat$beta.exposure[index] * -1
    # dat$beta.outcome[index] = dat$beta.outcome[index] * -1
    p = ggplot(data = dat, aes(x = beta.exposure, y = beta.outcome)) +
        geom_errorbar(aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome), colour = "grey", width = 0) +
        geom_errorbarh(aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), colour = "grey", height = 0) +
        geom_point(size = 3, shape = 16, alpha = 0.6) +
        geom_abline(data = res, aes(intercept = a, slope = b, colour = method), show.legend = TRUE) +
        labs(color = "Method", x = "GWAS effect on exposure", y = "GWAS effect on outcome") +
        theme_bw() +
        theme(legend.position = 'bottom', legend.direction = "vertical",
            legend.text = element_text(size = 12),
            panel.background = element_blank(),
            axis.text = element_text(size = 12, margin = margin(t = 10, r = 10, b = 10, l = 10)),  # axis tick
            axis.title = element_text(size = 14, margin = margin(t = 10, r = 10, b = 10, l = 10))) +  # axis title
        guides(color = guide_legend(title.position = "left", nrow = 1)) 
    print(p)
}


#================================================================
# temp
#================================================================
## compare degree
# e.g., compare_degree(links, 'UKB: MONOCYTE PERCENTAGE', 'IEU: MONOCYTE COUNT')
compare_degree = function(links, trait1, trait2) {
  n1 = length(unique(links %>% dplyr::filter(exposure == trait1) %>% pull(outcome)))
  n2 = length(unique(links %>% dplyr::filter(exposure == trait2) %>% pull(outcome)))
  trait_better = ifelse(length(n1) > length(n2), trait1, trait2)
  res = sprintf("%s: %d outcomes \n%s: %d outcomes\n%s is better\n", trait1, n1, trait2, n2, trait_better)
  cat(res)
  return(trait_better)
}