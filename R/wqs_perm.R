#' WQS permutation test
#' 
#' \code{wqs_perm} takes a `gwqs` object as an input and runs the permutation test (Day 
#' et al, 2022) to obtain an estimate for the p-value significance for the WQS coefficient.  
#' 
#' `wqs_perm` uses a `gwqs` object (from the gWQS package) as an input. To use `wqs_perm`,
#' we first need to run an initial WQS regression while setting `validation=0`. We will 
#' use this `gwqs` object as the model argument for the `wqs_perm` function. 
#' 
#' The argument `boots` is the number of bootstraps for the WQS run in each permutation 
#' test iteration. Note that we may elect a bootstrap count `boots` lower than that 
#' specified in the model object for the sake of efficiency. If `boots` is not specified, 
#' then we will use the same bootstrap count in the permutation test WQS runs as that 
#' specified in the model argument.
#'
#' The arguments `b1_pos` and `rs` should be consistent with the inputs chosen in the 
#' model object. The seed should ideally be consistent with the seed set in the 
#' model object, though this is not required.
#' 
#' For full details on how to use this function, please reference the vignette.
#'
#' @param model A \code{gwqs} object as generated from the \code{gWQS} package.  
#' @param niter Number of permutation test iterations. 
#' @param boots Number of bootstrap samples for each permutation test \code{wqs} run.  
#' If `boots` is not specified, then we will use the same bootstrap count in the 
#' permutation test WQS runs as that specified in the main WQS run.
#' @param b1_pos A logical value that indicates whether beta values should be positive 
#' or negative.
#' @param rs A logical value indicating whether random subset implementation should be 
#' performed. 
#' @param plan_strategy Evaluation strategy for the plan function. You can choose among 
#' "sequential", "transparent", "multisession", "multicore", "multiprocess", "cluster" 
#' and "remote." See gWQS documentation for full details. 
#' @param seed (optional) Random seed for the permutation test WQS reference run. 
#' This should be the same random seed as used for the main WQS run. 
#'
#' @return \code{wqs_perm} returns an object of class `wqs_perm`, which contains 
#' three sublists: 
#' 
#' \item{perm_test}{Contains three objects: (1) `pval`: p-value for the proportion of 
#' permuted WQS coefficient values greater than the reference value, (2) `testbeta1`: reference WQS coefficient beta1 value, 
#' (3) `betas`: Vector of beta values from each permutation test run.}
#' \item{gwqs_main}{Main gWQS object (same as model input)}
#' \item{gwqs_perm}{Permutation test reference gWQS object (NULL if same number of bootstraps
#' as main gWQS object)}
#' @import gWQS ggplot2 viridis cowplot
#' @export wqs_perm
#'
#' @examples
#' library(gWQS)
#'
#' # mixture names
#' PCBs <- names(wqs_data)[1:34]
#' 
#' # create reference wqs object with 1000 bootstraps
#' wqs_main <- gwqs(yLBX ~ wqs, mix_name = PCBs, data = wqs_data, q = 10, validation = 0,
#'                  b = 1000, b1_pos = T, plan_strategy = "multicore", family = "gaussian", seed = 16)
#' 
#' # run permutation test
#' perm_test_res <- wqs_perm(wqs_main, niter = 10, b1_pos = T)
#' 
#' # Note: The default value of niter = 200 is the recommended parameter values. 
#' # This example has a lower niter in order to serve as a shorter test run. 
#' 
#' @references
#' Day, D., Collett, B., Barrett, E., ... & Sathyanarayana, S. (2021). Phthalate mixtures 
#' in pregnancy, autistic traits, and adverse childhood behavioral outcomes. Environment 
#' International, 147, 106330.
#'
#' Loftus, C. T., Bush, N. R., Day, D.... & LeWinn, K. Z. (2021). Exposure to prenatal 
#' phthalate mixtures and neurodevelopment in the Conditions Affecting Neurocognitive 
#' Development and Learning in Early childhood (CANDLE) study. Environment international, 
#' 150, 106409.
#' 

wqs_perm <- function(model, niter = 200, boots = NULL, b1_pos = TRUE, rs = FALSE, 
                    plan_strategy = "multicore", seed = NULL) {
  
  if (class(model) == "gwqs") {
    if (model$family$family != "gaussian" | model$family$link != "identity"){
      stop("The permutation test is currently only set up to accomodate the Gaussian family with an 
           identity link.")
    }
  } else stop("'model' must be of class 'gwqs' (see gWQS package).")
  
  mm <- model$fit
  formchar <- as.character(formula(mm))
  
  if (!is.null(model$stratified) | grepl("wqs:", formchar[3], fixed = T))
  {
    # TODO: We should be able to accomodate stratified weights though we haven't tested that yet,
    # and I'm not sure it makes sense to have stratified weights without a WQS interaction term.
    stop("This permutation test is not yet set up to accomodate stratified weights or 
         WQS interaction terms.")
  }  
  
  cl = match.call()
  Data <- model$data[model$vindex, -which(names(model$data) %in% c("wqs", "wghts"))]
  yname <- as.character(formula(mm))[2]
  mix_name <- names(model$bres)[names(model$bres) %in% model$final_weights$mix_name]

  if (!is.null(model$qi)) {
    nq <- max(sapply(model$qi, length)) - 1
  } else {
    # this is for cases when there is no quantile transformation or it's already been
    # done in the data frame
    nq <- NULL
  }

  # reference WQS run 
  if (is.null(boots)){
    boots <- length(model$bindex)
  }
  
  if (boots == length(model$bindex)){
    perm_ref_wqs <- model
    ref_beta1 <- mm$coef[2]
  }
  
  else{
    perm_ref_wqs <- gwqs(formula = formula(mm), data = Data, mix_name = mix_name, 
                         q = nq, b = boots, rs = rs, validation = 0, plan_strategy = plan_strategy,
                         b1_pos = b1_pos, seed = seed)
    
    ref_beta1 <- perm_ref_wqs$fit$coef[2]
  }
  
  
  if (length(mm$coef) > 2) {
    # This is the permutation test algorithm when there are multiple independent variables in
    # the model
    lm_form <- formula(paste0(formchar[2], formchar[1], gsub("wqs + ", "", formchar[3], fixed = T)))
    fit.partial <- lm(lm_form, data = Data)
    partial.yhat <- predict(fit.partial)
    partial.resid <- resid(fit.partial)
    reorgmat <- matrix(NA, dim(Data)[1], niter)
    reorgmat <- apply(reorgmat, 2, function(x) partial.yhat + sample(partial.resid, replace = F))
  } else {
    # This is the permutation test algorithm when there is only one independent variable in
    # the model
    reorgmat <- matrix(NA, dim(Data)[1], niter)
    reorgmat <- apply(reorgmat, 2, function(x) sample(Data[, yname]))
  }
  
  getbetas <- function(x) {
    
    newDat <- Data
    newDat[, yname] <- x
    names(newDat) <- c(names(Data))
    formchar <- as.character(formula(mm))
    
    if (length(mm$coef) > 2) {
      form1 <- formula(paste0(formchar[2], formchar[1], formchar[3]))
    } else {
      form1 <- formula(paste0(formchar[2], formchar[1], "wqs"))
    }
    
    gwqs1 <- tryCatch({
      suppressWarnings(gwqs(formula = form1, data = newDat, mix_name = mix_name, 
                            q = nq, b = boots, rs = rs, validation = 0, plan_strategy = plan_strategy,
                            b1_pos = b1_pos))
    }, error = function(e) NULL, 
      warning = function(e) ifelse(rs == TRUE, message("WQSRS failed"), message("WQS failed")))
    
    if (is.null(gwqs1))
      lm1 <- NULL else lm1 <- gwqs1$fit
    if (is.null(lm1)) {
      retvec <- NA
    } else {
      retvec <- lm1$coef[2]
    }
    return(retvec)
  }
  
  pbapply::pboptions(type = "timer")
  betas <- pbapply::pbapply(reorgmat, 2, getbetas)
  
  if (any(is.na(betas))) {
    print(paste0(length(which(is.na(betas))), " failed model attempts"))
  }
  
  calculate_pval <- function(x, true, posb1 = b1_pos) {
    if (posb1) {
      length(which(x > true))/length(betas)
    } else {
      length(which(x < true))/length(betas)
    }
  }
  
  pval <- calculate_pval(betas, ref_beta1, b1_pos)
  
  perm_retlist <- list(pval = pval, testbeta1 = ref_beta1, betas = betas, call = cl)
  
  model$b1_pos <- b1_pos
  perm_ref_wqs$b1_pos <- b1_pos
    
  ret_ref_wqs <- ifelse(boots == length(model$bindex), NULL, perm_ref_wqs)
  
  results <- list(gwqs_main = model, 
                  gwqs_perm = ret_ref_wqs, 
                  perm_test = perm_retlist)
  
  class(results) <- "wqs_perm"
  
  results
}

#' @rawNamespace S3method(print, wqs_perm)
#' @rdname methods

print.wqs_perm <- function(x, ...){
  
  cat("Permutation test WQS coefficient p-value: \n", 
      x$perm_test$pval,
      "\n")

  main_sum <- summary(x$gwqs_main)

  print(main_sum)

}

#' @rawNamespace S3method(summary, wqs_perm)
#' @rdname methods

summary.wqs_perm <- function(x, ...){
  
  cat("Permutation test WQS coefficient p-value: \n", 
      x$perm_test$pval,
      "\n")
  
  main_sum <- summary(x$gwqs_main)

  print(main_sum)

}

#' Plotting method for wqsperm object 
#'
#' @param wqspermresults an object of class “wqs_perm”
#' @param FixedPalette If true, the heatmap color key for the mixture weights has 
#' ategorical cutoffs with the following categories: <0.1, 0.1 - <0.2, 0.2 - <0.3, 
#' and ≥0.3. If false, the heatmap color key is continuous and dependent on the 
#' weight values.
#' @param InclKey If true, a horizontal heatmap legend is included at the bottom 
#' of the full plot.
#' @param AltMixName Defaults to NULL. If not NULL, these are alternative names 
#' for the mixture components to be displayed on the heatmap y axis.
#' @param AltOutcomeName Defaults to NULL. If not NULL, this is an alternative name 
#' for the outcome to be displayed on the heatmap x axis.
#' @param ViridisPalette  This is the color palette to be used for the viridisLite 
#' package-based coloring of the heatmap, with possible values from “A” to “E”. 
#' Defaults to “D”.
#' @param StripTextSize This is the text size for the plot strip labels. Defaults to 14.
#' @param AxisTextSize.Y This is the text size for the y axis text. Defaults to 12.
#' @param AxisTextSize.X This is the text size for the x axis text. Defaults to 12.
#' @param LegendTextSize This is the text size for the legend text. Defaults to 14.
#' @param PvalLabelSize This is the geom_text size for the permutation test p-value 
#' label. Defaults to 5.
#' @param HeatMapTextSize This is the geom_text size for the mixture weight heatmap labels. 
#' Defaults to 5.
#'
#' @return Returns a list of plots. 
#' \item{Fullplot}{}
#' \item{CoefPlot}{}
#' \item{WtPlot}{}
#' \item{WtLegend}{}
#' 
#' @export
#'
#' @examples
plot.wqs_perm <- function(wqspermresults, FixedPalette=F, InclKey=F, AltMixName=NULL, 
                          AltOutcomeName=NULL, ViridisPalette="D", StripTextSize=14, 
                          AxisTextSize.Y=12, AxisTextSize.X=12, LegendTextSize=14, 
                          PvalLabelSize=5, HeatMapTextSize=5) {

  thisfit<-wqspermresults$gwqs_main$fit
  b1pos<-wqspermresults$gwqs_main$b1_pos
  if(b1pos) thisdir<-"Positive" else thisdir<-"Negative"
  if(!is.null(AltOutcomeName)) outname<-AltOutcomeName else 
    outname<-as.character(attr(thisfit$terms,"variables")[[2]])
  WQSResults<-data.frame(Outcome=outname,Direction=thisdir,Beta=thisfit$coef['wqs'],
                         LCI=suppressMessages(confint(thisfit)[2,1]),
                         UCI=suppressMessages(confint(thisfit)[2,2]),pval=summary(thisfit)$coef["wqs","Pr(>|t|)"],
                         PTp=wqspermresults$perm_test$pval)
  WQSResults$PTlabel<-paste0("PTp=",signif(WQSResults$PTp,3))
  WQSResults$FacetLabel<-"Coefficient"
  cirange<-WQSResults$UCI-WQSResults$LCI
  widercirange<-c(WQSResults$LCI-(WQSResults$LCI/10),WQSResults$UCI+(WQSResults$UCI/10))
  if(widercirange[1]<0&widercirange[2]>0){
    gg1<-ggplot(WQSResults,aes(x=Outcome,y=Beta))+geom_point(size=3)+theme_bw()+
      geom_errorbar(aes(ymin=LCI,ymax=UCI),size=1,width=0.75)+
      geom_hline(yintercept=0)+
      geom_text(aes(label=PTlabel,y=UCI+cirange/10),size=PvalLabelSize)+
      facet_grid(FacetLabel~Direction)+
      theme(strip.text=element_text(size=StripTextSize),
            axis.text.y=element_text(size=AxisTextSize.Y),
            axis.text.x=element_blank(),
            axis.title=element_blank(),
            axis.ticks.x=element_blank())
  } else {
    gg1<-ggplot(WQSResults,aes(x=Outcome,y=Beta))+geom_point(size=3)+theme_bw()+
      geom_errorbar(aes(ymin=LCI,ymax=UCI),size=1,width=0.75)+
      geom_text(aes(label=PTlabel,y=UCI+cirange/10),size=PvalLabelSize)+
      facet_grid(FacetLabel~Direction)+
      theme(strip.text=element_text(size=StripTextSize),
            axis.text.y=element_text(size=AxisTextSize.Y),
            axis.text.x=element_blank(),
            axis.title=element_blank(),
            axis.ticks.x=element_blank())
  }
  
  WQSwts<-wqspermresults$gwqs_main$final_weights[wqspermresults$gwqs_main$mix_name,]
  WQSwts$FacetLabel<-"Weights"
  WQSwts$Outcome<-WQSResults$Outcome
  WQSwts$Direction<-WQSResults$Direction
  WQSwts$mix_name<-factor(as.character(WQSwts$mix_name),levels=wqspermresults$gwqs_main$mix_name)
  if(!is.null(AltMixName)) levels(WQSwts$mix_name) <- AltMixName
  WQSwts$mix_name<-factor(WQSwts$mix_name,levels=rev(levels(WQSwts$mix_name)))
  names(WQSwts)[1:2]<-c("Exposure","Weight")
  if(FixedPalette){
    mypal<-viridis::viridis_pal(option=ViridisPalette)(4)
    WQSwts$Wt<-WQSwts$Weight
    WQSwts$Weight<-factor(ifelse(WQSwts$Wt<0.1,"<0.1",
                                 ifelse(WQSwts$Wt>=0.1&WQSwts$Wt<0.2,"0.1-0.2",
                                        ifelse(WQSwts$Wt>=0.2&WQSwts$Wt<0.3,"0.2-0.3",
                                               paste0("\u2265","0.3")))),
                          levels=c("<0.1","0.1-0.2","0.2-0.3",paste0("\u2265","0.3")))
    Virclr<-ifelse(WQSwts$Weight=="<0.1",mypal[1],
                   ifelse(WQSwts$Weight=="0.1-0.2",mypal[2],
                          ifelse(WQSwts$Weight=="0.2-0.3",mypal[3],
                                 ifelse(is.na(WQSwts$Weight)==T,"grey50",mypal[4]))))
    names(Virclr)<-as.character(WQSwts$Weight)
    legplot<-ggplot(data.frame(Weight=factor(levels(WQSwts$Weight),levels=levels(WQSwts$Weight))),
                    aes(x=1,y=Weight))+
      geom_tile(aes(fill=Weight))+scale_fill_manual(values=mypal)+
      theme(legend.position = "bottom",legend.title=element_text(size=14,face="bold"),
            legend.text=element_text(size=14))
    l1<-cowplot::get_legend(legplot)
    
    gg2<-ggplot(WQSwts,aes(x=Outcome,y=Exposure))+theme_classic()+
      geom_tile(aes(fill=Weight),alpha=0.7)+geom_text(aes(label=round(Wt,2)),size=HeatMapTextSize)+
      scale_fill_manual(values=Virclr)+
      facet_grid(FacetLabel~Direction)+
      theme(strip.text.x=element_blank(),
            strip.text.y=element_text(size=StripTextSize),
            axis.text.x=element_text(size=AxisTextSize.X),
            axis.text.y=element_text(size=AxisTextSize.Y),
            axis.title=element_blank(),
            strip.background.y=element_rect(fill = "grey85",colour = "grey20"),
            legend.position="bottom",
            legend.title = element_text(size=LegendTextSize,face="bold"),
            legend.text = element_text(size=LegendTextSize))
  } else {
    gg2<-ggplot(WQSwts,aes(x=Outcome,y=Exposure))+theme_classic()+
      geom_tile(aes(fill=Weight),alpha=0.7)+geom_text(aes(label=round(Weight,2)),size=HeatMapTextSize)+
      scale_fill_viridis_c(option=ViridisPalette)+
      facet_grid(FacetLabel~Direction)+
      theme(strip.text.x=element_blank(),
            strip.text.y=element_text(size=StripTextSize),
            axis.text.x=element_text(size=AxisTextSize.X),
            axis.text.y=element_text(size=AxisTextSize.Y),
            axis.title=element_blank(),
            strip.background.y=element_rect(fill = "grey85",colour = "grey20"),
            legend.position="bottom",
            legend.title = element_text(size=LegendTextSize,face="bold"),
            legend.text = element_text(size=LegendTextSize),
            legend.key.size = unit(0.4,units='in'))
    l1<-cowplot::get_legend(gg2)
  }
  
  if(InclKey){
    gg2<-gg2+theme(legend.position="none")
    fullplot<-cowplot::plot_grid(cowplot::plot_grid(gg1,gg2,ncol=1,align="v",rel_heights=c(0.4,0.6)),l1,
                                 ncol=1,rel_heights=c(1,0.1))
  } else {
    gg2<-gg2+theme(legend.position="none")
    fullplot<-cowplot::plot_grid(gg1,gg2,ncol=1,rel_heights=c(0.4,0.6),align="v")
  }
  
  return(list(FullPlot=fullplot,CoefPlot=gg1,WtPlot=gg2,WtLegend=l1))
  
}
