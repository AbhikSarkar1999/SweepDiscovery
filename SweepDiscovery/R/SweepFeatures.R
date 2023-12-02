#'@title SweepFeatures
#' @param path Datasets
#' @param chromosome Chromosome Number
#' @param start Start Number
#' @param finish Finish Number
#' @import stats utils randomForest PopGenome
#' @return
#' \itemize{
#'   \item df: Results
#' }
#' @export
#'
#' @examples
#' \donttest{
#' library("SweepDiscovery")
#' Path<-file.path(system.file("exdata/ExampleVCF.vcf.gz", package = "SweepDiscovery"))
#' pred<-SweepFeatures(path=Path,"1",20253,1976067)
#' }
#' @references
#' \itemize{
#'\item Pavlidis, P., Alachiotis, N. A survey of methods and tools to detect recent and strong positive selection. J of Biol Res-Thessaloniki 24, 7 (2017). https://doi.org/10.1186/s40709-017-0064-0
#' }

SweepFeatures<-function(path,chromosome,start,finish) {

  num <- floor((finish-start+1)/11)
  end_points <- seq(from = start, to = finish, by = num)
  end_points[length(end_points)+1] <- finish
  ChrSeg <- list()
  for (i in 1:11) {
    ChrSeg[[i]]<- readVCF(path,1000,chromosome,end_points[i], end_points[i+1])
  }

  df <- data.frame(matrix(nrow=1, ncol=0))
  ChrSeg[1:11]
  for (i in 1:11) {
    class<-ChrSeg[[i]]
    GENOME.class.slide <- sliding.window.transform(class,2000,1000, type=2)
    genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
      split <- strsplit(x," ")[[1]][c(1,3)]
      val <- mean(as.numeric(split))
      return(val)
    })
    # total length
    L<-length(GENOME.class.slide@region.names)
    #Diversity Statistics
    slide <- diversity.stats(GENOME.class.slide)
    #Calculation of Pi
    nucdiv <- slide@Pi
    nucdhapdiv <- slide@hap.diversity.within
    Hd<-sum(nucdhapdiv)
    Hd<-Hd/L
    pi_wini<-(0.0081*((Hd)**2))
    df <- cbind(df, pi_wini)
    colnames(df)[i] <- paste0("pi_win",i)
  }

  for (i in 1:11) {
    class<-ChrSeg[[i]]
    GENOME.class.slide <- sliding.window.transform(class,2000,1000, type=2)
    genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
      split <- strsplit(x," ")[[1]][c(1,3)]
      val <- mean(as.numeric(split))
      return(val)
    })
    # total length
    L<-length(GENOME.class.slide@region.names)
    #Diversity Statistics
    slide <- diversity.stats(GENOME.class.slide)
    #Calculation of Wattersons_theta
    head(nucdiv)
    slide<- neutrality.stats(GENOME.class.slide,detail=TRUE)
    get.neutrality(slide)
    theta_Watterson<-slide@theta_Watterson
    theta_Watterson<-na.omit(theta_Watterson)
    thetawat<-sum(theta_Watterson)
    thetaW_wini<-(thetawat/L)
    df <- cbind(df, thetaW_wini)
    colnames(df)[i+11] <- paste0("thetaW_win",i)
  }

  for (i in 1:11) {
    class<-ChrSeg[[i]]
    GENOME.class.slide <- sliding.window.transform(class,2000,1000, type=2)
    genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
      split <- strsplit(x," ")[[1]][c(1,3)]
      val <- mean(as.numeric(split))
      return(val)
    })
    # total length
    L<-length(GENOME.class.slide@region.names)
    #Diversity Statistics
    slide <- diversity.stats(GENOME.class.slide)
    #Calculation of Tajima.D Statistic
    head(nucdiv)
    slide<- neutrality.stats(GENOME.class.slide,detail=TRUE)
    get.neutrality(slide)
    Tajima.D<-slide@Tajima.D
    ActTaj<-na.omit(Tajima.D)
    Tajima.DSum<-sum(ActTaj)
    tajD_wini<-(Tajima.DSum/L)
    df <- cbind(df, tajD_wini)
    colnames(df)[i+22] <- paste0("tajD_win",i)
  }

  for (i in 1:11) {
    class<-ChrSeg[[i]]
    GENOME.class.slide <- sliding.window.transform(class,2000,1000, type=2)
    genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
      split <- strsplit(x," ")[[1]][c(1,3)]
      val <- mean(as.numeric(split))
      return(val)
    })
    # total length
    L<-length(GENOME.class.slide@region.names)
    #Diversity Statistics
    slide <- diversity.stats(GENOME.class.slide)
    #Calculation of Kelly.ZnS Statistic
    GENOME.class.slide<-linkage.stats(GENOME.class.slide)
    get.linkage(GENOME.class.slide)
    Kelly.Z_nS<-GENOME.class.slide@Kelly.Z_nS
    Kelly.Z_nS<-na.omit(Kelly.Z_nS)
    kelly.Z_nsSUM<-sum(Kelly.Z_nS)
    ZnS_wini<-(kelly.Z_nsSUM/L)
    df <- cbind(df, ZnS_wini)
    colnames(df)[i+33] <- paste0("ZnS_win",i)
  }

  for (i in 1:11) {
    class<-ChrSeg[[i]]
    GENOME.class.slide <- sliding.window.transform(class,2000,1000, type=2)
    genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
      split <- strsplit(x," ")[[1]][c(1,3)]
      val <- mean(as.numeric(split))
      return(val)
    })
    # total length
    L<-length(GENOME.class.slide@region.names)
    #Diversity Statistics
    slide <- diversity.stats(GENOME.class.slide)
    #Calculation of  Omega Statistic
    GENOME.class.slide<-linkage.stats(GENOME.class.slide)
    get.linkage(GENOME.class.slide)
    Kelly.Z_nS<-GENOME.class.slide@Kelly.Z_nS
    Kelly.Z_nS<-na.omit(Kelly.Z_nS)
    kelly.Z_nsSUM<-sum(Kelly.Z_nS)
    ZnS_wini<-(kelly.Z_nsSUM/L)
    k<-ZnS_wini
    k
    Omega_wini<-(-((log(1-(k/2)))/2))
    Omega_wini
    df <- cbind(df, Omega_wini)
    colnames(df)[i+44] <- paste0("Omega_win",i)
  }
  df
  return(df)}

#'@title SweepPrediction
#' @param Datapath Datasets
#' @param Chromosome Chromosome Number
#' @param Start Start Number
#' @param Finish Finish Number
#' @import stats utils randomForest PopGenome
#' @return
#' \itemize{
#'   \item Prediction: Results
#' }
#' @export
#'
#' @examples
#' \donttest{
#' library("SweepDiscovery")
#' Path<-file.path(system.file("exdata/ExampleVCF.vcf.gz", package = "SweepDiscovery"))
#' pred<-SweepPrediction(Datapath=Path,"1",20253,1976067)
#' }
#' @references
#' \itemize{
#'\item Pavlidis, P., Alachiotis, N. A survey of methods and tools to detect recent and strong positive selection. J of Biol Res-Thessaloniki 24, 7 (2017). https://doi.org/10.1186/s40709-017-0064-0
#' }

SweepPrediction<-function(Datapath,Chromosome,Start,Finish) {

  #Creating Windows/Partitions

  num <- floor((Finish-Start+1)/11)
  end_points <- seq(from = Start, to = Finish, by = num)
  end_points[length(end_points)+1] <- Finish

  #Listing 11 windows
  ChrSeg <- list()
  for (i in 1:11) {
    ChrSeg[[i]]<- readVCF(Datapath,1000,Chromosome,end_points[i], end_points[i+1])
  }

  #Calculating features for 11 windows

  df <- data.frame(matrix(nrow=1, ncol=0))
  ChrSeg[1:11]
  for (i in 1:11) {
    class<-ChrSeg[[i]]
    GENOME.class.slide <- sliding.window.transform(class,2000,1000, type=2)
    genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
      split <- strsplit(x," ")[[1]][c(1,3)]
      val <- mean(as.numeric(split))
      return(val)
    })
    # total length
    L<-length(GENOME.class.slide@region.names)
    #Diversity Statistics
    slide <- diversity.stats(GENOME.class.slide)

    #Calculation of Pi
    nucdiv <- slide@Pi
    nucdhapdiv <- slide@hap.diversity.within
    Hd<-sum(nucdhapdiv)
    Hd<-Hd/L
    pi_wini<-(0.0081*((Hd)**2))
    df <- cbind(df, pi_wini)
    colnames(df)[i] <- paste0("pi_win",i)
  }

  for (i in 1:11) {
    class<-ChrSeg[[i]]
    GENOME.class.slide <- sliding.window.transform(class,2000,1000, type=2)
    genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
      split <- strsplit(x," ")[[1]][c(1,3)]
      val <- mean(as.numeric(split))
      return(val)
    })
    # total length
    L<-length(GENOME.class.slide@region.names)
    #Diversity Statistics
    slide <- diversity.stats(GENOME.class.slide)

    #Calculation of Wattersons_theta
    head(nucdiv)
    slide<- neutrality.stats(GENOME.class.slide,detail=TRUE)
    get.neutrality(slide)
    theta_Watterson<-slide@theta_Watterson
    theta_Watterson<-na.omit(theta_Watterson)
    thetawat<-sum(theta_Watterson)
    thetaW_wini<-(thetawat/L)
    df <- cbind(df, thetaW_wini)
    colnames(df)[i+11] <- paste0("thetaW_win",i)
  }

  for (i in 1:11) {
    class<-ChrSeg[[i]]
    GENOME.class.slide <- sliding.window.transform(class,2000,1000, type=2)
    genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
      split <- strsplit(x," ")[[1]][c(1,3)]
      val <- mean(as.numeric(split))
      return(val)
    })
    # total length
    L<-length(GENOME.class.slide@region.names)
    #Diversity Statistics
    slide <- diversity.stats(GENOME.class.slide)

    #Calculation of Tajima's_D Statistic
    head(nucdiv)
    slide<- neutrality.stats(GENOME.class.slide,detail=TRUE)
    get.neutrality(slide)
    Tajima.D<-slide@Tajima.D
    ActTaj<-na.omit(Tajima.D)
    Tajima.DSum<-sum(ActTaj)
    tajD_wini<-(Tajima.DSum/L)
    df <- cbind(df, tajD_wini)
    colnames(df)[i+22] <- paste0("tajD_win",i)
  }

  for (i in 1:11) {
    class<-ChrSeg[[i]]
    GENOME.class.slide <- sliding.window.transform(class,2000,1000, type=2)
    genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
      split <- strsplit(x," ")[[1]][c(1,3)]
      val <- mean(as.numeric(split))
      return(val)
    })
    # total length
    L<-length(GENOME.class.slide@region.names)
    #Diversity Statistics
    slide <- diversity.stats(GENOME.class.slide)

    #Calculation of Kelly's_ZnS Statistic

    GENOME.class.slide<-linkage.stats(GENOME.class.slide)
    get.linkage(GENOME.class.slide)
    Kelly.Z_nS<-GENOME.class.slide@Kelly.Z_nS
    Kelly.Z_nS<-na.omit(Kelly.Z_nS)
    kelly.Z_nsSUM<-sum(Kelly.Z_nS)
    ZnS_wini<-(kelly.Z_nsSUM/L)
    df <- cbind(df, ZnS_wini)
    colnames(df)[i+33] <- paste0("ZnS_win",i)
  }

  for (i in 1:11) {
    class<-ChrSeg[[i]]
    GENOME.class.slide <- sliding.window.transform(class,2000,1000, type=2)
    genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
      split <- strsplit(x," ")[[1]][c(1,3)]
      val <- mean(as.numeric(split))
      return(val)
    })
    # total length
    L<-length(GENOME.class.slide@region.names)
    #Diversity Statistics
    slide <- diversity.stats(GENOME.class.slide)

    #Calculation of  Omega Statistic
    GENOME.class.slide<-linkage.stats(GENOME.class.slide)
    get.linkage(GENOME.class.slide)
    Kelly.Z_nS<-GENOME.class.slide@Kelly.Z_nS
    Kelly.Z_nS<-na.omit(Kelly.Z_nS)
    kelly.Z_nsSUM<-sum(Kelly.Z_nS)
    ZnS_wini<-(kelly.Z_nsSUM/L)
    k<-ZnS_wini
    Omega_wini<-(-((log(1-(k/2)))/2))
    Omega_wini
    df <- cbind(df, Omega_wini)
    colnames(df)[i+44] <- paste0("Omega_win",i)
  }
  #Fetching the Data Frame with features
  df
  #Fetching RDS File for Prediction
  utils::globalVariables("model")
  model<-readRDS("inst/exdata/rf.rds")

  #Predicting Selective Sweep Class
  Prediction <- predict(model, newdata = df)
  Prediction<-data.frame(Prediction)
  return(Prediction)
}

