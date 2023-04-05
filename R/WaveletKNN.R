
#' @title Wavelet Based K-Nearest Neighbor Model
#' @param ts Time Series Data
#' @param MLag Maximum Lags
#' @param split_ratio Training and Testing Split
#' @param wlevels Number of Wavelet Levels
#' @import caret dplyr caretForecast Metrics tseries stats wavelets
#' @return
#' \itemize{
#'   \item Lag: Lags used in model
#'   \item Parameters: Parameters of the model
#'   \item Train_actual: Actual train series
#'   \item Test_actual: Actual test series
#'   \item Train_fitted: Fitted train series
#'   \item Test_predicted: Predicted test series
#'   \item Accuracy: RMSE and MAPE of the model
#' }
#'
#' @export
#'
#' @examples
#' library("WaveletKNN")
#' data<- rnorm(100,100, 10)
#' WG<-WaveletKNN(ts=data)
#' @references
#' \itemize{
#'\item Aminghafari, M. and Poggi, J.M. 2012. Nonstationary time series forecasting using wavelets and kernel smoothing. Communications in Statistics-Theory and Methods, 41(3),485-499.

#' \item Paul, R.K. A and Anjoy, P. 2018. Modeling fractionally integrated maximum temperature series in India in presence of structural break. Theory and Applied Climatology 134, 241–249.

#' }
WaveletKNN<-function(ts,MLag=12,split_ratio=0.8,wlevels=3){
  TSModel<-NULL
  SigLags<-NULL
  SigLags<-function(Data,MLag){
    ts<-as.ts(na.omit(Data))
    adf1<-adf.test(na.omit(ts))
    if (adf1$p.value>0.05){
      ts<-ts
    } else {
      ts<-diff(ts)
    }
    adf2<-adf.test(ts)
    if (adf2$p.value>0.05){
      ts<-ts
    } else {
      ts<-diff(ts)
    }

    CorrRes<-NULL
    for (i in 1:MLag) {
      # i=1
      ts_y<-dplyr::lag(as.vector(ts), i)
      t<-cor.test(ts,ts_y)
      corr_res<-cbind(Corr=t$statistic,p_value=t$p.value)
      CorrRes<-rbind(CorrRes,corr_res)
    }
    rownames(CorrRes)<-seq(1:MLag)
    Sig_lags<-rownames(subset(CorrRes,CorrRes[,2]<=0.05))
    maxlag<-max(as.numeric(Sig_lags))
    return(list(Result=as.data.frame(CorrRes),SigLags=as.numeric(Sig_lags),MaxSigLag=maxlag))
  }
  ntest<-round(length(ts)*(1-split_ratio), digits = 0)
  Split1 <- caretForecast::split_ts(as.ts(ts), test_size = ntest)
  train_data1 <- Split1$train
  test_data1 <- Split1$test
  Wvlevels<-wlevels
  mraout <- wavelets::modwt(as.vector(ts), filter="haar", n.levels=Wvlevels)
  WaveletSeries <- cbind(do.call(cbind,mraout@W),mraout@V[[Wvlevels]])
  MaxL<-NULL
  model_par<-NULL
  ts_fitted<-NULL
  ts_foreast<-NULL

  TSModel<-function(ts,MLag,split_ratio,method){
    ntest<-round(length(ts)*(1-split_ratio), digits = 0)
    Split <- caretForecast::split_ts(as.ts(ts), test_size = ntest)
    train_data <- Split$train
    test_data <- Split$test
    maxl<-SigLags(Data=train_data,MLag = MLag)$MaxSigLag
    model <- caretForecast::ARml(train_data, cv= TRUE, max_lag = maxl, caret_method = method,
                                 verbose = FALSE)
    model_par<-as.data.frame(model$model$bestTune)
    ts_fitted<-as.vector(model$fitted)
    ts_foreast<-as.vector(forecast(model, h=ntest)$mean)
    return(list(Model=model,Maxl=maxl,Param=model_par,Train_actual=train_data,Test_actual=test_data, Train_fitted=ts_fitted,Test_predicted=ts_foreast))
  }

  for (j in 1:ncol(WaveletSeries)) {
    w<-as.ts(WaveletSeries[,j])
    model<-TSModel(ts=w,MLag=MLag,split_ratio=split_ratio,method="knn")
    MaxL<- rbind(MaxL,model$Maxl)
    model_par<-rbind(model_par,model$Param)
    ts_fitted<-cbind(ts_fitted,model$Train_fitted)
    ts_foreast<-cbind(ts_foreast,model$Test_predicted)
  }
  rownames(MaxL)<- c(paste0("W",seq(1:Wvlevels)),"V")
  rownames(model_par)<-c(paste0("W",seq(1:Wvlevels)),"V")

  trainf <- apply(ts_fitted,1,sum)
  testf <- apply(ts_foreast,1,sum)

  RMSE<-c(Train=Metrics::rmse(train_data1[-c(1,MLag)],trainf[-c(1:MLag)]),Test=Metrics::rmse(test_data1,testf))
  MAPE<-c(Train=Metrics::mape(train_data1[-c(1,MLag)],trainf[-c(1:MLag)]),Test=Metrics::mape(test_data1,testf))
  accuracy<-rbind(RMSE,MAPE)
  return(list(Lag=MaxL,Parameters=model_par,Train_actual=train_data1,Test_actual=test_data1,Train_fitted=trainf,Test_predicted=testf, Accuracy=accuracy))
}

