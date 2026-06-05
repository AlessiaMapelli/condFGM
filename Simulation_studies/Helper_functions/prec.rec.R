# Function to compute precision, recall, TPR, FPR, and F1 score for binary adjacency matrices.

# @param G.true  p x p binary adjacency matrix (ground truth)
# @param G.mat   p x p binary adjacency matrix (estimated)
# @param type    "AND" (edge requires mutual recognition) or "OR" (either direction suffices)
# @return        list with prec (precision), rec (recall), TPR, FPR, F1
prec.rec <- function(G.true, G.mat, type=c("AND","OR")){
  
  p <- nrow(G.true)
  TP <- 0; TN <- 0; FP <- 0; FN <- 0
  
  if(type=="AND"){
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        if(G.mat[i,j] == 1 & G.mat[j,i] == 1){
          if(G.true[i,j] == 1)
            TP <- TP + 1
          else
            FP <- FP + 1
        }else{
          if(G.true[i,j] == 1)
            FN <- FN + 1
          else
            TN <- TN + 1
        }
      }
    }
  }else if(type=="OR"){
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        if(G.mat[i,j] == 1 | G.mat[j,i] == 1){
          if(G.true[i,j] == 1)
            TP <- TP + 1
          else
            FP <- FP + 1
        }else{
          if(G.true[i,j] == 1)
            FN <- FN + 1
          else
            TN <- TN + 1
        }
      }
    }
  }
  if(TP+FP==0){prec <- 0
  }else{prec <- TP / (TP + FP)}

  if(TP+FN==0){ rec <- 0
  }else{rec <- TP / (TP + FN)}

  if(sum(G.mat) ==0 & sum(G.true)==0 ){
    TPR <- 0
    FPR <- 0
    F1 <- 1
  } else{
   TPR <- TP / (TP + FN)
   FPR <- FP / (FP + TN)
   F1 <- 2*TP/(2*TP+FP+FN)
  }
  return(list(prec=prec, rec=rec, TPR=TPR, FPR=FPR, F1=F1))

}
