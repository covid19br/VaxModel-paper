res_corr <- c()
for(i in 1:5){
  for(j in 1:11){
    val <- c(results[i,j],results[i,j+11],results[i,j+22],results[i,j+33],results[i,j+44])
    print(val)
    res_corr <- rbind(res_corr,val)
  }
}