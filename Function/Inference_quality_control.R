#Environment setup
library(ggplot2)

#Data input
path <- '../Output/'
sub_path <- c('0','L','H')

##########################################################################
#Collect info
for (m in 1:length(sub_path)) {
  fit_result_files <- dir(paste(path, sub_path[m], 'fit_result', sep = '/'))
  info <- array(0, dim = c(1,8))
  ##Collect in sequence
  for (i in 1:length(fit_result_files)) {
    file_info <- unlist(strsplit(fit_result_files[i], split = '_'))
    workbook <- paste(path, sub_path[m], 'fit_result', fit_result_files[i], sep = '/')
    fit_result <- read.csv(workbook)
    info <- rbind(info, cbind(file_info[1], file_info[2], fit_result$X, 
                              fit_result$RMSE_o, fit_result$RMSE_o/fit_result$Std_o, 
                              fit_result$Test.set.score, fit_result$ymax_test, fit_result$ymin_test))
  }
  
  ##Make input
  info <- info[-1,]
  colnames(info) <- c('target_otu', 'Theta', 'State', 'RMSE', 
                      'RMSE/STD', 'Test_score', 'max', 'min')
  info <- as.data.frame(info)
  for (i in 1:ncol(info)) {
    info[,i] <- as.numeric(as.character(info[,i]))
  }
  info$Theta <- as.character(info$Theta)
  
  ##Plot
  p <- 
  ggplot(info, aes(x=Theta, y=RMSE))+
    theme(panel.grid.major =element_line(size = 1),
          axis.title = element_text(size = 16, colour = "black"),
          axis.text = element_text(size = 14, colour = "black"),
          axis.line = element_line(linetype = "solid", size = 0.5),
          panel.background = element_rect(color = "black", size = 0.5),
          panel.grid.minor = element_blank()
    )+
    geom_boxplot(fill=c('#C52427', '#87AED6', '#7BAC53')[1])+
    scale_x_discrete(limits=c('0.1', '0.5', '1', '5', '10'))+
    ylim(0,2)
  
  print(p) # The plot can help you to determine the best theta to use. Commonly speaking, the smaller RMSE the better quality.
  
  p <- 
    ggplot(info, aes(x=Theta, y=`RMSE/STD`))+
    theme(panel.grid.major =element_line(size = 1),
          axis.title = element_text(size = 16, colour = "black"),
          axis.text = element_text(size = 14, colour = "black"),
          axis.line = element_line(linetype = "solid", size = 0.5),
          panel.background = element_rect(color = "black", size = 0.5),
          panel.grid.minor = element_blank()
    )+
    geom_boxplot(fill=c('#C52427', '#87AED6', '#7BAC53')[1])+
    scale_x_discrete(limits=c('0.1', '0.5', '1', '5', '10'))+
    ylim(0,2)+
    geom_hline(yintercept=1, linetype='dashed', color='red', size=1.5)
  
  print(p) # The plot can help you to determine the best theta to use. Commonly speaking, the RMSE/STD of a good inference should be lower than 1.
  
  ##Select best theta
  best_theta <- array(NA, dim = c(nrow(fit_result)*13, ncol(info)))
  colnames(best_theta) <- colnames(info)

  for (i in 0:12) {
    loc_otu <- which(info[,1] == i)
    focus <- info[loc_otu,]

    for (j in 1:nrow(fit_result)) {
      loc_state <- which(focus$State == (j-1))
      state_info <- focus[loc_state,]

      best_theta[i*nrow(fit_result)+j,] <- as.matrix(state_info[which.min(state_info[,4]),]) # state_info[,4] is used because we care about the RMSE of Jacobian inference.
    }
  }

  ##Collect the coefs via best theta
  coefs_file <- dir(paste(path, sub_path[m], 'coefs', sep = '/'))
  collected_coefs <- array(NA, dim = c(1, 13))

  for (i in 1:nrow(best_theta)) {
    # best_coefs <- paste(best_theta[i,1], 10, 'coefs.csv', sep = '_') # Simply use the Jacobian element at theta of 10.
    best_coefs <- paste(best_theta[i,1], best_theta[i,2], 'coefs.csv', sep = '_') # Use the Jacobian element which can give the minimum RMSE.
    workbook <- paste(path, sub_path[m], 'coefs', best_coefs, sep = '/')
    best_coefs <- read.csv(workbook, row.names = 1)

    colnames(collected_coefs) <- colnames(best_coefs)
    collected_coefs <- rbind(collected_coefs, best_coefs[as.numeric(best_theta[i,3])+1,])
  }
  collected_coefs <- collected_coefs[-1,]

  ##Output best coefs
  output_name <- paste(sub_path[m], 'best_coefs.csv', sep = '_')
  write.csv(collected_coefs,
            file = paste(path, output_name, sep = '/'))
}
