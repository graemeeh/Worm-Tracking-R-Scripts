# The below script makes plots and does t tests on outputted data from the tracker. 
# this is for calculating chemotaxis index over time (procedure from Margie et al. 2013) and doing stats on that. 
# the marked circles outputted by the tracker are used to figure out quadrants. so make sure plates are aligned on tracker (i.e. lines 
# are vertical or horizontal). 

library(ggplot2)
library(dplyr)
library(tidyr)
library(pracma)
library(purrr)
library(stringr)

folder = "./GE Octanol quadrants - Copy"
timelapse_interval_s = 10
controls = "tgMOR"
framesize = c(2456, 2052)

#these should be appropriately named .csv files containing worm centroid positions over time
position_files = list.files(path = folder, pattern = "*.avi_worms.csv", ignore.case = T, recursive=T, full.names = T) 
position_data = lapply(position_files, read.csv)
names(position_data) = position_files

#these should be appropriately named .csv files containing information about the label sizes/shapes
#so that the program can figure out where the quadrants are. 
label_files = list.files(path = folder,pattern = "*.avi_details.csv",ignore.case = T, recursive=T, full.names = T)
label_data = lapply(label_files, read.csv, blank.lines.skip = F)
names(label_data) = label_files

min_dist = function(x1, y1, lx2, ly2){
  d2 = c()
  for (i in 1:length(lx2)){
    d2 = append(d2, sqrt((x1 - lx2[i])**2 + (y1 - ly2[i])**2))}
  return(which.min(d2))
}

strainrenamer = function(st){
  for (i in 1:length(st)){
    if (is.na(st[i])){
      st[i] = NA
    }
    else if (st[[i]] == "XMN1408"){
      st[i] = "tgMOR"
    }
    else if (st[i] == "VG1062"){
      st[i] = "dop-3"
    }
    else if (st[i] == "VG1049"){
      st[i] = "cat-2"
    }
    else if (st[i] == "VG1038"){
      st[i] = "dat-1"
    }
    else if (st[i] == "VG998"){
      st[i] = "comt-4"
    }
    else if (st[i] == "VG997"){
      st[i] = "faah-1"
    }
    else if (st[i] == "VG996"){
      st[i] = "irk-2"
    }
    else if (st[i] == "N2"){
      st[i] = "WT"
    }
  }
  return (st)
}

reframe = function(l1){
  new_frames = c()
  frame_c = l1[1]
  for (i in 1:length(l1)){
    if (i == 1) {
      new_frames = append(new_frames, frame_c)
    }
    else if (l1[i] != l1[i-1]){
      frame_c = frame_c + 1
      new_frames = append(new_frames, frame_c)
    }
    else {
      new_frames = append(new_frames, frame_c)
    }
  }
  return(new_frames)
}

exclude_circle = function(x, y, circ_x, circ_y, rad){
  return (!((x-circ_x)^2 + (y-circ_y)^2 <= rad^2))
}

newdf = data.frame(time=numeric(), 
                   chemotaxis_index = numeric(), 
                   n = numeric(), 
                   strain = character(), 
                   stringsAsFactors = F)

for (i in position_files) {
  strain = strsplit(strsplit(strsplit(i, "/")[[1]][4], "_")[[1]][1], ".", fixed = TRUE)[[1]][1]
  print(strain)
  position = position_data[[i]]
  position$frame = reframe(position$frame)
  labels_file= gsub("worms", "details", i)
  labels = label_data[[labels_file]]
  labels = labels[labels[,"type"] == "circle"]
  center_circle = labels %>% slice(min_dist(framesize[1]/2, framesize[2]/2, labels[,'x1'], labels[,'y1']))
  
  # set cutoff lines for x and y
  x_line = center_circle[,"x1"][1]
  y_line = center_circle[,"y1"][1]
  
  # filter out center circle worm
  position$centre = exclude_circle(position$centroid_x, position$centroid_y, x_line, y_line, 200)
  position = position[position$centre == T,]
  # populate df_new with T or F, based on below conditions
  position$oct_quadrant = (ifelse(
    #bottom right corner 
    position[,"centroid_x"] > x_line & position[,"centroid_y"] > y_line, T,
    #bottom left corner
    ifelse(position[,"centroid_x"] < x_line & position[,"centroid_y"] > y_line, F,
           #top right corner
           ifelse(position[,"centroid_x"] > x_line & position[,"centroid_y"] < y_line, F,
                  #top left corner
                  ifelse(position[,"centroid_x"] < x_line & position[,"centroid_y"] < y_line, T, F)))))
  
  octanol = position[position$oct_quadrant == T,]
  ethanol = position[position$oct_quadrant == F,]
  df_wormsums = data.frame(frame = unique(position$frame), n = unique(position$frame), octanol = unique(position$frame, ethanol = unique(position$frame)))
  for (j in 1:length(df_wormsums$frame)){
    df_wormsums$n[j] = nrow(position[position$frame == j-1,])[[1]]
    df_wormsums$octanol[j] = nrow(octanol[octanol$frame == j-1,])[[1]]
    df_wormsums$ethanol[j] = nrow(ethanol[ethanol$frame == j-1,])[[1]]
  }
  df_wormsums$octanol = as.integer(df_wormsums$octanol)
  df_wormsums$ethanol = as.integer(df_wormsums$ethanol)
  df_wormsums$n = as.integer(df_wormsums$n)
  df_wormsums$time = df_wormsums$frame * timelapse_interval_s / 60
  df_wormsums$chemotaxis_index = (df_wormsums$octanol - df_wormsums$ethanol)/(df_wormsums$n)
  df_wormsums$file = rep(i, length(df_wormsums$frame))
  df_wormsums$strain = rep(strain, times = length(df_wormsums$frame), length.out = NA, each = 1)
  newdf = dplyr::bind_rows(df_wormsums, newdf)
  
}

newdf$strain = strainrenamer(newdf$strain)
# THis bit is to make both strains have the same length ie if one has worms appear earlier the timepoints from earlier are removed 
# necessary for AUC
print(unique(newdf$strain))
strain_cleanup = function(df){
  non_nan_newdf = na.omit(df)
  
  strainmins = c()
  for (i in unique(non_nan_newdf$file)){
    strainmins = append(strainmins, min(non_nan_newdf[non_nan_newdf$file == i,]$time))
  }
  strainmaxes = c()
  for (i in unique(non_nan_newdf$file)){
    strainmaxes = append(strainmaxes, max(non_nan_newdf[non_nan_newdf$file == i,]$time))
  }
  
  non_nan_newdf = non_nan_newdf[non_nan_newdf$time>max(strainmins, na.rm = T),]
  non_nan_newdf = non_nan_newdf[non_nan_newdf$time<min(strainmaxes, na.rm = T),]
  return (non_nan_newdf)
}



for (i in unique(newdf$strain)){
  if (i != controls){
    df = newdf[newdf$strain == i | newdf$strain == controls,]
    if (length(unique(df$strain)) > 1){
      
      
      folders = unique(sapply(df[df$strain == i,]$file, FUN = function(x) {strsplit(x, "/")[[1]][3]}))
      df <- df[sapply(df$file, FUN = function(x) { strsplit(x, "/")[[1]][3] }) %in% folders, ]
      non_nan_newdf = strain_cleanup(df)
      df <- non_nan_newdf %>%
        group_by(strain, time) %>%
        filter(sum(!is.na(chemotaxis_index)) > 1) %>%
        summarise(
          c_sd = sd(chemotaxis_index, na.rm = T),
          c_mean = mean(chemotaxis_index, na.rm = T),
          c_confint = mean(chemotaxis_index, na.rm = T) - t.test(chemotaxis_index)$conf.int[1]
        )
      df = df %>% arrange(factor(strain, levels = c('tgMOR', i))) %>% mutate(strain=factor(strain, levels = c('tgMOR', i)))
      length(df[df$strain == 'irk-2',]$strain)
      pi = ggplot(data = df, aes(x = time, y = c_mean, colour = strain, ymin = c_mean-c_confint, ymax = c_mean+c_confint)) + 
        geom_line() + 
        theme_classic() + 
        scale_color_brewer(palette = "Set2") + 
        geom_errorbar(alpha = 0.3) + 
        xlab("Time (min)") + 
        ylab("Chemotaxis Index")
      
      file_name <- paste0("my_plot_", i, ".png")
      ggsave(
        filename = file_name,
        plot = pi, 
        width = 5,
        height = 3.5, 
        dpi = 600)
    }
    last10mins = non_nan_newdf[non_nan_newdf$time > 2,]
    
    df.summary_auc <- non_nan_newdf %>%
      group_by(strain, file) %>%
      summarise(
        c_auc = trapz(time, chemotaxis_index)
      )
    df.summary_last10 <- last10mins %>%
      group_by(strain, file) %>%
      summarise(
        c_auc = trapz(time, chemotaxis_index)
      )
    
    print(i)
    print("full analysis, bonferroni corrected:")
    print(t.test(c_auc ~ strain, data = df.summary_auc)$p.value * (length(unique(newdf$strain))-2))
    print("last 10 mins, bonferroni corrected:")
    print(t.test(c_auc ~ strain, data = df.summary_last10)$p.value * (length(unique(newdf$strain))-2))
    
    }

  
  
}








