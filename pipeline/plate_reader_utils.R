## 


read_plate = function(plate, layout,random=F){
  plate_in = read.csv(plate, header=T)
  
  time = plate_in$Time
  # Only use the first 100 readings ... 
  #
  plate_data = plate_in[1:min(300,nrow(plate_in)),3:ncol(plate_in)]
  #layout = factor(layout, levels=layout)
  conditions = split(1:96,f=layout)
  conditions_data = list()
  i=1
  for(split_cond in conditions){
    conditions_data[[i]] = plate_data[,split_cond]
    
    i = i + 1
  }
  names(conditions_data) = names(conditions)
  conditions_data = lapply(conditions_data, function(x){if(class(x) == "numeric"){return(as.data.frame(x))}else{return(x)}})
  i = 1#rand
  for(split_cond in conditions){
    # print()
    colnames(conditions_data[[i]]) = colnames(plate_data)[split_cond]
    i =i +1
  }
  
  return(conditions_data)
}

read_layout = function(conditions, samples, random){
  cond_in = read.csv(conditions,header=F,stringsAsFactors = F)
  samples_in = read.csv(samples, header=F,stringsAsFactors = F,na.strings = "")
  mat_samples  = matrix(ncol=12,nrow=8)
  #  mat_samples[row,column] = 
  
  final_matrix = rep(NA,nrow(cond_in),ncol(cond_in))
  #print(cond_in)
  #print(samples_in)
  for(i in 1:nrow(cond_in)){
    for(j in 1:ncol(cond_in)){
      # print((i-1) * ncol(cond_in)  + j)
      final_matrix[(i - 1) * ncol(cond_in) + (j)] = paste(samples_in[i,j], cond_in[i,j],sep="_")
    }
  }
  #rint(final_matrix)
  #print(dim(final_matrix))
  final_matrix[grep("NA",final_matrix)] = NA
  #print(final_matrix)
  return(final_matrix)
}

convert_to_df_long_format = function(plate){
  aa = lapply(plate,function(x){ x$time = as.numeric(rownames(x)); melt(x,id.vars="time")})
  x = list()
  
  
  for(i in 1:length(aa)){
    x[[i]] = aa[[i]]
    x[[i]]$cond_sample = names(aa)[i]
    x[[i]]$sample = unlist(lapply(strsplit(x[[i]]$cond_sample,split = "_"),function(x){x[1]}))
    x[[i]]$cond = unlist(lapply(strsplit(x[[i]]$cond_sample,split = "_"),function(x){x[2]}))
    
    #x[[i]]$time = 1:nrow(x[[i]])
  }
  x2 = do.call("rbind",x)
  return(x2)
}

get_maximum_growth_rates = function(plate, lower, upper){
  plate_out_list = list()
  i = 1
  for(plate_in in plate){
    plate_out_list[[i]] = apply(plate_in,2, function(x){
      #print(x)
      x = x[x>lower & x <upper]
      #print(x)
      g = seq(length(x))
      #  print(x)
      summary_lm = (summary(lm(log(x)  ~    g)))
      return(summary_lm$coef[2,1])
    })
    i = i + 1
  }
  names(plate_out_list) = names(plate)
  return(plate_out_list)
}

get_geometric_mean_rates = function(plate, lower, upper){
  plate_out_list = list()
  i = 1
  for(plate_in in plate){
    plate_out_list[[i]] = apply(plate_in,2, function(x){
      value_one = x[which(x >= lower)[1]]
      value_two =x[which(x >= upper)[1]]
      
      return(log(value_two/value_one)/(abs(which(x >= lower)[1] - which(x >= upper)[1])))
      #return(summary_lm$coef[2,1])
    })
    i = i + 1
  }
  names(plate_out_list) = names(plate)
  return(plate_out_list)
}

get_geometric_mean_rates_spline = function(plate, lower, upper){
  plate_out_list = list()
  i = 1
  for(plate_in in plate){
    plate_out_list[[i]] = apply(plate_in,2, function(x){
      #print(x)
      if(sum(x > lower) == 0){
        return(0)
      }
      s_fun = splinefun(x,1:length(x))
      value_two = s_fun(upper)
      value_one = s_fun(lower)
      return(log(upper/lower)/(value_two-value_one))
      #return(summary_lm$coef[2,1])
    })
    i = i + 1
  }
  names(plate_out_list) = names(plate)
  return(plate_out_list)
}

process_phenotype_plates = function(plate, conditions, samples){
  
  #df_pheno = read.csv(mutant_pheno,stringsAsFactors = F)
  layout = read_layout(conditions, samples, random=random)
  plate = read_plate(plate,layout, random=random)
  max_growth_rates = get_geometric_mean_rates_spline(plate,lower=0.2,upper=.8) 
  switching = get_geometric_mean_rates_spline(plate, lower=0.8,upper=1.1)
  inflection = get_geometric_mean_rates_spline(plate, lower=1.1, upper=1.3)
  low = get_geometric_mean_rates_spline(plate,lower=0.15,upper=.3) 
  mid = get_geometric_mean_rates_spline(plate,lower=0.2,upper=.4) 
  mid2 = get_geometric_mean_rates_spline(plate,lower=0.4,upper=.6) 
  high = get_geometric_mean_rates_spline(plate,lower=0.6,upper=.8)
  df = stack(max_growth_rates)
  df_switch = stack(switching)
  df_inflect = stack(inflection)
  df_low = stack(low)
  df_mid = stack(mid)
  df_mid2 = stack(mid2)
  df_high = stack(high)
  df$low =  1/(log(2)/log(1+df_low$values)/4)
  df$mid =  1/(log(2)/log(1+df_mid$values)/4)
  df$high = 1/(log(2)/log(1+df_high$values)/4)
  df$mid2 = 1/(log(2)/log(1+df_mid2$values)/4)
  df$doubling = 1/(log(2)/log(1+df$values)/4)
  df$cond = unlist(lapply(strsplit(as.character(df$ind),"_"),function(x){x[2]}))
  df$sample = unlist(lapply(strsplit(as.character(df$ind),"_"),function(x){x[1]}))
  df$switching = df_switch$values
  df$switching = 1/(log(2)/log(1+df$switching)/4)
  df$inflection = df_inflect$values
  df$inflection = 1/(log(2)/log(1+df$inflection)/4)
  df$well = rownames(df)
  df_long = convert_to_df_long_format(plate)
  df_long$well = df_long$variable
  inflection_point_search = function(value){
    # print(value)
    first_idx = which(value > 0.2)[1]
    second_idx = which(value <= 1.2)[sum(value <=0.2)]
    #print(first_idx)
    #print(second_idx)
    first_delta_0.2 = abs(value[first_idx] - 0.2)
    second_delta_0.2 = abs(value[second_idx] -.2)
    
    avg = (second_delta_0.2 * first_idx + first_delta_0.2  * second_idx)/ (     (first_delta_0.2 + second_delta_0.2 ))
    
    return(avg)
  }
  inflection_point =df_long %>% group_by(well) %>% dplyr::summarize(inflection_point = inflection_point_search(value))
  df_long = merge(df_long, inflection_point,by="well")
  #df_long_m = df_long_m[order(df_long_m$time),]
  df_long$newtime = df_long$time - df_long$inflection_point + 10
  #df_long_m_s= split(df_long_m,df_long_m$well)
  return(list(df_long=df_long,df=df))
  #df_long$sample_id = as.numeric(unlist(lapply(strsplit(df_long$sample,split = "-"),function(x){x[[2]]})))
  
}
cor.cross <- function(x0, y0, i=0) {
  #
  # Sample autocorrelation at (integral) lag `i`:
  # Positive `i` compares future values of `x` to present values of `y`';
  # negative `i` compares past values of `x` to present values of `y`.
  #
  if (i < 0) {x<-y0; y<-x0; i<- -i}
  else {x<-x0; y<-y0}
  n <- length(x)
  cor(x[(i+1):n], y[1:(n-i)], use="complete.obs")
}
