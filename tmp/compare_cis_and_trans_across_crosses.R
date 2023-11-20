

cross_data$A$trans$hotspot_table$cross = "A"
cross_data$B$trans$hotspot_table$cross = "B"
cross_data$`3004`$trans$hotspot_table$cross = "3004"


A_hot = cross_data$A$trans$hotspot_bins_final #%>% filter(count > cross_data$A$trans$hotspot_threshold) %>% mutate(seqnames=chr,start=pos,width=1) %>% as_granges()
B_hot = cross_data$B$trans$hotspot_bins_final# %>% filter(count > cross_data$B$trans$hotspot_threshold) %>% mutate(seqnames=chr,start=pos,width=1) %>% as_granges()
C_hot = cross_data$`3004`$trans$hotspot_bins_final# %>% filter(count > cross_data$`3004`$trans$hotspot_threshold) %>% mutate(seqnames=chr,start=pos,width=1) %>% as_granges()

length(A_hot)
length(B_hot)
length(C_hot)
aa = A_hot %>% join_overlap_left(B_hot) %>% join_overlap_left(C_hot)
sum(is.na(aa$hotspot_pos.y)&is.na(aa$hotspot_pos))
aa = B_hot %>% join_overlap_left(A_hot) %>% join_overlap_left(C_hot)
sum(is.na(aa$hotspot_pos.y)&is.na(aa$hotspot_pos))
aa = C_hot %>% join_overlap_left(A_hot) %>% join_overlap_left(B_hot)
sum(is.na(aa$hotspot_pos.y)&is.na(aa$hotspot_pos))

aaa = c(A_hot,B_hot,C_hot) %>% reduce_ranges()
table(aaa$bin)

#A_hot$bin

#4
A_hot %>% join_overlap_inner(C_hot)
#2
B_hot %>% join_overlap_inner(C_hot)
#1

A_old_hot = cross_data$A_bulk$trans$hotspot_bins_final# %>% filter(count > cross_data$A_bulk$trans$hotspot_threshold) %>% mutate(seqnames=chr,start=pos,width=1) %>% as_granges()
a_old_hot_sum = cross_data$A_bulk$trans$hotspot_list %>% group_by(chrom, hotspot_pos) %>% summarise(chrom=chrom[1],bin=bin[1],pos=hotspot_pos[1])
a_hot_sum =  cross_data$A$trans$hotspot_list %>% group_by(chrom, hotspot_pos) %>% summarise(chrom=chrom[1],bin=bin[1],pos=hotspot_pos[1])


#cross_data$A$trans$hotspot_table$pos == cross_data$A_bulk$trans$hotspot_table

novel_hotspots_match = a[is.na(a$bin.y)]
a_hot_sum[a_hot_sum$bin %in%novel_hotspots_match$bin.x,]


a = A_hot %>% join_overlap_left(A_old_hot)sum(is.na(a$n.y))
cross_data$`3004`$trans$hotspot_list %>% group_by(chrom, hotspot_pos) %>% summarise(chrom=chrom[1],bin=bin[1],pos=hotspot_pos[1])


A_old


a_hot_sum %>% filter(pos %in% a_old_hot_sum$pos)


A_old_hot %>% join_overlap_left(A_hot)

# 5 /12
length(A_hot)
A_hot %>% join_overlap_inner(B_hot)
# 4/ 12
A_hot %>% join_overlap_inner(c(B_hot,C_hot))
# 2/12 
length(C_hot)
length(B_hot)
B_hot %>% join_overlap_inner(c(C_hot,A_hot))
# 1 /6 
B_hot %>% join_overlap_inner(c(A_hot,C_hot))
# 4 / 6

#B_hot

#length(C_hot)
table(cross_data$A_bulk$trans$hotspot_list$bin)
cross_data$A_bulk$trans$hotspot_list[cross_data$A_bulk$trans$hotspot_list$bin == "chrXIV:245104-490207",]

cross_data$A$trans$hotspot_table %>% mutate(cross_data$A$trans$combined_peaks)

join_overlap_inner(A_hot,A_old_hot)

library(regioneR)
lets_check_overlap = c()
join_overlap_inner(A_hot,C_hot)

aaa = (lapply(cross_data$A$trans$hotspot_enrichments_and_overlaps, function(x){x$eqtl_region}))
unlist(aaa[[1]])
#lapply(aaa,function(x){x[1]})
as.list(aaa)
str(aaa)
unlist(unlist(aaa))

A_hot = unlist(GRangesList(aaa))
A_old_hot = unlist(GRangesList(lapply(cross_data$A_bulk$trans$hotspot_enrichments_and_overlaps, function(x){x$eqtl_region})))

#A_old_hot = A_old_hot[-4]

A_hot$test ="h"
A_old_hot$test2= "h"
ab = A_hot %>% join_overlap_inner(A_old_hot)
lets_check_overlap = c()
for(i in 1:100){
  #print(i)
  aaa = randomizeRegions(A_hot,genome="sacCer3")
  lets_check_overlap = c(lets_check_overlap,(length(aaa %>% join_overlap_inner(A_old_hot))))
}
sum(lets_check_overlap > length(ab))

A_old_hot = unlist(GRangesList(lapply(cross_data$B$trans$hotspot_enrichments_and_overlaps, function(x){x$eqtl_region})))
ab = A_hot %>% join_overlap_inner(A_old_hot)
lets_check_overlap = c()
for(i in 1:100){
  #print(i)
  aaa = randomizeRegions(A_hot,genome="sacCer3")
  lets_check_overlap = c(lets_check_overlap,(length(aaa %>% join_overlap_inner(A_old_hot))))
}
sum(lets_check_overlap > length(ab))

A_old_hot = unlist(GRangesList(lapply(cross_data$`3004`$trans$hotspot_enrichments_and_overlaps, function(x){x$eqtl_region})))
A_old_hot$test2 = "test3"
ab = A_hot %>% join_overlap_inner(A_old_hot)

lets_check_overlap = c()
for(i in 1:100){
  print(i)
  aaa = randomizeRegions(A_hot,genome="sacCer3")
  lets_check_overlap = c(lets_check_overlap,(length(aaa %>% join_overlap_inner(A_old_hot))))
}
sum(lets_check_overlap > length(ab))

## ## 

A_hot
