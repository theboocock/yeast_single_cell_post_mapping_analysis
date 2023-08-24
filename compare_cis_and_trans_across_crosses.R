

cross_data$A$trans$hotspot_table$cross = "A"
cross_data$B$trans$hotspot_table$cross = "B"
cross_data$`3004`$trans$hotspot_table$cross = "3004"


A_hot = cross_data$A$trans$hotspot_bins_final #%>% filter(count > cross_data$A$trans$hotspot_threshold) %>% mutate(seqnames=chr,start=pos,width=1) %>% as_granges()
B_hot = cross_data$B$trans$hotspot_bins_final# %>% filter(count > cross_data$B$trans$hotspot_threshold) %>% mutate(seqnames=chr,start=pos,width=1) %>% as_granges()
C_hot = cross_data$`3004`$trans$hotspot_bins_final# %>% filter(count > cross_data$`3004`$trans$hotspot_threshold) %>% mutate(seqnames=chr,start=pos,width=1) %>% as_granges()



A_old_hot = cross_data$A_bulk$trans$hotspot_bins_final# %>% filter(count > cross_data$A_bulk$trans$hotspot_threshold) %>% mutate(seqnames=chr,start=pos,width=1) %>% as_granges()

a = A_hot %>% join_overlap_left(A_old_hot)
sum(!is.na(a$n.y))

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
