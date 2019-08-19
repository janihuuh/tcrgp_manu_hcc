
## Clonality assessment

message("Analyze clonality")

a <- getGini(viz_df, covariate = "hbv_reactive") %>% mutate(sample.type = "all")
b <- getGini(viz_df, covariate = "epitope_spec") %>% mutate(sample.type = "all")
c <- getGini(viz_df, covariate = "Patient")      %>% mutate(sample.type = "all")
d <- getGini(viz_df, covariate = "Sample.Type")  %>% mutate(sample.type = "all")

e <- getGini(subset(viz_df, Sample.Type == "TTC"), covariate = "hbv_reactive") %>% mutate(sample.type = "TTC")
f <- getGini(subset(viz_df, Sample.Type == "TTC"), covariate = "epitope_spec") %>% mutate(sample.type = "TTC")

g <- getGini(subset(viz_df, Sample.Type == "NTC"), covariate = "hbv_reactive") %>% mutate(sample.type = "NTC")
h <- getGini(subset(viz_df, Sample.Type == "NTC"), covariate = "epitope_spec") %>% mutate(sample.type = "NTC")

i <- getGini(subset(viz_df, Stage == "I/II"), covariate = "epitope_spec") %>% mutate(sample.type = "I/II")
j <- getGini(subset(viz_df, Stage == "IVB"), covariate = "epitope_spec")  %>% mutate(sample.type = "IVB")

k <- getGini(subset(viz_df, Sample.Type == "TTC" & Stage == "I/II"), covariate = "epitope_spec") %>% mutate(sample.type = "I/II tumor")
l <- getGini(subset(viz_df, Sample.Type == "TTC" & Stage == "IVB"), covariate = "epitope_spec")  %>% mutate(sample.type = "IVB tumor")


gini_df <- rbind(a,b,c,d,e,f,g,h,i,j,k,l)
write.table(gini_df, paste0(output_dir, "gini.txt"), sep = "\t", quote = F, row.names = F)
