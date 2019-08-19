
messsage("Calculate Fisher's exact...")

## Create a df that will be used for plotting as well
# viz_df  <- data.frame(umap_df, "cluster" = Idents(hcc_cd8), hcc_cd8@meta.data)
# viz_df  <- merge(viz_df, stage_df, by = "Patient", all.x = T)
# colnames(viz_df)  <- make.names(colnames(viz_df), unique=TRUE)
# write.table(viz_df, paste0(output_dir, "viz_df.txt"), sep = "\t", quote = F, row.names = F)



## Run Fisher's
fisherGreater  <- RunFisherEpitope(viz_df, alternative = 'greater')  %>% ModifyFisherOutput
fisherLess     <- RunFisherEpitope(viz_df, alternative = 'less') %>% ModifyFisherOutput
fisherTwosided <- RunFisherEpitope(viz_df, alternative = 'two.sided') %>% ModifyFisherOutput



## For manuscript, split for smaller groups
fisherGreaterAllPatients        <- fisherGreater %>% filter(patient_names == "P0205, P0508, P1116, P0322, P0407") %>% ModifyFisherOutput
fisherGreaterPerPatients        <- fisherGreater %>% filter(patient_names != "P0205, P0508, P1116, P0322, P0407") %>% ModifyFisherOutput

fisherLessAllPatients           <- fisherLess    %>% filter(patient_names == "P0205, P0508, P1116, P0322, P0407") %>% ModifyFisherOutput
fisherLessPerPatients           <- fisherLess    %>% filter(patient_names != "P0205, P0508, P1116, P0322, P0407") %>% ModifyFisherOutput

fisherGreaterAllPatientsTumor   <- fisherGreater %>% filter(patient_names == "P0205, P0508, P1116, P0322, P0407" & tissue_names == "TTC") %>% ModifyFisherOutput
fisherGreaterAllPatientsBlood   <- fisherGreater %>% filter(patient_names == "P0205, P0508, P1116, P0322, P0407" & tissue_names == "PTC") %>% ModifyFisherOutput

fisherGreaterEarlyPatients      <- fisherGreater %>% filter(patient_names == "P0205, P0508, P0407" & stage_names == "I/II") %>% ModifyFisherOutput
fisherGreaterEarlyPatientsTumor <- fisherGreater %>% filter(patient_names == "P0205, P0508, P0407" & stage_names == "I/II" & tissue_names == "TTC") %>% ModifyFisherOutput

fisherGreaterLatePatients       <- fisherGreater %>% filter(patient_names == "P1116, P0322" & stage_names == "IVB")  %>% ModifyFisherOutput
fisherGreaterLatePatientsTumor  <- fisherGreater %>% filter(patient_names == "P1116, P0322" & stage_names == "IVB" & tissue_names == "TTC")  %>% ModifyFisherOutput


## Write results
dir.create(paste0(output_dir, "fisher/"), showWarnings = F)

write.table(fisherGreater,  paste0(output_dir, "fisher/fisher_greater.txt"), sep = "\t", quote = F, row.names = F)
write.table(fisherLess,     paste0(output_dir, "fisher/fisher_less.txt"), sep = "\t", quote = F, row.names = F)
write.table(fisherTwosided, paste0(output_dir, "fisher/fisher_twosided.txt"), sep = "\t", quote = F, row.names = F)

write.table(fisherGreaterAllPatients,  paste0(output_dir, "fisher/fisherGreaterAllPatients.txt"), sep = "\t", quote = F, row.names = F)
write.table(fisherGreaterPerPatients,  paste0(output_dir, "fisher/fisherGreaterPerPatients.txt"), sep = "\t", quote = F, row.names = F)

write.table(fisherLessAllPatients,  paste0(output_dir, "fisher/fisherLessAllPatients.txt"), sep = "\t", quote = F, row.names = F)
write.table(fisherLessPerPatients,  paste0(output_dir, "fisher/fisherLessPerPatients.txt"), sep = "\t", quote = F, row.names = F)

write.table(fisherGreaterAllPatientsTumor,  paste0(output_dir, "fisher/fisherGreaterAllPatientsTumor.txt"), sep = "\t", quote = F, row.names = F)
write.table(fisherGreaterAllPatientsBlood,  paste0(output_dir, "fisher/fisherGreaterAllPatientsBlood.txt"), sep = "\t", quote = F, row.names = F)

write.table(fisherGreaterEarlyPatients,  paste0(output_dir, "fisher/fisherGreaterEarlyPatients.txt"), sep = "\t", quote = F, row.names = F)
write.table(fisherGreaterLatePatients,  paste0(output_dir, "fisher/fisherGreaterLatePatients.txt"), sep = "\t", quote = F, row.names = F)

write.table(fisherGreaterEarlyPatientsTumor,  paste0(output_dir, "fisher/fisherGreaterEarlyPatientsTumor.txt"), sep = "\t", quote = F, row.names = F)
write.table(fisherGreaterLatePatientsTumor,  paste0(output_dir, "fisher/fisherGreaterLatePatientsTumor.txt"), sep = "\t", quote = F, row.names = F)
