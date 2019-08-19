
message("Preprocess TCRGP file...")

n_epitopes <- 12

tcrgp    <- tcrgp[,2:c(n_epitopes+3)]
tcrab    <- merge(tcrab, tcrgp, by.x = "Cell Name", by.y = "Cell name")

stage_df <- data.frame(Patient = c("P0205", "P0508", "P1116", "P0322", "P0407"),
                       Stage   = c("I/II", "I/II", "IVB", "IVB", "I/II"))
tcrab    <- merge(tcrab, stage_df, by = "Patient", all.x = T)
colnames(tcrab)  <- make.names(colnames(tcrab), unique=TRUE)


## Add names to df
gene_names <- make.unique(counts$symbol)
gene_names[is.na(gene_names)] <- "empty"
row.names(counts) <- gene_names

## Take only cells that have TCRab
counts    <- counts[ ,colnames(counts) %in% tcrab$Cell.Name]
tcrab     <- tcrab[tcrab$Cell.Name %in% colnames(counts), ]
rownames(tcrab) <- tcrab$Cell.Name

## Add metadata about the predictions with TCRGP
tcrab$hbv_reactive <- "FALSE"
tcrab$hbv_reactive <- ifelse(tcrab$pol282_b > 0 | tcrab$pol387_b > 0 |  tcrab$core169_b > 0 |  tcrab$core195_b > 0, TRUE, FALSE)
tcrab$hbv_reactive[is.na(tcrab$hbv_reactive)] <- FALSE

tcrab$epitope_spec <- "None"
tcrab$epitope_spec[tcrab$pol282_b  > 0] = "pol282_b"
tcrab$epitope_spec[tcrab$pol387_b  > 0] = "pol387_b"
tcrab$epitope_spec[tcrab$core169_b > 0] = "core169_b"
tcrab$epitope_spec[tcrab$core195_b > 0] = "core195_b"

tcrab$viral_reactive <- "FALSE"
tcrab$viral_reactive <- ifelse(tcrab$IPSINVHHY_trb > 0 | tcrab$TPRVTGGGAM_trb > 0 |  tcrab$NLVPMVATV_trb > 0 |  tcrab$GLCTLVAML_trb > 0  |  tcrab$RAKFKQLL_trb > 0  |  tcrab$YVLDHLIVV_trb > 0  |  tcrab$GILGFVFTL_trb > 0 |  tcrab$PKYVKQNTLKLAT_trb > 0, TRUE, FALSE)
tcrab$viral_reactive[is.na(tcrab$viral_reactive)] <- FALSE

tcrab$viral_epitope_spec <- "None"
tcrab$viral_epitope_spec[tcrab$IPSINVHHY_trb  > 0]    = "IPSINVHHY"
tcrab$viral_epitope_spec[tcrab$TPRVTGGGAM_trb  > 0]   = "TPRVTGGGAM"
tcrab$viral_epitope_spec[tcrab$NLVPMVATV_trb > 0]     = "NLVPMVATV"
tcrab$viral_epitope_spec[tcrab$GLCTLVAML_trb > 0]     = "GLCTLVAML"
tcrab$viral_epitope_spec[tcrab$RAKFKQLL_trb > 0]      = "RAKFKQLL"
tcrab$viral_epitope_spec[tcrab$YVLDHLIVV_trb > 0]     = "YVLDHLIVV"
tcrab$viral_epitope_spec[tcrab$GILGFVFTL_trb > 0]     = "GILGFVFTL"
tcrab$viral_epitope_spec[tcrab$PKYVKQNTLKLAT_trb > 0] = "PKYVKQNTLKLAT"

tcrab$phenotype    <- ifelse(tcrab$T.cell.Cluster %in% c("C1_CD8-LEF1", "C2_CD8-CX3CR1", "C4_CD8-LAYN", "C5_CD8-GZMK"), "CD8", "CD4")
tcrab$phenotype[tcrab$phenotype == "C3_CD8-SLC4A10"] <- "Else"

