# ---- 0. Load Required Libraries ----
required_packages <- c(
  "readxl", "stringr", "tidyr", "httr", "jsonlite", "purrr",
  "tibble", "ggplot2", "forcats", "dplyr", "ggpubr", "BiocManager"
)

# Install and load required libraries
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}
install_and_load(required_packages)

# ---- 1. Load and Clean Data ----
clean_data <- function(filepath, sheet_name) {
  sarco <- read_excel(filepath, sheet = sheet_name) %>%
    dplyr::select(where(~ any(!is.na(.)))) %>%
    filter(if_any(everything(), ~ !is.na(.))) %>%
    mutate(
      Gene_Symbol = str_remove(str_extract(Gene_protein, "^[^_;]+"), ";?GN=+$"),
      Accession   = str_extract(Gene_protein, "\\b[A-Z][0-9][A-Z0-9]{3}[0-9]\\b")
    )
  
  sample_cols <- grep("^(EAA|PRE|POST|PPS)", names(sarco), value = TRUE)
  sarco <- sarco %>% mutate(detection_score = rowSums(across(all_of(sample_cols)), na.rm = TRUE))
  
  list(sarco = sarco, sample_cols = sample_cols)
}

data <- clean_data("13_MASTER_SEER proteome analysis (4-4-23).xlsx", "MASTER_SARCO_data")
sarco <- data$sarco
sample_cols <- data$sample_cols

# ---- 2. Query UniProt for Missing Gene Symbols ----
query_uniprot <- function(accessions) {
  url <- "https://rest.uniprot.org/uniprotkb/search"
  q <- paste(accessions, collapse = " OR ")
  r <- GET(url, query = list(query = q, fields = "accession,gene_primary", format = "json", size = 500))
  if (status_code(r) == 200) {
    j <- httr::content(r, as = "parsed", type = "application/json")
    tibble(
      Accession = sapply(j$results, function(x) x$primaryAccession),
      Gene_Symbol_API = sapply(j$results, function(x)
        if (!is.null(x$genes[[1]]$geneName$value)) x$genes[[1]]$geneName$value else NA_character_)
    )
  } else stop("UniProt query failed")
}

# Process UniProt queries
accessions_to_query <- sarco %>%
  group_by(Accession, Gene_Symbol) %>%
  slice_max(detection_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(is.na(Gene_Symbol) | Gene_Symbol == "") %>%
  pull(Accession) %>% unique() %>% na.omit()

uniprot_results <- accessions_to_query %>%
  split(ceiling(seq_along(.) / 100)) %>%
  map_dfr(query_uniprot) %>%
  distinct(Accession, .keep_all = TRUE)

# ---- 3. Merge and Clean Gene Symbols ----
sarco <- sarco %>%
  left_join(uniprot_results, by = "Accession") %>%
  mutate(
    Gene_Symbol = if_else(is.na(Gene_Symbol) | Gene_Symbol == "", Gene_Symbol_API, Gene_Symbol),
    Gene_Symbol = str_remove(Gene_Symbol, "^isoform_")
  ) %>%
  select(-Gene_Symbol_API)

# ---- 4. Remove Duplicates and Reshape Data ----
sarco_clean <- sarco %>%
  group_by(Accession, Gene_Symbol) %>%
  slice_max(detection_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(Gene_Symbol) %>%
  slice_max(detection_score, n = 1, with_ties = FALSE) %>%
  ungroup()

long_df <- sarco_clean %>%
  select(-matches("^T-test_|^DELTA_|^CORREL_|^Young_|^detection_")) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  filter(!str_detect(Sample, "_Avg|_StdDev")) %>%
  mutate(
    Timepoint = case_when(
      str_detect(Sample, "PRE") ~ "Pre",
      str_detect(Sample, "POST") ~ "Post",
      TRUE ~ NA_character_
    ),
    Group = case_when(
      str_detect(Sample, "EAA|PPS") ~ "Young",
      TRUE ~ "MA"
    )
  )

# Save cleaned data
write.csv(long_df, "Clean_Long2.csv", row.names = FALSE)
write.csv(sarco_clean, "Clean_Wide2.csv", row.names = FALSE)

# ---- 5. Define Categories ----
categories <- list(
  Oxidative_Damage = list(
    ROS_Generation = c("NOX1", "NOX2", "NOX3", "NOX4", "DUOX1", "DUOX2", "CYBA", "CYBB", "XDH", "AIFM1", "AIFM2", "MPO", "NOS2"),
    Damage_Products = c("HMOX1", "HMOX2", "ALDH2", "ALOX5", "ALOX12", "ALOX15"),
    Redox_Signaling = c("NFE2L2", "KEAP1", "NQO1", "BACH1", "MAFG", "CUL3", "SQSTM1")
  ),
  Antioxidant_Response = list(
    Superoxide_Dismutases = c("SOD1", "SOD2", "SOD3"),
    Catalases = c("CAT"),
    Glutathione_Peroxidases = c("GPX1", "GPX2", "GPX3", "GPX4"),
    Peroxiredoxins = c("PRDX1", "PRDX2", "PRDX3", "PRDX5", "PRDX6"),
    Thiol_System = c("TXN", "TXNRD1", "TXNRD2", "GSR", "GLRX", "GLRX2"),
    Glutathione_Synthesis = c("GCLC", "GCLM", "GSS"),
    Glutathione_Transferases = c("GSTA1", "GSTP1", "MGST1", "GSTM1"),
    NADPH_Regenerators = c("G6PD", "PGD", "IDH2", "NNT")
  ),
  ER_Stress_UPR = list(
    Sensors = c("HSPA5", "EIF2AK3", "ERN1", "ATF6"),
    Pro_apoptotic_Effectors = c("DDIT3", "BAX", "BAK1", "PMAIP1", "BBC3"),
    Adaptive_Effectors = c("XBP1s", "ATF4", "PPP1R15A", "HERPUD1", "CHAC1", "ASNS"),
    Chaperones = c("HSP90B1", "PDIA3", "PDIA4", "PDIA6", "CALR", "CANX", "ERP44", "UGGT1", "DNAJB9", "DNAJC3")
  ),
  Protein_Synthesis = list(
    mTORC1_Components = c("MTOR", "RPTOR", "MLST8", "PRAS40", "DEPTOR"),
    Upstream_Regulators = c("AKT1", "AKT2", "TSC1", "TSC2", "RHEB"),
    Downstream_Readouts = c("RPS6KB1", "RPS6KB2", "EIF4EBP1", "EIF4EBP2", "EIF4E", "EIF4G1", "ULK1", "LARP1")
  ),
  Proteostasis = list(
    Macroautophagy = c("BECN1", "ATG3", "ATG5", "ATG7", "ATG12", "ATG16L1", "MAP1LC3A", "MAP1LC3B", "SQSTM1", "WIPI2", "ATG14", "UVRAG", "LAMP2"),
    Chaperone_Mediated_Autophagy = c("LAMP2A", "HSPA8", "HSP90AA1"),
    Ubiquitin_Proteasome = c("PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7", "PSMC1", "PSMC2", "PSMC3", "PSMC4", "PSMC5", "PSMC6", "PSMD1", "PSMD2", "PSMD3", "PSMD4", "UBB", "UBC", "UBA1", "UBE2D1", "UBE2D2", "UBE2L3", "CUL1", "SKP1", "RBX1"),
    Calpain_System = c("CAPN1", "CAPN2", "CAPNS1", "CAST"),
    Apoptotic_Caspases = c("CASP3", "CASP7", "CASP8", "CASP9", "APAF1")
  ),
  Mitochondrial_Homeostasis = list(
    Fusion = c("MFN1", "MFN2", "OPA1", "OPA3", "MIRO1", "MIRO2"),
    Fission = c("DNM1L", "FIS1", "MFF", "MIEF1", "MIEF2", "MTP18"),
    Mitophagy = c("PINK1", "PRKN", "BNIP3", "BNIP3L", "FUNDC1", "BCL2L13", "CALCOCO2", "OPTN", "TBK1"),
    Biogenesis = c("PPARGC1A", "PPARGC1B", "NRF1", "NRF2", "ESRRA", "TFAM", "TFB1M", "TFB2M", "POLG")
  )
)

# ---- 6. Helper Functions ----
calculate_pathway_data <- function(category, data) {
  abundance_matrix <- sapply(groups, function(group) {
    subset_data <- data %>%
      filter(
        (group == "Young" & Group == "Young") |
          (group == "Pre" & Group == "MA" & Timepoint == "Pre") |
          (group == "Post" & Group == "MA" & Timepoint == "Post")
      )
    sapply(category, function(genes) mean(subset_data$Abundance[subset_data$Gene_Symbol %in% genes], na.rm = TRUE))
  })
  list(
    normalized = abundance_matrix / abundance_matrix[, "Young"],
    raw = abundance_matrix
  )
}

plot_bar <- function(data, pathway_name) {
  bar_data <- as.data.frame(t(data$raw)) %>%
    rownames_to_column("Group") %>%
    pivot_longer(-Group, names_to = "Module", values_to = "Abundance") %>%
    mutate(Group = factor(Group, levels = group_order))
  
  ggplot(bar_data, aes(x = Group, y = Abundance, fill = Group)) +
    geom_col(position = "dodge", color = "black") +
    facet_wrap(~ Module, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = colors) +
    labs(
      title = paste(pathway_name, "Pathway Abundance Bar Plot"),
      x = "Group",
      y = "Mean Abundance",
      fill = "Group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

# ---- 7. Calculate Significance Values ----
calculate_significance <- function(categories, data) {
  purrr::imap_dfr(categories, function(genes, module) {
    df <- dplyr::filter(data, Gene_Symbol %in% genes)
  # 1. Pre vs Post within MA
    df_pp <- dplyr::filter(df, Group == "MA", Timepoint %in% c("Pre","Post"))
    p_pp <- if (length(unique(df_pp$Timepoint))==2) {
      t.test(Abundance ~ Timepoint, df_pp)$p.value
    } else NA_real_
  # 2. Young vs Pre
    df_yp <- dplyr::filter(df, Group == "Young" | (Group=="MA" & Timepoint=="Pre"))
    p_yp <- if (length(unique(df_yp$Group))==2) {
      t.test(Abundance ~ Group, df_yp)$p.value
    } else NA_real_
  # 3. Young vs Post
    df_ypo <- dplyr::filter(df, Group == "Young" | (Group=="MA" & Timepoint=="Post"))
    p_ypo <- if (length(unique(df_ypo$Group))==2) {
      t.test(Abundance ~ Group, df_ypo)$p.value
    } else NA_real_
    
    tibble::tibble(
      Pathway_Module   = module,
      p_Pre_vs_Post    = p_pp,
      p_Young_vs_Pre   = p_yp,
      p_Young_vs_Post  = p_ypo
    )
  }) %>%
  # Benjamini–Hochberg FDR correction column of p-values
  dplyr::mutate(dplyr::across(dplyr::starts_with("p_"),
                              ~ p.adjust(., method = "BH"),
                              .names = "{.col}_adj"))
}
# Apply to each major pathway
signif_oxidative  <- calculate_significance(categories$Oxidative_Damage,long_df)
signif_antiox     <- calculate_significance(categories$Antioxidant_Response,long_df)
signif_er_up       <- calculate_significance(categories$ER_Stress_UPR,long_df)
signif_protsyn    <- calculate_significance(categories$Protein_Synthesis,long_df)
signif_proteolysis<- calculate_significance(categories$Proteolysis,long_df)
signif_mitoqc     <- calculate_significance(categories$Mitochondrial_Homeostasis,long_df)

# ---- 8. Generate Bar Plots ----
groups <- c("Young", "Pre", "Post")
colors <- c("#1b9e77", "#d95f02", "#7570b3")
group_order <- c("Young", "Pre", "Post")

data_oxidative_damage <- calculate_pathway_data(categories$Oxidative_Damage, long_df)
data_antioxidant_response <- calculate_pathway_data(categories$Antioxidant_Response, long_df)
data_er_stress_upr <- calculate_pathway_data(categories$ER_Stress_UPR, long_df)
data_protein_synthesis <- calculate_pathway_data(categories$Protein_Synthesis, long_df)
data_proteostasis <- calculate_pathway_data(categories$Proteostasis, long_df)