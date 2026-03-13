build_cluster_feature_map <- function(
  duckdb_parquet_path,
  output_path = NULL
) {
    
  # --- helpers ---
  .standardize_cols <- function(df, mapping, required) {
    # mapping: named list like list(std_name = original_name)
    nm <- names(df)
    for (std in required) {
      orig <- mapping[[std]]
      if (!std %in% nm) {
        if (!is.null(orig) && orig %in% nm) {
          names(df)[match(orig, names(df))] <- std
        } else {
          stop(sprintf("Required column '%s' not found (nor mapped) in data frame.", std), call. = FALSE)
        }
      }
    }
    df
  }

    # list of columns 
    cols = list(
    gene_protein   = list(genome_id = "genome_ids", protein_id = "protein_ids", Gene = "Gene"),
    struct         = list(struct    = "struct",     genome_id  = "genome_id",   value = "value"),
    domain         = list(protein_id= "AccNum",     domain_id  = "DB.ID"),
    cluster_members= list(cluster   = "cluster",     member     = "member")
  )
    
   # --- connect the duckdb of parquets ---  
    con <- DBI::dbConnect(duckdb::duckdb(), normalizePath(duckdb_parquet_path))
    
  # --- read & prep gene_protein ---
    gene_protein_table <- "genome_gene_protein"
    
  gp <- DBI::dbReadTable(con, gene_protein_table) |>
    tibble::as_tibble() |>
    dplyr::distinct()
    
  gp <- .standardize_cols(
    gp,
    mapping  = cols$gene_protein,
    required = c("genome_id", "protein_id", "Gene")
  ) 
  # --- read & prep struct ---
    struct_table <- "struct"
  st <- DBI::dbReadTable(con, struct_table) |>
    tibble::as_tibble()
  st <- .standardize_cols(
    st,
    mapping  = cols$struct,
    required = c("struct", "genome_id", "value")
  )
  st <- st |>
    dplyr::filter(.data$value == 1) |>
    dplyr::mutate(Gene = .data$struct) |>
    dplyr::mutate(.row_id = dplyr::row_number()) |>
    tidyr::separate_rows(.data$Gene, sep = "\\.") |>
    dplyr::group_by(.data$.row_id) |>
    dplyr::mutate(gene_order = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::select(-.data$.row_id) |>
    dplyr::inner_join(gp, by = c("Gene", "genome_id")) |>
    dplyr::distinct(struct, protein_id, gene_order)

  # --- read & prep domain ---
    domain_table <- "domain_names"
  dm <- DBI::dbReadTable(con, domain_table) |>
    tibble::as_tibble() |>
    dplyr::distinct()
  dm <- .standardize_cols(
    dm,
    mapping  = cols$domain,
    required = c("protein_id", "domain_id")
  )
  dm <- dm |>
    dplyr::distinct(.data$protein_id, .data$domain_id)

  # --- merge gpd ---
  gpd <- dm |>
    dplyr::full_join(gp, by = "protein_id") |>
    dplyr::full_join(st, by = "protein_id") |>
    dplyr::distinct(.data$protein_id, .data$domain_id, .data$Gene, .data$struct, .data$gene_order)

  # --- read & prep cluster members ---
    cluster_members_table <- "protein_members"
  cm <- DBI::dbReadTable(con, cluster_members_table) |>
    tibble::as_tibble() |>
    dplyr::distinct()
  cm <- .standardize_cols(
    cm,
    mapping  = cols$cluster_members,
    required = c("cluster", "member")
  )

  # --- long two-column mapping: cluster ~ feature ---
  cluster_feature <- gpd |>
    dplyr::select(-.data$gene_order) |>
    tidyr::pivot_longer(-.data$protein_id, names_to = "scale", values_to = "feature") |>
    dplyr::distinct(.data$protein_id, .data$feature) |>
    tidyr::drop_na(.data$protein_id, .data$feature) |>
    dplyr::mutate(protein_id = gsub("_len$", "", .data$protein_id)) |>
    dplyr::left_join(cm, by = dplyr::join_by("protein_id" == "member")) |>
    tidyr::pivot_longer(cols = c("protein_id", "feature"), names_to = "type", values_to = "feature") |>
    dplyr::distinct(.data$cluster, .data$feature) |>
    tidyr::drop_na(.data$cluster, .data$feature)
    
  out_dir <- if (is.null(output_path)) dirname(duckdb_parquet_path) else normalizePath(output_path)
    
cluster_feature_parquet <- file.path(out_dir, "cluster_feature.parquet")

  writeCompressedParquet <- function(df, path) {
    arrow::write_parquet(
      df,
      path,
      compression = "zstd",
      compression_level = 9,
      use_dictionary = TRUE
    )
  }
    
    cluster_feature |>  writeCompressedParquet(cluster_feature_parquet)
    
   DBI::dbExecute(con, sprintf("CREATE OR REPLACE VIEW gene_seqs AS SELECT * FROM read_parquet('%s')", cluster_feature_parquet))

}
