#' @title json_to_seurat
#'
#' @description Convert a JSON file into a Seurat object
#'
#' @param json_file
#'
#' @importFrom RcppSimdJson fload
#' @importFrom Matrix Matrix
#' @importFrom SeuratObject CreateAssayObject CreateDimReducObject as.Graph `DefaultAssay<-` `Idents<-`
#' @importFrom purrr map reduce pluck
#' @importFrom tibble enframe column_to_rownames as_tibble
#' @importFrom tidyr unnest
#' @importFrom dplyr left_join filter pull mutate across
#' @importFrom magrittr set_rownames set_colnames extract extract2 use_series
#' @importFrom rlang set_names
#' @importFrom janitor make_clean_names
#'
#' @return A Seurat object
#' @export
#'
#' @examples
json_to_seurat <- function(json_file, idents = NULL){
  obj <- fload(json_file)

  if (obj[["flavor"]] == "anndata"){
    if (!all(c("counts", "var", "obs") %in% names(obj))){
      stop("Either the 'counts', 'var', or 'obs' members are missing and this really isn't going to work without all three")
    }
  }

  obs <-
    map(
      .x = obj[["obs"]],
      .f = \(x) unnest(
        data = enframe(x),
        cols = value
      )
    ) |>
    reduce(
      .f = left_join,
      by = "name"
    ) |>
    set_colnames(
      value = c("cell", make_clean_names(names(obj[["obs"]])))
    ) |>
    column_to_rownames(var = "cell")

  var <-
    map(
      .x = obj[["var"]],
      .f = \(x) unnest(
        data = enframe(x),
        cols = value
      )
    ) |>
    reduce(
      .f = left_join,
      by = "name"
    ) |>
    set_colnames(
      value = c("gene", make_clean_names(names(obj[["var"]])))
    ) |>
    column_to_rownames(var = "gene")

  if (nrow(obj[["counts"]]) == nrow(var)){
    counts <-
      obj[["counts"]] |>
      Matrix(sparse=TRUE, dimnames = list(rownames(var), rownames(obs)))
  } else if (nrow(obj[["counts"]]) == nrow(obs)){
    counts <-
      obj[["counts"]] |>
      t() |>
      Matrix(sparse=TRUE, dimnames = list(rownames(var), rownames(obs)))
  }

  RNA <- CreateAssayObject(
    counts = counts
  )
  RNA@key <- "rna_"

  rna_var_features <-
    as_tibble(
      x = var,
      rownames = "gene"
    ) |>
    filter(highly_variable == TRUE) |>
    pull("gene")
  RNA@meta.features <- var
  RNA@var.features <- rna_var_features

  if ("scale_data" %in% names(obj)){
    if (nrow(obj[["scale_data"]]) == nrow(var)){
      scale_data <-
        obj[["scale_data"]] |>
        set_rownames(rownames(var)) %>%
        set_colnames(rownames(obs))
    } else if (nrow(obj[["scale_data"]]) == nrow(obs)){
      scale_data <-
        obj[["scale_data"]] |>
        t() |>
        set_rownames(rownames(var)) %>%
        set_colnames(rownames(obs))
    } else {
      message("scale_data dimensions do not match count data")
    }

    RNA@scale.data <- scale_data
  }

  seurat_obj <- CreateSeuratObject(
    counts = RNA,
    meta.data = obs,
  )

  reductions <-
    map(
      .x = names(obj[["obsm"]]),
      .f = \(x) map(
        .x = obj[["obsm"]][[x]],
        .f = \(y) unnest(
          enframe(y),
          cols = value
        )
      ) %>%
        reduce(
          .f = left_join,
          by="name"
        ) %>%
        set_colnames(
          value = c("cell", names(obj[["obsm"]][[x]]))
        )
    ) %>%
    set_names(names(obj[["obsm"]]))

  if ("paga" %in% names(obj)){
    paga_groups <-
      obs %>%
      extract2(
        obj %>%
          use_series("paga") %>%
          use_series("groups")
      ) %>%
      unique() %>%
      as.numeric() %>%
      sort()

    paga <- list(
      connectivities           = obj[["paga"]][["connectivities"]] |>
        extract(paga_groups+1, paga_groups+1) |>
        set_rownames(levels(paga_groups)) |>
        set_colnames(levels(paga_groups)),
      connectivities_tree      = obj[["paga"]][["connectivities_tree"]],
      group_name               = obj[["paga"]][["groups"]],
      groups                   = paga_groups,
      position                 = as_tibble(
        cbind(
          paga_groups,
          extract(obj[["paga"]][["pos"]], paga_groups+1,)
        ),
        .name_repair = ~make.names(c("group","x", "y"))
      ) |>
        mutate(
          across(
            x:y,
            .fns = as.numeric
            ),
        )
    )

    paga[["edges"]] <- tibble(
      group1 = paga[["groups"]][row(paga$connectivities)[upper.tri(paga$connectivities)]],
      group2 = paga[["groups"]][col(paga$connectivities)[upper.tri(paga$connectivities)]],
      weight = paga[["connectivities"]][upper.tri(paga$connectivities)] %>% as.numeric()
    ) %>%
      mutate(
        x1 = paga[["position"]][["x"]][match(.[["group1"]], paga[["position"]][["group"]])] |>  as.numeric(),
        y1 = paga[["position"]][["y"]][match(.[["group1"]], paga[["position"]][["group"]])] |> as.numeric(),
        x2 = paga[["position"]][["x"]][match(.[["group2"]], paga[["position"]][["group"]])] |> as.numeric(),
        y2 = paga[["position"]][["y"]][match(.[["group2"]], paga[["position"]][["group"]])] |> as.numeric()
      ) |>
      filter(weight >= 0.1)
  }

  if ("pca" %in% names(reductions)){
    PCA <- CreateDimReducObject(
      embeddings =
        reductions |>
        use_series("pca") |>
        column_to_rownames(var = "cell") %>%
        set_colnames(paste("PC", seq(ncol(.)), sep="_")) |>
        as.matrix(),
      loadings =
        obj |>
        use_series("varm") |>
        pluck(1) |>
        set_rownames(rownames(var)) %>%
        set_colnames(paste("PC", seq(ncol(.)), sep = "_")) |>
        as.matrix(),
      assay = "RNA",
      global = FALSE,
      key = "PC_"
    )

    seurat_obj[["pca"]] <- PCA

    reductions[["pca"]] <- NULL
  }

  for (x in names(reductions)){
    y <- CreateDimReducObject(
      embeddings =
        column_to_rownames(
          .data = reductions[[x]],
          var = "cell"
          ) |>
        as.matrix(),
      key = paste0(x, "_"),
      assay = "RNA"
    )
    seurat_obj[[x]] <- y
  }

  # It does not appear that the connectivities are used
  # connections <-
  #   obj[["obsp"]][["connectivities"]] %>%
  #   set_rownames(rownames(obs)) %>%
  #   set_colnames(rownames(obs)) %>%
  #   Matrix(sparse=TRUE)

  # SeuratDisk is supposed to set this to obj@graphs@ASSAY_nn
  # Looks like we need to save some of the information from the 'uns' member
  # after all
  if ("obsp" %in% names(obs)){
    if ("distances" %in% names(obs[["obsp"]])){
      distances <-
        obj[["obsp"]][["distances"]] |>
        set_rownames(rownames(obs)) |>
        set_colnames(rownames(obs)) |>
        Matrix(sparse=TRUE)
      seurat_obj@graphs[["RNA_nn"]] <- as.Graph(distances)
    }
  }

  DefaultAssay(seurat_obj) <- "RNA"

  if (!is.null(idents)){
    if (idents %in% colnames(var)){
      Idents(seurat_obj) <- idents
    } else {
      Idents(seurat_obj) <- "Converted"
    }
  } else {
    Idents(seurat_obj) <- "Converted"
  }
  seurat_obj
}
