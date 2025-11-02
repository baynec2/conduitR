#' Getting Cazyme infromation for a vector of cazyme ids
#'
#' @param uniprot_annotated_data
#'
#' @returns
#' @export
#'
#' @examples
get_cazyme_info = function(cazyme_ids){

# Reading in cazyme infomation
  url <-"https://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.07312019.fam-activities.txt"

  cazyme_db <- readr::read_delim(url,skip = 1,col_names = c("cazy_id",
                                                            "cazy_description"))

  # Annotating with cazyme class
  class_lookup <- tibble::tibble(
    code = c("GH","GT","PL","CE","AA","CBM"),
    cazy_class = c(
      "glycoside_hydrolase",
      "glycosyl_transferase",
      "polysaccharide_lyases",
      "carbohydrate_esterases",
      "auxiliary_activities",
      "carbohydrate_binding_module"
    )
  )

  # Extract family code (first 2 letters for most, but CBM is 3)
  cazyme_db <- cazyme_db |>
    dplyr::mutate(
      code = ifelse(stringr::str_detect(cazy_id, "^CBM"), "CBM", stringr::str_sub(cazy_id, 1, 2)),
      cazy_class = dplyr::left_join(tibble::tibble(code = code), class_lookup, by = "code")$cazy_class
    )

  xref_cazy = tibble::tibble(xref_cazy = cazyme_ids)

  # Join UniProt annotated data with CAZyme annotations
  annotated <- xref_cazy |>
    dplyr::inner_join(cazyme_db, by = c("xref_cazy" = "cazy_id"))

  # Split multiple CAZy families into separate rows
  uniprot_expanded <- uniprot_annotated_data |>
    tidyr::separate_rows(xref_cazy, sep = ";")

  # Join with CAZyme DB
  annotated <- uniprot_expanded |>
    dplyr::inner_join(cazyme_db, by = c("xref_cazy" = "cazy_id"))

  # Extracting Cazyme Information
  out <- annotated |>
    dplyr::select(protein_id,xref_cazy,cazy_class,cazy_description)

  return(out)
}
