#' Downloads data from the COAT dataportal
#'
#' This function is based on scripts provided in the COATnor repository:
#' https://github.com/COATnor/coat_data_portal_example_scripts/blob/master/download_data_from_coat_data_portal.R
#' 
#' @param COAT_key string. API key for the COAT dataportal. Should be saved as an environmental variable.
#' @param COATdataset.name string. Name of the dataset you want to download (without version specifier).
#' @param COATdataset.version integer. Version of the dataset to download.
#' @param forceDownload logical. If TRUE, downloads data anew even if a local version is present.
#' If FALSE (default), downloads data only if it does not already exist locally.
#'
#' @return a dataframe or list of dataframes containing data downloaded from the COAT dataportal
#' @export
#'
#' @examples


downloadData_COAT <- function(COAT_key, COAT_module, 
                              COATdataset.name, COATdataset.version,
                              forceDownload = FALSE){
  
  ## Write full dataset name
  COATdataset.name2 <- paste0(COATdataset.name, "_v", COATdataset.version)
  
  ## Set up local storage directory for downloaded data (quicker access while COAT database is slow)
  if(!dir.exists("Data")){
    dir.create("Data")
  }
  
  if(!dir.exists("Data/COATdb_local")){
    dir.create("Data/COATdb_local")
  }
  
  ## Check whether data has already been downloaded and saved locally
  data_present <- ifelse(file.exists(paste0("Data/COATdb_local/", COATdataset.name2, ".rds")), TRUE, FALSE)
  
  ## Optional: load local file, otherwise download
  
  if(data_present & !forceDownload){
    
    message("Reading local datafile:")
    print((paste0("Data/COATdb_local/", COATdataset.name2, ".rds")))
    dat <- readRDS((paste0("Data/COATdb_local/", COATdataset.name2, ".rds")))
                   
  }else{
    
    message("Downloading dataset from https://data.coat.no:")
    print(COATdataset.name2)
    
    #### ------------------------------------------------------------------------------------------------------------ ####
    #### DOWNLOAD DATA FROM THE COAT DATA PORTAL
    #### ------------------------------------------------------------------------------------------------------------ ####
    
    ## this script can be used to download datasets from the COAT data portal
    ## the data can either be loaded into R or can be saved to a computer
    
    ## the development version of the ckanr package has to be installed (remotes::install_github("ropensci/ckanr"))
    
    ## ---------------------------------- ##
    ## SETUP
    ## ---------------------------------- ##
    
    ## setup the connection to the data portal
    COAT_url <- "https://data.coat.no/"  # write here the url to the COAT data portal
    COAT_key <- COAT_key  # write here your API key if you are a registered user, continue without API key if you are not registered
    # the API can be found on you page on the COAT data portal (log in and click on your name in the upper right corner of the page)
    # The use of an API key allows the user to access also non-public data
    
    ckanr::ckanr_setup(url = COAT_url, key = COAT_key)  # set up the ckanr-API
    
    ## Source functions provided by COAT
    source("https://github.com/COATnor/data_management_scripts/blob/master/download_data_from_coat_data_portal.R?raw=TRUE")
    
    ## ---------------------------------- ##
    ## DOWNLOAD DATA
    ## ---------------------------------- ##
    
    
    ## Check that selected module is available
    if(!(COAT_module %in% ckanr::organization_list(as = "table")$name)){
      stop("Selected module is not available on the COAT data portal. Please check the spelling of the module name and make sure that the module is available on the COAT data portal.")
    }
    
    ## Check that selected dataset is available
    if(!(COATdataset.name2 %in% list_datasets(module = COAT_module, printContents = FALSE)$name)){
      #if(!(COATdataset.name2 %in% list_datasets(module = COAT_module)$name)){
      stop("Selected dataset is not available in the selected module on the COAT data portal. Please check the spelling of the dataset name and make sure that the dataset is available in the selected module on the COAT data portal.")
    }
    
    ## list the names of all data files of the selected dataset (optional)
    #filenames <- list_data_files(COATdataset.name2, printContents = FALSE)
    message("Available files:")
    filenames <- list_data_files(COATdataset.name2)
    
    
    ## download data
    coat_dat <- download_coat_data(name = COATdataset.name2, # write here name of the dataset (choose from the list above)
                                   filenames = filenames[!grepl("readme|ReadMe|aux|coordinates", filenames)], # names of the files that should be downloaded, e.g. all data files (without readme, aux and coordinate file)
                                   store = "session", #"session" for importing the data in the current R session, "disk" for saving the data on your computer
                                   out.dir = NA) # if store = "disk", you have to specify the path to a folder where the data should be stored
    
    ## Bind yearly tables of the same format together
    groups <- vapply(coat_dat, function(x){
      paste(colnames(x), collapse = "|")
    }, character(1))
    
    dat <- unname(lapply(
      split(coat_dat, factor(groups, levels = unique(groups))),
      function(x) do.call(rbind, x)
    ))
    
    ## Save downloaded dataset to folder
    saveRDS(dat, paste0("Data/COATdb_local/", COATdataset.name2, ".rds"))
    
  }
  
  
  ## Return data
  return(dat)
}



