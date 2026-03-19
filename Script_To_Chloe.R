############################################
## WORKFLOW TELEMTRY POLAR FOXES TO IPM 
############################################

## MAB 18032026


###########################
## 0. LOAD LIBRARIES
###########################

library(tidyverse)
library(sf)
library(trip)
library(raster)
library(sp)
library(dplyr)
library(lubridate)
library(tidyr)
library(geosphere)
library(ggplot2)
library(purrr)
library(trip)
library(patchwork)

############################
## 1. METADATA
############################
## Careful some PTT have been re-used on different ID
## Some ID have been instrumented several times
## Use the metadata with deployment date, 
## date of death confirmed by carcass or tag retrieval 
## date of death estimated by expert looking at telemetry plots 
## and add the estimated date of death with telemetry data
## in addition we have a last date alive (essentially last date of telemetry
## unless the collar malfunctionned)

setwd("D:/MarieHD/MArieHD_UiT/New_PC_NPI/Foxes")

metadata1<-read.csv("Metadata_Fox_corrected.csv",sep=",") %>% 
  mutate(Date_deploy=lubridate::dmy(paste(Day,Month,Year)))%>% 
  dplyr::select(Date_deploy,ID_Fox,PTT,LAT_N,LON_E,area,kl,Sex,Est.Age.CORR,Weight) %>% 
  mutate(ID_PTT=paste(ID_Fox,PTT,sep="-")) %>% 
  mutate(datetime_dep=lubridate::ymd_hm(paste(Date_deploy,kl,sep=" "))) %>% 
  mutate(year=lubridate::year(Date_deploy)) %>% 
  mutate(Fox_PTT_year=paste(ID_Fox,PTT,year,sep="-")) %>% 
  rename(Age=Est.Age.CORR)


metadata2<-read.csv(
  "Dead or alive_tag_04022026_Corrected FOX_ID_+corrected.csv") %>% 
  rename(ID_Fox=ID.Fox) %>% 
  mutate(Date_deploy=lubridate::dmy(paste(Day,Month,Year))) %>% 
  mutate(Date_death=lubridate::my(paste(Month.dead,Year.dead))) %>% 
  mutate(Date_death_est=lubridate::my(paste(Month.Est.dead,Year.Est.dead))) %>% 
  mutate(Last_date_alive=lubridate::my(paste(Month_Expert,Year_Expert))) %>% 
  dplyr::select(ID_Fox,PTT,Sex,Date_deploy,Date_death,Date_death_est,Last_date_alive) %>% 
  mutate(ID_PTT=paste(ID_Fox,PTT,sep="-"))

## sanity check
length(unique(metadata1$ID_Fox))==length(unique(metadata2$ID_Fox))## TRUE
length(unique(metadata1$ID_PTT))==length(unique(metadata2$ID_PTT))## TRUE
metadata1$ID_PTT==metadata2$ID_PTT## TRUE

## Join both metadata
metadata <- metadata1 %>%
  left_join(metadata2 %>% dplyr::select(ID_PTT,Date_deploy,Date_death,Date_death_est,
                                 Last_date_alive), by = "ID_PTT")

metadata$Date_deploy.x==metadata$Date_deploy.y

metadata<-metadata %>% 
  dplyr::select(-Date_deploy.x) %>% 
  rename(Date_deploy=Date_deploy.y)

rm(metadata1,metadata2)

###############################################
## 2. IMPORT MAPS
###############################################

## The polygons for Nordenskjoldland and Kongsfjorden and SV are hand drawn 
## with the siny app "DrawSavePolygon.R"

KF<-st_read("Kongsfjorden.shp")
NS<-st_read("Nordenskjoldland.shp")
SV<-st_read("Spitsbergen.shp")
KF_NS <- st_union(KF, NS)
# ── Remove KF and NS from SV ──────────────────────────────────────────────────
SV_clean <- st_difference(SV, st_union(KF_NS))

## Find instrumentation location
meta_sf <- st_as_sf(meta,
                    coords = c("LON_E", "LAT_N"),
                    crs    = 32633)  # UTM 33N

# ── Transform to WGS84 to match polygons ──────────────────────────────────────
meta_sf <- st_transform(meta_sf, crs = 4326)

# ── Spatial join ──────────────────────────────────────────────────────────────
meta_sf <- st_join(meta_sf, areas, join = st_within)

# ── Check ─────────────────────────────────────────────────────────────────────
table(meta_sf$area, useNA = "ifany")



####################################################
## 3. BATCH DOWNLOAD TELEMTRY DATA FROM NPI SERVER
####################################################
## One needs to be physically at NPI or to use a FTP

Platform_list<-list.files("//npdata/PRODUCTION/satellite-tracking/argos/product/program-11660/csv/")## 91 platforms

theDir<-c("//npdata/PRODUCTION/satellite-tracking/argos/product/program-11660/csv/")
setwd(theDir)

Foxes_Locs_raw<-read_csv(Platform_list[1]) %>% 
  mutate_all(as.character)
for ( i in 2:length(Platform_list)){
  print(i)
  t<-read_csv(Platform_list[i]) %>% mutate_all(as.character)
  Foxes_Locs_raw<-bind_rows(Foxes_Locs_raw,t)
}

setwd("G:/MArieHD_UiT/New_PC_NPI/Foxes")

Foxes_Locs_raw$measured<-lubridate::ymd_hms(Foxes_Locs_raw$measured)

dim(Foxes_Locs_raw)## 983796
length(unique(Foxes_Locs_raw$platform_id))##89


## Join with metadata
Foxes_Locs_raw <- map_dfr(unique(metadata$PTT), function(ptt) {
  
  meta  <- metadata %>% filter(PTT == ptt)
  locs  <- Foxes_Locs_raw %>% filter(platform == ptt)
  
  if (nrow(locs) == 0) return(NULL)
  
  if (nrow(meta) == 1) {
    # ── Single deployment: join all metadata columns directly ──────────────
    locs %>%
      bind_cols(meta %>% dplyr::select(-PTT) %>% slice(1))
    
  } else {
    # ── Multiple deployments: split by datetime and join correct metadata ──
    map_dfr(1:nrow(meta), function(j) {
      
      # Define time window for this deployment
      t_start <- meta$datetime_dep[j]
      t_end   <- if (j < nrow(meta)) meta$datetime_dep[j + 1] else as.POSIXct(Inf)
      
      locs %>%
        filter(measured >= t_start & measured < t_end) %>%
        bind_cols(meta %>% dplyr::select(-PTT) %>% slice(j))
    })
  }
}) 
 

dim(Foxes_Locs_raw)##983715
length(unique(Foxes_Locs_raw$ID_PTT))##94

## The alternative is to load the R.Data file "Foxes_Locs_raw2026-02-05.RData"
## the object is Foxes_Locs_raw, but not associated to the full metadata

#####################################################################
##save(Foxes_Locs_raw, file=paste("Foxes_Locs_raw",as.Date(Sys.time()),".RData",sep=""))
###################################################################


####################################################
## 4. CLEAN FORMAT TELEMETRY DATA
####################################################

## Each data point has a time stamp
## Not every data point has a lat and lon 
## Deployed is the deployment date from metadata

Foxes_Locs_raw_formatted<-Foxes_Locs_raw%>%
  dplyr::select(ID_Fox,ID_PTT,platform,year,Fox_PTT_year,Date_deploy,
                measured,latitude,longitude,lc,hdop,error_radius,
                semi_major,semi_minor,orientation,
                temperature,activity_3_days_ago,activity_today,activity_yesterday,
                LAT_N,LON_E,area,kl,Sex,Age,Weight,Date_death,Date_death_est,Last_date_alive) %>% 
  dplyr::filter(!is.na(measured)) %>% 
  group_by(ID_PTT) %>% 
  dplyr::arrange(measured,.by_group =TRUE) %>%
  dplyr::filter(measured>Date_deploy) %>% ## subset by tagging date
  ungroup() %>% 
  mutate_at(vars(ID_PTT),~as.character(.)) %>% 
  mutate_at(vars(semi_major),~as.numeric(.)) %>% 
  mutate_at(vars(semi_minor),~as.numeric(.)) %>% 
  mutate_at(vars(orientation),~as.numeric(.)) %>% 
  mutate_at(vars(temperature),~as.numeric(.)) %>% 
  mutate_at(vars(latitude),~as.numeric(.)) %>% 
  mutate_at(vars(longitude),~as.numeric(.)) %>% 
  mutate_at(vars(hdop),~as.numeric(.)) %>% 
  mutate_at(vars(error_radius),~as.numeric(.)) %>% 
  dplyr::filter(lc!="Z") %>% 
  dplyr::rename(date=measured) %>% 
  group_by(ID_PTT) %>% 
  dplyr::arrange(date,.by_group =TRUE) %>%
  mutate(dtime_sec=as.vector(difftime(date, lag(date),units="sec"))) %>% 
  mutate(date=trip::adjust.duplicateTimes(date,ID_PTT)) %>% 
  dplyr::arrange(date,.by_group =TRUE) %>%
  mutate(dtime_sec=as.vector(difftime(date, lag(date),units="sec"))) %>% 
  rename(area_capture=area)
  

## Sanity checks
dim(Foxes_Locs_raw_formatted)##966610
length(unique(Foxes_Locs_raw_formatted$platform))## 90; Fox ID 36 PTT 167930 has not transmitted
length(unique(Foxes_Locs_raw_formatted$ID_PTT))## 93
length(which(Foxes_Locs_raw_formatted$dtime_sec==0))## 0
length(which(is.na(Foxes_Locs_raw_formatted$latitude)))##0

################################################
## 5. ASSIGN A LOCATION TO EACH LAT/LON
################################################

## "LAT_N, LON_E" deployment locations in UTM 33N
## locations from the sender are in 4326

# ── Convert locations to sf object ───────────────────────────────────────────
Foxes_sf <- st_as_sf(Foxes_Locs_raw_formatted,
                     coords = c("longitude", "latitude"),
                     crs    = 4326)

# ── Add polygon name to each layer before joining ────────────────────────────
KF       <- KF       %>% mutate(area = "Kongsfjorden")
NS       <- NS       %>% mutate(area = "Nordenskjoldland")
SV_clean <- SV_clean %>% mutate(area = "Spitsbergen")

# ── Combine all polygons into one sf object ───────────────────────────────────
areas <- bind_rows(
  KF       %>% dplyr::select(area),
  NS       %>% dplyr::select(area),
  SV_clean %>% dplyr::select(area)
)
# ── Spatial join — assigns area name to each location ────────────────────────
Foxes_sf <- st_join(Foxes_sf, areas, join = st_within)

Foxes_Locs_raw_formatted$area <- Foxes_sf$area
  
rm(Foxes_sf)

###########################################################################
## save(Foxes_Locs_raw_formatted, file=paste("Foxes_Locs_raw_formatted",as.Date(Sys.time()),".RData",sep=""))
###########################################################################


############################################################################
## 6. CREATE MONTHLY SUMMARY DATA
#############################################################################
## The monthly speed is calculated based on locations with lc==3
## and speed > 10kmh is removed because likely due to large ARGOS errors
## Temp > 40C is not considered as likely the result of a wrong sensor reading
## The area refers to teh area the fox is in during that month only based on lc=3
## if the fox were between several areas they are collated
## the area_capture is the area where the fox was captured


summarize_fox_data <- function(data,output_dir = "fox_boxplots") {
  
  
  fox_data <-data%>%
    dplyr::arrange(ID_PTT,date) %>%
    mutate(month_year = paste0(lubridate::month(date, label = TRUE, abbr = TRUE), "-", lubridate::year(date))) %>% 
    mutate(activity_today = as.numeric(activity_today)) %>%
    mutate(
      lag_lat = lag(latitude),
      lag_lon = lag(longitude),
      lag_time = lag(date)
    ) %>%
    mutate(
      dist_m = distHaversine(cbind(longitude, latitude), cbind(lag_lon, lag_lat)),
      time_diff = as.numeric(difftime(date, lag_time, units = "secs"))
    ) %>%
    filter(!is.na(dist_m), !is.na(time_diff), time_diff > 0) %>%
    mutate(
      # Calculate speed only for lc = 3
      speed_kmh = ifelse(lc == 3, (dist_m / time_diff) * 3.6, NA),
      speed_kmh = ifelse(speed_kmh > 10, NA, speed_kmh) # Remove speeds > 10 km/h
    ) %>%
    mutate(
      # Flag rows where the time difference to the next point is greater than 30 days
      is_before_gap = ifelse(lead(time_diff, default = 0) > 30 * 24 * 60 * 60, TRUE, FALSE))
  
  # Summarize data by ID_PTT and month-year
  summary_table <- fox_data %>%
    group_by(ID_PTT, month_year) %>%
    summarise(
      avg_speed = mean(speed_kmh, na.rm = TRUE),      # Average speed
      avg_temperature = mean(ifelse(temperature < 40, temperature, NA), na.rm = TRUE), # Average temperature
      avg_activity = mean(activity_today, na.rm = TRUE), # Average activity
      num_data_points = n(),                        # Number of data points
      num_lc_3 = sum(lc == "3", na.rm = TRUE),# Number of locations with lc = 3
      ID_Fox          = first(ID_Fox),
      Fox_PTT_year    = first(Fox_PTT_year),
      Sex             = first(Sex),
      Age             = first(Age),
      Weight          = first(Weight),
      area = paste(unique(na.omit(area[lc == "3"])), collapse = " / "),  # lc==3 only
      Date_deploy     = first(Date_deploy),
      area_capture=first(area_capture),
      Date_death      = first(Date_death),
      Date_death_est  = first(Date_death_est),
      Last_date_alive = first(Last_date_alive),
      LAT_N           = first(LAT_N),
      LON_E           = first(LON_E),
    ) %>%
    ungroup() %>% 
    # Convert month_year to a proper date format for sorting
    mutate(area            = na_if(area, ""),
      month_year_date = parse_date_time(month_year, orders = "my")) %>%
    # Arrange by Fox_PTT and month_year_date
    arrange(ID_PTT, month_year_date) %>%
    # Drop the temporary sorting column
    dplyr::select(-month_year_date)
  
  complete_table <- summary_table %>%
    mutate(
      month_year_date = as.Date(parse_date_time(month_year, orders = "my")) # Convert month_year to Date format
    ) %>%
    group_by(ID_PTT) %>%
    complete(
      month_year_date = seq.Date(
        from = min(month_year_date, na.rm = TRUE),
        to = max(month_year_date, na.rm = TRUE),
        by = "month"
      )
    ) %>%
    ungroup() %>%
    # Recreate the month_year column for display
    mutate(month_year = paste0(month(month_year_date, label = TRUE, abbr = TRUE), "-", year(month_year_date))) %>%
    # Arrange by Fox_PTT and month_year_date
    arrange(ID_PTT, month_year_date) %>%
    # Drop the temporary sorting column
    dplyr::select(-month_year_date)
  
  # --- Generate and Save Combined Boxplots ---
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Loop through each Fox_PTT to create combined plots
  unique_foxes <- unique(complete_table$ID_PTT)
  
  for (fox in unique_foxes) {
    # Filter data for the current Fox_PTT
    data <- fox_data %>%
      filter(ID_PTT == fox)
    
    # Ensure months without data are included
    data <- data %>%
      mutate(
        month_year = paste0(lubridate::month(date, label = TRUE, abbr = TRUE), "-", lubridate::year(date))
      ) %>%
      complete(
        month_year = unique(complete_table$month_year[complete_table$ID_PTT == fox]),
        fill = list(activity_today = NA, temperature = NA, speed_kmh = NA)
      ) %>%
      mutate(
        month_year = factor(month_year, levels = unique(complete_table$month_year[complete_table$ID_PTT == fox])) # Order months
      )
    
    
    # Create boxplots for activity, temperature, and speed
    p1 <- ggplot(data, aes(x = month_year, y = activity_today)) +
      geom_boxplot(outlier.color = "red", fill = "lightblue") +
      labs(title = "Activity", x = "Month-Year", y = "Activity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    p2 <- ggplot(data %>% filter(temperature < 40), aes(x = month_year, y = temperature)) +
      geom_boxplot(outlier.color = "red", fill = "lightgreen") +
      labs(title = "Temperature", x = "Month-Year", y = "Temperature (°C)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    p3 <- ggplot(data, aes(x = month_year, y = speed_kmh)) +
      geom_boxplot(outlier.color = "red", fill = "lightpink") +
      labs(title = "Speed", x = "Month-Year", y = "Speed (km/h)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Combine the three plots into one figure
    combined_plot <- (p1 / p2 / p3) + 
      plot_annotation(title = paste("Monthly Boxplots for Fox:", fox))
    
    # Save the combined plot as a PNG file
    ggsave(filename = file.path(output_dir, paste0(fox, "_monthly_boxplots.png")), plot = combined_plot, width = 12, height = 10)
  }
  
  return(complete_table)
}

summary_table <- summarize_fox_data(Foxes_Locs_raw_formatted)

dim(summary_table)
length(unique(summary_table$ID_PTT))##93


############################################################################
## 7. ADD DEAD OR ALIVE STATUS
#############################################################################
## dead= temp<0 and activity <10
## When we have NAs sometimes the sensors do not work and then assigned as collar malfunction
## alive = temp>0 and activity >25
## maybe = everything else

result <- assign_alive_status(summary_table)

assign_alive_status <- function(summary_table, temp_dead = 0, activity_dead = 10, temp_alive = 0, activity_alive = 25) {
  # Ensure the input table has the required columns
  required_cols <- c("avg_temperature", "avg_activity")
  if (!all(required_cols %in% names(summary_table))) {
    stop("The input table must contain the following columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Assign alive status based on criteria
  summary_table <- summary_table %>%
    mutate(
      alive = case_when(
        # Dead: Temperature < temp_dead AND Activity < activity_dead
        avg_temperature < temp_dead & avg_activity < activity_dead ~ "N",
        
        # Dead: When we have NAs in temperature or activity
        is.na(avg_temperature) | is.na(avg_activity) ~ "Possible collar malfunction",
        
        # Alive: Temperature > temp_alive AND Activity > activity_alive
        avg_temperature > temp_alive & avg_activity > activity_alive ~ "Y",
        
        # Maybe: Everything else
        TRUE ~ "Maybe"
      )
    )
  
  return(summary_table)
}

summary_table$alive<-result$alive
rm(result)

