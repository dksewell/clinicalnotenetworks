pacman::p_load(lubridate,
               dplyr,
               magrittr,
               tidygraph,
               RSQLite,
               DBI)

pat_id = "36257634596"
pat_id = "10525016660"
pat_id = "78340080308"
db =
  dbConnect(SQLite(),
            "E:/R01-SMART_cancer_care/smart_cancer_care.db")
site = "UCDavis"
include_med_students = FALSE
verbose = TRUE
only_earliest_read = TRUE

rm(list = setdiff(ls(),
                  c("db","site","include_med_students","verbose")))


pat_ids =
  list.files("E:/R01-SMART_cancer_care/patient_networks/ucdavis/bp_network_data/")
pat_ids %<>%
  str_replace("bp-","") %>%
  str_replace(".RDS","")


# Find patient with notes having their first entry as view ----------------

first_entry_is_view =
  db %>%
  tbl(paste0("access_logs_",site)) %>%
  filter(NOTE_TYPE != "Telephone Encounter") %>%
  group_by(NEW_DATA_OBFUS_ID) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  filter(EVENT_ACTION == "View") %>%
  collect()


# Find a patient with med students modifying ------------------------------



pb = txtProgressBar(0,length(pat_ids),style=3)
for(id in 1:length(pat_ids)){
  id0 = pat_ids[id]
  alogs =
    db %>%
    tbl(paste0("access_logs_",site)) %>%
    filter(PAT_OBFUS_ID == id0, # PAT_OBFUS_ID needs to be of same format (chr vs. numeric) as patient level data
           NOTE_TYPE != "Telephone Encounter") %>%
    collect() %>%
    mutate(date_time = dmy_hms(ACCESS_TIME_CHAR)) %>% # Convert access_time_char to date
    arrange(date_time) %>%
    mutate(ACCESS_USER_OBFUS_ID = as.character(ACCESS_USER_OBFUS_ID), # This will help merge with note ids later.
           log_event_number = 1:n()) # This is to help with several operations later

  med_student_authored_notes =
    alogs %>%
    filter(ACCESS_USER_PROV_TYPE == ".STUDENT: MEDICAL",
           EVENT_ACTION == "Modify")

  setTxtProgressBar(pb,id)
  if(nrow(med_student_authored_notes) > 0) stop("Found one!")
}


pat_id = "10525016660"


# Find patients with multiple authorship ----------------------------------


pat_ids =
  list.files("E:/R01-SMART_cancer_care/patient_networks/ucdavis/bp_network_data/")
pat_ids %<>%
  str_replace("bp-","") %>%
  str_replace(".RDS","")


complex_notes_ids =
  db %>%
  tbl(paste0("access_logs_",site)) %>%
  filter(EVENT_ACTION == "Modify",
         NOTE_TYPE != "Telephone Encounter",
         ACCESS_USER_PROV_TYPE != ".STUDENT: MEDICAL") %>%
  select(PAT_OBFUS_ID,
         ACCESS_USER_PROV_TYPE,
         NEW_DATA_OBFUS_ID) %>%
  distinct() %>%
  collect() %>%
  group_by(NEW_DATA_OBFUS_ID) %>%
  summarize(NEW_DATA_OBFUS_ID = NEW_DATA_OBFUS_ID[1],
            PAT_OBFUS_ID = PAT_OBFUS_ID[1],
            n = n()) %>%
  filter(n > 1)
complex_notes_ids %>%
  arrange(desc(n))

pat_id = "78340080308"

