## This is the code to run script_RASo_class in a shiny enviornment

######## INSTALLATION ALL PACKAGES AND CONNECTIONS
library(shiny)
require("RPostgreSQL")
library(RMySQL)
require(stringr)
library(DT)


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  
  # Application title
  tags$br(),
  fluidRow(
    column(5, tags$img(height = 80, src = "Cobranding_IGTP_Hospital_ENG.png"),offset = 9)
    )
  ,
  titlePanel(h1("RASopathy-related variant classification", align = "center")),
  tags$br(),
  tags$br(),
  
  # Sidebar with a slider input for number of bins 
  # Numeric Input with variant identifier in Pandora 
  sidebarLayout(
    sidebarPanel(     
      helpText("Please, indicate the variant's identifier at Pandora
               [ctl + shift + j]"),
      
      textInput(inputId = "id", label = "Specify variant ID in Pandora"),
      numericInput(inputId = "denovo_noconfirmed", label = "Number of the novo cases reported, paternity non-confirmed", value = 0, min = 0),
      numericInput(inputId = "denovo_confirmed", label = "Number of the novo cases reported, paternity confirmed", value = 0, min = 0),
      numericInput(inputId = "cosegregation", label = "Number of the cosegregated families reported", value = 0, min = 0),
      selectInput(inputId = "Functional_evidence", label = "Are functional studies demostrating variant pathogenicity?", 
                  choices = list("There is no evidence" = 0, 
                                 "There is Functional_evidence" = 1), selected = 0),
      selectInput(inputId = "PPAT_evidence", label = "Are relevant references demostrating variant pathogenicity?", 
                         choices = list("There is no evidence" = 0, 
                                        "There is PPAT_evidence" = 1), selected = 0),
      selectInput(inputId = "PPOL_evidence", label = "Are relevant references demostrating variant neutrality?", 
                        choices = list("There is no evidence" = 0, 
                             "There is PPOL_evidence" = 1), selected = 0),
      actionButton("go", "Search")
  ), 
    

    mainPanel(
      tags$h4(tags$style('h4 {color:darkblue;}')),
      fluidRow( 
                column(9, tags$h4("Variant to classify"), offset = 4),
                column(9, tableOutput("Variant"), offset = 0)
     ),
     
      tags$hr(),
     
   fluidRow( 
     column(9, tags$h4("Criteria used to classify this variant"), offset = 3),
     column(9, tableOutput("AutoClass"), offset = 4)
    
   ),
      
      tags$hr(),
 
     fluidRow( 
     column(9, tags$h4("Variant classification following ACMG guidelines"), offset = 2),
     column(9, tableOutput("FinalClass"), offset = 2)
     )
    )
  )
)

# Define server logic required to connect to DB and run all functions to classifiy RASopathy variant
server <- function(input, output, session) {
  
    ## conection to UCSC
  my_connection <- dbConnect(
    MySQL(),
    user="genome",
    dbname="hg38",
    host="genome-euro-mysql.soe.ucsc.edu",
    port=3306
  )
  
  ## connection to Pandora DB
  get_db_parameters <- function(db_conf) {
    params <- read.table("params.csv", sep = ",", stringsAsFactors = FALSE)
    return(list(user = params$V2[1],
                password = params$V2[2],
                dbname = params$V2[3],
                host = params$V2[4],
                port = params$V2[5]))
  }
  
  ## connects to the NGS BD using a config file
  ## param db_conf: a full path to a config file (CSV)
  
  db_connect_postgres <- function(db_conf) {
    drv <- dbDriver("PostgreSQL")
    
    db_conf <- get_db_parameters(db_conf)
    
    con <- dbConnect(drv,
                     user = db_conf[['user']],
                     password = db_conf[['password']],
                     dbname = db_conf[['dbname']],
                     host = db_conf[['host']],
                     port = db_conf[['port']])
    
    return(con)
  }
  
  con <- db_connect_postgres('db_pandora.conf')
  
  ####### Stablished criteria are different depending on GOF genes and LOF genes. 
  
  #LOF Genes: NF1 & SPRED1. BP1 % PVS1 criteria are only applicable to LOF genes
  #GOF Genes: BRAF, CBL, HRAS, KRAS, LZTR1, NF1, NRAS, MAP2K1, MAP2K2, PTPN11, RAF1, RIT1, SHOC2, SOS1, SOS2, SPRED1 i RASA1,
  
  ##### Calling RASopathies gene information
  domain_groupRAF <- read.csv("files/domini_grupRAF.csv")
  domain_groupRAS <- read.csv("files/domini_grupRAS.csv")
  domain_groupSOS <- read.csv("files/domini_grupSOS.csv")
  domain_groupMAPK <- read.csv("files/domini_grupMAPK.csv")
  Transcripts_RASos <- read.csv("files/Transcripts_RASos.csv")
      
  ####### AUTOMATIC CRITERIA FUNCTION
  
      
 Automatic_criteria_AMCG<- function(id = input$id, con = con, denovo_noconfirmed = input$denovo_noconfirmed,
                                         denovo_confirmed = input$denovo_confirmed, cosegregation = input$cosegregation,
                                         PPAT_evidence = as.numeric(as.character(input$PPAT_evidence)), 
                                         PPOL_evidence = as.numeric(as.character(input$PPOL_evidence)),
                                    Functional_evidence = as.numeric(as.character(input$Functional_evidence))){
   
   if (is.null(input$id) | is.null (input$denovo_noconfirmed) | is.null (input$denovo_confirmed) |
       is.null(input$cosegregation) | is.null (input$PPAT_evidence) |
       is.null (input$PPOL_evidence) | is.null(input$Functional_evidence) | input$go == 0) {
     
     ## nothing to do here
     return(NULL)
   }
  
   
        #dataframe to fill with all criteria (each criteria in one line)
        criteria<-data.frame(criteria=c(rep(NA, 37)), row.names = c("PS1", "PS2_veryStrong", "PS2", "PS3", "PS4_strong", "PS4_moderate", "PS4_supporting", "PM5_strong",  "PM6_veryStrong", "PM6","PM1", "PM2", "PM4",  "PP1_strong", "PP1_moderate", "PP1_supporting", "PP2", "PP3","PP5", "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7", "BS1_supporting", "PM6_strong", "PM5", "PM5_supporting", "BP8", "PVS1"))
 

        group_position<-NULL
        num<-NULL
        ### BP1 & PVS1 cirteria. 
        #BP1: Missense variant in a gene for which primarily truncating variants are known to cause the disease
        #PVS1: NUll variant. Only applicable to NF1 and SPRED1
        if (is.null(id)) {
          
          ## nothing to do here
          return(NULL)}
        
        else{
        
        bp1.1<-paste('SELECT "start", "end", "chr", "validatedEffect", "effect", "symbol"
                     FROM "VA_VariantsInTranscripts" AS a
                     LEFT OUTER JOIN "UniqueVariantsInGenome" AS b
                     ON a."uniqueVariantId"=b."uniqueVariantId"
                     LEFT OUTER JOIN "VA_VariantNomenclature" AS c
                     ON c."uniqueVariantId"=b."uniqueVariantId"
                     LEFT OUTER JOIN "Genes" AS d
                     ON c."transcriptId"=d."mainTranscriptId"
                     WHERE a."uniqueVariantId"=', id, 'AND a."isMainTranscript"=\'TRUE\' AND d."symbol"<> \'NA\' ORDER BY a."date" DESC LIMIT 1') #obtenim les coordenades, el cromosoma, l'efecte i el gen de la nostra variant. 
        
    
        
        bp1<-dbGetQuery(con, bp1.1)
        }
        gen<-bp1[1,"symbol"]
        
        if (!is.na(bp1$validatedEffect)){
          #criteria["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&bp1[1,"validatedEffect"]=="missense (non-synonymous)"]<-1 #If the variant is a validated missense in NF1 or SPRED1 gene, BP1 = 1
          #criteria["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&bp1[1,"validatedEffect"]!="missense (non-synonymous)"]<-0 #If the variant is not a validated missense in NF1 or SPRED1 gene, BP1 = 0 in criteria matrix
          criteria["BP1",1][(bp1[1,"symbol"]!="NF1"&bp1[1,"symbol"]!="SPRED1")&(bp1[1,"validatedEffect"]=="nonsense (stop gain)"|bp1[1,"validatedEffect"]=="nonsense/splice"|bp1[1,"validatedEffect"]=="frameshift")]<-1 #If the variant is a validated null variant in any RASo gene (not NF1 nor SPRED1 genes), BP1 = 1
          criteria["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&(bp1[1,"validatedEffect"]=="nonsense (stop gain)"|bp1[1,"validatedEffect"]=="nonsense/splice"|bp1[1,"validatedEffect"]=="frameshift")]<-0 #If the variant is a validated null variant in any NF1 nor SPRED1 genes, BP1 = 0
          criteria["BP1",1][bp1[1,"validatedEffect"]=="missense (non-synonymous)"]<-0 #if variant is a validated missense, the value is 0
          criteria["PVS1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&(bp1[1,"validatedEffect"]=="nonsense (stop gain)"|bp1[1,"validatedEffect"]=="nonsense/splice"|bp1[1,"validatedEffect"]=="stopgain"|bp1[1,"validatedEffect"]=="frameshift")]<-1 #If the variant is a validated null variant in NF1 or SPRED1 gene, PVS1 = 1
          
        }
        
        else if (is.na(bp1$validatedEffect)){
          #criteria["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&bp1[1,"effect"]=="nonsynonymous SNV"]<-1 #si el gen es NF1 o SPRED1 i ?s una missense, s'assigna el valor 1 a la matriu criteria BP1
          #criteria["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&bp1[1,"effect"]!="nonsynonymous SNV"]<-0#si el gen es NF1 o SPRED1 i no  ?s una missense, s'assigna el valor 0 a la matriu criteria BP1
          criteria["BP1",1][(gen!="NF1"&gen!="SPRED1")&(bp1[1,"effect"]=="stopgain"|bp1[1,"effect"]=="stoploss"|bp1[1,"effect"]=="frameshift deletion")]<-1 #If the variant is a null variant in any RASo gene (not NF1 nor SPRED1 genes), BP1 = 1
          criteria["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&(bp1[1,"effect"]=="stopgain"|bp1[1,"effect"]=="stoploss"|bp1[1,"effect"]=="frameshift deletion")]<-0 #If the variant is a null variant in any NF1 nor SPRED1 genes, BP1 = 0
          criteria["BP1",1][bp1[1,"effect"]=="nonsynonymous SNV"]<-0 #if variant is a missense, the value is 0
          criteria["PVS1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&(bp1[1,"effect"]=="stopgain"|bp1[1,"effect"]=="stoploss"|bp1[1,"effect"]=="frameshift deletion")]<-1 #If the variant is a null variant in NF1 or SPRED1 gene, PVS1 = 1
        }
        
        criteria["PVS1",1][is.na(criteria["PVS1",1])]<-0 #we assign PVS1 = 0 if none of the previous conditions has been acomplished.
        
        
        ###PS1: Same amino acid change as a previously established pathogenic variant regardless of nucleotide change
        
        if (is.null(id)) {
          
          ## nothing to do here
          return(NULL)}
        
        else{
        ps1_aa<-paste('SELECT "cDNAAnnotation", "proteinAnnotation", "transcriptId"
                      FROM "VA_VariantsInTranscripts"
                      WHERE "isMainTranscript"=\'TRUE\' AND "uniqueVariantId"=', id, 'ORDER BY date DESC LIMIT 1')
       
        
        ps1_a<-dbGetQuery(con, ps1_aa) 
        }
        cdna<-unlist(strsplit(ps1_a$cDNAAnnotation,">")) #it separates the text string when it encounters ">"
        cdna<-paste(cdna[1], "%", sep="") #add "%" to the first element of the string, to search for any pattern that starts with cdna[1] followed by anything.
        
        if (is.null(id)) {
          
          ## nothing to do here
          return(NULL)}
        
        else{
        join_ps1<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification"
                         FROM "VA_VariantsInTranscripts"
                         LEFT OUTER JOIN "VA_CustomClassification"
                         ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                         WHERE "VA_VariantsInTranscripts"."transcriptId"=',ps1_a$transcriptId,'AND "VA_VariantsInTranscripts"."cDNAAnnotation" LIKE \'',cdna,'\' AND "VA_VariantsInTranscripts"."proteinAnnotation"=\'', ps1_a$proteinAnnotation, '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
        ps1.1<-dbGetQuery(con, join_ps1)
        }
        ps1<-sum(ps1.1$classification=="pat") #counts rows containing word PAT of synonymous variants
        criteria["PS1", 1][ps1>=1]<-1 
        
        
        if (gen=="SOS1"|gen=="SOS2"|gen=="BRAF"|gen=="RAF1"|gen=="MAP2K1"|gen=="MAP2K2"|gen=="HRAS"|gen=="NRAS"|gen=="KRAS"){
          if (ps1==0 & !is.na(ps1_a$proteinAnnotation)&!is.na(bp1[1,"effect"])&(bp1[1,"effect"]!="synonymous SNV"|bp1[1,"effect"]=="nonsynonymous SNV")){
            prot<-unlist(strsplit(ps1_a$proteinAnnotation,"p.")) # it separates the aa number from "p."
            proti<-prot[2] # it takes the second element, which is the informative one
            first_aa<-str_sub(proti,1,1) # the first aa is extracted
            num<-as.integer(str_sub(proti,2,str_length(proti)-1)) # it is taken from position 2 to -1, that is, the position of aa
            aa_canvi<-str_sub(proti, str_length(proti), str_length(proti))
            group_position<-c()
            if (gen=="SOS1"){ #if it is SOS1 gene
              functional_domain<-domain_groupSOS[num>=domain_groupSOS$SOS1.start&num<=domain_groupSOS$SOS1.end,] # it checks if the variant is within a known domain
              relative_position<-num-functional_domain$SOS1.start # it establishes the relative position within that domain
              group_position<-relative_position+functional_domain$SOS2.start # it looks at what position is equivalent to the other domain of the same group of genes, in that case, SOS2
            }
            
            else if (gen=="SOS2"){
              functional_domain<-domain_groupSOS[num>=domain_groupSOS$SOS2.start&num<=domain_groupSOS$SOS2.end,]
              relative_position<-num-functional_domain$SOS2.start
              group_position<-relative_position+functional_domain$SOS1.start
            }
            
            else if(gen=="BRAF"){
              functional_domain<-domain_groupRAF[num>=domain_groupRAF$BRAF.start&num<=domain_groupRAF$BRAF.end,]
              relative_position<-num-functional_domain$BRAF.start
              group_position<-relative_position+functional_domain$RAF1.start
            }
            else if (gen=="RAF1"){
              functional_domain<-domain_groupRAF[num>=domain_groupRAF$RAF1.start&num<=domain_groupRAF$RAF1.end,]
              relative_position<-num-functional_domain$RAF1.start
              group_position<-relative_position+functional_domain$BRAF.start
            }
            else if (gen=="MAP2K1"){
              functional_domain<-domain_groupMAPK[num>=domain_groupMAPK$MAP2K1.start&num<=domain_groupMAPK$MAP2K1.end,]
              relative_position<-num-functional_domain$MAP2K1.start
              group_position<-relative_position+functional_domain$MAP2K2.start}
            
            else if (gen=="MAP2K2"){
              functional_domain<-domain_groupMAPK[num>=domain_groupMAPK$MAP2K2.start&num<=domain_groupMAPK$MAP2K2.end,]
              relative_position<-num-functional_domain$MAP2K2.start
              group_position<-relative_position+functional_domain$MAP2K1.start
            }
            else if (gen=="NRAS"){
              functional_domain<-domain_groupRAS[num>=domain_groupRAS$NRAS.start&num<=domain_groupRAS$NRAS.end,]
              relative_position<-num-functional_domain$NRAS.start
              group_position<-c(relative_position+functional_domain$HRAS.start,relative_position+functional_domain$KRAS.start )
            }
            else if (gen=="HRAS"){
              functional_domain<-domain_groupRAS[num>=domain_groupRAS$HRAS.start&num<=domain_groupRAS$HRAS.end,]
              relative_position<-num-functional_domain$HRAS.start
              if (nrow(functional_domain)>0){
                if (functional_domain$NRAS.start==0){
                  group_position<-relative_position+functional_domain$KRAS.start
                }
                else if (functional_domain$NRAS.start!=0){
                  group_position<-c(relative_position+functional_domain$KRAS.start,relative_position+functional_domain$NRAS.start )
                }
              }}
            else if (gen=="KRAS"){
              functional_domain<-domain_groupRAS[num>=domain_groupRAS$KRAS.start&num<=domain_groupRAS$KRAS.end,]
              relative_position<-num-functional_domain$KRAS.start
              if (nrow(functional_domain)>0){
                if (functional_domain$NRAS.start==0){
                  group_position<-relative_position+functional_domain$HRAS.start
                }
                else if (functional_domain$NRAS.start!=0){
                  group_position<-c(relative_position+functional_domain$HRAS.start,relative_position+functional_domain$NRAS.start )
                }}}
            
            else if(length(group_position)!=0){
              group_position2<-paste("p.",first_aa, group_position, aa_canvi, sep="")
              join_ps1<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification"
                               FROM "VA_VariantsInTranscripts"
                               LEFT OUTER JOIN "VA_CustomClassification"
                               ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                               WHERE "VA_VariantsInTranscripts"."transcriptId"=',Transcripts_RASos[gen,1],'AND  "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', group_position2[1], '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,' ORDER BY "VA_CustomClassification".date LIMIT 1', sep="")
              ps1_analogs_groups<-dbGetQuery(con,join_ps1)
              
              # it compares now if the aa changed is the same as what we are looking for
              #analogs_groups<-str_sub(pm5_analogs_groups$proteinAnnotation, str_length(pm5_analogs_groups$proteinAnnotation), str_length(pm5_analogs_groups$proteinAnnotation))
              #loop to detect if analog groups and aa_change are equal and pathogenic PM1 =1, and also to see if aa is different and pathogenic, then PM5 = 1
              criteria["PS1",1][ps1_analogs_groups$classification=="pat"|ps1_analogs_groups$classification=="ppat"]<-1 # it assigns 1 if classification is PAT or PPAT
              if (length(group_position)==2){ # it looks if = 2 to consider group of genes 2 that are 3 genes.
                join_ps1<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification"
                                 FROM "VA_VariantsInTranscripts"
                                 LEFT OUTER JOIN "VA_CustomClassification"
                                 ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                                 WHERE "VA_VariantsInTranscripts"."transcriptId"=',Transcripts_RASos[gen,2],'AND  "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', group_position2[2], '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
                ps1_analogs_groups2<-dbGetQuery(con,join_ps1)
                criteria["PS1",1][ps1_analogs_groups2$classification=="pat"|ps1_analogs_groups2$classification=="ppat"]<-1 # it assigns 1 if classification is PAT or PPAT
              }
            }
          }}
        
        
        criteria["PS1", 1][is.na(criteria["PS1", 1])]<-0 # it assigns 0 if none of the above is PAT or PPAT.
        
        ###pm5 i bp8
        #PM5: Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before
        #BP8
        pm5<-0  
        bp8<-0
        if (is.null(id)) {
          
          ## nothing to do here
          return(NULL)}
        
        else{
        join_pm5<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification", "effect"
                         FROM "VA_VariantsInTranscripts"
                         LEFT OUTER JOIN "VA_CustomClassification"
                         ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                         WHERE "VA_VariantsInTranscripts"."transcriptId"=', ps1_a$transcriptId,'AND "VA_VariantsInTranscripts"."cDNAAnnotation" LIKE \'',cdna,'\' AND "VA_VariantsInTranscripts"."proteinAnnotation"<>\'', ps1_a$proteinAnnotation, '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
        
      
        
        pm5.1<-dbGetQuery(con, join_pm5)
    }
        if(!is.na(ps1_a$proteinAnnotation)& bp1$effect!="synonymous SNV"){
          prot_sense<-ps1_a$proteinAnnotation
          if(str_detect(prot_sense,"fs")==F){ #when it is different a frameshift
            if(str_detect(prot_sense,"_")==F){
              prot_cerca<-str_sub(prot_sense,1,str_length(prot_sense)-1)
              prot_cerca<-paste("\'",prot_cerca, "%","\'", sep="")
              
              if (is.null(id)) {
                
                ## nothing to do here
                return(NULL)}
              
              else{
              join_pm5b<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification", "effect"
                                FROM "VA_VariantsInTranscripts"
                                LEFT OUTER JOIN "VA_CustomClassification"
                                ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                                WHERE "VA_VariantsInTranscripts"."transcriptId"=', ps1_a$transcriptId,'AND "VA_VariantsInTranscripts"."proteinAnnotation" LIKE ',prot_cerca,'AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
              
             
              pm5.1b<-dbGetQuery(con, join_pm5b)
              }
              duplicatsb<-duplicated(pm5.1b$uniqueVariantId)
              pm5.2b<-NULL
              for (i in 1:nrow(pm5.1b)){ # it removes duplicate entries
                if (nrow(pm5.1b>0)){
                  if (duplicatsb[i]==FALSE){
                    pm5.3b<-pm5.1b[i,]
                    pm5.2b<-rbind(pm5.2b,pm5.3b)
                  }
                  else if (duplicatsb[i]==TRUE& pm5.1b$date[i]>pm5.2b$date[pm5.2b$uniqueVariantId==pm5.1b$uniqueVariantId[i]]){
                    pm5.2b$date[pm5.2b$uniqueVariantId==pm5.1b$uniqueVariantId[i]]<-pm5.1b$date[i]
                    pm5.2b$classification[pm5.2b$uniqueVariantId==pm5.1b$uniqueVariantId[i]]<-pm5.1b$classification[i]
                  }}
              }
              
              pm5<-sum(pm5.2b$classification=="pat", pm5.2b$classification=="ppat")
              pm5[is.na(pm5)]<-0
              bp8<-sum(pm5.2b$classification=="pol", pm5.2b$classification=="ppol")
            }}}
        pm5.2<-NULL
        duplicats<-duplicated(pm5.1$uniqueVariantId)
        
        
        for (i in 1:nrow(pm5.1)){ # it removes duplicate entries
          if (nrow(pm5.1>0)){
            if (duplicats[i]==FALSE){
              pm5.3<-pm5.1[i,]
              pm5.2<-rbind(pm5.2,pm5.3)
            }
            else if (duplicats[i]==TRUE& pm5.1$date[i]>pm5.2$date[pm5.2$uniqueVariantId==pm5.1$uniqueVariantId[i]]){
              pm5.2$date[pm5.2$uniqueVariantId==pm5.1$uniqueVariantId[i]]<-pm5.1$date[i]
              pm5.2$classification[pm5.2$uniqueVariantId==pm5.1$uniqueVariantId[i]]<-pm5.1$classification[i]
            }}
        }
        
        pm5<-sum(pm5,pm5.2$classification=="pat",pm5.2$classification=="ppat") #it counts rows which are pat / ppat having different amino acid
        bp8<-sum(bp8,pm5.2$classification=="pol",pm5.2$classification=="ppol") #it counts rows which are pol / ppol having different amino acid
        criteria[c("PM5_strong","PM5", "PM5_supporting"), 1][pm5>=2]<-c(1,0,0) #if any row is assigned, PM5 = 1 in the criteria array
        criteria[c("PM5_strong","PM5", "PM5_supporting"), 1][pm5==1]<-c(0,1,0)
        criteria[c("PM5_strong","PM5"), 1][pm5==0]<-c(0,0)
        criteria["BP8", 1][bp8>=1]<-1
        
        if (is.null(id)) {
          
          ## nothing to do here
          return(NULL)}
        
        else{
        
        if(length(group_position)!=0&pm5==0){
          group_position2<-paste("p.",first_aa, group_position,"_", sep="")
          join_pm5<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification", "effect"
                           FROM "VA_VariantsInTranscripts"
                           LEFT OUTER JOIN "VA_CustomClassification"
                           ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                           WHERE "VA_VariantsInTranscripts"."transcriptId"=',Transcripts_RASos[gen,1], 'AND  "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', group_position2[1], '\'AND "VA_VariantsInTranscripts"."proteinAnnotation"<>\'', ps1_a$proteinAnnotation, '\' AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,' ORDER BY "VA_CustomClassification".date', sep="")
      
          }
        }
          if (is.null(id)) {
            
            ## nothing to do here
            return(NULL)}
          
          else{
          pm5_analogs_groups<-dbGetQuery(con,join_pm5)
          #loop to detect if analog groups and aa_change are equal and pathogenic PM1 = 1, and also to see if aa is different and pathogenic, then PM5 = 1
          pm5<-sum(pm5_analogs_groups$classification=="pat"|pm5_analogs_groups$classification=="ppat")
          bp8<-sum(pm5_analogs_groups$classification=="pol"|pm5_analogs_groups$classification=="ppol")
          criteria["PM5_supporting",1][pm5>=1]<-1 # PM5 = 1 if the variant is PAT or PPAT
          criteria["BP8", 1][bp8>=1]<-1
          if (length(group_position)==2){ # it looks if = 2 to consider group of genes 2 that are 3 genes.
            join_pm5<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification"
                             FROM "VA_VariantsInTranscripts"
                             LEFT OUTER JOIN "VA_CustomClassification"
                             ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                             WHERE "VA_VariantsInTranscripts"."transcriptId"=',Transcripts_RASos[gen,2],'AND "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', group_position2[1], '\'AND  "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', group_position2[2], '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
          
            }
            
            pm5_analogs_groups2<-dbGetQuery(con,join_pm5)
            criteria["PM5_supporting",1][pm5_analogs_groups2$classification=="pat"|pm5_analogs_groups2$classification=="ppat"]<-1 #assignem 1 si ?s pat o ppat
            bp8<-sum(pm5_analogs_groups2$classification=="pol"|pm5_analogs_groups2$classification=="ppol")
            bp8<-sum(pm5_analogs_groups$classification=="pol"|pm5_analogs_groups$classification=="ppol")
          }
        
        
        pm5[is.null(pm5)]<-0
        criteria["PM5_supporting", 1][pm5==0]<-0 #si no hi ha cap fila no es compleix
        criteria["BP8", 1][is.na(criteria["BP8", 1])]<-0
        
        ## PM1: Located in a mutational hot spot and/or critical and well-established functional domain (e.g., active site of an enzyme) without benign variation
        
        if (!is.na(ps1_a$proteinAnnotation)& (gen=="NF1"|gen=="SPRED1")){
          prot<-unlist(strsplit(ps1_a$proteinAnnotation,"p.")) # it separates the aa number from "p."
          proti<-prot[2] #it takes the second element, which is the informative one
          first_aa<-str_sub(proti,1,1) # the first aa is extracted
          if(str_detect(proti,"fs")==F){ # when it is different a frameshift
            if(str_detect(proti,"_")==F){
              num<-as.integer(str_sub(proti,2,str_length(proti)-1))
            }
            else if(str_detect(proti,"_")==T){ #it takes inframe positions
              sep<-unlist(strsplit(proti, "_"))
              num<-sep[1]
            }
            
          }
          else{ #if it is a frameshift
            num<-as.integer(str_sub(proti,2,str_length(proti)-2))
          }
        }
        
        #hot spots and domains of RASo genes
        BRAF<-c(238:286,439:474,531,459:474,594:627)
        HRAS<-c(12,13,14,58,59,60,61,62,63)
        KRAS<-c(12,13,14,58,59,60,61,62,63)
        MAP2K1<-c(43:61,124:134)
        MAP2K2<-c(47:65,128:138)
        PTPN11<-c(308,4,7,8,9,58:63,69:77,247,251,255,256,258,261,265,278:281,284)
        RAF1<-c(251:266)
        SHOC2<-c(2)
        SOS1<-c(269,552,420:500)
        SPRED1<-c(6:123, 233-285, 334:442)
        NF1<-c(543:909,1095:1176, 1170:1218, 1198:1530, 1510:1570, 1560:1696,1713:1816, 2619:2719)
        
        hot_spots<-list(BRAF,HRAS, KRAS,MAP2K1,MAP2K2, PTPN11, RAF1,SHOC2, SOS1,SPRED1, NF1)
        names(hot_spots)<-c("BRAF", "HRAS", "KRAS", "MAP2K1", "MAP2K2", "PTPN11", "RAF1", "SHOC2", "SOS1", "SPRED1", "NF1")
        suma<-0
        if (!is.na(ps1_a$proteinAnnotation)){
          hot_sp<-hot_spots[[gen]]==num
          suma<-sum(hot_sp==T)
        }
        criteria["PM1",1][suma==1]<-1 # PM1 = 1 if it is PAT or PPAT
        criteria["PM1",1][suma==0|is.na(suma)]<-0
        
        ##BA1 , BS1 i PM2
        #BA1: Allele frequency is above 5% in Exome Sequencing Project, 1000 Genomes, or ExAC.
        #BS1: Allele frequency is greater than expected for these disorders
        #PM2: Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes, or ExAC
       
         if (is.null(id)) {
          
          ## nothing to do here
          return(NULL)}
        
        else{
        freq<-paste('SELECT "PopFreqMax" FROM "VA_Frequencies" WHERE "uniqueVariantId"=', id,'ORDER BY "date" DESC LIMIT 1') #agafem la informaci? de frequencia m?s actualitzada de la variant d'interes
        BA1<- dbGetQuery(con, freq)
        criteria["BA1",1][BA1>=0.0005]<-1 # BA1 = 1 if the poblacional frequency is greater than 0.0005
        criteria["BA1",1][BA1<0.0005|is.na(BA1)]<-0 # BA1 = 0 if the poblacional frequency is lower than 0.0005
        criteria["BS1",1][BA1>=0.00025&BA1<0.0005]<-1 # BS1 = 1 if the poblacional frequency is greater than 0.00025
        criteria["BS1",1][BA1<0.00025|is.na(BA1)|BA1>=0.0005]<-0 # BS1 = 1 if the poblacional frequency is lower than 0.00025
        criteria["BS1_supporting",1][!is.na(BA1)]<-0
        criteria["PM2", 1][BA1==0|is.na(BA1)]<-1
        criteria["PM2", 1][BA1!=0&!is.na(BA1)]<-0
        
        }
        
        
        ##PM4 i BP3
        #PM4: Protein length changes as a result of in- frame deletions/insertions in a nonrepeat region or stop-loss variants
        #BP3: In-frame deletions/insertions in a repetitive region without a known function
        
        criteria["PM4", 1][gen!="NF1"|gen!="SPRED1"]<-0
        criteria["BP3", 1][gen!="NF1"|gen!="SPRED1"]<-0
       
        if (is.null(id)) {
          
          ## nothing to do here
          return(NULL)}
        
        else {
        
         if (gen=="NF1"|gen=="SPRED1"){
          pm4.1<-paste('SELECT b."start", b."end", b."chr", "validatedEffect", "strand"
                   FROM "VA_VariantsInTranscripts" AS a
                   LEFT OUTER JOIN "UniqueVariantsInGenome" AS b
                   ON a."uniqueVariantId"=b."uniqueVariantId"
                   LEFT OUTER JOIN "VA_VariantNomenclature" AS c
                   ON c."uniqueVariantId"=b."uniqueVariantId"
                   LEFT OUTER JOIN "TranscriptsPosition" as d
                   ON d."transcriptId"=c."transcriptId"
                   WHERE a."uniqueVariantId"=', id, 'AND a."isMainTranscript"=\'TRUE\' ORDER BY c."date" DESC LIMIT 1') #it gets information about the coordinates of start, end, effect and the string in which the gene is.
          pm4<-dbGetQuery(con, pm4.1)
          criteria["PM4", 1][pm4$validatedEffect!="frameshift"]<-0 # PM4 = 0 if variant is not a frameshift
          criteria["BP3", 1][!pm4$validatedEffect %in% "in frame"]<-0 #BP3 = 0 if it is in frames
          if(!is.na(pm4$validatedEffect)){
            result<-NULL
          }
          else if(!is.na(pm4$validatedEffect)){
            if (pm4$validatedEffect=="frameshift"|pm4$validatedEffect %in% "in frame"){ # if it equals a frameshift or inframe
              library(RMySQL)
              query<-paste('SELECT * FROM rmsk WHERE genoStart< ',pm4$start, ' AND genoEnd>',pm4$end,' AND genoName=\'chr',pm4$chr,'\' ', sep="") #agafem informacio del UCSC RepeatMasker
              pm4.2<- dbGetQuery(my_connection, query)
              result<-nrow(pm4.2$strand[pm4.2$strand==pm4$strand]) # it looks they are in the same DNA strand
            }
            criteria["PM4", 1][pm4$validatedEffect=="frameshift"&is.null(result)]<-1 # frameshift and there are not results in a none repetitive region, PM4 = 1
            criteria["PM4", 1][pm4$validatedEffect=="frameshift"&!is.null(result)]<-0 # frameshift and there are a result in a repetitive region, PM4 = 0
            criteria["BP3", 1][pm4$validatedEffect %in% "in frame" & !is.null (result)]<-1 # in frame and there are a result in a repetitive region, BP3 = 1
            criteria["BP3", 1][pm4$validatedEffect %in% "in frame" & is.null (result)]<-0 # in frame and and there are not results in a none repetitive region, BP3 = 0
          }
        }
        
        }
        
        ##PP2: Missense variant in a gene that has a low rate of benign missense variation and in which missense variants are a common mechanism of variation
        
        if (is.null(id)) {
          
          ## nothing to do here
          return(NULL)}
        
        else {
        pp2.1<- paste('SELECT "effect","validatedEffect"
                      FROM "VA_VariantsInTranscripts"
                      LEFT OUTER JOIN "VA_VariantNomenclature"
                      ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_VariantNomenclature"."uniqueVariantId"
                      WHERE "VA_VariantsInTranscripts"."uniqueVariantId"=', id, 'AND "isMainTranscript"=\'TRUE\' ORDER BY "VA_VariantsInTranscripts"."date" DESC LIMIT 1') # AND "validatedEffect"<> \'<NA>\' li diem que validated effect no sigui NA
        
        }
        
        pp2<-dbGetQuery(con, pp2.1)
        efecte_validat<-pp2$validatedEffect[!is.na(pp2$validatedEffect)]
        efecte_novalidat<-pp2$effect[!is.na(pp2$effect)]
        criteria["PP2",1][length(efecte_validat)==0&length(efecte_novalidat)==0]<-0
        if( length(efecte_validat)>0){
          if(efecte_validat[1]=="missense (non-synonymous)"){
            criteria["PP2", 1]<-1
          }
          else if (efecte_validat[1]!="missense (non-synonymous)"){
            criteria["PP2", 1]<-0
          }}
        else if (length(efecte_validat)==0 &length(efecte_novalidat)>0){ #it considers the ones that have no value in validated effect and look at the effect field
          if (efecte_novalidat[1]=="nonsynonymous SNV"){
            criteria["PP2", 1]<-1
          }
          else if(efecte_novalidat[1]!="nonsynonymous SNV"){
            criteria["PP2", 1]<-0
          }}
        
        ##PP3 & BP4
        #PP3: Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact,etc)
        #BP4: Multiple lines of computational evidence suggest no impact on gene or gene product (conservation, evolutionary, splicing, etc)
        
        
        if (is.null(id)) {
          
          ## nothing to do here
          return(NULL)}
        
        else{
        pp3.1<- paste('SELECT "SIFT_pred", "Polyphen2_HVAR_pred","MutationTaster_pred","PROVEAN_pred", "CADD_phred" 
                      FROM "VA_VariantsInTranscripts"
                      LEFT OUTER JOIN "VA_InSilicoPathogenicity" AS a
                      ON "VA_VariantsInTranscripts"."uniqueVariantId"= a."uniqueVariantId"
                      WHERE "VA_VariantsInTranscripts"."uniqueVariantId"=', id, 'AND
                      "isMainTranscript"=\'TRUE\' ORDER BY a."date" DESC LIMIT 1')
      
        }
        
        pp3<-dbGetQuery(con, pp3.1)
        pp3$CADD_phred[pp3$CADD_phred>19]<-"D"
        pp3$CADD_phred[pp3$CADD_phred<10]<-"N"
        pp3$Polyphen2_HVAR_pred[pp3$Polyphen2_HVAR_pred=="B"]<-"N"
        pp3$Polyphen2_HVAR_pred[pp3$Polyphen2_HVAR_pred=="P"]<-"D"
        pp3$MutationTaster_pred[pp3$MutationTaster_pred=="A"]<-"D"
        pp3$MutationTaster_pred[pp3$MutationTaster_pred=="P"]<-"N"
        pp3$SIFT_pred[pp3$SIFT_pred=="T"]<-"N"
        pp3<-t(pp3)
        table<-table(pp3[,1])
        criteria["PP3",1][table["D"]==5]<-1
        criteria["PP3",1][table["D"]<=4|table["N"]==5]<-0
        criteria["BP4",1][table["N"]==5]<-1
        criteria["BP4",1][table["N"]<=4|table["D"]==5]<-0
        
        
        ##BS2 & PS4
        #BS2: Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked (hemizygous) disorder with full penetrance expected at an early age
        #PS4: The prevalence of the variant in affected individuals is significantly increased compared with the prevalence in controls
        
        if (is.null(id)) {
          
          ## nothing to do here
          return(NULL)}
        
        else{
        bs2.1<-paste('SELECT "pedigreeProgenyId", "startDate",f."name"
                     FROM "UniqueVariants" AS a
                     LEFT OUTER JOIN "Variants" AS b
                     ON b."uniqueVariantId"=a."id"
                     LEFT OUTER JOIN "Samples" AS c
                     ON c."id"=b."sampleId"
                     LEFT OUTER JOIN "Individual" AS d
                     ON d."id"=c."individualId"
                     LEFT OUTER JOIN "Analysis" AS e
                     ON b."analysisId"=e."id"
                     LEFT OUTER JOIN "GeneCascade" AS f
                     ON f."id" = e."geneCascadeId"
                     WHERE a."id"=', id,'')
        
        
        }
        
        bs2<-dbGetQuery(con, bs2.1)
        independent_control<-sum(bs2$name!="NF1_SPRED1"&bs2$name!="Rasopathies"|is.na(bs2$name))
        affected_independently<-sum(bs2$name=="NF1_SPRED1"&!is.na(bs2$name)|bs2$name=="Rasopathies"&!is.na(bs2$name))
        
        criteria[c("BS2","PS4_strong", "PS4_moderate", "PS4_supporting"),1][independent_control>=3]<-c(1,0,0,0)
        criteria[c("PS4_strong", "PS4_moderate", "PS4_supporting"),1][affected_independently==0]<-c(0,0,0)
        if (independent_control==0){
          criteria[c("BS2","PS4_strong","PS4_moderate", "PS4_supporting"),1][affected_independently>=6]<-c(0,1,0,0)
          criteria[c("BS2","PS4_strong","PS4_moderate", "PS4_supporting"),1][affected_independently==4|affected_independently==5]<-c(0,0,1,0)
          criteria[c("BS2","PS4_strong","PS4_moderate", "PS4_supporting"),1][affected_independently==3|affected_independently==2]<-c(0,0,0,1)
          criteria[c("BS2","PS4_strong","PS4_moderate", "PS4_supporting"),1][affected_independently==0]<-c(0,0,0,0)
        }
        else if (independent_control<3&independent_control>0){
        }
        
        
        ## Addatpation PM2 & BS1:
        #PM2:
        #Bs1:
        criteria["PM2",1][criteria["PM2",1]==1&independent_control>100]<-0 # independently that Pandora says variant is not described in population databases, if the variant is detected in more than 100 samples analyzed in Pandora and not showing RASopathy-like phentoype, and then used as control individuals, PM2 = 0
        inhousefreq<-paste('SELECT "freqs","inHouseFrequency","inHouseNumber" FROM "VA_InHouseFrequencies" WHERE "uniqueVariantId"=', id,'ORDER BY "date" DESC LIMIT 1') # it takes the most updated frequency information of the variant of interest
        inhouse<-dbGetQuery(con,inhousefreq)
        
        if (is.null(inhousefreq)) {
          
          ## nothing to do here
          return(NULL)
        }
        else {
        if (is.na(BA1$PopFreqMax)){
          if (is.null(independent_control)) {
            
            ## nothing to do here
            return(NULL)
          }
          criteria["BS1",1][independent_control>=1000&inhouse$inHouseFrequency[1]>=0.30]<-1 #independently that Pandora says variant is not described in population databases, if the variant is detected in more than 1000 samples analyzed in Pandora and not showing RASopathy-like phentoype, and then used as control individuals or in house-frequency >0,30, BS1 = 1
          criteria["BS1_supporting",1][independent_control>500&independent_control<1000&inhouse$inHouseFrequency[1]>0.10&inhouse$inHouseFrequency[1]<0.30]<-1
          criteria["BS1_supporting",1][is.na(criteria["BS1_supporting",1])]<-0
        } 
        }
        
        ###PS2 - de novo cases reported
        if (is.null(denovo_confirmed | is.null (input$denovo_confirmed))) {
          
          return(NULL)}
        else {
        criteria[c("PS2_veryStrong", "PS2"),1][denovo_confirmed>=2]<-c(1,0)
        criteria[c("PS2_veryStrong", "PS2"),1][denovo_confirmed==1 & denovo_noconfirmed>=2]<-c(1,0)
        criteria[c("PS2_veryStrong", "PS2"),1][denovo_confirmed==1 & denovo_noconfirmed<2]<-c(0,1)
        criteria[c("PS2_veryStrong", "PS2"),1][denovo_confirmed==0]<-c(0,0)
        
        }
        
        ###PM6  - de novo cases reported no confirmed
      
          if (is.null(denovo_confirmed | is.null (input$denovo_confirmed))) {
          
          return(NULL)}
        else {
          
          criteria[c("PM6_veryStrong","PM6_strong","PM6"),1][denovo_confirmed==0 & denovo_noconfirmed>3]<-c(1,0,0)
          criteria[c("PM6_veryStrong","PM6_strong","PM6"),1][denovo_confirmed==0 & denovo_noconfirmed==3 | denovo_noconfirmed==2]<-c(0,1,0)
          criteria[c("PM6_veryStrong","PM6_strong","PM6"),1][denovo_confirmed==0 & denovo_noconfirmed==1]<-c(0,0,1)
          criteria[c("PM6_veryStrong","PM6_strong","PM6"),1][denovo_confirmed==0 & denovo_noconfirmed==0]<-c(0,0,0)
        }
        
        ###PP1 - cosegregation cases reported
        if (is.null(cosegregation)) {
          
          ## nothing to do here
          return(NULL)}
          else {
        criteria[c("PP1_strong", "PP1_moderate", "PP1_supporting"),1][cosegregation>=7]<-c(1,0,0)
        criteria[c("PP1_strong", "PP1_moderate", "PP1_supporting"),1][cosegregation==5|cosegregation==6]<-c(0,1,0)
        criteria[c("PP1_strong", "PP1_moderate", "PP1_supporting"),1][cosegregation==3|cosegregation==4]<-c(0,0,1)
        criteria[c("PP1_strong", "PP1_moderate", "PP1_supporting"),1][cosegregation<3]<-c(0,0,0)
                }
      
        
        ###PS3 - Functional studies
        
        if (is.null(Functional_evidence)) {
          
          ## nothing to do here
          return(NULL)}
        else {
        criteria[c("PS3"),1]<-c(Functional_evidence[1])
              }
        
        ###PP5, BP6
        if (is.null(PPAT_evidence) | is.null (input$PPOL_evidence)) {
          
          ## nothing to do here
          return(NULL)}
        else {
        criteria[c("PP5", "BP6"),1]<-c(PPAT_evidence[1], PPOL_evidence[2])
        
        }
        
    ACMG<-data.frame(AMGC=c("PS1", "PS2_veryStrong", "PS2", "PS3", "PS4_strong", "PS4_moderate", "PS4_supporting", "PM5_strong",  "PM6_veryStrong", "PM6","PM1", "PM2", "PM4",  "PP1_strong", "PP1_moderate", "PP1_supporting", "PP2", "PP3","PP5", "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7", "BS1_supporting", "PM6_strong", "PM5", "PM5_supporting", "BP8", "PVS1"))
    if (is.null(criteria)) {
      
      ## nothing to do here
      return(NULL)}
    else {
      criteria <- cbind.data.frame(ACMG, criteria) 
      criteria <- criteria[criteria$criteria == 1,]
    criteria_filtered <- data.frame(criteria[!is.na(criteria$criteria),])
      
   return(criteria_filtered)
    }
    
      }

  ### FINAL CLASSIFICATION

  Final_classificationB<- function(id = input$id, con = con, denovo_noconfirmed = input$denovo_noconfirmed,
                                     denovo_confirmed = input$denovo_confirmed, cosegregation = input$cosegregation,
                                     PPAT_evidence = as.numeric(as.character(input$PPAT_evidence)), 
                                     PPOL_evidence = as.numeric(as.character(input$PPOL_evidence)),
                                   Functional_evidence = as.numeric(as.character(input$Functional_evidence))){

    if (is.null(input$id) | is.null (input$denovo_noconfirmed) | is.null (input$denovo_confirmed) |
        is.null(input$cosegregation) | is.null (input$PPAT_evidence) |
        is.null (input$PPOL_evidence) | is.null(input$Functional_evidence) | input$go == 0) {
      
      ## nothing to do here
      return(NULL)
    }
    
    #dataframe to fill with all criteria (each criteria in one line)
    criteria<-data.frame(criteria=c(rep(NA, 37)), row.names = c("PS1", "PS2_veryStrong", "PS2", "PS3", "PS4_strong", "PS4_moderate", "PS4_supporting", "PM5_strong",  "PM6_veryStrong", "PM6","PM1", "PM2", "PM4",  "PP1_strong", "PP1_moderate", "PP1_supporting", "PP2", "PP3","PP5", "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7", "BS1_supporting", "PM6_strong", "PM5", "PM5_supporting", "BP8", "PVS1"))
    
    
    group_position<-NULL
    num<-NULL
    ### BP1 & PVS1 cirteria. 
    #BP1: Missense variant in a gene for which primarily truncating variants are known to cause the disease
    #PVS1: NUll variant. Only applicable to NF1 and SPRED1
    if (is.null(id)) {
      
      ## nothing to do here
      return(NULL)}
    
    else{
      
      bp1.1<-paste('SELECT "start", "end", "chr", "validatedEffect", "effect", "symbol"
                   FROM "VA_VariantsInTranscripts" AS a
                   LEFT OUTER JOIN "UniqueVariantsInGenome" AS b
                   ON a."uniqueVariantId"=b."uniqueVariantId"
                   LEFT OUTER JOIN "VA_VariantNomenclature" AS c
                   ON c."uniqueVariantId"=b."uniqueVariantId"
                   LEFT OUTER JOIN "Genes" AS d
                   ON c."transcriptId"=d."mainTranscriptId"
                   WHERE a."uniqueVariantId"=', id, 'AND a."isMainTranscript"=\'TRUE\' AND d."symbol"<> \'NA\' ORDER BY a."date" DESC LIMIT 1') #obtenim les coordenades, el cromosoma, l'efecte i el gen de la nostra variant. 
      
      
      
      bp1<-dbGetQuery(con, bp1.1)
    }
    gen<-bp1[1,"symbol"]
    
    if (!is.na(bp1$validatedEffect)){
      #criteria["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&bp1[1,"validatedEffect"]=="missense (non-synonymous)"]<-1 #If the variant is a validated missense in NF1 or SPRED1 gene, BP1 = 1
      #criteria["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&bp1[1,"validatedEffect"]!="missense (non-synonymous)"]<-0 #If the variant is not a validated missense in NF1 or SPRED1 gene, BP1 = 0 in criteria matrix
      criteria["BP1",1][(bp1[1,"symbol"]!="NF1"&bp1[1,"symbol"]!="SPRED1")&(bp1[1,"validatedEffect"]=="nonsense (stop gain)"|bp1[1,"validatedEffect"]=="nonsense/splice"|bp1[1,"validatedEffect"]=="frameshift")]<-1 #If the variant is a validated null variant in any RASo gene (not NF1 nor SPRED1 genes), BP1 = 1
      criteria["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&(bp1[1,"validatedEffect"]=="nonsense (stop gain)"|bp1[1,"validatedEffect"]=="nonsense/splice"|bp1[1,"validatedEffect"]=="frameshift")]<-0 #If the variant is a validated null variant in any NF1 nor SPRED1 genes, BP1 = 0
      criteria["BP1",1][bp1[1,"validatedEffect"]=="missense (non-synonymous)"]<-0 #if variant is a validated missense, the value is 0
      criteria["PVS1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&(bp1[1,"validatedEffect"]=="nonsense (stop gain)"|bp1[1,"validatedEffect"]=="nonsense/splice"|bp1[1,"validatedEffect"]=="stopgain"|bp1[1,"validatedEffect"]=="frameshift")]<-1 #If the variant is a validated null variant in NF1 or SPRED1 gene, PVS1 = 1
      
    }
    
    else if (is.na(bp1$validatedEffect)){
      #criteria["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&bp1[1,"effect"]=="nonsynonymous SNV"]<-1 #si el gen es NF1 o SPRED1 i ?s una missense, s'assigna el valor 1 a la matriu criteria BP1
      #criteria["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&bp1[1,"effect"]!="nonsynonymous SNV"]<-0#si el gen es NF1 o SPRED1 i no  ?s una missense, s'assigna el valor 0 a la matriu criteria BP1
      criteria["BP1",1][(gen!="NF1"&gen!="SPRED1")&(bp1[1,"effect"]=="stopgain"|bp1[1,"effect"]=="stoploss"|bp1[1,"effect"]=="frameshift deletion")]<-1 #If the variant is a null variant in any RASo gene (not NF1 nor SPRED1 genes), BP1 = 1
      criteria["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&(bp1[1,"effect"]=="stopgain"|bp1[1,"effect"]=="stoploss"|bp1[1,"effect"]=="frameshift deletion")]<-0 #If the variant is a null variant in any NF1 nor SPRED1 genes, BP1 = 0
      criteria["BP1",1][bp1[1,"effect"]=="nonsynonymous SNV"]<-0 #if variant is a missense, the value is 0
      criteria["PVS1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&(bp1[1,"effect"]=="stopgain"|bp1[1,"effect"]=="stoploss"|bp1[1,"effect"]=="frameshift deletion")]<-1 #If the variant is a null variant in NF1 or SPRED1 gene, PVS1 = 1
    }
    
    criteria["PVS1",1][is.na(criteria["PVS1",1])]<-0 #we assign PVS1 = 0 if none of the previous conditions has been acomplished.
    
    ###PS1: Same amino acid change as a previously established pathogenic variant regardless of nucleotide change
    
    
    if (is.null(id)) {
      
      ## nothing to do here
      return(NULL)}
    
    else{
      ps1_aa<-paste('SELECT "cDNAAnnotation", "proteinAnnotation", "transcriptId"
                    FROM "VA_VariantsInTranscripts"
                    WHERE "isMainTranscript"=\'TRUE\' AND "uniqueVariantId"=', id, 'ORDER BY date DESC LIMIT 1')
      
      
      ps1_a<-dbGetQuery(con, ps1_aa) 
    }
    cdna<-unlist(strsplit(ps1_a$cDNAAnnotation,">")) #it separates the text string when it encounters ">"
    cdna<-paste(cdna[1], "%", sep="") #add "%" to the first element of the string, to search for any pattern that starts with cdna[1] followed by anything.
    
    if (is.null(id)) {
      
      ## nothing to do here
      return(NULL)}
    
    else{
      join_ps1<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification"
                       FROM "VA_VariantsInTranscripts"
                       LEFT OUTER JOIN "VA_CustomClassification"
                       ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                       WHERE "VA_VariantsInTranscripts"."transcriptId"=',ps1_a$transcriptId,'AND "VA_VariantsInTranscripts"."cDNAAnnotation" LIKE \'',cdna,'\' AND "VA_VariantsInTranscripts"."proteinAnnotation"=\'', ps1_a$proteinAnnotation, '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
      ps1.1<-dbGetQuery(con, join_ps1)
    }
    ps1<-sum(ps1.1$classification=="pat") #counts rows containing word PAT of synonymous variants
    criteria["PS1", 1][ps1>=1]<-1 
    
    
    if (gen=="SOS1"|gen=="SOS2"|gen=="BRAF"|gen=="RAF1"|gen=="MAP2K1"|gen=="MAP2K2"|gen=="HRAS"|gen=="NRAS"|gen=="KRAS"){
      if (ps1==0 & !is.na(ps1_a$proteinAnnotation)&!is.na(bp1[1,"effect"])&(bp1[1,"effect"]!="synonymous SNV"|bp1[1,"effect"]=="nonsynonymous SNV")){
        prot<-unlist(strsplit(ps1_a$proteinAnnotation,"p.")) # it separates the aa number from "p."
        proti<-prot[2] # it takes the second element, which is the informative one
        first_aa<-str_sub(proti,1,1) # the first aa is extracted
        num<-as.integer(str_sub(proti,2,str_length(proti)-1)) # it is taken from position 2 to -1, that is, the position of aa
        aa_canvi<-str_sub(proti, str_length(proti), str_length(proti))
        group_position<-c()
        if (gen=="SOS1"){ #if it is SOS1 gene
          functional_domain<-domain_groupSOS[num>=domain_groupSOS$SOS1.start&num<=domain_groupSOS$SOS1.end,] # it checks if the variant is within a known domain
          relative_position<-num-functional_domain$SOS1.start # it establishes the relative position within that domain
          group_position<-relative_position+functional_domain$SOS2.start # it looks at what position is equivalent to the other domain of the same group of genes, in that case, SOS2
        }
        
        else if (gen=="SOS2"){
          functional_domain<-domain_groupSOS[num>=domain_groupSOS$SOS2.start&num<=domain_groupSOS$SOS2.end,]
          relative_position<-num-functional_domain$SOS2.start
          group_position<-relative_position+functional_domain$SOS1.start
        }
        
        else if(gen=="BRAF"){
          functional_domain<-domain_groupRAF[num>=domain_groupRAF$BRAF.start&num<=domain_groupRAF$BRAF.end,]
          relative_position<-num-functional_domain$BRAF.start
          group_position<-relative_position+functional_domain$RAF1.start
        }
        else if (gen=="RAF1"){
          functional_domain<-domain_groupRAF[num>=domain_groupRAF$RAF1.start&num<=domain_groupRAF$RAF1.end,]
          relative_position<-num-functional_domain$RAF1.start
          group_position<-relative_position+functional_domain$BRAF.start
        }
        else if (gen=="MAP2K1"){
          functional_domain<-domain_groupMAPK[num>=domain_groupMAPK$MAP2K1.start&num<=domain_groupMAPK$MAP2K1.end,]
          relative_position<-num-functional_domain$MAP2K1.start
          group_position<-relative_position+functional_domain$MAP2K2.start}
        
        else if (gen=="MAP2K2"){
          functional_domain<-domain_groupMAPK[num>=domain_groupMAPK$MAP2K2.start&num<=domain_groupMAPK$MAP2K2.end,]
          relative_position<-num-functional_domain$MAP2K2.start
          group_position<-relative_position+functional_domain$MAP2K1.start
        }
        else if (gen=="NRAS"){
          functional_domain<-domain_groupRAS[num>=domain_groupRAS$NRAS.start&num<=domain_groupRAS$NRAS.end,]
          relative_position<-num-functional_domain$NRAS.start
          group_position<-c(relative_position+functional_domain$HRAS.start,relative_position+functional_domain$KRAS.start )
        }
        else if (gen=="HRAS"){
          functional_domain<-domain_groupRAS[num>=domain_groupRAS$HRAS.start&num<=domain_groupRAS$HRAS.end,]
          relative_position<-num-functional_domain$HRAS.start
          if (nrow(functional_domain)>0){
            if (functional_domain$NRAS.start==0){
              group_position<-relative_position+functional_domain$KRAS.start
            }
            else if (functional_domain$NRAS.start!=0){
              group_position<-c(relative_position+functional_domain$KRAS.start,relative_position+functional_domain$NRAS.start )
            }
          }}
        else if (gen=="KRAS"){
          functional_domain<-domain_groupRAS[num>=domain_groupRAS$KRAS.start&num<=domain_groupRAS$KRAS.end,]
          relative_position<-num-functional_domain$KRAS.start
          if (nrow(functional_domain)>0){
            if (functional_domain$NRAS.start==0){
              group_position<-relative_position+functional_domain$HRAS.start
            }
            else if (functional_domain$NRAS.start!=0){
              group_position<-c(relative_position+functional_domain$HRAS.start,relative_position+functional_domain$NRAS.start )
            }}}
        
        else if(length(group_position)!=0){
          group_position2<-paste("p.",first_aa, group_position, aa_canvi, sep="")
          join_ps1<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification"
                           FROM "VA_VariantsInTranscripts"
                           LEFT OUTER JOIN "VA_CustomClassification"
                           ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                           WHERE "VA_VariantsInTranscripts"."transcriptId"=',Transcripts_RASos[gen,1],'AND  "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', group_position2[1], '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,' ORDER BY "VA_CustomClassification".date LIMIT 1', sep="")
          ps1_analogs_groups<-dbGetQuery(con,join_ps1)
          
          # it compares now if the aa changed is the same as what we are looking for
          #analogs_groups<-str_sub(pm5_analogs_groups$proteinAnnotation, str_length(pm5_analogs_groups$proteinAnnotation), str_length(pm5_analogs_groups$proteinAnnotation))
          #loop to detect if analog groups and aa_change are equal and pathogenic PM1 =1, and also to see if aa is different and pathogenic, then PM5 = 1
          criteria["PS1",1][ps1_analogs_groups$classification=="pat"|ps1_analogs_groups$classification=="ppat"]<-1 # it assigns 1 if classification is PAT or PPAT
          if (length(group_position)==2){ # it looks if = 2 to consider group of genes 2 that are 3 genes.
            join_ps1<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification"
                             FROM "VA_VariantsInTranscripts"
                             LEFT OUTER JOIN "VA_CustomClassification"
                             ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                             WHERE "VA_VariantsInTranscripts"."transcriptId"=',Transcripts_RASos[gen,2],'AND  "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', group_position2[2], '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
            ps1_analogs_groups2<-dbGetQuery(con,join_ps1)
            criteria["PS1",1][ps1_analogs_groups2$classification=="pat"|ps1_analogs_groups2$classification=="ppat"]<-1 # it assigns 1 if classification is PAT or PPAT
          }
        }
      }}
    
    
    criteria["PS1", 1][is.na(criteria["PS1", 1])]<-0 # it assigns 0 if none of the above is PAT or PPAT.
    
    ###pm5 i bp8
    #PM5: Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before
    #BP8
    pm5<-0  
    bp8<-0
    if (is.null(id)) {
      
      ## nothing to do here
      return(NULL)}
    
    else{
      join_pm5<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification", "effect"
                       FROM "VA_VariantsInTranscripts"
                       LEFT OUTER JOIN "VA_CustomClassification"
                       ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                       WHERE "VA_VariantsInTranscripts"."transcriptId"=', ps1_a$transcriptId,'AND "VA_VariantsInTranscripts"."cDNAAnnotation" LIKE \'',cdna,'\' AND "VA_VariantsInTranscripts"."proteinAnnotation"<>\'', ps1_a$proteinAnnotation, '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
      
      
      
      pm5.1<-dbGetQuery(con, join_pm5)
    }
    if(!is.na(ps1_a$proteinAnnotation)& bp1$effect!="synonymous SNV"){
      prot_sense<-ps1_a$proteinAnnotation
      if(str_detect(prot_sense,"fs")==F){ #when it is different a frameshift
        if(str_detect(prot_sense,"_")==F){
          prot_cerca<-str_sub(prot_sense,1,str_length(prot_sense)-1)
          prot_cerca<-paste("\'",prot_cerca, "%","\'", sep="")
          
          if (is.null(id)) {
            
            ## nothing to do here
            return(NULL)}
          
          else{
            join_pm5b<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification", "effect"
                              FROM "VA_VariantsInTranscripts"
                              LEFT OUTER JOIN "VA_CustomClassification"
                              ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                              WHERE "VA_VariantsInTranscripts"."transcriptId"=', ps1_a$transcriptId,'AND "VA_VariantsInTranscripts"."proteinAnnotation" LIKE ',prot_cerca,'AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
            
            
            pm5.1b<-dbGetQuery(con, join_pm5b)
          }
          duplicatsb<-duplicated(pm5.1b$uniqueVariantId)
          pm5.2b<-NULL
          for (i in 1:nrow(pm5.1b)){ # it removes duplicate entries
            if (nrow(pm5.1b>0)){
              if (duplicatsb[i]==FALSE){
                pm5.3b<-pm5.1b[i,]
                pm5.2b<-rbind(pm5.2b,pm5.3b)
              }
              else if (duplicatsb[i]==TRUE& pm5.1b$date[i]>pm5.2b$date[pm5.2b$uniqueVariantId==pm5.1b$uniqueVariantId[i]]){
                pm5.2b$date[pm5.2b$uniqueVariantId==pm5.1b$uniqueVariantId[i]]<-pm5.1b$date[i]
                pm5.2b$classification[pm5.2b$uniqueVariantId==pm5.1b$uniqueVariantId[i]]<-pm5.1b$classification[i]
              }}
          }
          
          pm5<-sum(pm5.2b$classification=="pat", pm5.2b$classification=="ppat")
          pm5[is.na(pm5)]<-0
          bp8<-sum(pm5.2b$classification=="pol", pm5.2b$classification=="ppol")
        }}}
    pm5.2<-NULL
    duplicats<-duplicated(pm5.1$uniqueVariantId)
    
    
    for (i in 1:nrow(pm5.1)){ # it removes duplicate entries
      if (nrow(pm5.1>0)){
        if (duplicats[i]==FALSE){
          pm5.3<-pm5.1[i,]
          pm5.2<-rbind(pm5.2,pm5.3)
        }
        else if (duplicats[i]==TRUE& pm5.1$date[i]>pm5.2$date[pm5.2$uniqueVariantId==pm5.1$uniqueVariantId[i]]){
          pm5.2$date[pm5.2$uniqueVariantId==pm5.1$uniqueVariantId[i]]<-pm5.1$date[i]
          pm5.2$classification[pm5.2$uniqueVariantId==pm5.1$uniqueVariantId[i]]<-pm5.1$classification[i]
        }}
    }
    
    pm5<-sum(pm5,pm5.2$classification=="pat",pm5.2$classification=="ppat") #it counts rows which are pat / ppat having different amino acid
    bp8<-sum(bp8,pm5.2$classification=="pol",pm5.2$classification=="ppol") #it counts rows which are pol / ppol having different amino acid
    criteria[c("PM5_strong","PM5", "PM5_supporting"), 1][pm5>=2]<-c(1,0,0) #if any row is assigned, PM5 = 1 in the criteria array
    criteria[c("PM5_strong","PM5", "PM5_supporting"), 1][pm5==1]<-c(0,1,0)
    criteria[c("PM5_strong","PM5"), 1][pm5==0]<-c(0,0)
    criteria["BP8", 1][bp8>=1]<-1
    
    if (is.null(id)) {
      
      ## nothing to do here
      return(NULL)}
    
    else{
      
      if(length(group_position)!=0&pm5==0){
        group_position2<-paste("p.",first_aa, group_position,"_", sep="")
        join_pm5<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification", "effect"
                         FROM "VA_VariantsInTranscripts"
                         LEFT OUTER JOIN "VA_CustomClassification"
                         ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                         WHERE "VA_VariantsInTranscripts"."transcriptId"=',Transcripts_RASos[gen,1], 'AND  "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', group_position2[1], '\'AND "VA_VariantsInTranscripts"."proteinAnnotation"<>\'', ps1_a$proteinAnnotation, '\' AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,' ORDER BY "VA_CustomClassification".date', sep="")
        
      }
    }
    if (is.null(id)) {
      
      ## nothing to do here
      return(NULL)}
    
    else{
      pm5_analogs_groups<-dbGetQuery(con,join_pm5)
      #loop to detect if analog groups and aa_change are equal and pathogenic PM1 = 1, and also to see if aa is different and pathogenic, then PM5 = 1
      pm5<-sum(pm5_analogs_groups$classification=="pat"|pm5_analogs_groups$classification=="ppat")
      bp8<-sum(pm5_analogs_groups$classification=="pol"|pm5_analogs_groups$classification=="ppol")
      criteria["PM5_supporting",1][pm5>=1]<-1 # PM5 = 1 if the variant is PAT or PPAT
      criteria["BP8", 1][bp8>=1]<-1
      if (length(group_position)==2){ # it looks if = 2 to consider group of genes 2 that are 3 genes.
        join_pm5<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification"
                         FROM "VA_VariantsInTranscripts"
                         LEFT OUTER JOIN "VA_CustomClassification"
                         ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                         WHERE "VA_VariantsInTranscripts"."transcriptId"=',Transcripts_RASos[gen,2],'AND "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', group_position2[1], '\'AND  "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', group_position2[2], '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
        
      }
      
      pm5_analogs_groups2<-dbGetQuery(con,join_pm5)
      criteria["PM5_supporting",1][pm5_analogs_groups2$classification=="pat"|pm5_analogs_groups2$classification=="ppat"]<-1 #assignem 1 si ?s pat o ppat
      bp8<-sum(pm5_analogs_groups2$classification=="pol"|pm5_analogs_groups2$classification=="ppol")
      bp8<-sum(pm5_analogs_groups$classification=="pol"|pm5_analogs_groups$classification=="ppol")
    }
    
    
    pm5[is.null(pm5)]<-0
    criteria["PM5_supporting", 1][pm5==0]<-0 #si no hi ha cap fila no es compleix
    criteria["BP8", 1][is.na(criteria["BP8", 1])]<-0
    
    ## PM1: Located in a mutational hot spot and/or critical and well-established functional domain (e.g., active site of an enzyme) without benign variation
    
    if (!is.na(ps1_a$proteinAnnotation)& (gen=="NF1"|gen=="SPRED1")){
      prot<-unlist(strsplit(ps1_a$proteinAnnotation,"p.")) # it separates the aa number from "p."
      proti<-prot[2] #it takes the second element, which is the informative one
      first_aa<-str_sub(proti,1,1) # the first aa is extracted
      if(str_detect(proti,"fs")==F){ # when it is different a frameshift
        if(str_detect(proti,"_")==F){
          num<-as.integer(str_sub(proti,2,str_length(proti)-1))
        }
        else if(str_detect(proti,"_")==T){ #it takes inframe positions
          sep<-unlist(strsplit(proti, "_"))
          num<-sep[1]
        }
        
      }
      else{ #if it is a frameshift
        num<-as.integer(str_sub(proti,2,str_length(proti)-2))
      }
    }
    
    #hot spots and domains of RASo genes
    BRAF<-c(238:286,439:474,531,459:474,594:627)
    HRAS<-c(12,13,14,58,59,60,61,62,63)
    KRAS<-c(12,13,14,58,59,60,61,62,63)
    MAP2K1<-c(43:61,124:134)
    MAP2K2<-c(47:65,128:138)
    PTPN11<-c(308,4,7,8,9,58:63,69:77,247,251,255,256,258,261,265,278:281,284)
    RAF1<-c(251:266)
    SHOC2<-c(2)
    SOS1<-c(269,552,420:500)
    SPRED1<-c(6:123, 233-285, 334:442)
    NF1<-c(543:909,1095:1176, 1170:1218, 1198:1530, 1510:1570, 1560:1696,1713:1816, 2619:2719)
    
    hot_spots<-list(BRAF,HRAS, KRAS,MAP2K1,MAP2K2, PTPN11, RAF1,SHOC2, SOS1,SPRED1, NF1)
    names(hot_spots)<-c("BRAF", "HRAS", "KRAS", "MAP2K1", "MAP2K2", "PTPN11", "RAF1", "SHOC2", "SOS1", "SPRED1", "NF1")
    suma<-0
    if (!is.na(ps1_a$proteinAnnotation)){
      hot_sp<-hot_spots[[gen]]==num
      suma<-sum(hot_sp==T)
    }
    criteria["PM1",1][suma==1]<-1 # PM1 = 1 if it is PAT or PPAT
    criteria["PM1",1][suma==0|is.na(suma)]<-0
    
    ##BA1 , BS1 i PM2
    #BA1: Allele frequency is above 5% in Exome Sequencing Project, 1000 Genomes, or ExAC.
    #BS1: Allele frequency is greater than expected for these disorders
    #PM2: Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes, or ExAC
    
    
    if (is.null(id)) {
      
      ## nothing to do here
      return(NULL)}
    
    else{
      freq<-paste('SELECT "PopFreqMax" FROM "VA_Frequencies" WHERE "uniqueVariantId"=', id,'ORDER BY "date" DESC LIMIT 1') #agafem la informaci? de frequencia m?s actualitzada de la variant d'interes
      BA1<- dbGetQuery(con, freq)
      criteria["BA1",1][BA1>=0.0005]<-1 # BA1 = 1 if the poblacional frequency is greater than 0.0005
      criteria["BA1",1][BA1<0.0005|is.na(BA1)]<-0 # BA1 = 0 if the poblacional frequency is lower than 0.0005
      criteria["BS1",1][BA1>=0.00025&BA1<0.0005]<-1 # BS1 = 1 if the poblacional frequency is greater than 0.00025
      criteria["BS1",1][BA1<0.00025|is.na(BA1)|BA1>=0.0005]<-0 # BS1 = 1 if the poblacional frequency is lower than 0.00025
      criteria["BS1_supporting",1][!is.na(BA1)]<-0
      criteria["PM2", 1][BA1==0|is.na(BA1)]<-1
      criteria["PM2", 1][BA1!=0&!is.na(BA1)]<-0
      
    }
    
    
    ##PM4 i BP3
    #PM4: Protein length changes as a result of in- frame deletions/insertions in a nonrepeat region or stop-loss variants
    #BP3: In-frame deletions/insertions in a repetitive region without a known function
    
    criteria["PM4", 1][gen!="NF1"|gen!="SPRED1"]<-0
    criteria["BP3", 1][gen!="NF1"|gen!="SPRED1"]<-0
    
    if (is.null(id)) {
      
      ## nothing to do here
      return(NULL)}
    
    else {
      
      if (gen=="NF1"|gen=="SPRED1"){
        pm4.1<-paste('SELECT b."start", b."end", b."chr", "validatedEffect", "strand"
                     FROM "VA_VariantsInTranscripts" AS a
                     LEFT OUTER JOIN "UniqueVariantsInGenome" AS b
                     ON a."uniqueVariantId"=b."uniqueVariantId"
                     LEFT OUTER JOIN "VA_VariantNomenclature" AS c
                     ON c."uniqueVariantId"=b."uniqueVariantId"
                     LEFT OUTER JOIN "TranscriptsPosition" as d
                     ON d."transcriptId"=c."transcriptId"
                     WHERE a."uniqueVariantId"=', id, 'AND a."isMainTranscript"=\'TRUE\' ORDER BY c."date" DESC LIMIT 1') #it gets information about the coordinates of start, end, effect and the string in which the gene is.
        pm4<-dbGetQuery(con, pm4.1)
        criteria["PM4", 1][pm4$validatedEffect!="frameshift"]<-0 # PM4 = 0 if variant is not a frameshift
        criteria["BP3", 1][!pm4$validatedEffect %in% "in frame"]<-0 #BP3 = 0 if it is in frames
        if(!is.na(pm4$validatedEffect)){
          result<-NULL
        }
        else if(!is.na(pm4$validatedEffect)){
          if (pm4$validatedEffect=="frameshift"|pm4$validatedEffect %in% "in frame"){ # if it equals a frameshift or inframe
            library(RMySQL)
            query<-paste('SELECT * FROM rmsk WHERE genoStart< ',pm4$start, ' AND genoEnd>',pm4$end,' AND genoName=\'chr',pm4$chr,'\' ', sep="") #agafem informacio del UCSC RepeatMasker
            pm4.2<- dbGetQuery(my_connection, query)
            result<-nrow(pm4.2$strand[pm4.2$strand==pm4$strand]) # it looks they are in the same DNA strand
          }
          criteria["PM4", 1][pm4$validatedEffect=="frameshift"&is.null(result)]<-1 # frameshift and there are not results in a none repetitive region, PM4 = 1
          criteria["PM4", 1][pm4$validatedEffect=="frameshift"&!is.null(result)]<-0 # frameshift and there are a result in a repetitive region, PM4 = 0
          criteria["BP3", 1][pm4$validatedEffect %in% "in frame" & !is.null (result)]<-1 # in frame and there are a result in a repetitive region, BP3 = 1
          criteria["BP3", 1][pm4$validatedEffect %in% "in frame" & is.null (result)]<-0 # in frame and and there are not results in a none repetitive region, BP3 = 0
        }
      }
      
    }
    
    ##PP2: Missense variant in a gene that has a low rate of benign missense variation and in which missense variants are a common mechanism of variation
    
    if (is.null(id)) {
      
      ## nothing to do here
      return(NULL)}
    
    else {
      pp2.1<- paste('SELECT "effect","validatedEffect"
                    FROM "VA_VariantsInTranscripts"
                    LEFT OUTER JOIN "VA_VariantNomenclature"
                    ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_VariantNomenclature"."uniqueVariantId"
                    WHERE "VA_VariantsInTranscripts"."uniqueVariantId"=', id, 'AND "isMainTranscript"=\'TRUE\' ORDER BY "VA_VariantsInTranscripts"."date" DESC LIMIT 1') # AND "validatedEffect"<> \'<NA>\' li diem que validated effect no sigui NA
      
    }
    
    pp2<-dbGetQuery(con, pp2.1)
    efecte_validat<-pp2$validatedEffect[!is.na(pp2$validatedEffect)]
    efecte_novalidat<-pp2$effect[!is.na(pp2$effect)]
    criteria["PP2",1][length(efecte_validat)==0&length(efecte_novalidat)==0]<-0
    if( length(efecte_validat)>0){
      if(efecte_validat[1]=="missense (non-synonymous)"){
        criteria["PP2", 1]<-1
      }
      else if (efecte_validat[1]!="missense (non-synonymous)"){
        criteria["PP2", 1]<-0
      }}
    else if (length(efecte_validat)==0 &length(efecte_novalidat)>0){ #it considers the ones that have no value in validated effect and look at the effect field
      if (efecte_novalidat[1]=="nonsynonymous SNV"){
        criteria["PP2", 1]<-1
      }
      else if(efecte_novalidat[1]!="nonsynonymous SNV"){
        criteria["PP2", 1]<-0
      }}
    
    
    ##PP3 & BP4
    #PP3: Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact,etc)
    #BP4: Multiple lines of computational evidence suggest no impact on gene or gene product (conservation, evolutionary, splicing, etc)
    
    
    if (is.null(id)) {
      
      ## nothing to do here
      return(NULL)}
    
    else{
      pp3.1<- paste('SELECT "SIFT_pred", "Polyphen2_HVAR_pred","MutationTaster_pred","PROVEAN_pred", "CADD_phred" 
                    FROM "VA_VariantsInTranscripts"
                    LEFT OUTER JOIN "VA_InSilicoPathogenicity" AS a
                    ON "VA_VariantsInTranscripts"."uniqueVariantId"= a."uniqueVariantId"
                    WHERE "VA_VariantsInTranscripts"."uniqueVariantId"=', id, 'AND
                    "isMainTranscript"=\'TRUE\' ORDER BY a."date" DESC LIMIT 1')
      
    }
    
    pp3<-dbGetQuery(con, pp3.1)
    pp3$CADD_phred[pp3$CADD_phred>19]<-"D"
    pp3$CADD_phred[pp3$CADD_phred<10]<-"N"
    pp3$Polyphen2_HVAR_pred[pp3$Polyphen2_HVAR_pred=="B"]<-"N"
    pp3$Polyphen2_HVAR_pred[pp3$Polyphen2_HVAR_pred=="P"]<-"D"
    pp3$MutationTaster_pred[pp3$MutationTaster_pred=="A"]<-"D"
    pp3$MutationTaster_pred[pp3$MutationTaster_pred=="P"]<-"N"
    pp3$SIFT_pred[pp3$SIFT_pred=="T"]<-"N"
    pp3<-t(pp3)
    table<-table(pp3[,1])
    criteria["PP3",1][table["D"]==5]<-1
    criteria["PP3",1][table["D"]<=4|table["N"]==5]<-0
    criteria["BP4",1][table["N"]==5]<-1
    criteria["BP4",1][table["N"]<=4|table["D"]==5]<-0
    
    
    ##BS2 & PS4
    #BS2: Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked (hemizygous) disorder with full penetrance expected at an early age
    #PS4: The prevalence of the variant in affected individuals is significantly increased compared with the prevalence in controls
    
    #PS4: The prevalence of the variant in affected individuals is significantly increased compared with the prevalence in controls
    
    if (is.null(id)) {
      
      ## nothing to do here
      return(NULL)}
    
    else{
      bs2.1<-paste('SELECT "pedigreeProgenyId", "startDate",f."name"
                   FROM "UniqueVariants" AS a
                   LEFT OUTER JOIN "Variants" AS b
                   ON b."uniqueVariantId"=a."id"
                   LEFT OUTER JOIN "Samples" AS c
                   ON c."id"=b."sampleId"
                   LEFT OUTER JOIN "Individual" AS d
                   ON d."id"=c."individualId"
                   LEFT OUTER JOIN "Analysis" AS e
                   ON b."analysisId"=e."id"
                   LEFT OUTER JOIN "GeneCascade" AS f
                   ON f."id" = e."geneCascadeId"
                   WHERE a."id"=', id,'')
      
      
    }
    
    bs2<-dbGetQuery(con, bs2.1)
    independent_control<-sum(bs2$name!="NF1_SPRED1"&bs2$name!="Rasopathies"|is.na(bs2$name))
    affected_independently<-sum(bs2$name=="NF1_SPRED1"&!is.na(bs2$name)|bs2$name=="Rasopathies"&!is.na(bs2$name))
    
    criteria[c("BS2","PS4_strong", "PS4_moderate", "PS4_supporting"),1][independent_control>=3]<-c(1,0,0,0)
    criteria[c("PS4_strong", "PS4_moderate", "PS4_supporting"),1][affected_independently==0]<-c(0,0,0)
    if (independent_control==0){
      criteria[c("BS2","PS4_strong","PS4_moderate", "PS4_supporting"),1][affected_independently>=6]<-c(0,1,0,0)
      criteria[c("BS2","PS4_strong","PS4_moderate", "PS4_supporting"),1][affected_independently==4|affected_independently==5]<-c(0,0,1,0)
      criteria[c("BS2","PS4_strong","PS4_moderate", "PS4_supporting"),1][affected_independently==3|affected_independently==2]<-c(0,0,0,1)
      criteria[c("BS2","PS4_strong","PS4_moderate", "PS4_supporting"),1][affected_independently==0]<-c(0,0,0,0)
    }
    else if (independent_control<3&independent_control>0){
    }
    
    
    ## Addatpation PM2 & BS1:
    #PM2:
    #Bs1:
    criteria["PM2",1][criteria["PM2",1]==1&independent_control>100]<-0 # independently that Pandora says variant is not described in population databases, if the variant is detected in more than 100 samples analyzed in Pandora and not showing RASopathy-like phentoype, and then used as control individuals, PM2 = 0
    inhousefreq<-paste('SELECT "freqs","inHouseFrequency","inHouseNumber" FROM "VA_InHouseFrequencies" WHERE "uniqueVariantId"=', id,'ORDER BY "date" DESC LIMIT 1') # it takes the most updated frequency information of the variant of interest
    inhouse<-dbGetQuery(con,inhousefreq)
    
    if (is.null(inhousefreq)) {
      
      ## nothing to do here
      return(NULL)
    }
    else {
      if (is.na(BA1$PopFreqMax)){
        if (is.null(independent_control)) {
          
          ## nothing to do here
          return(NULL)
        }
        criteria["BS1",1][independent_control>=1000&inhouse$inHouseFrequency[1]>=0.30]<-1 #independently that Pandora says variant is not described in population databases, if the variant is detected in more than 1000 samples analyzed in Pandora and not showing RASopathy-like phentoype, and then used as control individuals or in house-frequency >0,30, BS1 = 1
        criteria["BS1_supporting",1][independent_control>500&independent_control<1000&inhouse$inHouseFrequency[1]>0.10&inhouse$inHouseFrequency[1]<0.30]<-1
        criteria["BS1_supporting",1][is.na(criteria["BS1_supporting",1])]<-0
      } 
    }
    
    ###PS2 - de novo cases reported
    if (is.null(denovo_confirmed | is.null (input$denovo_confirmed))) {
      
      return(NULL)}
    else {
      criteria[c("PS2_veryStrong", "PS2"),1][denovo_confirmed>=2]<-c(1,0)
      criteria[c("PS2_veryStrong", "PS2"),1][denovo_confirmed==1 & denovo_noconfirmed>=2]<-c(1,0)
      criteria[c("PS2_veryStrong", "PS2"),1][denovo_confirmed==1 & denovo_noconfirmed<2]<-c(0,1)
      criteria[c("PS2_veryStrong", "PS2"),1][denovo_confirmed==0]<-c(0,0)
      
    }
    
    ###PM6  - de novo cases reported no confirmed
    
    if (is.null(denovo_confirmed | is.null (input$denovo_confirmed))) {
      
      return(NULL)}
    else {
    
    criteria[c("PM6_veryStrong","PM6_strong","PM6"),1][denovo_confirmed==0 & denovo_noconfirmed>3]<-c(1,0,0)
    criteria[c("PM6_veryStrong","PM6_strong","PM6"),1][denovo_confirmed==0 & denovo_noconfirmed==3 | denovo_noconfirmed==2]<-c(0,1,0)
    criteria[c("PM6_veryStrong","PM6_strong","PM6"),1][denovo_confirmed==0 & denovo_noconfirmed==1]<-c(0,0,1)
    criteria[c("PM6_veryStrong","PM6_strong","PM6"),1][denovo_confirmed==0 & denovo_noconfirmed==0]<-c(0,0,0)
    }
    
    ###PP1 - cosegregation cases reported
    if (is.null(cosegregation)) {
      
      ## nothing to do here
      return(NULL)}
    else {
      criteria[c("PP1_strong", "PP1_moderate", "PP1_supporting"),1][cosegregation>=7]<-c(1,0,0)
      criteria[c("PP1_strong", "PP1_moderate", "PP1_supporting"),1][cosegregation==5|cosegregation==6]<-c(0,1,0)
      criteria[c("PP1_strong", "PP1_moderate", "PP1_supporting"),1][cosegregation==3|cosegregation==4]<-c(0,0,1)
      criteria[c("PP1_strong", "PP1_moderate", "PP1_supporting"),1][cosegregation<3]<-c(0,0,0)
    }
    
    
    ###PS3 - Functional studies
    
    if (is.null(Functional_evidence)) {
      
      ## nothing to do here
      return(NULL)}
    else {
      criteria[c("PS3"),1]<-c(Functional_evidence[1])
    }
    
    ###PP5, BP6
    if (is.null(PPAT_evidence) | is.null (input$PPOL_evidence)) {
      
      ## nothing to do here
      return(NULL)}
    else {
      criteria[c("PP5", "BP6"),1]<-c(PPAT_evidence[1], PPOL_evidence[2])
      
    }
    
    ## final classification
    
    #dataframe to fill with all criteria (each criteria in one line)
    
    if (is.null(criteria) ) {
      
      ## nothing to do here
      return(NULL)}
    else {
      haz.cero.na=function(x){
      ifelse(is.na(x),0,x)}
    criteria=data.frame(sapply(criteria,haz.cero.na))
          }
    
    suma_criteria<- data.frame(verystrong_pat=sum(criteria[2,1], criteria[9,1], criteria[37,1]),
                               strong_pat=sum(criteria[1,1], criteria[3:5,1], criteria[8,1], criteria[14,1], criteria[33,1]),
                               moderate_pat=sum(criteria[6,1],criteria[10:13,1],criteria[15,1],criteria[34,1]), 
                               supporting_pat=sum(criteria[16:19,1], criteria[7,1], criteria[35,1]), 
                               very_strong_benign=sum(criteria[20,1]), 
                               strong_benign=sum(criteria[21:24,1]), 
                               support_benign=sum(criteria[25:32,1],criteria[36,1]))
    
    classification<- data.frame(PAT=0, Likely_PAT=0, Neutral=0, Likely_Neutral=0, VUS=0)
    
    if (is.null(suma_criteria) ) {
      
      ## nothing to do here
      return(NULL)}
    else {
     
    classification$PAT[  suma_criteria$verystrong_pat>=1 & suma_criteria$supporting_pat>=2|
                         suma_criteria$verystrong_pat>=1 & suma_criteria$moderate_pat>=2|
                         suma_criteria$verystrong_pat>=1 & suma_criteria$moderate_pat==1 & suma_criteria$supporting_pat==1|
                         suma_criteria$verystrong_pat>=1 & suma_criteria$supporting_pat>=2|
                         suma_criteria$verystrong_pat>=1 & suma_criteria$strong_pat>=1|
                         suma_criteria$strong_pat>=2|
                         suma_criteria$strong_pat==1 & suma_criteria$moderate_pat>=3|
                         suma_criteria$strong_pat==1 & suma_criteria$moderate_pat==2 & suma_criteria$supporting_pat>=2|
                         suma_criteria$strong_pat==1 & suma_criteria$moderate_pat==1 & suma_criteria$supporting_pat>=4]<-1
    
    classification$Likely_PAT[  suma_criteria$verystrong_pat==1 & suma_criteria$moderate_pat==1|
                                suma_criteria$strong_pat==1 & suma_criteria$moderate_pat>=1|
                                suma_criteria$strong_pat==1 & suma_criteria$supporting_pat>=2|
                                suma_criteria$moderate_pat>=3|
                                suma_criteria$moderate_pat==2 & suma_criteria$supporting_pat>=2|
                                suma_criteria$moderate_pat==1 & suma_criteria$supporting_pat>=4]<-1
    
    classification$Neutral[ suma_criteria$very_strong_benign==1|
                             suma_criteria$strong_benign>=2]<-1
    
    classification$Likely_Neutral[  suma_criteria$strong_benign==1 & suma_criteria$support_benign==1|
                                    suma_criteria$support_benign>=2]<-1
   
     classification$VUS[ classification$PAT==0 & classification$Likely_PAT== 0 & classification$Neutral==0 & classification$Likely_Neutral==0|
                         classification$PAT == 1 & classification$Neutral==1 | classification$PAT==1 & classification$Likely_Neutral==1 |
                        classification$Likely_PAT==1 & classification$Neutral==1 | classification$Likely_PAT==1 & classification$Likely_Neutral==1]<-1
    
    return(classification)
  }
  }
  
  
  #### VARIANT CONFIRMATION
     
      variant_reactive <- eventReactive(c(input$go, input$id),
                                       
                                       if (input$go == 0){
                                         return(NULL)
                                       }
                                       else {
                              (dbGetQuery(con, (bp1.1<-paste('SELECT "cDNAAnnotation", "proteinAnnotation", "symbol","validatedEffect", "effect"
                              FROM "VA_VariantsInTranscripts" AS a
                              LEFT OUTER JOIN "UniqueVariantsInGenome" AS b
                              ON a."uniqueVariantId"=b."uniqueVariantId"
                              LEFT OUTER JOIN "VA_VariantNomenclature" AS c
                              ON c."uniqueVariantId"=b."uniqueVariantId"
                              LEFT OUTER JOIN "Genes" AS d
                              ON c."transcriptId"=d."mainTranscriptId"
                              WHERE a."uniqueVariantId"=', input$id, 'AND 
                              a."isMainTranscript"=\'TRUE\' AND d."symbol"<> \'NA\' ORDER BY a."date" DESC LIMIT 1'))))})
      
      
      output$Variant <- renderTable(expr = variant_reactive(), rownames = TRUE, bordered = FALSE)
      
  
  #### AUTOMATIC VARIANT CLASSIFICATION
  
        
      AutomClass_reactive <-eventReactive(c(input$go, input$id, input$denovo_noconfirmed, input$denovo_confirmed, input$cosegregation,
                                            input$PPAT_evidence,input$PPOL_evidence, input$Functional_evidence),
                                          if (input$go == 0 | input$id == 0){
                                            return(NULL)
                                          }
                                          else 
                                          {Automatic_criteria_AMCG(id = input$id, con = con, denovo_noconfirmed = input$denovo_noconfirmed,
                                                              denovo_confirmed = input$denovo_confirmed, cosegregation = input$cosegregation,
                                                              PPAT_evidence = as.numeric(as.character(input$PPAT_evidence)), 
                                                              PPOL_evidence = as.numeric(as.character(input$PPOL_evidence)),
                                                              Functional_evidence = as.numeric(as.character(input$Functional_evidence)))})
      
      
      criteria <-  renderTable(if (input$go == 0 | input$id == 0){
        return(NULL)
      }
      else  
      {expr = AutomClass_reactive()},rownames = FALSE, bordered = FALSE)
      output$AutoClass <- criteria


  #### FINAL VARIANT CLASSIFICATION
     
     FinalClass_reactive <-eventReactive(c(input$go, input$id, input$denovo_noconfirmed, input$denovo_confirmed, input$cosegregation,
                                           input$PPAT_evidence,input$PPOL_evidence, input$Functional_evidence),
                                         if (input$go == 0){
                                           return(NULL)
                                         }
                                         else 
                                         {Final_classificationB(id = input$id, con = con, denovo_noconfirmed = input$denovo_noconfirmed,
                                                           denovo_confirmed = input$denovo_confirmed, cosegregation = input$cosegregation,
                                                           PPAT_evidence = as.numeric(as.character(input$PPAT_evidence)), 
                                                           PPOL_evidence = as.numeric(as.character(input$PPOL_evidence)),
                                                           Functional_evidence = as.numeric(as.character(input$Functional_evidence)))})

  
                                 
     output$FinalClass <- renderTable(expr = FinalClass_reactive(),rownames = TRUE, bordered = FALSE)

    ## session end clean up
    cancel.onSessionEnded <- session$onSessionEnded(function() {
        RPostgreSQL::dbDisconnect(con)
        RMySQL::dbDisconnect(my_connection)
  })
    
     
}
     

# Run the application 
shinyApp(ui = ui, server = server)
