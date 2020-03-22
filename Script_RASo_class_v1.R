

##Pandora
#install.packages("RPostgreSQL")
require("RPostgreSQL")
drv <- dbDriver("PostgreSQL")  # loads the PostgreSQL driver
#con <- dbConnect(drv, dbname = "ngsdiagnosticsdevel", host = "172.19.26.99", port = 5432, user = "pandora", password =  "oreo no more")
con <- dbConnect(drv, dbname = "ngsdiagnostics", host = "172.19.26.99", port = 5432, user = "pandora", password =  "oreo no more")


##UCSC
#install.packages("RMySQL")
library(RMySQL)
my_connection <- dbConnect(
  MySQL(),
  user="genome",
  dbname="hg38",
  host="genome-euro-mysql.soe.ucsc.edu",
  port=3306
)

#install.packages("stringr")
require(stringr)

### creem la funció myFunction
##gens possibles : BRAF, CBL, HRAS, KRAS, LZTR1, NF1, NRAS, MAP2K1, MAP2K2, PTPN11, RAF1, RIT1, SHOC2, SOS1, SOS2, SPRED1 i RASA1,

myFunction<- function(id){
  #fem un dataframe on cada linia es un criteri
  criteris<-data.frame(criteris=c(rep(NA, 37)), row.names = c("PS1", "PS2_veryStrong", "PS2", "PS3", "PS4_strong", "PS4_moderate", "PS4_supporting", "PM5_strong",  "PM6_veryStrong", "PM6","PM1", "PM2", "PM4",  "PP1_strong", "PP1_moderate", "PP1_supporting", "PP2", "PP3","PP5", "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7", "BS1_supporting", "PM6_strong", "PM5", "PM5_supporting", "BP8", "PVS1"))
  #establim els dominis per grups
  domini_grup2<-data.frame(row.names<-c("RBD1", "RBD2"),BRAF.start=c(157,201), BRAF.end=c(200,261), RAF1.start=c(58,105), RAF1.end=c(101,185), row.names=T)
  domini_grup3<-data.frame(row.names<-c("Effector region", "GTP", "GTP1", "GTP2"),HRAS.start=c(32,10,57,116), HRAS.end=c(40,17,61,119), KRAS.start=c(32,10,59,116), KRAS.end=c(38,18,60,119),NRAS.start=c(0,10,57,116), NRAS.end=c(0,18,61,119),row.names=T)
  domini_grup4<-data.frame(row.names<-c("NES", "Negative regulatory region", "P-loop", "ATP binding","Catalytic loop", "Activation loop", "Protein rich domain", "Versatile docking"),MAP2K1.start=c(32,44,74,143,192,208,262,362), MAP2K1.end=c(44,51,82,146,195,233,307,396), MAP2K2.start=c(36,48,78,147,196,212,266,370), MAP2K2.end=c(48,55,96,150,199,237,311,400), row.names=T)
  domini_grup5<-data.frame(row.names<-c("DH", "PH", "N-terminal Ras-GEF", "Ras-GEF"),SOS1.start=c(200,444,597,780), SOS1.end=c(390,548,741,1019), SOS2.start=c(198,442,595,778), SOS2.end=c(388,546,739,1017), row.names=T)
  
  posicio_grup<-NULL
  num<-NULL
  ###bp1 i pvs1 (només SPRED1 i NF1)
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
  gen<-bp1[1,"symbol"]
  
  if (!is.na(bp1$validatedEffect)){
    #criteris["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&bp1[1,"validatedEffect"]=="missense (non-synonymous)"]<-1 #si el gen es NF1 o SPRED1 i és una missense, s'assigna el valor 1 a la matriu criteris BP1
    #criteris["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&bp1[1,"validatedEffect"]!="missense (non-synonymous)"]<-0#si el gen es NF1 o SPRED1 i no  ?s una missense, s'assigna el valor 0 a la matriu criteris BP1
    criteris["BP1",1][(bp1[1,"symbol"]!="NF1"&bp1[1,"symbol"]!="SPRED1")&(bp1[1,"validatedEffect"]=="nonsense (stop gain)"|bp1[1,"validatedEffect"]=="nonsense/splice"|bp1[1,"validatedEffect"]=="frameshift")]<-1 #si no ?s NF1 o SPRED1, i es una trucating variant s'assigna un 1.
    criteris["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&(bp1[1,"validatedEffect"]=="nonsense (stop gain)"|bp1[1,"validatedEffect"]=="nonsense/splice"|bp1[1,"validatedEffect"]=="frameshift")]<-0 #si no ?s NF1 o SPRED1, i es una trucating variant s'assigna un 1.
    criteris["BP1",1][bp1[1,"validatedEffect"]=="missense (non-synonymous)"]<-0 #si es una missense s'assigna un 0
    criteris["PVS1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&(bp1[1,"validatedEffect"]=="nonsense (stop gain)"|bp1[1,"validatedEffect"]=="nonsense/splice"|bp1[1,"validatedEffect"]=="stopgain"|bp1[1,"validatedEffect"]=="frameshift")]<-1 #si no ?s NF1 o SPRED1, i es una trucating variant s'assigna un 1.
    
  }
  
  if (is.na(bp1$validatedEffect)){
    #criteris["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&bp1[1,"effect"]=="nonsynonymous SNV"]<-1 #si el gen es NF1 o SPRED1 i és una missense, s'assigna el valor 1 a la matriu criteris BP1
    #criteris["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&bp1[1,"effect"]!="nonsynonymous SNV"]<-0#si el gen es NF1 o SPRED1 i no  ?s una missense, s'assigna el valor 0 a la matriu criteris BP1
    criteris["BP1",1][(gen!="NF1"&gen!="SPRED1")&(bp1[1,"effect"]=="stopgain"|bp1[1,"effect"]=="stoploss"|bp1[1,"effect"]=="frameshift deletion")]<-1 #si no és NF1 o SPRED1, i es una trucating variant s'assigna un 1.
    criteris["BP1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&(bp1[1,"effect"]=="stopgain"|bp1[1,"effect"]=="stoploss"|bp1[1,"effect"]=="frameshift deletion")]<-0
    criteris["BP1",1][bp1[1,"effect"]=="nonsynonymous SNV"]<-0 #si es una missense s'assigna un 0
    criteris["PVS1",1][(bp1[1,"symbol"]=="NF1"|bp1[1,"symbol"]=="SPRED1")&(bp1[1,"effect"]=="stopgain"|bp1[1,"effect"]=="stoploss"|bp1[1,"effect"]=="frameshift deletion")]<-1
    
  }
  
  criteris["PVS1",1][is.na(criteris["PVS1",1])]<-0 #assignem 0 a PVS1 si no s'ha complert cap de les condicions anteriors
  
  
  ###ps1
  ps1_aa<-paste('SELECT "cDNAAnnotation", "proteinAnnotation", "transcriptId"
        FROM "VA_VariantsInTranscripts"
        WHERE "isMainTranscript"=\'TRUE\' AND "uniqueVariantId"=', id, 'ORDER BY date DESC LIMIT 1')
  ps1_a<-dbGetQuery(con, ps1_aa) 
  cdna<-unlist(strsplit(ps1_a$cDNAAnnotation,">")) #separa la cadena de text quan troba >
  cdna<-paste(cdna[1], "%", sep="") #afegeix % al primer element de la cadena, per tal de que cerqui qualsevol patró que comenci com cdna[1] seguit del que sigui
  join_ps1<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification"
  FROM "VA_VariantsInTranscripts"
  LEFT OUTER JOIN "VA_CustomClassification"
  ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
  WHERE "VA_VariantsInTranscripts"."transcriptId"=',ps1_a$transcriptId,'AND "VA_VariantsInTranscripts"."cDNAAnnotation" LIKE \'',cdna,'\' AND "VA_VariantsInTranscripts"."proteinAnnotation"=\'', ps1_a$proteinAnnotation, '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
  ps1.1<-dbGetQuery(con, join_ps1)
  ps1<-sum(ps1.1$classification=="pat") #compta les files que contenen la paraula pat de les variants sinònimes
  criteris["PS1", 1][ps1>=1]<-1 
  transcripts_id<-data.frame(posicio=c(283,181,229,67,30,34,240,7,7),posicio2=c(rep(NA,6),235,240,235), row.names = c("SOS1", "SOS2", "MAP2K1", "MAP2K2", "BRAF", "RAF1","HRAS","NRAS","KRAS"))
  
  if (gen=="SOS1"|gen=="SOS2"|gen=="BRAF"|gen=="RAF1"|gen=="MAP2K1"|gen=="MAP2K2"|gen=="HRAS"|gen=="NRAS"|gen=="KRAS"){
    if (ps1==0 & !is.na(ps1_a$proteinAnnotation)&!is.na(bp1[1,"effect"])&(bp1[1,"effect"]!="synonymous SNV"|bp1[1,"effect"]=="nonsynonymous SNV")){
      prot<-unlist(strsplit(ps1_a$proteinAnnotation,"p.")) #es separa pel p.
      proti<-prot[2]#s’agafa el segon element, que es on hi ha la informacio
      primer_aa<-str_sub(proti,1,1) #s’extreu el primer aa
      num<-as.integer(str_sub(proti,2,str_length(proti)-1)) #s’agafa des de la posició 2 a la -1, es a dir la posició del aa
      aa_canvi<-str_sub(proti, str_length(proti), str_length(proti))
      posicio_grup<-c()
      if (gen=="SOS1"){ #si el gen es SOS1
        domini_funcional<-domini_grup5[num>=domini_grup5$SOS1.start&num<=domini_grup5$SOS1.end,] #comprovem si es troba dins del domini
        posicio_relativa<-num-domini_funcional$SOS1.start #establim la posicio relavtiva dins daquell domini
        posicio_grup<-posicio_relativa+domini_funcional$SOS2.start #mirem a quina posicio equival de l'altre domini
      }
      
      if (gen=="SOS2"){
        domini_funcional<-domini_grup5[num>=domini_grup5$SOS2.start&num<=domini_grup5$SOS2.end,]
        posicio_relativa<-num-domini_funcional$SOS2.start
        posicio_grup<-posicio_relativa+domini_funcional$SOS1.start
      }
      
      if(gen=="BRAF"){
        domini_funcional<-domini_grup2[num>=domini_grup2$BRAF.start&num<=domini_grup2$BRAF.end,]
        posicio_relativa<-num-domini_funcional$BRAF.start
        posicio_grup<-posicio_relativa+domini_funcional$RAF1.start
      }
      if (gen=="RAF1"){
        domini_funcional<-domini_grup2[num>=domini_grup2$RAF1.start&num<=domini_grup2$RAF1.end,]
        posicio_relativa<-num-domini_funcional$RAF1.start
        posicio_grup<-posicio_relativa+domini_funcional$BRAF.start
      }
      if (gen=="MAP2K1"){
        domini_funcional<-domini_grup4[num>=domini_grup4$MAP2K1.start&num<=domini_grup4$MAP2K1.end,]
        posicio_relativa<-num-domini_funcional$MAP2K1.start
        posicio_grup<-posicio_relativa+domini_funcional$MAP2K2.start}
      
      
      if (gen=="MAP2K2"){
        domini_funcional<-domini_grup4[num>=domini_grup4$MAP2K2.start&num<=domini_grup4$MAP2K2.end,]
        posicio_relativa<-num-domini_funcional$MAP2K2.start
        posicio_grup<-posicio_relativa+domini_funcional$MAP2K1.start
      }
      if (gen=="NRAS"){
        domini_funcional<-domini_grup3[num>=domini_grup3$NRAS.start&num<=domini_grup3$NRAS.end,]
        posicio_relativa<-num-domini_funcional$NRAS.start
        posicio_grup<-c(posicio_relativa+domini_funcional$HRAS.start,posicio_relativa+domini_funcional$KRAS.start )
      }
      if (gen=="HRAS"){
        domini_funcional<-domini_grup3[num>=domini_grup3$HRAS.start&num<=domini_grup3$HRAS.end,]
        posicio_relativa<-num-domini_funcional$HRAS.start
        if (nrow(domini_funcional)>0){
          if (domini_funcional$NRAS.start==0){
            posicio_grup<-posicio_relativa+domini_funcional$KRAS.start
          }
          if (domini_funcional$NRAS.start!=0){
            posicio_grup<-c(posicio_relativa+domini_funcional$KRAS.start,posicio_relativa+domini_funcional$NRAS.start )
          }
        }}
      if (gen=="KRAS"){
        domini_funcional<-domini_grup3[num>=domini_grup3$KRAS.start&num<=domini_grup3$KRAS.end,]
        posicio_relativa<-num-domini_funcional$KRAS.start
        if (nrow(domini_funcional)>0){
          if (domini_funcional$NRAS.start==0){
            posicio_grup<-posicio_relativa+domini_funcional$HRAS.start
          }
          if (domini_funcional$NRAS.start!=0){
            posicio_grup<-c(posicio_relativa+domini_funcional$HRAS.start,posicio_relativa+domini_funcional$NRAS.start )
          }}}
      
      if(length(posicio_grup)!=0){
        posicio_grup2<-paste("p.",primer_aa, posicio_grup, aa_canvi, sep="")
        join_ps1<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification"
  FROM "VA_VariantsInTranscripts"
                   LEFT OUTER JOIN "VA_CustomClassification"
                   ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                   WHERE "VA_VariantsInTranscripts"."transcriptId"=',transcripts_id[gen,1],'AND  "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', posicio_grup2[1], '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,' ORDER BY "VA_CustomClassification".date LIMIT 1', sep="")
        ps1_grups_analegs<-dbGetQuery(con,join_ps1)
        #ara comparar si aa canviat és el mateix que en la que nosaltres busquem
        #grups_analegs<-str_sub(pm5_grups_analegs$proteinAnnotation, str_length(pm5_grups_analegs$proteinAnnotation), str_length(pm5_grups_analegs$proteinAnnotation))
        #fer un bucle perque si grups analegs i aa_canvi son iguals i es patogenica es posi PM1 i fer un bucle perque si són aa diferent i es pato es posi PM5
        criteris["PS1",1][ps1_grups_analegs$classification=="pat"|ps1_grups_analegs$classification=="ppat"]<-1 #assignem 1 si és pat o ppat
        if (length(posicio_grup)==2){ #mirem si és igual a dos perque tingui en compte grup de gens 2 que són 3 gens.
          join_ps1<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification"
                     FROM "VA_VariantsInTranscripts"
                     LEFT OUTER JOIN "VA_CustomClassification"
                     ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                     WHERE "VA_VariantsInTranscripts"."transcriptId"=',transcripts_id[gen,2],'AND  "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', posicio_grup2[2], '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
          ps1_grups_analegs2<-dbGetQuery(con,join_ps1)
          criteris["PS1",1][ps1_grups_analegs2$classification=="pat"|ps1_grups_analegs2$classification=="ppat"]<-1 #assignem 1 si és pat o ppat
        }
      }
    }}
  
  
  criteris["PS1", 1][is.na(criteris["PS1", 1])]<-0 #assignem 0 si no es pat o ppat cap dels anteriors.
  
  
  
  ###pm5 i bp8
  pm5<-0  
  bp8<-0
  join_pm5<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification", "effect"
  FROM "VA_VariantsInTranscripts"
                   LEFT OUTER JOIN "VA_CustomClassification"
                    ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                    WHERE "VA_VariantsInTranscripts"."transcriptId"=', ps1_a$transcriptId,'AND "VA_VariantsInTranscripts"."cDNAAnnotation" LIKE \'',cdna,'\' AND "VA_VariantsInTranscripts"."proteinAnnotation"<>\'', ps1_a$proteinAnnotation, '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
  pm5.1<-dbGetQuery(con, join_pm5)
  if(!is.na(ps1_a$proteinAnnotation)& bp1$effect!="synonymous SNV"){
    prot_sense<-ps1_a$proteinAnnotation
    if(str_detect(prot_sense,"fs")==F){ #quan no es un frameshift
      if(str_detect(prot_sense,"_")==F){
        prot_cerca<-str_sub(prot_sense,1,str_length(prot_sense)-1)
        prot_cerca<-paste("\'",prot_cerca, "%","\'", sep="")
        join_pm5b<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification", "effect"
  FROM "VA_VariantsInTranscripts"
                   LEFT OUTER JOIN "VA_CustomClassification"
                   ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                   WHERE "VA_VariantsInTranscripts"."transcriptId"=', ps1_a$transcriptId,'AND "VA_VariantsInTranscripts"."proteinAnnotation" LIKE ',prot_cerca,'AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
        
        pm5.1b<-dbGetQuery(con, join_pm5b)
        duplicatsb<-duplicated(pm5.1b$uniqueVariantId)
        pm5.2b<-NULL
        for (i in 1:nrow(pm5.1b)){ #eliminem les entrades duplicades
          if (nrow(pm5.1b>0)){
            if (duplicatsb[i]==FALSE){
              pm5.3b<-pm5.1b[i,]
              pm5.2b<-rbind(pm5.2b,pm5.3b)
            }
            if (duplicatsb[i]==TRUE& pm5.1b$date[i]>pm5.2b$date[pm5.2b$uniqueVariantId==pm5.1b$uniqueVariantId[i]]){
              pm5.2b$date[pm5.2b$uniqueVariantId==pm5.1b$uniqueVariantId[i]]<-pm5.1b$date[i]
              pm5.2b$classification[pm5.2b$uniqueVariantId==pm5.1b$uniqueVariantId[i]]<-pm5.1b$classification[i]
            }}
        }
        
        pm5<-sum(pm5.2b$classification=="pat", pm5.2b$classification=="ppat")
        pm5[is.na(pm5)]<-0
        bp8<-sum(pm5.2b$classification=="pol", pm5.2b$classification=="ppol")
      }}}
  # valor_unic<-unique(pm5.1$uniqueVariantId)
  pm5.2<-NULL
  duplicats<-duplicated(pm5.1$uniqueVariantId)
  
  
  for (i in 1:nrow(pm5.1)){ #eliminem les entrades duplicades
    if (nrow(pm5.1>0)){
      if (duplicats[i]==FALSE){
        pm5.3<-pm5.1[i,]
        pm5.2<-rbind(pm5.2,pm5.3)
      }
      if (duplicats[i]==TRUE& pm5.1$date[i]>pm5.2$date[pm5.2$uniqueVariantId==pm5.1$uniqueVariantId[i]]){
        pm5.2$date[pm5.2$uniqueVariantId==pm5.1$uniqueVariantId[i]]<-pm5.1$date[i]
        pm5.2$classification[pm5.2$uniqueVariantId==pm5.1$uniqueVariantId[i]]<-pm5.1$classification[i]
      }}
  }
  
  
  pm5<-sum(pm5,pm5.2$classification=="pat",pm5.2$classification=="ppat")#compta les files que són pat/ppat que tenen aminoacid diferent
  bp8<-sum(bp8,pm5.2$classification=="pol",pm5.2$classification=="ppol")#compta les files que són pol/ppol que tenen aminoacid diferent
  criteris[c("PM5_strong","PM5", "PM5_supporting"), 1][pm5>=2]<-c(1,0,0) #si alguna fila es pat li assigna 1 a la matriu de criteris fila pm5
  criteris[c("PM5_strong","PM5", "PM5_supporting"), 1][pm5==1]<-c(0,1,0)
  criteris[c("PM5_strong","PM5"), 1][pm5==0]<-c(0,0)
  criteris["BP8", 1][bp8>=1]<-1
  
  if(length(posicio_grup)!=0&pm5==0){
    posicio_grup2<-paste("p.",primer_aa, posicio_grup,"_", sep="")
    join_pm5<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification", "effect"
  FROM "VA_VariantsInTranscripts"
                   LEFT OUTER JOIN "VA_CustomClassification"
                   ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                   WHERE "VA_VariantsInTranscripts"."transcriptId"=',transcripts_id[gen,1], 'AND  "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', posicio_grup2[1], '\'AND "VA_VariantsInTranscripts"."proteinAnnotation"<>\'', ps1_a$proteinAnnotation, '\' AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,' ORDER BY "VA_CustomClassification".date', sep="")
    pm5_grups_analegs<-dbGetQuery(con,join_pm5)
    #fer un bucle perque si grups analegs i aa_canvi son iguals i es patogenica es posi PM1 i fer un bucle perque si són aa diferent i es pato es posi PM5
    pm5<-sum(pm5_grups_analegs$classification=="pat"|pm5_grups_analegs$classification=="ppat")
    bp8<-sum(pm5_grups_analegs$classification=="pol"|pm5_grups_analegs$classification=="ppol")
    criteris["PM5_supporting",1][pm5>=1]<-1 #assignem 1 si és pat o ppat
    criteris["BP8", 1][bp8>=1]<-1
    if (length(posicio_grup)==2){ #mirem si és igual a dos perque tingui en compte grup de gens 2 que són 3 gens.
      join_pm5<- paste('SELECT  "VA_VariantsInTranscripts"."uniqueVariantId", "cDNAAnnotation", "proteinAnnotation", "transcriptId",  "VA_CustomClassification"."date", "classification"
                     FROM "VA_VariantsInTranscripts"
                     LEFT OUTER JOIN "VA_CustomClassification"
                     ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_CustomClassification"."uniqueVariantId"
                     WHERE "VA_VariantsInTranscripts"."transcriptId"=',transcripts_id[gen,2],'AND "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', posicio_grup2[1], '\'AND  "VA_VariantsInTranscripts"."proteinAnnotation"LIKE\'', posicio_grup2[2], '\'AND "VA_VariantsInTranscripts"."isMainTranscript"=\'TRUE\' AND "VA_VariantsInTranscripts"."uniqueVariantId" <>', id,'', sep="")
      pm5_grups_analegs2<-dbGetQuery(con,join_pm5)
      criteris["PM5_supporting",1][pm5_grups_analegs2$classification=="pat"|pm5_grups_analegs2$classification=="ppat"]<-1 #assignem 1 si és pat o ppat
      bp8<-sum(pm5_grups_analegs2$classification=="pol"|pm5_grups_analegs2$classification=="ppol")
      bp8<-sum(pm5_grups_analegs$classification=="pol"|pm5_grups_analegs$classification=="ppol")
    }
  }
  
  
  
  
  pm5[is.null(pm5)]<-0
  criteris["PM5_supporting", 1][pm5==0]<-0 #si no hi ha cap fila no es compleix
  criteris["BP8", 1][is.na(criteris["BP8", 1])]<-0
  
  ##pm1
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
  if (!is.na(ps1_a$proteinAnnotation)& (gen=="NF1"|gen=="SPRED1")){
    prot<-unlist(strsplit(ps1_a$proteinAnnotation,"p.")) #es separa pel p.
    proti<-prot[2]#s’agafa el segon element, que es on hi ha la informacio
    primer_aa<-str_sub(proti,1,1) #s’extreu el primer aa
    if(str_detect(proti,"fs")==F){ #quan no es un frameshift
      if(str_detect(proti,"_")==F){
        num<-as.integer(str_sub(proti,2,str_length(proti)-1))
      }
      if(str_detect(proti,"_")==T){ #per agafar les posicions de les inframes
        sep<-unlist(strsplit(proti, "_"))
        num<-sep[1]
      }
      
    }
    else{ #si es una frameshift
      num<-as.integer(str_sub(proti,2,str_length(proti)-2))
    }
  }
  
  
  hot_spots<-list(BRAF,HRAS, KRAS,MAP2K1,MAP2K2, PTPN11, RAF1,SHOC2, SOS1,SPRED1, NF1)
  names(hot_spots)<-c("BRAF", "HRAS", "KRAS", "MAP2K1", "MAP2K2", "PTPN11", "RAF1", "SHOC2", "SOS1", "SPRED1", "NF1")
  suma<-0
  if (!is.na(ps1_a$proteinAnnotation)){
    hot_sp<-hot_spots[[gen]]==num
    suma<-sum(hot_sp==T)
  }
  criteris["PM1",1][suma==1]<-1 #assignem 1 si és pat o ppat
  criteris["PM1",1][suma==0|is.na(suma)]<-0
  
  
  ##BA1 , BS1 i PM2
  freq<-paste('SELECT "PopFreqMax" FROM "VA_Frequencies" WHERE "uniqueVariantId"=', id,'ORDER BY "date" DESC LIMIT 1') #agafem la informaci? de frequencia m?s actualitzada de la variant d'interes
  BA1<- dbGetQuery(con, freq)
  criteris["BA1",1][BA1>=0.0005]<-1 #si es major de 0.0005 la frequencia assignem 1 a BA1
  criteris["BA1",1][BA1<0.0005|is.na(BA1)]<-0 #si ?s menor de 0.0005 la frequencia assignem 0 a BA1
  criteris["BS1",1][BA1>=0.00025&BA1<0.0005]<-1 #si ?s major de 0.00025 la frequencia assignem 1 a BS1
  criteris["BS1",1][BA1<0.00025|is.na(BA1)|BA1>=0.0005]<-0 #si ?s menor de 0.00025 la frequencia assignem 0 a BS1
  criteris["BS1_supporting",1][!is.na(BA1)]<-0
  criteris["PM2", 1][BA1==0|is.na(BA1)]<-1
  criteris["PM2", 1][BA1!=0&!is.na(BA1)]<-0
  
  ##pm4 i bp3
  
  criteris["PM4", 1][gen!="NF1"|gen!="SPRED1"]<-0
  criteris["BP3", 1][gen!="NF1"|gen!="SPRED1"]<-0
  if (gen=="NF1"|gen=="SPRED1"){
    pm4.1<-paste('SELECT b."start", b."end", b."chr", "validatedEffect", "strand"
   FROM "VA_VariantsInTranscripts" AS a
   LEFT OUTER JOIN "UniqueVariantsInGenome" AS b
   ON a."uniqueVariantId"=b."uniqueVariantId"
   LEFT OUTER JOIN "VA_VariantNomenclature" AS c
   ON c."uniqueVariantId"=b."uniqueVariantId"
   LEFT OUTER JOIN "TranscriptsPosition" as d
   ON d."transcriptId"=c."transcriptId"
   WHERE a."uniqueVariantId"=', id, 'AND a."isMainTranscript"=\'TRUE\' ORDER BY c."date" DESC LIMIT 1') #agafa informacio sobre la coordenada d'inici, final, efecte i la cadena en que est? el gen.
    pm4<-dbGetQuery(con, pm4.1)
    criteris["PM4", 1][pm4$validatedEffect!="frameshift"]<-0 #si no es una frameshift ja assigna 0 a PM4
    criteris["BP3", 1][!pm4$validatedEffect %in% "in frame"]<-0 #si no es una in frame ja assigna 0 a BP3
    if(!is.na(pm4$validatedEffect)){
      resultat<-NULL
    }
    if(!is.na(pm4$validatedEffect)){
      if (pm4$validatedEffect=="frameshift"|pm4$validatedEffect %in% "in frame"){ #si es igual a frameshift o conte la paraula in frame
        library(RMySQL)
        query<-paste('SELECT * FROM rmsk WHERE genoStart< ',pm4$start, ' AND genoEnd>',pm4$end,' AND genoName=\'chr',pm4$chr,'\' ', sep="") #agafem informacio del UCSC RepeatMasker
        pm4.2<- dbGetQuery(my_connection, query)
        resultat<-nrow(pm4.2$strand[pm4.2$strand==pm4$strand]) #mirem que estiguin en la mateixa cadena de DNA
      }
      criteris["PM4", 1][pm4$validatedEffect=="frameshift"&is.null(resultat)]<-1 #si es frameshift i no hi ha cap resultat (no esta en regio repetitiva) assigna un 1 a PM4
      criteris["PM4", 1][pm4$validatedEffect=="frameshift"&!is.null(resultat)]<-0 #si es frameshift i hi ha algun resultat (esta en regio repetitiva) assigna un 0 a PM4
      criteris["BP3", 1][pm4$validatedEffect %in% "in frame" & !is.null (resultat)]<-1 #si conte in frame el resulat i hi ha alguna resultat (esta en regio repetitiav) assiga un 1 a BP3
      criteris["BP3", 1][pm4$validatedEffect %in% "in frame" & is.null (resultat)]<-0#si cont? in frame el resulat i no  hi ha cap resultat (no esta en regio repetitiva) assiga un 0 a BP3
    }
  }
  
  ##PP2
  
  pp2.1<- paste('SELECT "effect","validatedEffect"
              FROM "VA_VariantsInTranscripts"
               LEFT OUTER JOIN "VA_VariantNomenclature"
               ON "VA_VariantsInTranscripts"."uniqueVariantId"="VA_VariantNomenclature"."uniqueVariantId"
               WHERE "VA_VariantsInTranscripts"."uniqueVariantId"=', id, 'AND "isMainTranscript"=\'TRUE\' ORDER BY "VA_VariantsInTranscripts"."date" DESC LIMIT 1') # AND "validatedEffect"<> \'<NA>\' li diem que validated effect no sigui NA
  pp2<-dbGetQuery(con, pp2.1)
  efecte_validat<-pp2$validatedEffect[!is.na(pp2$validatedEffect)]
  efecte_novalidat<-pp2$effect[!is.na(pp2$effect)]
  #criteris["PP2",1][gen=="NF1"|gen=="SPRED1"]<-0 #si els gens son NF1 o SPRED1 ja assigna 0 a PP2
  #if(gen!="NF1"& gen!="SPRED1"){ 
  criteris["PP2",1][length(efecte_validat)==0&length(efecte_novalidat)==0]<-0
  if( length(efecte_validat)>0){
    if(efecte_validat[1]=="missense (non-synonymous)"){
      criteris["PP2", 1]<-1
    }
    if (efecte_validat[1]!="missense (non-synonymous)"){
      criteris["PP2", 1]<-0
    }}
  if (length(efecte_validat)==0 &length(efecte_novalidat)>0){ #te en compte les que no tenen valor en validated effect i mira el camp effect
    if (efecte_novalidat[1]=="nonsynonymous SNV"){
      criteris["PP2", 1]<-1
    }
    if(efecte_novalidat[1]!="nonsynonymous SNV"){
      criteris["PP2", 1]<-0
    }}
  
  #} 
  
  
  #pp3 i bp4
  
  pp3.1<- paste('SELECT "SIFT_pred", "Polyphen2_HVAR_pred","MutationTaster_pred","PROVEAN_pred", "CADD_phred" 
              FROM "VA_VariantsInTranscripts"
               LEFT OUTER JOIN "VA_InSilicoPathogenicity" AS a
               ON "VA_VariantsInTranscripts"."uniqueVariantId"= a."uniqueVariantId"
               WHERE "VA_VariantsInTranscripts"."uniqueVariantId"=', id, 'AND "isMainTranscript"=\'TRUE\' ORDER BY a."date" DESC LIMIT 1')
  pp3<-dbGetQuery(con, pp3.1)
  pp3$CADD_phred[pp3$CADD_phred>19]<-"D"
  pp3$CADD_phred[pp3$CADD_phred<10]<-"N"
  pp3$Polyphen2_HVAR_pred[pp3$Polyphen2_HVAR_pred=="B"]<-"N"
  pp3$Polyphen2_HVAR_pred[pp3$Polyphen2_HVAR_pred=="P"]<-"D"
  pp3$MutationTaster_pred[pp3$MutationTaster_pred=="A"]<-"D"
  pp3$MutationTaster_pred[pp3$MutationTaster_pred=="P"]<-"N"
  pp3$SIFT_pred[pp3$SIFT_pred=="T"]<-"N"
  pp3<-t(pp3)
  taula<-table(pp3[,1])
  criteris["PP3",1][taula["D"]==5]<-1
  criteris["PP3",1][taula["D"]<=4|taula["N"]==5]<-0
  criteris["BP4",1][taula["N"]==5]<-1
  criteris["BP4",1][taula["N"]<=4|taula["D"]==5]<-0
  
  
  #bs2, ps4
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
  bs2<-dbGetQuery(con, bs2.1)
  independent_control<-sum(bs2$name!="NF1_SPRED1"&bs2$name!="Rasopathies"|is.na(bs2$name))
  independent_afectats<-sum(bs2$name=="NF1_SPRED1"&!is.na(bs2$name)|bs2$name=="Rasopathies"&!is.na(bs2$name))
  
  criteris[c("BS2","PS4_strong", "PS4_moderate", "PS4_supporting"),1][independent_control>=3]<-c(1,0,0,0)
  criteris[c("PS4_strong", "PS4_moderate", "PS4_supporting"),1][independent_afectats==0]<-c(0,0,0)
  if (independent_control==0){
    criteris[c("BS2","PS4_strong","PS4_moderate", "PS4_supporting"),1][independent_afectats>=5]<-c(0,1,0,0)
    criteris[c("BS2","PS4_strong","PS4_moderate", "PS4_supporting"),1][independent_afectats==3|independent_afectats==4]<-c(0,0,1,0)
    criteris[c("BS2","PS4_strong","PS4_moderate", "PS4_supporting"),1][independent_afectats==2|independent_afectats==1]<-c(0,0,0,1)
    criteris[c("BS2","PS4_strong","PS4_moderate", "PS4_supporting"),1][independent_afectats==0]<-c(0,0,0,0)
  }
  if (independent_control<3&independent_control>0){
    print (paste("bs2: S'han trobat", independent_control, "individus sans amb la variant"))
    print (paste("ps4: S'han trobat", independent_afectats, "individus afectes amb la variant"))
  }
  
  
  ##adaptacio pm2 i BS1
  criteris["PM2",1][criteris["PM2",1]==1&independent_control>100]<-0 #per molt que segons Pandora no estigui descrit en bases de dades poblacionals, si en la mostra de pandora hi ha més de 100 individus control amb la variant PM2 es considera 0
  inhousefreq<-paste('SELECT "freqs","inHouseFrequency","inHouseNumber" FROM "VA_InHouseFrequencies" WHERE "uniqueVariantId"=', id,'ORDER BY "date" DESC LIMIT 1') #agafem la informacio de frequencia m?s actualitzada de la variant d'interes
  inhouse<-dbGetQuery(con,inhousefreq)
  if (is.na(BA1$PopFreqMax)){
    criteris["BS1",1][independent_control>=1000&inhouse$inHouseFrequency[1]>=0.30]<-1 #per molt que segons Pandora no estigui descrit en bases de dades poblacionals, si en la mostra de pandora hi ha més de 1000 individus control amb la variant i la frequencia inhouse es superior a 0.30 considerem 1 a BS1
    criteris["BS1_supporting",1][independent_control>500&independent_control<1000&inhouse$inHouseFrequency[1]>0.10&inhouse$inHouseFrequency[1]<0.30]<-1
    criteris["BS1_supporting",1][is.na(criteris["BS1_supporting",1])]<-0
  } 
  return(criteris)
}

criteris<-myFunction(id)


####CRITERIS MANUALS

manual<-function(criteris){
  
  #fem que ens pregunti els criteris manuals
  print(paste("En quantes de novo s'ha trobat?"))
  denovo<-as.integer(scan(n=1, what="integer"))
  if (denovo!=0){
    print ("En quantes d'aquestes s'ha estudiat als progenitors? Primer entra el nombre que tenen pare i mare confirmat i després el nombre que tenen només 1 dels dos confirmat")
    denovo_confirmats<-as.integer(scan(n=2, what="integer"))
  }
  if (denovo==0){
    denovo_confirmats=c(0,0)
  }
  
  
  denovo_noconfirmats=denovo-denovo_confirmats[1]
  criteris["PS2_veryStrong",1][denovo_confirmats[1]<2|denovo_confirmats[2]<3]<-0
  criteris[c("PS2_veryStrong","PS2"),1][denovo_confirmats[1]>=2]<-c(1,0)
  criteris[c("PS2_veryStrong", "PS2"),1][denovo_confirmats[2]>=3]<-c(1,0)
  criteris[c("PS2_veryStrong","PS2"),1][denovo_confirmats[2]==2&denovo_confirmats[1]==1]<-c(1,0)
  if (criteris["PS2_veryStrong",1]==0){
    criteris["PS2",1][denovo_confirmats[1]==1]<-1#si hi ha 1 evidencia amb pare i mare confirmats es compleix PS2
    criteris["PS2",1][denovo_confirmats[1]==0]<-0#si no hi ha evidencia no es compleix
  }
  criteris[c("PM6_veryStrong","PM6"),1][denovo_noconfirmats>=4]<-c(1,0)#es compleix PM6_veryStrong quan hi ha almenys 4 de novo no confirmades
  criteris["PM6_veryStrong",1][denovo_noconfirmats<4]<-0
  criteris["PM6_strong",1][denovo_noconfirmats==2|denovo_noconfirmats==3]<-1
  criteris["PM6_strong",1][denovo_noconfirmats!=2&denovo_noconfirmats!=3]<-0
  criteris["PM6",1][denovo_noconfirmats==1]<-1
  criteris["PM6",1][denovo_noconfirmats!=1]<-0
  
  
  print(paste("Hi ha cosegregació de la variant amb la malaltia? Indica els casos en que cosegrega"))
  cosegregacio<-as.integer(scan(n=1, what="integer"))
  criteris[c("PP1_strong", "PP1_moderate", "PP1_supporting"),1][cosegregacio>=7]<-c(1,0,0)
  criteris[c("PP1_strong", "PP1_moderate", "PP1_supporting"),1][cosegregacio==5|cosegregacio==6]<-c(0,1,0)
  criteris[c("PP1_strong", "PP1_moderate", "PP1_supporting"),1][cosegregacio==3|cosegregacio==4]<-c(0,0,1)
  criteris[c("PP1_strong", "PP1_moderate", "PP1_supporting"),1][cosegregacio<3]<-c(0,0,0)
  
  print("Hi ha evidència en alguna base de dades fiable de que la variant és patogènica o benigne? Primer indica s'hi ha de patogènica, on 1 es que  n'hi ha i 0 si no n'hi ha i després el mateix però amb benigne")
  evidencia<-as.integer(scan(n=2, what="integer"))
  criteris[c("PP5", "BP6"),1]<-c(evidencia[1], evidencia[2])
  
  
  for(i in 1:nrow(criteris)){
    #si en aquellafila criteri es NA
    while(is.na(criteris[i,]))
    {
      #ensenya a l'usuari quin criteri correspon 
      print(paste("compleix el criteri:",row.names(criteris)[i]))
      #demana que s'introdueixi
      input<-scan(n=1,what="integer")
      if (input>1){ #si és més gran que 1 surt el missatge d'error
        print("input no valid, només posar 0 (si no compleix) o 1 (si compleix)")
      }
      else{
        #escribe lo que te diga el usuario en la tabla
        criteris[i,]<-as.integer(input)}
    }
  }
  
  return(criteris)
}

criteris_manuals<-manual(criteris)

###CLASSIFICACIO FINAL

classificacio_final<-function(criteris_manuals){
  suma_criteris<- data.frame(verystrong_pat=sum(criteris[2,1], criteris[9,1], criteris[37,1]),strong_pat=sum(criteris[1,1], criteris[3:5,1],criteris[8,1], criteris[33,1]),moderate_pat=sum(criteris[6,1],criteris[10:13,1],criteris[15,1],criteris[34,1]), supporting_pat=sum(criteris[16:19,1], criteris[7,1], criteris[35,1]), very_strong_benign=sum(criteris[20,1]), strong_benign=sum(criteris[21:24,1]), support_benign=sum(criteris[25:32,1],criteris[36,1]))
  
  classificacio<- data.frame(pat=0, likely_pat=0, ben=0, likely_ben=0, vsd=0)
  
  classificacio$pat[suma_criteris$verystrong_pat>=1&suma_criteris$strong_pat>=1|
                      suma_criteris$verystrong_pat>=1&suma_criteris$moderate_pat>=2|
                      suma_criteris$verystrong_pat>=1&suma_criteris$moderate_pat==1&suma_criteris$supporting_pat==1|
                      suma_criteris$verystrong_pat>=1&suma_criteris$supporting_pat>=2|
                      suma_criteris$strong_pat>=2|
                      suma_criteris$strong_pat==1&suma_criteris$moderate_pat>=3|
                      suma_criteris$strong_pat==1&suma_criteris$moderate_pat==2&suma_criteris$supporting_pat>=2|
                      suma_criteris$strong_pat==1&suma_criteris$moderate==1&suma_criteris$supporting_pat>=4]<-1
  
  classificacio$likely_pat[suma_criteris$verystrong_pat==1&suma_criteris$moderate_pat==1|
                             suma_criteris$strong_pat==1&suma_criteris$moderate_pat==1|
                             suma_criteris$strong_pat==1&suma_criteris$moderate_pat==2|
                             suma_criteris$strong_pat==1&suma_criteris$supporting_pat>=2|
                             suma_criteris$moderate_pat>=3|
                             suma_criteris$moderate_pat==2&suma_criteris$supporting_pat>=2|
                             suma_criteris$moderate_pat==1&suma_criteris$supporting_pat>=4]<-1
  
  classificacio$ben[suma_criteris$very_strong_benign==1|
                      suma_criteris$strong_benign>=2]<-1
  
  classificacio$likely_ben[suma_criteris$strong_benign&suma_criteris$support_benign==1|
                             suma_criteris$support_benign>=2|
                             suma_criteris$strong_benign>=1&criteris["BS1",1]==1&suma_criteris$very_strong_benign==0&suma_criteris$supporting_pat==0&suma_criteris$moderate_pat==0&suma_criteris$strong_pat==0&suma_criteris$verystrong_pat==0]<-1
  classificacio$vsd[classificacio$pat==0&classificacio$likely_pat==0&classificacio$ben==0&classificacio$likely_ben==0|classificacio$pat==1&classificacio$ben==1|classificacio$pat==1&classificacio$likely_ben==1|classificacio$likely_pat==1&classificacio$ben==1|classificacio$likely_pat==1&classificacio$likely_ben==1]<-1
  
  
  if (classificacio$vsd==1){
    print("Variant de Significat Desconegut")
    print(suma_criteris)
  }
  if(classificacio$vsd==0){
    if (classificacio$pat==1){
      print("Variant Patogènica")
      print(suma_criteris)
    }
    if (classificacio$pat==0&classificacio$likely_pat==1){
      print ("Variant Probablement Patogènica")
      print(suma_criteris)
    }
    if (classificacio$ben==1){
      print ("Variant Benigne")
      print(suma_criteris)
    }
    if(classificacio$ben==0&classificacio$likely_ben==1){
      print ("Variant Probablement Benigne")
      print(suma_criteris)
      
    }
  }
}

classificacio(criteris_manuals)
