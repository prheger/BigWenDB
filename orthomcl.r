#    orthomcl.r
#    the script performs those steps of the original OrthoMCL pipeline (Li et 
#    al., 2003) which use a MySQL database, within R
#
#    Copyright (C) 2019  Wen Zheng, Peter Heger
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#    Contact: Peter Heger, University of Cologne, Institute for Genetics, 
#    Zuelpicher Strasse 47a, 50674 Koeln, Germany
#    Email: peter.heger@uni-koeln.de, spacepeter@gmx.de


#!/usr/bin/Rscript

# Author: Wen Zheng

#orthomclpairs.r simi.txt.update
.libPaths("/home/wzheng/R/x86_64-unknown-linux-gnu-library/3.1") 
#library(RH2,lib.loc="/home/wzheng/R/x86_64-unknown-linux-gnu-library/3.1")
require(sqldf)

args<-commandArgs(TRUE)
data = read.table(args[1],sep="\t",header=F)
colnames(data) <- c("QUERY_ID","SUBJECT_ID","QUERY_TAXON_ID","SUBJECT_TAXON_ID","EVALUE_MANT","EVALUE_EXP","PERCENT_IDENTITY","PERCENT_MATCH")
min_evalue_exp = min(data[which(data[,"EVALUE_MANT"]!=0),"EVALUE_EXP"])
print(min_evalue_exp)

data[intersect(which(data[,"EVALUE_MANT"]==0),which(data[,"EVALUE_EXP"]==0)),"EVALUE_EXP"]=min_evalue_exp-1
interTaxonMatchView <- sqldf('SELECT QUERY_ID,SUBJECT_ID,SUBJECT_TAXON_ID,EVALUE_MANT,EVALUE_EXP FROM data WHERE QUERY_TAXON_ID != SUBJECT_TAXON_ID')

BestQueryTaxonScore <- sqldf('SELECT im.QUERY_ID, im.SUBJECT_TAXON_ID, low_EXP.EVALUE_EXP, min(im.EVALUE_MANT) as EVALUE_MANT from interTaxonMatchView im,
     (SELECT QUERY_ID, SUBJECT_TAXON_ID, min(EVALUE_EXP) as EVALUE_EXP
      from interTaxonMatchView
              group by QUERY_ID, SUBJECT_TAXON_ID) low_EXP
               where im.QUERY_ID = low_EXP.QUERY_ID
                and im.SUBJECT_TAXON_ID = low_EXP.SUBJECT_TAXON_ID
                and im.EVALUE_EXP = low_EXP.EVALUE_EXP group by im.QUERY_ID, im.SUBJECT_TAXON_ID, low_EXP.EVALUE_EXP')

################################################################################
############################### Orthologs  #####################################
################################################################################

BestHit2 <- sqldf('SELECT S.QUERY_ID, S.SUBJECT_ID,
                  S.QUERY_TAXON_ID, S.SUBJECT_TAXON_ID,
                  S.EVALUE_EXP, S.EVALUE_MANT
                  FROM data S, BestQueryTaxonScore CUTOFF
                  WHERE S.QUERY_ID = CUTOFF.QUERY_ID
                  AND S.SUBJECT_TAXON_ID = CUTOFF.SUBJECT_TAXON_ID
                  AND S.QUERY_TAXON_ID != S.SUBJECT_TAXON_ID
                  AND S.EVALUE_EXP <= -5
                  AND S.PERCENT_MATCH >= 50
                  AND (S.EVALUE_MANT < 0.01
                  OR S.EVALUE_EXP = CUTOFF.EVALUE_EXP
                  AND S.EVALUE_MANT = CUTOFF.EVALUE_MANT)')

OrthologTemp <- sqldf('SELECT BH1.QUERY_ID AS SEQUENCE_ID_A, BH1.SUBJECT_ID AS SEQUENCE_ID_B,
                      BH1.QUERY_TAXON_ID AS TAXON_ID_A, BH1.SUBJECT_TAXON_ID AS TAXON_ID_B,
                      CASE -- DONNOT TRY TO CALCULATE LOG(0) -- USE RIGGED EXPONENTS OF SIMSEQ
                      WHEN BH1.EVALUE_MANT < 0.01 OR BH2.EVALUE_MANT < 0.01
                      THEN (BH1.EVALUE_EXP + BH2.EVALUE_EXP) / -2
                      ELSE  -- SCORE = ( -log10(EVALUE1) - log10(EVALUE2) ) / 2
                      (log10(BH1.EVALUE_MANT * BH2.EVALUE_MANT)
                      + BH1.EVALUE_EXP + BH2.EVALUE_EXP) / -2
                      END AS UNNORMALIZED_SCORE
                      FROM BestHit2 BH1, BestHit2 BH2
                      WHERE BH1.QUERY_ID < BH1.SUBJECT_ID
                      AND BH1.QUERY_ID = BH2.SUBJECT_ID
                      AND BH1.SUBJECT_ID = BH2.QUERY_ID')

orthologTaxonSub <- function(a){
  if(a == ""){
    table_name <- OrthologTemp
  }else{
    table_name <- CoOrthologTemp
  }
  
  table <- sqldf('SELECT CASE
                 WHEN TAXON_ID_A < TAXON_ID_B
                 THEN TAXON_ID_A
                 ELSE TAXON_ID_B
                 END AS SMALLER_TAX_ID,
                 CASE
                 WHEN TAXON_ID_A < TAXON_ID_B
                 THEN TAXON_ID_B
                 ELSE TAXON_ID_A
                 END AS BIGGER_TAX_ID,
                 UNNORMALIZED_SCORE
                 FROM table_name')
  return(table)
}

normalizeOrthologsSub <- function(a){
  if(a == ""){
    table_name <- OrthologTaxon
    table_name2 <- paste0(a,"OrthologTable")
    table_name3 <- OrthologTemp
  }else{
    table_name <- CoOrthologTaxon
    table_name2 <- paste0(a,"OrthologTable")
    table_name3 <- CoOrthologTemp
  }
  
  table_ave_score <- sqldf('SELECT SMALLER_TAX_ID, BIGGER_TAX_ID, AVG(UNNORMALIZED_SCORE) AVG_SCORE
                           FROM table_name
                           GROUP BY SMALLER_TAX_ID, BIGGER_TAX_ID')
 orthologtable <- sqldf('SELECT OT.SEQUENCE_ID_A AS SEQUENCE_ID_A, OT.SEQUENCE_ID_B AS SEQUENCE_ID_B, OT.TAXON_ID_A AS TAXON_ID_A, OT.TAXON_ID_B AS TAXON_ID_B, OT.UNNORMALIZED_SCORE AS UNNORMALIZED_SCORE, OT.UNNORMALIZED_SCORE/A.AVG_SCORE AS NORMALIZED_SCORE
                        FROM table_name3 OT, table_ave_score A
                        WHERE MIN(OT.TAXON_ID_A, OT.TAXON_ID_B) = A.SMALLER_TAX_ID
                        AND MAX(OT.TAXON_ID_A, OT.TAXON_ID_B) = A.BIGGER_TAX_ID')
  return(list(table_ave_score,orthologtable))
}

OrthologTaxon <- orthologTaxonSub("")
OrthologAvgScore <- normalizeOrthologsSub("")[[1]]
OrthologTable <- normalizeOrthologsSub("")[[2]]
fn <- paste0(args[1],".orthol")
write.table(OrthologTable,file=fn,sep="\t",col.names=T, row.names=F)

################################################################################
############################### InParalogs  ####################################
################################################################################
BestInterTaxonScore <- sqldf('SELECT IM.QUERY_ID AS QUERY_ID, LOW_EXP.EVALUE_EXP AS EVALUE_EXP, MIN(IM.EVALUE_MANT) AS EVALUE_MANT
                             FROM BestQueryTaxonScore IM,
                             (SELECT QUERY_ID, MIN(EVALUE_EXP) AS EVALUE_EXP
                             FROM BestQueryTaxonScore
                             GROUP BY QUERY_ID) LOW_EXP
                             WHERE IM.QUERY_ID = LOW_EXP.QUERY_ID
                             AND IM.EVALUE_EXP = LOW_EXP.EVALUE_EXP
                             GROUP BY IM.QUERY_ID, LOW_EXP.EVALUE_EXP')

#here set WHERETAXONFILTER="",it can be other options, please see the original perl code. :P
WHERETAXONFILTER=""
UniqSimSeqsQueryId <- sqldf('SELECT DISTINCT S.QUERY_ID AS QUERY_ID FROM data S')

# $ANDTAXONFILTER=""
BetterHit <- sqldf('SELECT S.QUERY_ID AS QUERY_ID, S.SUBJECT_ID AS SUBJECT_ID,
                   S.QUERY_TAXON_ID AS TAXON_ID,
                   S.EVALUE_EXP AS EVALUE_EXP, S.EVALUE_MANT AS EVALUE_MANT
                   FROM data S, BestInterTaxonScore BIS
                   WHERE S.QUERY_ID != S.SUBJECT_ID
                   AND S.QUERY_TAXON_ID = S.SUBJECT_TAXON_ID
                   AND S.QUERY_ID = BIS.QUERY_ID
                   AND S.EVALUE_EXP <= -5
                   AND S.PERCENT_MATCH >= 50
                   AND (S.EVALUE_MANT < 0.001
                   OR S.EVALUE_EXP < BIS.EVALUE_EXP
                   OR (S.EVALUE_EXP = BIS.EVALUE_EXP AND S.EVALUE_MANT <= BIS.EVALUE_MANT))
                   -- . . . OR SIMILARITY FOR A PROTEIN WITH NO BESTINTERTAXONSCORE
                   --       (I.E. AN INTRATAXON MATCH FOR A PROTEIN WITH NO INTERTAXON
                   --        MATCH IN THE DATABASE)
                   UNION
                   SELECT S.QUERY_ID, S.SUBJECT_ID, S.QUERY_TAXON_ID AS TAXON_ID, S.EVALUE_EXP, S.EVALUE_MANT
                   FROM data S
                   WHERE S.QUERY_TAXON_ID = S.SUBJECT_TAXON_ID
                   AND S.EVALUE_EXP <= -5
                   AND S.PERCENT_MATCH >= 50
                   AND S.QUERY_ID IN 
                   (SELECT DISTINCT UST.QUERY_ID
                   FROM UniqSimSeqsQueryId UST
                   LEFT OUTER JOIN BestInterTaxonScore BIS ON BIS.QUERY_ID = UST.QUERY_ID
                   WHERE BIS.QUERY_ID IS NULL)')
fn <- paste0(args[1],".betterhit")
write.table(BetterHit,file=fn,sep="\t",col.names=T, row.names=F)

InParalogTemp <-  sqldf('SELECT BH1.QUERY_ID AS SEQUENCE_ID_A, BH1.SUBJECT_ID AS SEQUENCE_ID_B,
                        BH1.TAXON_ID AS TAXON_ID,
                        CASE -- DONNOT TRY TO CALCULATE LOG(0) -- USE RIGGED EXPONENTS OF SIMSEQ
                        WHEN BH1.EVALUE_MANT < 0.01 OR BH2.EVALUE_MANT < 0.01
                        THEN (BH1.EVALUE_EXP + BH2.EVALUE_EXP) / -2
                        ELSE  -- SCORE = ( -log10(EVALUE1) - log10(EVALUE2) ) / 2
                        (log10(BH1.EVALUE_MANT * BH2.EVALUE_MANT)
                        + BH1.EVALUE_EXP + BH2.EVALUE_EXP) / -2
                        END AS UNNORMALIZED_SCORE
                        FROM BetterHit BH1, BetterHit BH2
                        WHERE BH1.QUERY_ID < BH1.SUBJECT_ID
                        AND BH1.QUERY_ID = BH2.SUBJECT_ID
                        AND BH1.SUBJECT_ID = BH2.QUERY_ID')

InParalogTaxonAvg <- sqldf('SELECT AVG(I.UNNORMALIZED_SCORE) AVERAGE, I.TAXON_ID AS TAXON_ID FROM InParalogTemp I
                           GROUP BY I.TAXON_ID')

OrthologUniqueId <- sqldf('SELECT DISTINCT(SEQUENCE_ID) FROM (
                          SELECT SEQUENCE_ID_A AS SEQUENCE_ID FROM OrthologTable
                          UNION
                          SELECT SEQUENCE_ID_B AS SEQUENCE_ID FROM OrthologTable) I')

InplgOrthTaxonAvg <- sqldf('SELECT AVG(I.UNNORMALIZED_SCORE) AVERAGE, I.TAXON_ID AS TAXON_ID
                           FROM InParalogTemp I
                           WHERE I.SEQUENCE_ID_A IN
                           (SELECT SEQUENCE_ID FROM OrthologUniqueId)
                           OR I.SEQUENCE_ID_B IN
                           (SELECT SEQUENCE_ID FROM OrthologUniqueId)
                           GROUP BY I.TAXON_ID')

InParalogAvgScore <- sqldf('SELECT CASE
                           WHEN ORTH_I.AVERAGE IS NULL
                           THEN ALL_I.AVERAGE
                           ELSE ORTH_I.AVERAGE
                           END AS AVG_SCORE,
                           ALL_I.TAXON_ID
                           FROM InParalogTaxonAvg ALL_I LEFT OUTER JOIN InplgOrthTaxonAvg ORTH_I
                           ON ALL_I.TAXON_ID = ORTH_I.TAXON_ID')

InParalogTable <- sqldf('SELECT IT.SEQUENCE_ID_A AS SEQUENCE_ID_A, IT.SEQUENCE_ID_B AS SEQUENCE_ID_B, IT.TAXON_ID AS TAXON_ID, IT.UNNORMALIZED_SCORE AS UNNORMALIZED_SCORE, IT.UNNORMALIZED_SCORE/A.AVG_SCORE AS NORMALIZED_SCORE
                        FROM InParalogTemp IT, InParalogAvgScore A
                        WHERE IT.TAXON_ID = A.TAXON_ID')

fn <- paste0(args[1],".inpara")
write.table(InParalogTable,file=fn,sep="\t",col.names=T, row.names=F)

################################################################################
############################### CoOrthologs  ###################################
################################################################################

InParalog2Way <- sqldf('SELECT SEQUENCE_ID_A, SEQUENCE_ID_B FROM InParalogTable
                       UNION
                       SELECT SEQUENCE_ID_B AS SEQUENCE_ID_A, SEQUENCE_ID_A AS SEQUENCE_ID_B FROM InParalogTable')

Ortholog2Way <- sqldf('SELECT SEQUENCE_ID_A, SEQUENCE_ID_B FROM OrthologTable
                      UNION
                      SELECT SEQUENCE_ID_B AS SEQUENCE_ID_A, SEQUENCE_ID_A AS SEQUENCE_ID_B FROM OrthologTable')

InplgOrthoInplg <- sqldf('SELECT IP1.SEQUENCE_ID_A AS SEQUENCE_ID_A, IP2.SEQUENCE_ID_B AS SEQUENCE_ID_B
                         FROM  Ortholog2Way O, InParalog2Way IP2, InParalog2Way IP1
                         WHERE IP1.SEQUENCE_ID_B = O.SEQUENCE_ID_A
                         AND O.SEQUENCE_ID_B = IP2.SEQUENCE_ID_A')

InParalogOrtholog <- sqldf('SELECT IP.SEQUENCE_ID_A AS SEQUENCE_ID_A, O.SEQUENCE_ID_B AS SEQUENCE_ID_B
                           FROM InParalog2Way IP, Ortholog2Way O
                           WHERE IP.SEQUENCE_ID_B = O.SEQUENCE_ID_A')

CoOrthologCandidate <- sqldf('SELECT DISTINCT
                             MIN(SEQUENCE_ID_A, SEQUENCE_ID_B) AS SEQUENCE_ID_A,
                             MAX(SEQUENCE_ID_A, SEQUENCE_ID_B) AS SEQUENCE_ID_B
                             FROM (SELECT SEQUENCE_ID_A, SEQUENCE_ID_B FROM InplgOrthoInplg
                             UNION
                             SELECT SEQUENCE_ID_A, SEQUENCE_ID_B FROM InParalogOrtholog) T')

CoOrthNotOrtholog <- sqldf('SELECT CC.SEQUENCE_ID_A AS SEQUENCE_ID_A, CC.SEQUENCE_ID_B AS SEQUENCE_ID_B
                           FROM CoOrthologCandidate CC
                           LEFT OUTER JOIN OrthologTable O
                           ON CC.SEQUENCE_ID_A = O.SEQUENCE_ID_A
                           AND CC.SEQUENCE_ID_B = O.SEQUENCE_ID_B
                           WHERE O.SEQUENCE_ID_A IS NULL')

CoOrthologTemp <- sqldf('SELECT CANDIDATE.SEQUENCE_ID_A AS SEQUENCE_ID_A, CANDIDATE.SEQUENCE_ID_B AS SEQUENCE_ID_B,
                        AB.QUERY_TAXON_ID AS TAXON_ID_A, AB.SUBJECT_TAXON_ID AS TAXON_ID_B,
                        CASE  -- IN CASE OF 0 EVALUE, USE RIGGED EXPONENT
                        WHEN AB.EVALUE_MANT < 0.00001 OR BA.EVALUE_MANT < 0.00001
                        THEN (AB.EVALUE_EXP + BA.EVALUE_EXP) / -2
                        ELSE -- SCORE = ( -log10(EVALUE1) - log10(EVALUE2) ) / 2
                        (log10(AB.EVALUE_MANT * BA.EVALUE_MANT)
                        + AB.EVALUE_EXP + BA.EVALUE_EXP) / -2
                        END AS UNNORMALIZED_SCORE
                        FROM data AB, data BA, CoOrthNotOrtholog CANDIDATE
                        WHERE AB.QUERY_ID = CANDIDATE.SEQUENCE_ID_A
                        AND AB.SUBJECT_ID = CANDIDATE.SEQUENCE_ID_B
                        AND AB.EVALUE_EXP <= -5
                        AND AB.PERCENT_MATCH >= 50
                        AND BA.QUERY_ID = CANDIDATE.SEQUENCE_ID_B
                        AND BA.SUBJECT_ID = CANDIDATE.SEQUENCE_ID_A
                        AND BA.EVALUE_EXP <= -5
                        AND BA.PERCENT_MATCH >= 50')

CoOrthologTaxon <- orthologTaxonSub("Co")
CoOrthologAvgScore <- normalizeOrthologsSub("Co")[[1]]
CoOrthologTable <- normalizeOrthologsSub("Co")[[2]]
fn <- paste0(args[1],".coorthol")
write.table(CoOrthologTable,file=fn,sep="\t",col.names=T, row.names=F)
