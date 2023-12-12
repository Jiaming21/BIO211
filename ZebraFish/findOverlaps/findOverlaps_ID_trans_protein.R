trans <- read.csv('/Users/jiaming/Desktop/BIO211_pre/ZebraFish/Transcriptomics/Datas/PvD_trans.csv', header = T, stringsAsFactors = F)
trans_id <- trans$id[trans$symbol %in% c("agr2","aldh1a2","aldh1l1","bgnb","bhmt","c4","capn12","fabp1b.2","fgfbp2a","gyg1b","hapln1a","hgd","hsd17b7","mbpa","muc5.2","myh11b","paplnb","ppp1r1c","si:dkey-65b12.6","si:dkeyp-93a5.3","vcanb","vim","zgc:136930","zgc:172244","and1","and2","anxa1c","bco1","ca2","krt93","krt94","tgfbi")]

pros <- read.csv('/Users/jiaming/Desktop/BIO211_pre/ZebraFish/Proteomics/Datas/proteins.csv', header = T, stringsAsFactors = F)
pros_Uniprot <- pros$Uniprot[pros$symbol %in% c("agr2","aldh1a2","aldh1l1","bgnb","bhmt","c4","capn12","fabp1b.2","fgfbp2a","gyg1b","hapln1a","hgd","hsd17b7","mbpa","muc5.2","myh11b","paplnb","ppp1r1c","si:dkey-65b12.6","si:dkeyp-93a5.3","vcanb","vim","zgc:136930","zgc:172244","and1","and2","anxa1c","bco1","ca2","krt93","krt94","tgfbi")]

write.csv(trans_id, '/Users/jiaming/Desktop/BIO211_pre/ZebraFish/Transcriptomics/Datas/trans_id.csv')
write.csv(pros_Uniprot, '/Users/jiaming/Desktop/BIO211_pre/ZebraFish/Transcriptomics/Datas/pros_Uniprot.csv')







