getwd()
write.table(DNase, file = "data/DNase.tsv", sep = "\t", row.names = FALSE)

write.table(ChickWeight, file = "data/ChickWeight.csv", sep = ",", row.names = FALSE)

write.table(chickwts, file = "data/ChickWeight_feed.txt", sep = " ", row.names = FALSE, quote = FALSE)
# export chickwts data as xlsx file:
library(writexl)
write_xlsx(chickwts, path = "data/ChickWeight_feed.xlsx")

