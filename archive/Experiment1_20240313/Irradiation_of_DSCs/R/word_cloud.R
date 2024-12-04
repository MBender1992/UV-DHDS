install.packages(c("tm", "wordcloud","SnowballC","magrittr"))
library(tm)
library(wordcloud2)
library(SnowballC)
library(magrittr)
library(dplyr)

text <- readLines('Data/wc_trial.txt')
docs <- Corpus(VectorSource(text))
inspect(docs)
docs <- docs %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
docs <- tm_map(docs, content_transformer(tolower))
docs <- tm_map(docs, removeWords, stopwords("english")) ## remove stopwords is the important step
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
x <- names(words)
Encoding(x) <- 'latin1'
names(words) <- x
df <- data.frame(word = names(words),freq=words)
length(df$word)
wordcloud2(data=df, size=0.3, color='random-dark',shape = 'circle',minSize = 0.1)
#set.seed(1234) # for reproducibility 
#wordcloud(words = df$word, freq = df$freq, min.freq = 2,max.words=200,random.order=TRUE,rot.per=0.0,colors=brewer.pal(8, "Dark2"))

run1 <- c(240, 240, 220, 150, 225, 144, 240, 210, 174, 185, 157, 300, 270, 220, 270, 220, 260, 180, 270, 280, 270, 250, 260, 200, 300, 290)
run2 <- c(250, 290, 258, 180, 144, 220, 260, 140, 121, 158, 130, 200, 98, 125, 181, 153, 119, 95, 210, 157, 173, 168, 155, 114, 109, 125)

median(run2)


