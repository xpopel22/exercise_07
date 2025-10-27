#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("Biostrings")

#library("Biostrings")

scoreDNA <- readDNAStringSet("seq_score.fasta")
motifDNA <- readDNAStringSet("seq_motif.fasta")

Score <- function(start_idx, sequences, motif_length){
  len_seq <- length(sequences[[1]])
  for (i in 1:length(start_idx)){
    if ((start_idx[i] + motif_length) > len_seq){
      #print("prilis dlouhy motif pro index")
      return(-1)
    }
  }
  align_matrix <- list()
  for (i in 1:length(start_idx)){
    seq <- DNAString(sequences[[i]][start_idx[i]:(start_idx[i]+motif_length)])
    if (length(align_matrix) == 0){
      align_matrix <- list(seq)
    }else{
      align_matrix <- append(align_matrix, list(seq))
    }
  }
  align_matrix <- DNAStringSet(align_matrix)
  con_mat <- consensusMatrix(align_matrix, baseOnly=TRUE)
  score <- 0
  for (i in 1:(motif_length+1)){
    score <- score + max(con_mat[,i])
  }
  con_str <- consensusString(con_mat)
  return(score)
}

ScoreModif <- function(start_idx, num_seq, sequences, motif_length){
  len_seq <- length(sequences[[1]])
  for (i in 1:num_seq){
    if ((start_idx[i] + motif_length) > len_seq){
      #print("prilis dlouhy motif pro index")
      return(-1)
    }
  }
  align_matrix <- list()
  for (i in 1:num_seq){
    seq <- DNAString(sequences[[i]][start_idx[i]:(start_idx[i]+motif_length)])
    if (length(align_matrix) == 0){
      align_matrix <- list(seq)
    }else{
      align_matrix <- append(align_matrix, list(seq))
    }
  }
  align_matrix <- DNAStringSet(align_matrix)
  con_mat <- consensusMatrix(align_matrix, baseOnly=TRUE)
  score <- 0
  for (i in 1:(motif_length+1)){
    score <- score + max(con_mat[,i])
  }
  con_str <- consensusString(con_mat)
  return(score)
}

NextLeaf <- function(s, t, k){
  for (i in t:-1:1){
    if (s[i] < k){
      s[i] <- s[i] + 1
      return(s)
    }
    s[i] <- 1
  }
  return(s)
}

BFMotifSearch <- function(DNA, t, n, l){
  s <- rep(1, t)
  bestScore <- Score(s, DNA, l)
  while (TRUE){
    s <- NextLeaf(s, t, n - l + 1)
    tmp_score <- Score(s, DNA, l)
    if (tmp_score > bestScore){
      bestScore <- tmp_score
      bestMotif <- s
    }
    if (all(s == rep(1, t))){
      return(bestMotif)
    }
  }
}

NextVertex <- function(s, i, t, k){
  if (i < t){
    s[(i+1)] <- 1
    return(s, i + 1)
  }else{
    for (j in t:-1:1){
      if (s[j] < k){
        s[j] <- s[j] + 1
        return(s, j)
      }
    }
  }
  return(s, 0)
}

ByPass <- function(s, i, t, k){
  for (j in i:-1:1){
    if (s[j] < k){
      s[j] <- s[j] + 1
      return(s, j)
    }
  }
  return(c(s, 0))
}

BBMotifSearch(DNA, t, n, l){
  s <-  rep(1, t)
  bestScore <- 0
  i <- 1
  while (i > 0){
    if (i < t){
      optimisticScore <- ScoreModif(s, i, DNA, l) + (t - i) * l
      if (optimisticScore < bestScore){
        si <- ByPass(s, i, t, (n-l+1))
      }
      
    }
  }
}

motif_length <- 6
starting_indexes <- c(6, 4, 5, 3)
align_matrix <- Score(starting_indexes, scoreDNA, motif_length)
BFMotifSearch(motifDNA, length(motifDNA[,1]), length(motifDNA[[1]]), 3)
