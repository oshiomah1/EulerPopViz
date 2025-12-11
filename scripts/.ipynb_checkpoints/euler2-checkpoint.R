
# List of required packages
packages <- c(
  "eulerr",
  "UpSetR",
  "jsonlite",
  "gridExtra",
  "cowplot",
  "magick",
  "grid"
)

# Install missing packages
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load packages
lapply(packages, library, character.only = TRUE)

# Define all population selections in a named list
pop_sets <- list(
  "African-American + Proxies" = c("European", "Eastern_Bantu", "West_Central_African", "African-American","Indigenous_American"),
  "Afro-South American + Proxies" = c("European", "Eastern_Bantu", "West_Central_African", "Afro-South_American","Indigenous_American"),
  "Afro-Caribbean + Proxies" = c("European", "Eastern_Bantu", "West_Central_African", "Afro-Caribbean","Indigenous_American"),
  "KS Components" = c("European", "Khoe-San", "Eastern_Bantu", "Asian"),
  "East Africa / Back to Africa" = c("European", "Eastern_Bantu", "Back-to-African", "Fulani", "Northeast_African"),
  "Major Continental Africa" = c("Khoe-San", "Eastern_Bantu", "West_Central_African", "Northeast_African", "Rainforest_Hunter-Gatherer")
)

# Define your label replacements ONCE
label_map <- c(
  "European" = "European",
  "Eastern_Bantu" = "Eastern Bantu",
  "West_Central_African" = "West Central African",
  "African-American" = "African-American",
  "Indigenous_American" = "Indigenous American",
  "Afro-South_American" = "Afro-South American",
  "Afro-Caribbean" = "Afro Caribbean",
  "Asian" = "Asian",
  "Khoe-San" = "Khoe-San",
  "Back-to-African" = "Back to African",
  "Fulani" = "Fulani",
  "Northeast_African" = "Northeast African",
  "Rainforest_Hunter-Gatherer" = "Rainforest Hunter-Gatherer"
)




sum_snps <- function(pop_combo, pops, counts, cutoff) {
  pop_tf <- rep(FALSE, length(pops))
  pop_tf[pop_combo] <- TRUE
  populations <- pops[pop_combo]
  if (length(pop_combo) > 1) {
    shared_common_snps <- which(rowSums(counts[,populations]>=cutoff)==length(populations))
  } else {
    shared_common_snps <- which(counts[,populations]>=cutoff)
  }
  snp_counts <- data.frame(t(c(pop_tf, sum(counts[shared_common_snps, ncol(counts)]))))
  colnames(snp_counts) <- c(pops, "shared_common_snps")
  return(snp_counts)
}

generate_euler_plot <- function(pairwise_file, poplist, selected_pops=c(), cutoff=3, common_pop="") {
  data <- read.table(
    file=pairwise_file,
    col.names=c("geovar_code", "counts"),
    colClasses=c("character", "numeric")
  )
  pops <- read.table(file=poplist)$V1
  if (length(selected_pops) < 1) {
    selected_pops <- pops
  }
  
  split.pops <- sapply(data[,1], function(a){strsplit(as.character(a),"")[[1]]})
  split.pops <- apply(split.pops,2,as.numeric)
  split.pops <- as.data.frame(t(split.pops))
  row.names(split.pops) <- 1:nrow(split.pops)
  colnames(split.pops) <- c(as.character(pops))
  data <- cbind(split.pops, data)
  data <- data[,c(selected_pops,"geovar_code","counts")]
  
  if (common_pop!=""){
    data <- data[which(data[common_pop]>=cutoff),]
  }
  
  combo_func <- Map(combn, list(1:length(selected_pops)), seq_along(1:length(selected_pops)), simplify=FALSE)
  combinations <- unlist(combo_func, recursive=FALSE)
  
  shared_common_snps <- do.call(rbind,lapply(combinations, FUN=sum_snps, pops=selected_pops, counts=data, cutoff=cutoff))
  shared_common_snps[nrow(shared_common_snps),"unique_snps"] <- shared_common_snps[nrow(shared_common_snps),"shared_common_snps"]
  # print("here are the common snps")
  #  print(shared_common_snps)
  
  for (i in (nrow(shared_common_snps)-1):1) {
    shared_common_snps[i,"unique_snps"] <- shared_common_snps[i,"shared_common_snps"] - sum(
      shared_common_snps[
        apply(
          shared_common_snps[,1:(ncol(shared_common_snps)-2)] - shared_common_snps[rep(i,nrow(shared_common_snps)),1:(ncol(shared_common_snps)-2)], 
          MARGIN=1, 
          FUN=function(x){!any(x<0)}
        ),"unique_snps"
      ], na.rm=TRUE
    )
  }
  
  euler_data <- shared_common_snps$unique_snps
  print(shared_common_snps)
  #names(euler_data) <- apply(shared_common_snps, 1, function(x){paste(head(names(x),-2)[which(x==1)],collapse="&")}) #this is where code breaks, edit this so it applyies to n - ith column
  names(euler_data) <- apply(shared_common_snps[selected_pops], 1, function(x){paste(names(x)[which(x==1)],collapse="&")}) #new fix
  #print(names(euler_data))
  return(list("sets"=euler_data,"euler"=eulerr::euler(euler_data, shape = "ellipse"),"size"=sum(shared_common_snps$unique_snps)))
}


plot_euler_and_upset <- function(name, pops) {
  raw <- generate_euler_plot(
    "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/caapa_chr21.txt",
    "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/pop_list.txt",
    cutoff = 3,
    selected_pops = pops
  )
  
  # ✅ Euler Plot WITH pretty labels
  plot(
    raw$euler,
    labels = label_map,
    fills = list(alpha = 0.6),
    edges = TRUE,
    main = paste("Euler Plot -", name)
  )
  
  # ✅ UpSet Plot WITHOUT renaming labels (raw set names preserved)
  UpSetR::upset(
    fromExpression(raw$sets),
    order.by = "freq",
    mainbar.y.label = "Number of Common Variants",
    sets.x.label = "Set Size",
    #main = paste("UpSet Plot -", name),
    text.scale = c(2, 1.5, 1.5, 1.2, 1.5, 2)
  )
}

# UpSetR::upset(
#   fromExpression(raw_euler_1$sets),
#   order.by = "freq",
#   show.numbers = "no",
#   mainbar.y.label = "Number of Common Variants",
#   sets.x.label = "Set Size",
#   text.scale = c(2, 1.5, 1.5, 1.2, 1.5, 2)
# )
plot_euler_and_upset("African-American + Proxies", pop_sets[["African-American + Proxies"]])
plot_euler_and_upset("Afro-South American + Proxies", pop_sets[["Afro-South American + Proxies"]])
plot_euler_and_upset("Afro-Caribbean + Proxies", pop_sets[["Afro-Caribbean + Proxies"]])
plot_euler_and_upset("KS Components", pop_sets[["KS Components"]])
plot_euler_and_upset("East Africa / Back to Africa", pop_sets[["East Africa / Back to Africa"]])
plot_euler_and_upset("Major Continental Africa", pop_sets[["Major Continental Africa"]])




