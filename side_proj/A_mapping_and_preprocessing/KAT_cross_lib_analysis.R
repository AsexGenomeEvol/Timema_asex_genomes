sp <- '4_Tte'
ind <- 'Tte_01'

path <- paste('data', sp, 'lib_cross_test', ind, sep = '/')
stat_files <- dir(path, pattern = 'stats')

libs <- data.frame(lib1 = as.character(sapply(stat_files, substr, 6, 33)),
                   lib2 = as.character(sapply(stat_files, substr, 55, 82)),
                   stringsAsFactors=FALSE)

lib_names <- unique(c(as.character(libs$lib1), as.character(libs$lib2)))
dim <- length(lib_names)
stats <- lapply(stat_files, function(x)(readLines(paste(path, x, sep = '/'))))

get_value_from_line <- function(line, col = 5, sep = ' '){
    as.numeric(strsplit(line, sep)[[1]][col])
}

# extract_numbers <- function(stat_file, lib = 1){
#     total_kmers <- get_value_from_line(stat_file[5 + lib], 5)
#     shared_kmers <- get_value_from_line(stat_file[21 + lib], 9)
#     get_value_from_line(stat_file[38], 5)
# }

libs$dist <- sapply(stats, function(x) get_value_from_line(x[38], 5))
libs$eucl_dist <- sapply(stats, function(x) get_value_from_line(x[35], 5))

heat_matrix <- matrix(0, ncol = dim, nrow = dim)
colnames(heat_matrix) <- lib_names
rownames(heat_matrix) <- lib_names

for(i_1 in 1:dim){
    for(i_2 in i_1:dim){
        if(i_1 == i_2){
            next
        }
        lib1 <- lib_names[i_1]
        lib2 <- lib_names[i_2]
        distance <- libs[libs$lib1 == lib1 & libs$lib2 == lib2, 'eucl_dist']
        heat_matrix[i_1, i_2] <- distance
        heat_matrix[i_2, i_1] <- distance
    }
}

heatmap.2(heat_matrix)

library(plotly)
p <- plot_ly(x = lib_names, y = lib_names,
             z = heat_matrix, type = "heatmap")

htmlwidgets::saveWidget(p, "Tte_ReSeq_libs.html")
