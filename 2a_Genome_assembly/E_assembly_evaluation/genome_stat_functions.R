# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# SET OF FUNCTIONS FOR WORK WITH THE GENOMIC STAT FILES I COMPUTE USING THE PYTHON SCRIPT # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

require(AsexStats)

# prepare data frame (matrix -> na columns -> convert to frame without any rows)
LXcols <- paste0('L', seq(10,80, by = 10))
NXcols <- paste0('N', seq(10,80, by = 10))
variables <- c('file', 'dir', 'soft', 'kmer', 'kc', 'fasteris', 'fewdata', 'mse', 'pse', 'mpe', 'corrected',
               'total_sum', 'num_of_records', LXcols , NXcols, 'LG50', 'NG50')

# function for reading kmer and kc
# for k and kc the construct is the same, try to find kmer / kc and
# if you found it convert it to number
getk <- function(dir_name, pattern, start){
  k <- dir_name[grep(pattern, dir_name)]
  k <- ifelse(identical(k, character(0)),
              NA, as.numeric(substr(k, start, nchar(k))))
  return(k)
}

# key \y value line type parser
get_second <- function(stat_file, parameter){
  # grepl will find line with parameter
  # ssplit will cut the line by tab
  # conver the second evelement to number and return as the value
  return(as.numeric(ssplit(stat_file[grepl(parameter, stat_file)], '\t')[2]))
}

parse_stadard <- function(one_row, stat_file){
  # this construction is again bit ugly, but I have onliners that extract L[1-8]0 and Ns
  line <- stat_file[3:10]
  one_row[1, NXcols] <- unlist(lapply(strsplit(line, split = '\t'),
                                      function(x){return(as.numeric(x[2]))}))
  one_row[1, LXcols] <- unlist(lapply(strsplit(line, split = ':'),
                                      function(x){return(as.numeric(ssplit(x[2], "\t")[1]))}))

  one_row[1, 'LG50'] <- as.numeric(ssplit(stat_file[11], '\t')[2])
  one_row[1, 'NG50'] <- as.numeric(ssplit(stat_file[11], '\t')[3])

  return(one_row)
}

parse_old <- function(one_row, stat_file){
  for(parameter in c('N50', 'NG50', 'L50', 'LG50')){
    if(any(grepl(parameter, stat_file))){
      one_row[1, parameter] <- get_second(stat_file, parameter)
    }
  }
  return(one_row)
}


# gets a filename, reads it
make_dataframe_row <- function(filename, one_row){
        # read file & prepare table from template
        stat_file <- readLines(filename)
        one_row <- one_row[F,]

        line <- ssplit(filename, "/")
        one_row[1, 'file'] <- line[length(line)]

        dir_name <- line[4]
        one_row[1, 'dir'] <- dir_name
        one_row[1, 'fasteris'] <- grepl('fasteris',dir_name)
        one_row[1, 'fewdata'] <- grepl('fewdata',dir_name)
        ## ugly construct: sapply will return FF TF or FT if both se were used, pe only or none
        # which will return 1, 1 2 or 1 3 in these three cases
        # 3 - max of the value will quantify "how many" se reads were used (0 none, 1 pe, 2 both)
        one_row[1, 'mse'] <- !any(sapply(c('nomse','nose'), grepl, dir_name))
        one_row[1, 'pse'] <- !grepl('nose',dir_name)
        # was pe lib from mate pairs used? (flag nompe check just opposite therefore !)
        one_row[1, 'mpe'] <- !grepl('nompe',dir_name)
        one_row[1, 'corrected'] <- grepl('corrected',dir_name)

        # the software name, k and kc will be extracted using split dir name
        dir_name <- ssplit(dir_name, '_')
        one_row[1, 'soft'] <- ifelse(any(grepl('BESST', dir_name)),'BESST',dir_name[1])
        one_row[1, 'kmer'] <- getk(dir_name, 'k[1-9]', 2)
        one_row[1, 'kc'] <- getk(dir_name, 'kc', 3)

        one_row[1, 'total_sum'] <- get_second(stat_file, 'Total')
        one_row[1, 'num_of_records'] <- get_second(stat_file, 'records')

        # recognise OLD formating
        if(any(grepl('L20', stat_file))){
                 return(parse_stadard(one_row, stat_file))
        } else {
                 return(parse_old(one_row, stat_file))
        }

        return(one_row)
}
