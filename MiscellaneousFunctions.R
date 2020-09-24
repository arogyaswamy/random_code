title.case = function(inputstring){
  outputstring = 
    paste0(
      paste0(
        toupper(
          gsub("^(.).*", 
               "\\1",
               strsplit(
                 tolower(inputstring), split = " ")[[1]])),
        gsub("^(.)",
             "",
             strsplit(
               tolower(inputstring), split = " ")[[1]])),
      collapse = " ")
  return(outputstring)
}

sentence.case = function(inputstring){
  outputstring = 
    paste0(
      toupper(
        gsub("^(.).*", "\\1", inputstring)),
     
      tolower(
        gsub("^."    ,  ""  , inputstring)),
     
      collapse = " ")
  return(outputstring)
}