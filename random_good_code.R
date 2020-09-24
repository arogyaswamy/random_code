### Good implementation of switch to rename entire columns; can be useful when renaming factors with one of the apply functions:
for(i in 2:3){for(j in 1:18){gradientprofiles[j,i] = switch(gradientprofiles[j,i], "1" = "Uniform", "2" = "Increasing", "3" = "Decreasing")}}

