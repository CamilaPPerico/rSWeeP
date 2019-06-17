aa2int <- function(xfas){
  i=1
  for(i in 1:length(xfas)){
    if(xfas[i] == "a"||xfas[i] == "A"){
      xfas[i] <- 1 
    }
    if(xfas[i] == "r"||xfas[i] == "R"){
      xfas[i] <- 2 
    }
    if(xfas[i] == "n"||xfas[i] == "N"){
      xfas[i] <- 3 
    }
    if(xfas[i] == "d"||xfas[i] == "D"){
      xfas[i] <- 4 
    }
    if(xfas[i] == "c"||xfas[i] == "C"){
      xfas[i] <- 5 
    }
    if(xfas[i] == "q"||xfas[i] == "Q"){
      xfas[i] <- 6 
    }
    if(xfas[i] == "e"||xfas[i] == "E"){
      xfas[i] <- 7 
    }
    if(xfas[i] == "g"||xfas[i] == "G"){
      xfas[i] <- 8 
    }
    if(xfas[i] == "h"||xfas[i] == "H"){
      xfas[i] <- 9 
    }
    if(xfas[i] == "i"||xfas[i] == "I"){
      xfas[i] <- 10 
    }
    if(xfas[i] == "l"||xfas[i] == "L"){
      xfas[i] <- 11 
    }
    if(xfas[i] == "k"||xfas[i] == "K"){
      xfas[i] <- 12
    }
    if(xfas[i] == "m"||xfas[i] == "M"){
      xfas[i] <- 13
    }
    if(xfas[i] == "f"||xfas[i] == "F"){
      xfas[i] <- 14 
    }
    if(xfas[i] == "p"||xfas[i] == "P"){
      xfas[i] <- 15 
    }
    if(xfas[i] == "s"||xfas[i] == "S"){
      xfas[i] <- 16 
    }
    if(xfas[i] == "t"||xfas[i] == "T"){
      xfas[i] <- 17 
    }
    if(xfas[i] == "w"||xfas[i] == "W"){
      xfas[i] <- 18 
    }
    if(xfas[i] == "y"||xfas[i] == "Y"){
      xfas[i] <- 19 
    }
    if(xfas[i] == "v"||xfas[i] == "V"){
      xfas[i] <- 20 
    }
    if(xfas[i] == "b"||xfas[i] == "B"){
      xfas[i] <- 21 
    }
    if(xfas[i] == "z"||xfas[i] == "Z"){
      xfas[i] <- 22 
    }
    if(xfas[i] == "x"||xfas[i] == "X"){
      xfas[i] <- 23 
    }
    if(xfas[i] == "*"){
      xfas[i] <- 24 
    }
    if(xfas[i] == "-"){
      xfas[i] <- 25 
    }
  }
  xfas
}