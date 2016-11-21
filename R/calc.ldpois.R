calc.ldpois <-
function(y,lambda) {
  -lambda + y*log(lambda) - lgamma(y+1)
 }
