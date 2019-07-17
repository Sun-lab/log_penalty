
Create_Updating_Order = function(ChromosomesOrder = c(1:22), updateNum = 20){
  
  O0 = ChromosomesOrder 
  updating = list()
  numOr = 1
  updating[[numOr]] = O0
  numOr = numOr + 1
  
  while(numOr <= updateNum){
    O1 = sample(ChromosomesOrder)
    
    if(sum(abs(O1 - O0))!=0){
      updating[[numOr]] = O1
      numOr = numOr + 1
      O0 = O1
    }
  }    
  return(updating)
}
  
  








