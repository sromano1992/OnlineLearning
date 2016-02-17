kendallTau<-function (obs1,obs2){  
  #Kendall tau
  comparison =  names(obs1)==names(obs2);
  concordants = length(comparison[comparison==TRUE]);
  discordants = length(comparison[comparison==FALSE]);
  n = length(tmp);
  tau = (concordants-discordants)/(n/2*(n-1));
}
