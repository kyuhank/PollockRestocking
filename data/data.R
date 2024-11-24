
#########################
## make input data obj ##
#########################

MakeInputObj=function (nProj=14,               ## to 2035 (roughly two gen)
                       nRelease=0,
                       nReleasePeriod=10,      ## from 2015 to 2024
                       nogariL05=27,
                       nogariL95=30,
                       AdultFraction=0.5,
                       L1=16.69,               ## from Kim and Hyun, 2018
                       L2=47.49,               ## from Kim and Hyun, 2018
                       rho=0.88,               ## from Kim and Hyun, 2018
                       s1=1.05,                ## from Kim and Hyun, 2018
                       s2=3.22,                ## from Kim and Hyun, 2018
                       bins=13:57,
                       A=9,
                       sigR=0.6,
                       sigL05=1.6,
                       sigL95=1.6,
                       gamma=1,
                       ThresCollapse=0.1,
                       FindBrps=0,
                       steepPoint=0.2,
                       RdevCor=0,
                       L05devCor=0,
                       L95devCor=0,
                       JuvMaxVul=1,
                       AduMaxVul=1,
                       L05_low=13,
                       L05_up=15,
                       L95_low=24,
                       L95_up=27,
                       PmixL=5,
                       NatExponent=-1,
                       JuvCollapseYear=52,
                       AduCollapseYear=62,
                       nplatoons=5,
                       LrefForM=37,               ## approx. length at 95% maturity
                       Rprop=c(3.1, 23.7, 46.4, 23.7, 3.1)*0.01
                       ) {
  
  ### catch data from FAO and Stats Korea websites
  
  Catch=structure(c(10149.25373, 12089.55224, 11194.02985, 10597.01493, 
                    11194.02985, 12686.56716, 18358.20896, 17462.68657, 13283.58209, 
                    27462.68657, 30298.50746, 38208.95522, 38805.97015, 20895.52239, 
                    14925.37313, 13283.58209, 27313.43284, 21940.29851, 20149.25373, 
                    26268.65672, 20597.01493, 16865.67164, 28208.95522, 9552.238806, 
                    13418, 11241, 40492, 42628, 64512, 59862, 88102, 122851, 104318, 
                    79872, 96384, 165837, 137657, 85909, 106678, 84545, 79373, 33719, 
                    16240, 24217, 26325, 20220, 14545, 16610, 10748, 9165, 8270, 
                    7283, 6232, 1392, 766, 207, 215, 242, 64, 25, 60, 35, 0, 1, 1, 
                    1, 1, 1, 2, 3, 6, 1, 9, 0, 0, 0, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, 4540, 5546, 18584, 11295, 11716, 28112, 
                    50283, 38414, 29642, 39906, 46496, 46890, 20162, 13348, 15786, 
                    9798, 10104, 9504, 9043, 7605, 6903, 4445, 6373, 6232, 1392, 
                    766, 207, 215, 242, 64, 25, 60, 35, 0, 1, 1, 1, 1, 1, 2, 3, 6, 
                    1, 9, 0, 0, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, 55322, 82556, 104267, 93023, 68156, 68272, 115554, 99243, 
                    56267, 66772, 38049, 32483, 13557, 2892, 8431, 16527, 10116, 
                    5041, 7567, 3143, 2262, 3825, 910, 0, 0, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA), dim = c(76L, 3L), dimnames = list(c("1946", "1947", "1948", 
                                                             "1949", "1950", "1951", "1952", "1953", "1954", "1955", "1956", 
                                                             "1957", "1958", "1959", "1960", "1961", "1962", "1963", "1964", 
                                                             "1965", "1966", "1967", "1968", "1969", "1970", "1971", "1972", 
                                                             "1973", "1974", "1975", "1976", "1977", "1978", "1979", "1980", 
                                                             "1981", "1982", "1983", "1984", "1985", "1986", "1987", "1988", 
                                                             "1989", "1990", "1991", "1992", "1993", "1994", "1995", "1996", 
                                                             "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", 
                                                             "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", 
                                                             "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", 
                                                             "2021"), c("Total", "Adult", "Juv")))
  
  
  
  ## catch Separation
  Totals=Catch[1:29,1]
  Catch[1:29,2]=Totals*AdultFraction
  Catch[1:29,3]=Totals*(1-AdultFraction)
  
  row.names(Catch)=1946:2021
  
  Catch[is.na(Catch)]=0
  
  
  Release=rep(0, dim(Catch)[1]+nProj)
  Release[70:(70+nReleasePeriod-1)]=rep(nRelease,nReleasePeriod)
  
  
  #######################################################
  ####### Length dist (from Kim and Hyun 2018) ##########
  #######################################################
  
  
  a=1:A
  
  width=(bins[2]-bins[1])
  
  pMeanLengths=matrix(0, nrow=nplatoons, ncol=A)
  pSDs=matrix(0, nrow=nplatoons, ncol=A)
  
  SDs=c()
  MeanLengths=c()
  sigBetween=c()
  sigWithin=c()
  
    for (i in 1:A){
      SDs[i]=s1+(s2-s1)*(i-1)/(A-1)
      MeanLengths[i]=(L1+(L2-L1)*(1-rho^(i-1))/(1-rho^(A-1)))
      
      sigBetween[i]=(PmixL^2+1)^(-1/2)*SDs[i]
      sigWithin[i]=PmixL*(PmixL^2+1)^(-1/2)*SDs[i]
      
      pSDs[,i]=rep(sigWithin[i], nplatoons)
      pMeanLengths[,i]=MeanLengths[i]+sigBetween[i]*seq(-(nplatoons-1)/2, (nplatoons-1)/2, 1)
  }
  
  
  LengDist=array(0, dim=c(nplatoons, A, length(bins)))
  
  for (p in 1:nplatoons) {
    for (i in 1:A){
      LengDist[p,i,]=pnorm(bins+width/2, pMeanLengths[p,i], pSDs[p,i])-pnorm(bins-width/2, pMeanLengths[p,i], pSDs[p,i])
    }
    
    LengDist[p,,]=LengDist[p,,]/rowSums(LengDist[p,,])
    colnames(LengDist[p,,])<-bins
    rownames(LengDist[p,,])<-1:A
  }
  
    
  #######################################################
  ####### Length-weight (from Kim and Hyun 2018) ########
  #######################################################
  
  WL=0.015*bins^2.728*(1e-6)
  
  
  WeightAge1=0.015*MeanLengths[1]^2.728*(1e-6)
  WeightAge2=0.015*MeanLengths[2]^2.728*(1e-6)
  WeightAge3=0.015*MeanLengths[3]^2.728*(1e-6)
  
  aveWeightAge1Age2Age3=(WeightAge1+WeightAge2+WeightAge3)/3
  
  
  #######################################################
  ####### Fork-Total (from Kim and Hyun 2018) ###########
  #######################################################
  
  TL=1:60
  FL=-0.424+0.959*TL
  
  ##################################################
  ########## Maturity (from kooka 2012) ############
  ##################################################
  
  Matsteep=1.2
  Mat50=34.8
  
  mats=1/(1+exp(-Matsteep*(bins-Mat50)))
  
  
  ##################################################
  ########## Fecundity (from kooka 2012) ###########
  ##################################################
  
  fBeta=3.72
  fAlpha=1.60*10^(-1)
  
  
  Fecundity=fAlpha*bins^(fBeta)
  
  
  #################################################
  ########## Selectivity ##########################
  #################################################
  
  Nogari=1-Sel(bins, nogariL05, nogariL95)

   
  Inputs=list("bins"=bins,
              "nbins"=length(bins),
              "JYt"=Catch[,"Juv"],
              "AYt"=Catch[,"Adult"],
              "A"=A,
              "nCatchYears"=dim(Catch)[1],
              "nProj"=nProj,
              "nyears"=dim(Catch)[1]+nProj,
              "LengthWeight"=WL,
              "Maturity"=mats,
              "Fecundity"=Fecundity,
              "FemaleProp"=0.5,
              "LengthDist"=LengDist,
              #"Vul"=Vul,
              "Nogari"=Nogari,
              "Release"=Release,
              "gamma"=gamma,
              "sigR"=sigR,
              "JuvMaxVul"=JuvMaxVul,
              "AduMaxVul"=AduMaxVul,
              "JuvCollapseYear"=JuvCollapseYear,
              "AduCollapseYear"=AduCollapseYear,
              #"L05"=L05,
              "ThresCollapse"=ThresCollapse,
              "VulRec"=rep(1, dim(Catch)[1]+nProj-1),
              "FindBrps"=FindBrps,
              "steepPoint"=steepPoint,
              "RdevCor"=RdevCor,
              "L05devCor"=L05devCor,
              "L95devCor"=L95devCor,
              "sigL05"=sigL05,
              "sigL95"=sigL95,
              "L05_low"=L05_low,
              "L05_up"=L05_up,
              "L95_low"=L95_low,
              "L95_up"=L95_up,
              "LrefForM"=LrefForM,
              "Rprop"=Rprop,
              "NatExponent"=NatExponent,
              "nplatoons"=nplatoons)
  
  
  return(Inputs)
  
}

