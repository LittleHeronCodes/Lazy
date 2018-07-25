## Script courtesy of Dr. HS. Lee
HyperG_test = function(setA,setB,Uni){
	
	if( length(setA) >=3 && length(setB) >=3){
		#Uni=unique(c(Uni,setA,setB))
		setA=intersect(setA, Uni)
		setB=intersect(setB, Uni)
		U = length(Uni)									#항아리 전체! (전체 probe)
		A = length(setA)									# 흰공 : DEG_A
		B = length(setB)									# 뽑힌 공  : DEG_B
		interAB = intersect(setA,setB)						# 뽑힌 공 중에서 흰공 ! : 공통 DEG
		Acomp = setdiff(Uni,setA)
		
		p = phyper(length(interAB)-1,  A, length(Acomp) , B , lower.tail = FALSE, log.p =FALSE)
		return(p) 											# -log 로 변환된 p value값으로 리턴!!! 
	
	}else{
		return(NA)
	}
}


