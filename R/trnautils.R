newSJD1 = function(inMatrix,i1=1:ncol(inMatrix),i2=1:ncol(inMatrix),verbose=TRUE){

  docnames = colnames(as.matrix(inMatrix))
  tf = as.matrix(inMatrix)
  
  numDocs1 = length(i1)
  numDocs2 = length(i2)
  tfn =  t(tf)/colSums(tf)
  D.SJ = matrix(0, numDocs1, numDocs2)
  rownames(D.SJ) = colnames(as.matrix(inMatrix))[i1]
  colnames(D.SJ) = colnames(as.matrix(inMatrix))[i2]
  entropies = apply(tfn,1,entropy)
  for(i in 1:numDocs1) {
     for(j in 1:numDocs2) { 
        D.SJ[i,j] = entropy(0.5*(tfn[i1[i],]+tfn[i2[j],]))-0.5*(entropies[i1[i]]+entropies[i2[j]])
  		#D.JensenShannon.new(tfn[i1[i],], tfn[i2[j],])
     }
     if(verbose){print(i)}
  }
  return(D.SJ)
}

entropy = function(p){
	p = p[p!=0]
	-sum(p*log(p))
}

rc = function(x){
	sapply(strsplit(x,""),function(y){
		temp = y
		temp[y=="A"] = "T"
		temp[y=="T"] = "A"
		temp[y=="G"] = "C"
		temp[y=="C"] = "G"
		paste(rev(temp),collapse="")
	})
}

grassgreens9 = c("#FBFFF7","#EBF7DE","#DBEFC6","#CAE19E",
	"#AED66B","#92C642","#71B521","#519C08","#306B08")
purples9 = c("#FBF7FF","#EBDEF7","#DBC6EF","#CA9EE1","#AE6BD6",
	"#9242C6","#7121B5","#51089C","#30086B")
golds9 = c("#FFFBF7","#F7EBDE","#EFDBC6","#E1CA9E","#D6AE6B",
	"#C69242","#B57121","#9C5108","#6B3008")
reds9 = c("#FFF7FB","#F7DEEB","#EFC6DB","#E19ECA","#D66BAE",
	"#C64292","#B52171","#9C0851","#6B0830")
seagreens9 = c("#F7FFFB","#DEF7EB","#C6EFDB","#9EE1CA","#6BD6AE",
	"#42C692","#21B571","#089C51","#086B30")

maxDiffColors = c("#E7298A", "#008080", "#800000", "#35D6EA", "#000080", "#FFA500", "#800080",
		"#38DA10", "#F0E44F", "#008000", "#0428CD", "#E41A1C", "#808000", "#F781BF",
		"#A3A3A3", "#E3E3E3")

color.groups = function(z, col = tim.colors(50),
zlim = c(min(z, na.rm = TRUE),
max(z, na.rm = TRUE)),breaks = NULL){
    if(is.null(breaks)){
	    breaks <- seq(zlim[1], zlim[2], length = length(col) + 1)
    }
    zgroups <- cut(z, breaks = breaks, labels = FALSE,
    include.lowest = TRUE)
    return(col[zgroups])
}

padrange = function(x,amount=0.06){range(x,na.rm=TRUE)+diff(range(x,na.rm=TRUE))*amount*c(-1,1)}


scatterBox = function(x,y,nbreaks=10,myjitter=1,ptcols=blues9[3],ylim=range(y),pch=1,
		      horizontal=FALSE,...){
	temp = 1:length(unique(x[!is.na(x)]))
	names(temp) = as.character(sort(unique(x[!is.na(x)])))
	
	simplifyLife = tapply(y,x,function(x){return(x)})
	jitterAmount = sapply(simplifyLife,function(z){
		zlim = c(min(z, na.rm = TRUE),max(z, na.rm = TRUE))
		breaks <- seq(zlim[1], zlim[2], length = nbreaks + 1)
		if(zlim[1]==zlim[2]){
		zgroups = rep(1,length(z))
		}else{
		zgroups <- cut(z, breaks = breaks, labels = FALSE,
			include.lowest = TRUE)
		}
		groupKey = table(zgroups)
		groupKey[as.character(zgroups)]/max(groupKey)
	
	})
	if(length(ptcols)==1){
		pointcols = ptcols
	}else{
		pointcols = tapply(ptcols,x,function(x){return(x)})
	}
	
	if(!horizontal){
	plot(jitter(temp[as.character(rep(names(simplifyLife),sapply(simplifyLife,length)))],myjitter*unlist(jitterAmount)),unlist(simplifyLife),
		xaxt='n',xlim=c(1,length(temp))+c(-0.5,0.5),col=unlist(pointcols),xlab="",ylab="",ylim=ylim,pch=pch,bty="n")
	}else{
	plot(unlist(simplifyLife),jitter(temp[as.character(rep(names(simplifyLife),sapply(simplifyLife,length)))],myjitter*unlist(jitterAmount)),
		yaxt='n',ylim=c(1,length(temp))+c(-0.5,0.5),col=unlist(pointcols),xlab="",ylab="",xlim=ylim,pch=pch,bty="n")
	}
	boxplot(y~x,at=temp,add=TRUE,boxwex=0.5,horizontal=horizontal,...)
}

gppred = function(x, y, l = diff(range(x))/1.5, tau = 0.1){
	ymean = mean(y)
	
	calcKern = exp(-as.matrix(dist(x))^2/(2*l^2))
	temp = chol(calcKern + diag(rep(tau^2,length(x))))
	Kinv = chol2inv(temp)
	findkstar = function(newx){exp(-outer(x,newx,"-")^2/(2*l^2))}
	nx = seq(min(x)-1,max(x)+1,length.out=200)
	kstar = findkstar(nx)
	
	calcMean = t(kstar) %*% Kinv %*% (y-ymean) + ymean
	calcSig = sqrt(1 + tau^2 - diag(t(kstar) %*% Kinv %*% kstar))
	
	list(x = nx, mu = calcMean, sig = calcSig)
}

