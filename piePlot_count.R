#!/usr/bin/env Rscript



vectCol = function(vin){

  vout = c()
	#print(vin)
  for (classLSR in vin){
	  	#print(classLSR)
    if (classLSR == "P"){
      vout = append(vout, "#f78e04")
    }
	else if (classLSR == "cycle.P"){
      vout = append(vout, "#D2691E")
    }
    else if (classLSR == "B"){
      vout = append(vout, "#c1fff5")
    }
	else if (classLSR == "cycle.B"){
      vout = append(vout, "#FFFF00")
    }
    else if (classLSR == "F"){
      vout = append(vout, "#faff7f")
    }
    else if (classLSR == "cycle.F"){
      vout = append(vout, "#227700")
    }
    else if (classLSR == "Cl"){
      vout = append(vout, "#00680a")
    }
    else if (classLSR == "cycle.Cl"){
      vout = append(vout, "#77FF00")
    }
    else if (classLSR == "Br"){
      vout = append(vout, "#ef8e0e")
    }
    else if (classLSR == "cycle.Br"){
      vout = append(vout, "#00BBFF")
    }
    else if (classLSR == "Be"){
      vout = append(vout, "#ff19fb")
    }
    else if (classLSR == "cycle.Be"){
      vout = append(vout, "#003377")
    }
    else if (classLSR == "NO2"){
      vout = append(vout, "#6d0305")
    }
    else if (classLSR == "cycle.NO2"){
      vout = append(vout, "#FF1493")
    }
    else if (classLSR == "SO2"){
      vout = append(vout, "#847e02")
    }
	else if (classLSR == "cycle.SO2"){
      vout = append(vout, "#9900FF")
    }
    else if (classLSR == "S"){
      vout = append(vout, "#fff200")
    }
	else if (classLSR == "cycle.S"){
      vout = append(vout, "#00FFFF")
    }
    else if (classLSR == "CON"){
      vout = append(vout, "#00ff37")
    }
	else if (classLSR == "cycle.CON"){
      vout = append(vout, "#CCCCCC")
    }
    else if (classLSR == "COO"){
      vout = append(vout, "#ff0008")
    }
	else if (classLSR == "cycle.COO"){
      vout = append(vout, "#008080")
    }
    else if (classLSR == "onlyC"){
      vout = append(vout, "#383a38")
    }
	else if (classLSR == "cycle.onlyC"){
      vout = append(vout, "#00FF00")
    }
    else if (classLSR == "C.O"){
      vout = append(vout, "#965b4c")
    }
	else if (classLSR == "cycle.C.O"){
      vout = append(vout, "#F0F0F0")
    }
    else if (classLSR == "C.N"){
      vout = append(vout, "#040d89")
    }
	else if (classLSR == "cycle.C.N"){
      vout = append(vout, "#CCEEFF")
    }
    else if (classLSR == "C.O.N"){
      vout = append(vout, "#c489b9")
    }
	else if (classLSR == "cycle.C.O.N"){
      vout = append(vout, "#FF3EFF")
    }
    else if (classLSR == "other"){
      vout = append(vout, "#afafa8")
    }
	else if (classLSR == "cycle.other"){
      vout = append(vout, "#808000")
    }
    #else{
    #  vout = append(vout, "#ffffff")
    #}
  }

  names(vout) = vin
  vunique = unique(vout)
  return(list(vout,vunique))
}

pieType = function (d, path_out){
	
	#print (d[0,])
	leg = NULL
	title = NULL
	data =NULL
	for (l in names (d)){
		#print(l)
		leg = append (leg, paste (l, ", ", d[l], sep = ""))
		title = append (title, paste (l))
		data = append (data, paste (d[l]))
	}
	colors = vectCol(title)[[2]]
	print(colors)
	par (lwd = 1000)
	png(filename=paste(path_out,".png",sep = ""),1000,1000)
	pie(as.double(d), col = colors,label = '',radius=0.55)
	#pie(as.double(d), col = colors, label = leg, lwd = 10, cex = 1.5)
	legend("bottom", ncol=5,legend = leg, cex = 1.45,fill=colors)
	dev.off()

}

# MAIN #
args <- commandArgs(TRUE)
p_filin = args[1]


d = read.table (p_filin, header = TRUE, sep = "\t")
#ldorded = orderanddividedtable (d)
#write.table(ldorded[[1]], file = paste(ptable, "_nocycle", sep = ""), sep = "\t")
pieType  (d, p_filin)

