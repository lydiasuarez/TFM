#Cargamos las librerías

library(devtools)
library(hdf5r)
library(loomR)


# Definimos variables

mypath= "../tissues/"
tissues<-c("amygdala","cerebellum","cortex1","enteric","cortex2","cortex3","drg","medulla","hippocampus","hypothalamus","midbrainventral","midbraindorsal","olfactory","pons","spinalcord","striatumdorsal","striatumventral","sympathetic","thalamus")


#Conectamos los archivos con el environment

connectLoomFiles  <-  function(directorio){
  listloom <- dir(path = directorio,pattern = "*.loom")
  for (i in listloom){
    lfile <- gsub(".loom","",i)
    assign(lfile, loomR::connect(filename = paste0(directorio,i), mode = "r+"), .GlobalEnv)
  }
}

connectLoomFiles(mypath)



# Script

getActiveGenesFromTissue <- function(datasets, tissues, methods,visibility){
  
  # Creamos el dataframe. Lo rellenamos todos con FALSE.
  
  lista <- list("BrainRegion","CellType")
  data <- data.frame(matrix(ncol = length(lista), nrow= 2000)) # Creamos el dataframe
  colnames(data) <- lista
  data[,] <- FALSE
  if (datasets == "Zheisel"){
    
  # Iteramos sobre nuestros datos de transcriptómica. Por cada región cerebral, y, dentro de esta, sobre cada tipo celular.
    
    for (dataset in datasets){
      contador = 0 # Contador que define la fila en la que vamos a introducir los datos
      for (brainregion in tissues){ # Iteramos sobre cada tipo de región cerebral
        new_brainregion <- get(brainregion) # Seleccionamos el objeto en el environment con dicho nombre
        all_celltypes <- data.frame(table(new_brainregion$col.attrs$Class[]))
        each_celltypes <- all_celltypes[all_celltypes$Freq > 100, 1] # Filtro para escoger solo los tipos de celula que tengan mas de 100 muestras
        mymatrix = new_brainregion$matrix[,]
        colnames(mymatrix) = new_brainregion$row.attrs$Gene[]
        rownames(mymatrix) = new_brainregion$col.attrs$CellID[]
        
        for (x in each_celltypes){ # Iteramos sobre cada tipo de celula
          expr.data = mymatrix[new_brainregion$col.attrs$Class[] == x,]
          genes = colnames(expr.data)
          for (method in methods){
            
            
  # Método Variationcoefficient(CV)            
            
            
            if (method == "variationcoefficient"){
              contador = contador + 1 # Empieza una nueva fila
              data$BrainRegion[contador] <- brainregion #AÃ±adimos los datos a las filas
              data$CellType[contador] <- x
              cv = function (x, na.rm = FALSE)
              {
                sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm)
              }
              cat(dataset,"-",brainregion,"-",x,"-",method)
              cvs = apply(expr.data,2,cv,na.rm=T)
              genes_cvs = names(na.omit(cvs))
              cat(length(genes_cvs),"\n")
              for( gen in genes_cvs){
                data[contador,gen] <- TRUE
              }
            }
            
  # Método Rawexpression
            
            
            if (method == "rawexpression"){ 
              contador = contador + 1
              data$BrainRegion[contador] <- brainregion
              data$CellType[contador] <- x
              cat(dataset,"-",brainregion,"-",x,"-",method,"-",visibility,"\n")
              mask = colSums(expr.data > 0.5)  > (visibility * nrow(expr.data))
              genes_activos = colnames(expr.data)[mask]
              for(gen in genes_activos){
                data[contador,gen] <- TRUE
              }
            }
          }
        }
      }
    }
  }

  # Guardamos la matriz con los datos  
  
  data <- data[rowSums(is.na(data)) != ncol(data),]
  assign(paste0("ActiveGenes_",methods,visibility),data)
  save(list = paste0("ActiveGenes_",methods,visibility),file=paste0("ActiveGenes_",methods,visibility))
}

# Lanzamos el script
getActiveGenesFromTissue("Zheisel",tissues,"variationcoefficient",NA)
getActiveGenesFromTissue("Zheisel",tissues,"rawexpression",5)
