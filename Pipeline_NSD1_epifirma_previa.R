###### LIBRERÍAS ######

# Instalamos Bioconductor y el resto de paquetes necesarios para el desarrollo del pipeline:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Paquetes CRAN:
install.packages(c("MatchIt", "httr", "httr2", "openssl", "knitr", "stringr", "ggplot2",
                   "lmerTest", "caTools", "e1071", "readxl", "openxlsx", "pROC",
                   "dplyr", "ggfortify", "qqman"))

# Paquetes Bioconductor: 
BiocManager::install(c("rtracklayer", "AnnotationHub", "bumphunter", "limma", "minfi", "minfiData",
                       "minfiDataEPIC", "Gviz", "DMRcate", "FlowSorted.Blood.EPIC",
                       "GenomicRanges", "epimutacions", "methylclock", "methylclockData",
                       "wateRmelon", "ENmix", "missMethyl",
                       "IlluminaHumanMethylationEPICmanifest",
                       "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
                       "IlluminaHumanMethylationEPICv2manifest"))

# Activamos las librerías:
library(MatchIt) # para hacer matcheo de edad y sexo cuando hacemos el train 
library(rtracklayer) # ?
library(AnnotationHub) # estimacion de exo y edad
library(httr)
library(httr2)
library(openssl)
library(bumphunter) # para DMRs 
library(knitr)
library(limma) # estadistica
library(minfi) # estadistica
# library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest) # v1 (no usamos)
library(RColorBrewer) #colores
library(missMethyl)# no se usa (para cargar datos de metilacion)
library(minfiDataEPIC) # estadistica
library(Gviz) # visualizacion de datos
library(DMRcate) # para calculo de dmr
library(stringr)
library(ggplot2)
library(FlowSorted.Blood.EPIC) # Calculo de tipo celular
library(lmerTest) # estadistica
library(caTools) # ?
library(e1071) # machine leraning
library(readxl)
library(openxlsx)
library(GenomicRanges)
library(epimutacions) # calculo DMRs
library(methylclock) # edad
library(methylclockData) # edad
library(wateRmelon) # sexo
library(ENmix) # eliminacion de sufijos y DMRs
library(pROC) # estadistica
library(dplyr)
library(ggfortify) # PCA
library(qqman)# grafico
#library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
#library(IlluminaHumanMethylationEPICv2manifest)


###### ANOTACIÓN ###### 

# Obtenemos la anotación de los arrays Illumina EPIC v2 para la referencia del genoma hg38,
# la cual contiene información sobre las sondas CpG, sus posiciones en el genoma, genes asociados, etc. 
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(annEPIC) # Mostramos las primeras filas de la anotación para inspeccionarla.

# Extraemos los nombres de las sondas CpG de la anotación y las almacenamos en  el objeto "CpGs_EPICv2".
CpGs_EPICv2 <- rownames(annEPIC)
write.table(CpGs_EPICv2, file = "CpGs_EPICv2.txt", sep = "\t") # Guardamos los nombres de las sondas en un archivo


###### DIRECTORIO DE DATOS ###### 

dataDirectory <- "D:/ALEJANDRO" # Definimos la ruta donde están los datos de metilación.
a <- read.metharray.sheet(dataDirectory, pattern="SampleSheetOGS.csv") # Creamos un dataframe llamado "a"
# La función read.metharray.sheet() busca y carga el archivo SampleSheetOGS.csv, en este caso, dentro 
# del directorio especificado anteriormente. 

#Comprobamos de que los archivos de metilación (.idat) existen:
file.exists(paste0(a$Basename, "_Red.idat"))


#### LOAD INPUT ####

# La función read.metharray.exp(targets=a)  carga los datos de metilación desde los archivos .idat.
# Dentro del objeto A está la información de las muestras, incluyendo la columna Basename con las rutas de los archivos .idat.
# Lo que hacemos con el siguiente comando es cargar los datos de metilación de los archivos idat de las muestras
# contenidas en el objeto "a" y los volcamos en un nuevo objeto "rgSet", el cual es de clase RGChannelSet.

rgSet <- read.metharray.exp(targets=a) 
rgSet # Vemos que se han importado 307 muestras en total
which(!a$Sample_Name %in% rgSet$Sample_Name) # Comprobamos si hay muestras en el objeto "a" que no 
# estén en el objet rgset. La devolcuión ha de ser integer(0) para garantizar que se han copiado todas las muestras.

# Asignamos información de anotación al objeto rgSet:
#annotation(rgSet)["array"] = "IlluminaHumanMethylationEPICv2"
#annotation(rgSet)["annotation"] = "20a1.hg38"


############# FILTRADO DE CALIDAD #############

# Calculamos los p-valores de detección para cada sonda y muestra, los cuales indican la confianza en la detección de cada sonda.
detP <- detectionP(rgSet) # Generamos una matriz con los p valores de [Número de sondas] x [Número de muestras]

keep <- colMeans(detP) < 0.05 # colMeans(detP) calcula el promedio de los p-valores de todas las sondas para cada muestra.
# Filtramos las muestras con demasiadas sondas mal detectadas (si el promedio de p-valores es mayor o igual a 0.05, se eliminan).
# "keep" es un vector lógico con TRUE para las muestras a conservar y FALSE para las muestras a eliminar.

table(keep) # Comprobamos cuantas muestras presentan un p value inferior a 0.05 (TRUE): 
# Podemos ver que hay dos muestras que no presentan un p value adecuado para su posterior análisis

# Identificamos que muestras son las que presentan mala calidad (a través de Basename)
muestras_mala_calidad <- colnames(rgSet)[!keep]
print(muestras_mala_calidad)

rgSet <- rgSet[, keep, drop=FALSE] # Eliminamos las muestras con mala calidad de rgSet.
rgSet # Comprobamos que ahora hay 305 muestras en total (2 menos que antes)

# Filtramos las mismas muestras en detP, asegurando que detP y rgSet contengan exactamente las mismas muestras.
detP <- detP[,keep]
dim(detP) # Vemos que tambien obtenemos 305 muestras en el matriz de p valores. 


############# NORMALIZACIÓN #############

# Comprobamos la densidad de los datos antes de la normalización: 
densityPlot(rgSet, sampGroups=a$Slide,
            main="Prenormalizacion", legend=FALSE)

# Realizamos la normalización: 
mSetSq <- preprocessFunnorm(rgSet)

# Volvemos a representar la densidad, una vez se han normalizado los datos: 
densityPlot(getBeta(mSetSq), sampGroups=a$Slide,
            main="Postnormalizacion", legend=FALSE)

# Aseguramos que detP y mSetSq tengan el mismo número y orden de CpGs antes de seguir con los filtros.
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

############# PREPROCESADO DE LOS DATOS  ###############

# A continuación, eliminamos aquellas sondas que tengan un p value superior a 0.01. 
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep) # Podemos ver que hay 102800 sondas con p values superiores a 0.01
mSetSqFlt <- mSetSq[keep,] # Nos quedamos solo con las sondas con p value inferior a 0.01

# Quitamos las sondas ubicadas en los chrX y chrY
keep <- !(featureNames(mSetSqFlt) %in% annEPIC$Name[annEPIC$chr %in% 
                                                      c("chrX","chrY")])
table(keep) # Podemos ver que 19708 caen en los cromosomas sexuales
mSetSqFlt <- mSetSqFlt[keep,] # eliminamos dichas sondas. 

# Quitamos las sondas que coincidan con SNPs conocidos
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

# Quitamos xReactiveProbes
xReactiveProbes <- read.csv(file=paste("D:/ALEJANDRO",
                                       "13059_2016_1066_MOESM1_ESM_2.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep) # Podemos identificar 2828 sondas xReactiveProbes
mSetSqFlt <- mSetSqFlt[keep,]    #### OJOOOOOOO!!!! no tiene que salir TRUE todas

########### CÁLCULO DE BETA VALUES Y M VALUES ################

mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)
bVals <- rm.cgsuffix(bVals)

##### USAR EPIFIRMA PREVIA #############

# Leemos la epifirma, para ello tenemos que tener decargado el conjunto de sondas de dicha epifirma:
archivo_epifirma <- "D:/ALEJANDRO/epifirma_Sotos.txt"
cpgs_epifirma <- read_lines(archivo_epifirma)

# Filtramos bVals por los CpGs de la epifirma:
bVals_epi <- bVals[rownames(bVals) %in% cpgs_epifirma, , drop = FALSE]

# Nos aseguramos de que las muestras estén bien nombradas:
colnames(bVals_epi) <- a$Sample_Name

# Creamos la matriz de entrada para el PCA y SVM:
data_new <- t(bVals_epi)   

###### PCA ########
pca <- prcomp(data_new, center = TRUE, scale. = TRUE)

autoplot(pca, data = a, colour = "Affected") +
  labs(title = "PCA: Epifirma bibliografía")

####### SVM #######

# Seleccionar solo muestras control sano y muestras con variantes patogénicas en NSD1 para entrenar
idx_train   <- which(a$Affected %in% c("Control", "NSD1"))
data_train  <- data_new[idx_train, ]             
labels_train <- factor(a$Affected[idx_train], levels = c("Control", "NSD1"))

svm_df <- data.frame(Affected = labels_train, data_train)

# Entrenamos SVM lineal con probabilidades:
set.seed(42)
svm_model <- svm(Affected ~ ., data = svm_df,
                 kernel = "linear", probability = TRUE)

# Usamos el mismo modelo (svm_model) para predecir en el conjunto completo
pred_full  <- predict(svm_model, newdata = data_new, probability = TRUE)
probs_full <- attr(pred_full, "probabilities")[, "NSD1"]  # probabilidad de NSD1

# Construimos el dataframe para realizar la representación gráfica de SVM:
svmplot_full <- data.frame(
  Sample    = rownames(data_new),
  Affected  = a$Affected,
  P_NSD1    = probs_full

# Creamos la matriz de confusión 
conf_mat <- table(Predicho = pred_full[idx_train],
                  Real     = labels_train)

# Gráfico SVM
ggplot(svmplot_full, aes(x = Affected, y = P_NSD1, color = Affected)) +
  geom_jitter(width = 0.2, size = 2) +
  labs(x = "", y = "Probabilidad de NSD1") +
  theme_bw() +
  ylim(0, 1) +
  ggtitle("Probabilidades de pertenecer a NSD1")
