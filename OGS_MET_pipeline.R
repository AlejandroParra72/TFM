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
# library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
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
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)


###### ANOTACIÓN ###### 

# Obtenemos la anotación de los arrays Illumina EPIC v2 para la referencia del genoma hg38,
# la cual contiene información sobre las sondas CpG, sus posiciones en el genoma, genes asociados, etc. 
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
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
annotation(rgSet)["array"] = "IlluminaHumanMethylationEPICv2"
annotation(rgSet)["annotation"] = "20a1.hg38"


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
                                       "xReactiveProbes_V2.csv",
                                       sep="/"), sep = ";", stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep) # Podemos identificar 2828 sondas xReactiveProbes
mSetSqFlt <- mSetSqFlt[keep,]

########### CÁLCULO DE BETA VALUES Y M VALUES ################

mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

# Encontramos las muestras que aparecen tanto en rgSet como en el objeto a, ya que hay dos muestras
# que no han pasado el control de calidad (filtrado) 
common_samples <- intersect(rgSet$Sample_Name, a$Sample_Name)

# Conservamos solo las muestras comunes entre ambos objetos
a_filtered <- a[a$Sample_Name %in% common_samples, ]
a <- a_filtered
rownames(a) <- NULL  # Reiniciamos los índices

# Cambiamos los nombres de las columnas de mVals y bVals para que coincidan con los Sample_Name de "a"
colnames(mVals) <- a$Sample_Name
colnames(bVals) <- a$Sample_Name

# Eliminamos los sufijos que incluye EPICv2 en el nombre de las sondas
bv <- rm.cgsuffix(bVals)
mv <- rm.cgsuffix(mVals)


######### SEX ESTIMATION ############

bVals_flt <- getBeta(mSetSq)
bv_flt <- rm.cgsuffix(bVals_flt)
colnames(bv_flt) <- a$Sample_Name


estSex <- estimateSex(bv_flt, do_plot=T)
a$estSex <- estSex$predicted_sex

######### AGE ESTIMATION ############

estage <- DNAmAge(bv_flt, cell.count=TRUE, localHub=F) # se debe realizar con los bvals sin filtrar
a$estAge <- estage$Horvath


### Clasificacion de muestras  ###


a$Affected <- "Control"
a$Affected[c(290:305)] <- "Control_validation"

a$Affected[c(2,3,6:9,13,15)] <- "CDKN1C"
a$Affected[c(14)] <- "CDKN1C_validation"
a$Affected[c(1,4,5,10:12)] <- "CDKN1C_VUS"

a$Affected[c(16:21)] <- "CHD3"

a$Affected[c(25,27)] <- "FBN1"
a$Affected[c(22:24,26,28,29)] <- "FBN1"

a$Affected[c(30:36)] <- "GPC3"

a$Affected[c(37:42)] <- "NFIX"

a$Affected[c(45:55)] <- "NSD1"
a$Affected[c(57,59,60)] <- "NSD1_Validation"
a$Affected[c(44,43,56,58,61)] <- "NSD1_VUS"

a$Affected[c(62:67)] <- "PTEN"

a$Affected[c(69:73)] <- "RNF125"
a$Affected[c(68,74)] <- "RNF125_VUS"

a$Affected[c(75:79)] <- "SETD2"


################ PCA ###############

# Definimos la cohorte de entrenamiento (controles vs condición)
traininit <- a[c(45:55, 80:289),] 
rownames(traininit) <- NULL # Reiniciamos los índicies 

# matchit() realiza el emparejamiento, la selección de subconjuntos y la subclasificación con el objetivo
# de crear grupos de tratamiento y control equilibrados en las covariables incluidas (en este caso, Sex + estAge)

m.out0 <- matchit(as.factor(Affected) ~ Sex + estAge , data = traininit, 
                  method = "optimal", distance = "glm", ratio = 3) #aqui tenemos que cambiar el Ratio segun vayamos viendo los resultados en el summary

summary(m.out0) # Podemos observar que se han emparejado 8 controles con nuestras 8 muestras afectas

plot(m.out0, type = "jitter", interactive = FALSE)

train <- match.data(m.out0) # Extraemos las muestras que han sido emparejadas en el dataframe "train"
rownames(train) <- NULL # Reiniciamos los índices

# A continuiación, tenemos que repetir todo el proceso pero con nuestra cohorte de entrenamiento: 
rgSet <- read.metharray.exp(targets=train)
rgSet 

annotation(rgSet)["array"] = "IlluminaHumanMethylationEPICv2"
annotation(rgSet)["annotation"] = "20a1.hg38"

detP <- detectionP(rgSet)

keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
train <- train[keep,]
rownames(train) <- NULL

detP <- detP[,keep]
dim(detP)

mSetSq <- preprocessFunnorm(rgSet)

detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)

mSetSqFlt <- mSetSq[keep,]

keep <- !(featureNames(mSetSqFlt) %in% annEPIC$Name[annEPIC$chr %in% 
                                                      c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]

mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

xReactiveProbes <- read.csv(file=paste("D:/ALEJANDRO",
                                       "xReactiveProbes_V2.csv",
                                       sep="/"), sep = ";",stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)

mSetSqFlt <- mSetSqFlt[keep,]

mvtrain <- getM(mSetSqFlt) 
# mvtrain <- rm.cgsuffix(mvtrain)
colnames(mvtrain) <- train$Sample_Name

bvtrain <- getBeta(mSetSqFlt)
bvtrain <- rm.cgsuffix(bvtrain)
colnames(bvtrain) <- train$Sample_Name

dim(mvtrain)

############# CELL TYPE COMPOSITION ###############

# Estimamos la composición celular de nuestras muestras a partir de los datos de metilación

IDOLOptimizedCpGsBloodv2<- IDOLOptimizedCpGs[which(IDOLOptimizedCpGs%in%rownames(bvtrain))]
identical(rownames(IDOLOptimizedCpGs.compTable[IDOLOptimizedCpGsBloodv2,]), IDOLOptimizedCpGsBloodv2)
propEPIC <- projectCellType_CP(
  bvtrain[IDOLOptimizedCpGsBloodv2, ],
  IDOLOptimizedCpGs.compTable[IDOLOptimizedCpGsBloodv2,],
  contrastWBC = NULL, nonnegative = TRUE,
  lessThanOne = FALSE)

cells <- as.data.frame(propEPIC)

############# MODELO ESTADÍSTICO ############

# Creamos un modelo estadístico para poder comparar las dos condiciones ( afectados vs controles)

# Comenzamos definiendo las variables:
diagnose <- factor(train$Affected)
age <- as.numeric(train$estAge)
sex <- as.factor(train$Sex)
slide <- as.factor(train$Slide)
cd4t <- cells$CD4T
cd8t <- cells$CD8T
mono <- cells$Mono
NK <- cells$NK
bcell <- cells$Bcell
neu <- cells$Neu

# COnstruimos el modelo de diseño para el análisis de metilación diferencial.
design <- model.matrix(~0+diagnose+age+sex+cd4t+cd8t+neu+NK, data=train)

# Ajustamos un modelo lineal para cada CpG
fit <- lmFit(mvtrain, design)

# Definimos los grupos (contrastes)
contMatrix <- makeContrasts(diagnoseNSD1-diagnoseControl, levels=design)
contMatrix

# Aplicamos la comparación de interés ( en este caso CDKN1C vs controles)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2) # ajustamos la varianza para mejorar la detección de diferencias significativas

# Mostramos el resultado del contraste, obteniendo que CpGs están diferencialmente metiladas entre afectados y controles
summary(decideTests(fit2))



#########################################


# Extraemos anotaciones de los CpGs presentes en mvtrain desde annEPIC
annEPICSub <- annEPIC[match(rownames(mvtrain),annEPIC$Name),
                      c(1:4,12:19,24:ncol(annEPIC))]

# Extraemos los CpGs diferencialmente metilados (DMPs) a partir del modelo ajustado (fit2) y usa la anotación correspondiente.
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=annEPICSub) 
DLGtop <- DMPs

bvtrain <- getBeta(mSetSqFlt)
colnames(bvtrain) <- train$Sample_Name

# Calculamos el promedio de los valores beta en los grupos de control (beta_norm)
# y afectados (beta_dis) y obtenemos su diferencia (Delta_beta):
beta_norm <- rowMeans(bvtrain[,design[,1]==1])
beta_dis <- rowMeans(bvtrain[,design[,2]==1])
Delta_beta <- beta_dis - beta_norm

# Combinamos las métricas estadísticas de los DMPs con los valores de delta beta en una sola tabla.
ewasdb <- merge(DLGtop, Delta_beta, by = "row.names")

## THRESHOLD ##

# Definimos como significativos los CpGs que tienen un valor ajustado de p < 0.05 y un cambio absoluto mayor de 0.1.
ewasdb$sigDM <- ewasdb$adj.P.Val < 0.05 & abs(ewasdb$y) > 0.1
table(ewasdb$sigDM)

# Filtramos  la tabla para quedarte solo con esos CpGs
ewasdbfilt <- filter(ewasdb, sigDM == TRUE)
rownames(ewasdbfilt) <- ewasdbfilt$Row.names
ewasdbfilt <- ewasdbfilt[,-1]


## INTERACTION ##


ewasdb$logP <- (log(ewasdb$adj.P.Val))*-1 # Transformamos el adj.P.Val en una escala logarítmica negativa (-log(p))
ewasdb$score <- ewasdb$logP*abs(ewasdb$y) # Calculamos un score combinado que multiplica la significancia (logP) por la magnitud del cambio (|delta beta|).
ewasdb<-ewasdb %>% arrange(desc(score)) # Ordenamos  la tabla desde el CpG con mayor score hasta el menor.
ewasdbfilt <- ewasdb[c(1:2000),] # Seleccionamos los 2000 CpGs con mayor score, es decir, los más relevantes según ambos criterios (efecto + significancia).
rownames(ewasdbfilt) <- ewasdbfilt$Row.names
ewasdbfilt <- ewasdbfilt[,-1]

# Unimos la tabla de los top 2000 CpGs (ewasdbfilt) con sus correspondientes valores beta normalizados (bvtrain). La fusión la hacemos por el nombre de sonda. 
Ftmv <- merge(ewasdbfilt,bvtrain, by= "row.names")
rownames(Ftmv) <- Ftmv$Row.names

### ROC CURVES ###

Ftmv <- Ftmv[c(43:ncol(Ftmv))] # Eliminamos las primeras columnas que corresponden a la anotación y estadísticos de los CpGs. Dejamos solo las columnas que contienen valores beta de las muestras.
colnames(Ftmv) <- train$Sample_Name

diag <- train$Affected # Creamos diag, que  contiene la clase de cada muestra en el conjunto de entrenamiento. Esto será la variable dependiente en las curvas ROC
res<-{} # # Se inicializa un objeto vacío donde se almacenarán los resultados de la curvas ROC

# A continuación itermaos fila por fila en Ftmv. Para cada CpG, se calcula el area bajo la curva (AUC) usando la función auc() de pROC
for (i in 1:nrow(Ftmv)) {
  res1<-summary(auc(diag, as.numeric(Ftmv[i,])) )   
  res<-rbind(res,res1)}

row.names(res)<-row.names(Ftmv) # Asociamos a cada fila su CpG correspondiente.
res

# Guardamos el objeto res en un nuevo objeto:
stxROC <- res
head(stxROC)

# Aseguramos que stxROC esté en formato dataframe
stxROC <- as.data.frame(stxROC)

# Nos quedamos solo con las dos primeras columnas de stxROC, que se corresponden con el mínimo y el 1º cuartil
stxROC <- stxROC[,1:2]
head(stxROC)

# Renombramos la columna Min, como AUC
names(stxROC)[names(stxROC) == "Min."] <- "AUC"
head(stxROC)

# Ordenamos los CpGs de mayor a menor AUC para priorizar aquellos que mejor discriminan entre clases.
stxROC<-stxROC %>% arrange(desc(AUC))
head(stxROC)

# Seleccionamos las CpGs con AUC mayor a 0.95 (alta discriminación entre clases)
stx2Froc <- subset(stxROC, AUC > 0.95) 

# Unimos el dataframe de CpGs seleccionados (stx2Froc) con sus valores beta (Ftmv), usando los nombres de fila (CpGs) como clave.
# Eliminamos la columna Row.names y cualquier metadato adicional, dejando solo la matriz de valores beta.
roc_MV <- merge(stx2Froc, Ftmv, by = "row.names")
rownames(roc_MV) <- roc_MV[,1]
roc_MV <- roc_MV[c(4:ncol(roc_MV))]

mbft <- t(roc_MV) # Transponemos la matriz
# mbft será la matriz que usaremos de ahora en adelante.

### CORRELATION ###

# Calculamos la matriz de correlación de Pearson entre todas las CpGs seleccionadas
cor_matrix <- cor(mbft)
cor_matrix

# Eliminamos los duplicados, poniendo 0 en la parte superior y diagonal de la matriz (así evitamos eliminar ambos CpGs correlacionados (solo se elimina uno de cada par))
cor_matrix_rm <- cor_matrix                 
cor_matrix_rm[upper.tri(cor_matrix_rm)] <- 0
diag(cor_matrix_rm) <- 0
cor_matrix_rm

# Eliminamos las CpGs redundantes (r>0.9):
data_new <- mbft[ , !apply(cor_matrix_rm,
                           2,
                           function(x) any(x > 0.9))]
head(data_new)                               
dim(data_new) # Identificamos cuantas CpGs nos quedan para entrenar el modelo

### PCA PLOT OF SIGNATURE ###

# Ejecutamos el PCA con las muestras de la cohorte de entrenamiento, así podremos 
# observar como se separan controles vs afectos en el PCA
mv.pca1 <- prcomp(data_new, center = TRUE, scale. = TRUE)

# Realizamos la representación gráfica del PCA: 
mv.pca.plot1 <- autoplot(mv.pca1,
                         data = train,
                         colour = 'Affected')
mv.pca.plot1

#Ahora realizamos el PCA, pero añadiendo el resto de muestras de nuestra hoja de trabajo: 

# Primero, nos aseguramos de que as columnas de bVals estén correctamente nombradas con el campo Sample_Name del metadata (a).
colnames(bVals) <- a$Sample_Name

# Transponemos data_new para  obtener CpGs en filas, lo necesario para hacer un merge con bVals (que tiene CpGs como filas).
dnt <- as.data.frame(t(data_new))
dn_mv <- merge(dnt[,1:2],bVals,by="row.names")
rownames(dn_mv) <- dn_mv[,1]
dn_mv[,1:3] <- NULL

dn_mvt <- t(dn_mv) # Transponemos de nuevo, ahora tenemos muestras en filas y CpGs en columnas.
rownames(dn_mvt) <- a$Sample_Name

# Ahora, ya podemos realizar el PCA con todas las muestras:
mv.pca1 <- prcomp(dn_mvt, center = TRUE, scale. = TRUE)
mv.pca.plot1 <- autoplot(mv.pca1,
                         data = a,
                         colour = 'Affected')
mv.pca.plot1


### TRAINING MLC ###

mbft <- t(Ftmv) # Transponemos Ftmv, la matriz de CpGs seleccionados y filtrados, para que las muestras estén en filas y las variables (CpGs) en columnas.

# Nos quedamos con las columnas 1 y 13 columnas de train para el modelo (muestra y condición):
id <- train[,c(1,13)] 
rownames(id) <- NULL
rownames(data_new) <- NULL

# Fusionamos los metadatos (id) con la matriz de CpGs seleccionados y no redundantes (data_new), usando el nombre de muestra como clave.
trainR_id <- merge(id, data_new, by = "row.names")
rownames(trainR_id) <- trainR_id$Sample_Name
trainR_id$Row.names <- NULL
trainR_id$Sample_Name <- NULL

# Aseguramos que la variable Affected sea un factor categórico, como requiere svm() para clasificación.
trainR_id$Affected <- as.factor(as.character(trainR_id$Affected))

#######################
# El siguiente bloque solo se utiliza cuando haya CpGs diferentes entre ambso datasets

# setdiff(colnames(trainR_id), colnames(testR_id)) #Muy util para eliminar las CpGs que no están en ambos datasets
trainR_id1 <- trainR_id[, !colnames(trainR_id) %in% c("cg26227310_TC21")]
# 
trainR_id <- trainR_id1
#######################

# Realizamos una búsqueda (grid search) de los hiperparámetros del modelo SVM:
     # gamma: controla la influencia de cada punto de entrenamiento.
     # cost: penaliza los errores de clasificación.

tuned_parameters <- tune.svm(Affected~., data = trainR_id, gamma = 10^(-5:-1), cost = 10^(-3:1), type="C", kernel="linear")
summary(tuned_parameters )

# Extraemos gamma y cost: 
gamma <- tuned_parameters$best.parameters[,1]
cost <- tuned_parameters$best.parameters[,2]

# Entrenamos el modelo SVM con los parámetros optimizados sobre el conjunto de entrenamiento (trainR_id): 
# Activamos la opción probability = TRUE para obtener las probabilidades de pertenencia a cada clase, no solo la predicción directa.
model <- svm(Affected ~ ., kernel="linear", cost =cost , gamma =gamma ,data=trainR_id, type="C", probability = TRUE )


### TESTING MLC###

# Transponemos los valores beta para que las muestras estén en filas, como espera el modelo.
mValstestT <- t(bVals)
rownames(mValstestT) <- a$Sample_Name

# Extraemos la información necesaria (Sample_Name y Affected) de todas las muestras, no solo las del conjunto de entrenamiento.
id2 <- a[,c(1,13)]
rownames(id2) <- id2$Sample_Name

# Fusionamos las etiquetas con los datos de CpGs para construir el conjunto de test.
testR_id <- merge(id2, mValstestT,  by = "row.names")

rownames(testR_id) <- testR_id$Sample_Name
testR_id$Row.names <- NULL
testR_id$Sample_Name <- NULL
testR_id$Affected <- as.factor(as.character(testR_id$Affected))


# Aplicamos el modelo SVM al conjunto de test, excluyendo la columna Affected. De este modo, 
# obtenemos tanto las predicciones como las probabilidades asociadas.
y_pred <- predict(model, testR_id[-1], probability=TRUE)
y_pred

# Generamos una matriz de confusión que muestra el rendimiento del modelo: clases predichas vs clases reales.
table(y_pred, testR_id$Affected)

### SVM PLOT ###

# Extraemos las probabilidades de predicción para cada clase desde el objeto y_pred
yp <- attr(y_pred, "probabilities")
yp1 <- as.data.frame(yp)

# Unimos las probabilidades con la información de las muestras (Sample_Name, Affected), creando el dataframe svmplot que usaremos para generar el gráfico del SVM
svmplot <- merge(id2, yp1,  by.x = "Sample_Name" , by.y = "row.names")

# Realizamos el gráfico del SVM:
p <- ggplot(svmplot, aes(Affected, NSD1)) ## Cambiar el segundo argumento de aes según lo que quieras visualizar

p + geom_jitter(width = 0.2, aes(colour = Affected)) +labs(x="", y="Probability", 
                                                           col="")+theme(legend.title=element_blank()) +  
  theme_bw() + ylim(0,1)

