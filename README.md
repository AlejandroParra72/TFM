- **pipeline_NSD1.R**
  Pipeline completo para identificar la epifirma de pacientes con **síndrome de Sotos (NSD1)**. Incluye pasos de preprocesamiento, filtrado, detección de DMPs, selección de CpGs, PCA y entrenamiento/evaluación de un SVM.

- **pipeline_CDKN1C.R**  
  Pipeline completo para identificar la epifirma de pacientes con **síndrome de Beckwith-Wiedemann (CDKN1C)**. Incluye pasos de preprocesamiento, filtrado, detección de DMPs, selección de CpGs, PCA y entrenamiento/evaluación de un SVM..

- **pipeline_NSD1_epifirma_previa.R**  
  Script  que permite cargar y aplicar una lista de CpGs ya publicada en la literatura para NSD1. Lee un fichero de texto con la epifirma externa, filtra la matriz de β-values y ejecuta PCA y SVM para evaluar el poder discriminativo de esa firma en tu conjunto de muestras.
