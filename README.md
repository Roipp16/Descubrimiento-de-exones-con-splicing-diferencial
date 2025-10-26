# Descubrimiento-de-exones-con-splicing-diferencial
En este repositorio se presenta el código usado para analizar unos datos RNASeq, descubriendo los exones diferencialmente expresados y finalmente el splicing diferencial

Tenemos, realmente, dos análisis que paren de diferentes inputs, se detalla a continuación: 

1) Fichero deBAM_a_Exon.R
   Input: Ficheros BAM dentro de una carpeta.
   Output: Tablas comparativas para cada grupo XXXXX
   Código escrito en lenguaje R.

2) Script_Enrichment.R
   Input: Tablas resultado del paso anterior
   Output: Gráficos de enriquecimiento de los resultados

1) Fichero Script_differentialPSI
   Input: Ficheros FASTQ dentro de una carpeta con los nombes adientes a cada una de las muestras.
   Output:
             - Tabla de expresión diferencial estadística para cada una de las comparativas.
             - Tabla de filtrado de resultados para las comparativas indicadas.
   Código escrito en lenguaje bash.

2) Lenght plot and filtering
  Input: Tablas resultado extraídas del paso anterior y la tabla con todos los eventos.
  Output: Gráfico representativo de longitud de exones y eventos que superan los filtros mandados.
 
   
