## Tutorial para construir una red de coexpresión génica con la herramienta WGCNA (Weighted-Correlation Network Analysis)
PhD STU.Cynthia Soto Cardinualt; cyntsc10@gmail.com; linkedin.com/in/cynthiacardinault; Centro de Investigación Científica de Yucátan

**Package:** WGCNA <br>
**Versión:** 1.70-3 <br>
**Date:** 2021-02-2017 <br>

**Nombre de la herramienta:** Weighted Correlation Network Analysis <br>
**Autor:** Peter Langfelder <Peter.Langfelder@gmail.com> and Steve Horvath <SHorvath@mednet.ucla.edu> with contributions by Chaochao Cai, Jun Dong, Jeremy Miller, Lin Song, Andy Yip, and Bin Zhang <br>
**Maintainer:** Peter Langfelder <Peter.Langfelder@gmail.com> <br>

*keywords: transcriptómica; RNA-Seq; redes_ponderadas; WGCNA* <br>

<hr/>

#### Resumen
La rápida aceptación y beneficios que han traído las tecnologías de secuenciación de próxima generación, particularmente RNASeq, ha acelerado el despliegue de métodos de análisis disponibles para ómicas favoreciendo el descubrimiento de muchos aspectos biológicos relevantes. <br> 
En el análisis de datos de coexpresión génica cuando se desea explorar o descubrir patrones génicos de respuesta coordinada ante determinada condición de estudio, se puede recurrir a la utilización de métodos de Aprendizaje Automático (Machine Learning/ML) no supervisados, asumiendosé que se tiene datos de alta dimensión no etiquetados, y suficientes muestras. De entre los diversos enfoques de ML, las técnicas de Agrupamiento (Clustering) han demostrado ser útiles en el análisis de este tipo de datos, y de entre los muchos algoritmos existentes, el Agrupamiento Jerárquico (hCLust) ha sido uno de los más efectivos. <br> 
En este mini-curso, utilizando transcriptomas RNA-Seq de plantas infectadas por hongos, se muestra como construir una red de coexpresión para identificar patrones génicos que responden a la  infección causada por hongos. Se utiliza la herramienta Weighted Correlation Network Analysis (WGCNA) (Langfelder & Horvath, 2008) disponible para R a través de Bioconductor.<br>

#### Objetivos de aprendizaje:
<ol>
  <li>Conocer los conceptos generales de redes y del clustering jerárquico</li>
  <li>Conocer las principales funciones de la herramienta WGCNA</li>
  <li>Construir una red con signo con la herramienta WGCNA
    <ol>
      <li>Preparar el entorno en R</li>
      <li>Descargar los datos a utilizar en el análisis</li>
      <li>Calcular e interpretar las propiedades topológicas de la red</li>
      <li>Definir un umbrales de corte para la construcción de una red con signo (met=p)</li>
        <ol>
          <li>Contruir una red estandar basada en dissimilitudes</li>
          <li>Contruir una red fusionada a partir de la red estandar</li>
        </ol>
      <li>Vincular la red creada con meta-datos para identificar patrones de interés</li>
      <li>Aprender a identificar módulos altamente correlacionados, y sus propiedades</li>
      <li>Aprender a validar la importancia de un módulo a partir de sus métricas de significancia (Student asymptotic p-value) y membresía al módulo</li>
      <li>Aprender a guardar los objetos de datos para su post-análisis</li>
      <li>Aprender a descargar los módulos para su post-análisis</li>
    </ol>
  </li>
</ol> 

<hr/>

#### Prerequisitos:
La versión actual del paquete WGCNA solo funcionará con la versión R 3.0.0 y superior. Si tiene una versión anterior de R, actualice su R. <br><br>
El paquete WGCNA requiere la instalación de los siguientes paquetes: stats, grDevices, utils, matrixStats (0.8.1 or higher), Hmisc, splines, foreach, doParallel, fastcluster, dynamicTreeCut, survival, parallel, preprocessCore, GO.db, impute, and AnnotationDbi. Si su sistema no los tiene instalados, puedes seguir las instrucciones dadas en **[CoexpressionNetwork/Rpackages/WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/)** o también pueden adaptar las primeras líneas de código proporcionadas en el **script 1_***. <br>

#### Depencias: R (+3.0), dynamicTreeCut (+1.62), fastcluster <br>
Bibliotecas a importar: stats, grDevices, utils, matrixStats (>= 0.8.1), Hmisc, impute, splines, foreach, doParallel, preprocessCore, survival <br><br>

Si desea primero checar si ya tiene instalados algunos de estos paquetes puede ejecutar la siguiente línea en el prompt de R: <br>

<div style="background-color: rgb(50, 50, 50);">

``python
a<-installed.packages(); 
packages<-a[,1];
is.element("WGCNA", packages);  
``

</div>    
    
    
<hr/>

#### Material incluido:
<ul>
  <li>Metadatos</li>
    <ul>
      <li>data/athal_GO_terms.txt </li>
      <li>data/geneGO_metadata.csv </li>
      <li>data/Traits_Infected_hpi.csv       # datos de vinculación externa por hora de infección </li>
      <li>data/Traits_Infected_sample.csv    # datos de vinculación externa por muestra </li>
    </ul>
  </li>
</ul>                        

<ul>
  <li>Datasets</li>
    <ul>
      <li>data/matrix_E_infected_random25p.csv       # subset aleatorizado al 25% </li>
      <li>data/matrix_E_infected_random50p.csv       # subset aleatorizado al 50% </li>
      <li>data/matrix_E_infected_simulated5000.csv   # subset aleatorizado con una muestra de 5000 genes en los valores extremos altos y bajos. </li>
    </ul>
  </li>
</ul>            

<ul>
  <li>Código</li>
    <ul>
      <li>1_Infected_SignedNtw_E.R <br>
      Step 1/4: análisis de topología libre de escala para la identificación del umbral de corte dinámico, cálculo de la matriz de superposición topológica (TOM), cálculo de la matriz de adyacencias y de la matriz de disimilitudes, e identificación de módulos genéticos para red estándar y fusionada con signo (+) por el método de pearson.
      </li>
      <li>02_Infected_GS_MM_24hpi.R <br>
      Step 2/4: cuantificación de asociaciones módulo-tratamiento(rasgo) y visualización de correlaciones con sus valores-p, análisis intra-módulos utilizando análisis de significancia genética (GS) y membresía al módulo (MM). 
      </li>
      <li>03_Infected_Plot_Ntw_Eigengenes.R <br>
      Step 3/4: análisis de la red eigengene por tratamiento/rasgo, visualización topológica de los eigengenes de la red y visualización de la red (dendrogram con heatmap vinculada al tratamiento).
      </li>
      <li>04_Infected_Enrichment_analysis.R <br>
      Step 4/4: análisis de enriquecimiento de los eigengenes con un archivo de anotación personalizado, exportación de datos de enriquecimiento en formato VisANT y Cytoscape.
      </li>
    </ul>
  </li>
</ul>            

#### Otros objetos de datos pre-procesados:
Athal_Infected_*.RData           # objetos ya procesados para agilizar la demostración. <br> 
**NOTAS:**
Si tiene la capacidad de computar en paralelo, en linux puede consultar el número de núcleos con el comando $nproc o $lscpu -e; en windows puede consultarlo en la ficha rendimiento del administrador de tareas. <br> <br> 

Última modificación: 19 de abril 2022
