# Matriz de validación — PolyAIA

Estados: Pendiente, Validado técnicamente, Validado, Regresión.
Una fila solo pasa a Validado con confirmación de la persona usuaria.
La validación se hace con datos reales en la máquina de la persona usuaria.

## CH-20261117-01 — Memoria del proyecto
| ID | Caso | Pasos | Resultado esperado | Estado | Evidencia |
| --- | --- | --- | --- | --- | --- |
| VAL-CH-20261117-01-01 | La memoria no entra en el paquete | `devtools::build()` o `R CMD build .` | El tar.gz no contiene `AGENTS.md`, `CLAUDE.md` ni `docs/`; sin NOTE nueva | Pendiente | |

## Línea base — funciones críticas
| ID | Caso | Pasos | Resultado esperado | Estado | Evidencia |
| --- | --- | --- | --- | --- | --- |
| VAL-BASE-01 | Documentación e instalación | `devtools::document(); devtools::install()` | Sin errores; `NAMESPACE` con 33 exports | Pendiente | |
| VAL-BASE-02 | Carga y anotación | `UploadSce()` (con `run_azimuth = TRUE`) | Objeto con genes anotados, QC y CellType/Subpopulation | Pendiente | |
| VAL-BASE-03 | QC por tipo celular | `QCPlot()`, `QCFilter()` | Umbrales por tipo celular y gráficos correctos | Pendiente | |
| VAL-BASE-04 | Integración | `SeuratPipeline()` en RNA (y en polyA con nombres propios) | Reducciones `pca`/`harmony`/`umap_harmony`; la 2.ª corrida no sobrescribe la 1.ª | Pendiente | |
| VAL-BASE-05 | Visualización celular | `UmapPlot()`, `DotPlot()`, `ClusterTree()` | Gráficos con el estilo actual | Pendiente | |
| VAL-BASE-06 | Assay polyA | `UploadPolyAAssays()` → `FilterPeaks()` → `RequantifyPolyA()` → `Peaksdb()` / `PolyAPipeline()` | Assay `polyA` anotado; filtro por condición correcto | Pendiente | |
| VAL-BASE-07 | APA diferencial | `DEPsMatrix()` con `Pasta = TRUE`, `pair_rule = "max"` | Tabla RED con padj sin NA espurios; RED_LRT/RED_QL disponibles | Pendiente | |
| VAL-BASE-08 | Filtro → QCDensity | `QCDensity()` sobre el filtro de `DEPsMatrix` | Funciona sin cambios adicionales | Pendiente | |
| VAL-BASE-09 | DEGs | `DEGsMatrix()` con covariables | Resultados edgeR por tipo celular | Pendiente | |
| VAL-BASE-10 | Coverage | `PolyAPlot()` con overlay y `show_stats` | Título/subtítulo/leyenda según el estilo acordado | Pendiente | |
| VAL-BASE-11 | Bulk (WIP) | `DEPsMatrixBulk()` con conteos de polyApipe_bulk | RED/DPAU calculados | Pendiente | |
| VAL-BASE-12 | Ejemplo | Knit de `examples/PolyAIA_example.Rmd` | Se ejecuta de punta a punta | Pendiente | |
