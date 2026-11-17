# Plan maestro — PolyAIA

## Propósito y visión
Una librería R sólida que sea un pipeline integral para el análisis de
poliadenilación en single-cell: de los datos 10x y los picos polyA crudos hasta
el uso diferencial de sitios (RED/DPAU) y figuras listas para publicar.
Crece según las necesidades que surgen al analizar datos reales (cerebro AD/NS,
ratón 3xTg/B6 con alcohol, organoides infectados, bulk 3' RNA-seq).

## Principios no negociables
1. Cambios pequeños y reversibles en las funciones del paquete; GitHub
   (`test_PolyAIA`) es el respaldo. Los mocks (`../mock_DE_QCDensity/`) ya no son
   el paso previo.
2. Ningún cambio rompe a las funciones que consumen su salida.
3. Documentación roxygen al día con cada función exportada (con ejemplos).
4. La persona usuaria decide qué se valida, qué se commitea y qué se revierte.

## Historia resumida
- 2026-07-29 — v0.1.0 inicial: pipeline de tres etapas (single cell, assay polyA, APA diferencial).
- 2026-09 — `modified_residuals` en DEPsMatrix; `DEGsMatrix`; filtro de picos de
  DEPsMatrix alineado con DEGsMatrix y usable por QCDensity.
- 2026-10 (inicio) — QCPlot/DotPlot; gran rediseño de `PolyAPlot` (overlay, sitios p/d, RED/FDR en el subtítulo).
- 2026-10-13 a 11-02 — se prueba APAdeg (APAdegsMatrix/CompareDEGs) y se elimina.
- 2026-10-19 — `AzimuthAnnotation()` y `AddQCMetrics()` se separan de `UploadSce`; `UploadData` exportada.
- 2026-10-27 — covariables en `DEGsMatrix` (edgeR directo).
- 2026-11-09 — RED_QL e interacción en `DEPsMatrix`; `DEPsMatrixBulk` (WIP) para bulk polyApipe.
- 2026-11-11/13 — ejemplos en todos los Rd; `FilterPeaks` filtra por condición y divide los gráficos.
- 2026-11-17 — se adopta esta memoria del proyecto.
Detalle de las funciones: doc del proyecto "Changes in the app" y `git log`.

## Fases
### Fase 0 — Versión inicial 0.1.0
**Estado:** completada (2026-07-29)

### Fase 1 — Evolución continua en `test_PolyAIA`
**Estado:** en curso
**Objetivo:** mejorar y ampliar funciones según el análisis de datos reales.
**Trabajo:** el que vaya surgiendo (ver ESTADO.md).
**Criterio de salida:** (por definir; p. ej. fusionar a `main`, publicación).

### Fases siguientes
(por definir)

## Protocolo de trabajo
1. Revisar estado → 2. Cambio pequeño y reversible (backup/ si se revierte) →
3. Revisión propia del código y de las funciones consumidoras →
4. Validación de la persona usuaria con datos reales → 5. Registrar →
6. Proponer mensajes de commit (los ejecuta la persona usuaria).
