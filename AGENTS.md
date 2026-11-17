# PolyAIA: memoria operativa para sesiones futuras

Antes de modificar nada, lee `docs/ESTADO.md` y las últimas entradas de
`docs/CHANGELOG.md`. Para cambios de alcance o dirección, lee `docs/PLAN.md`.

## Qué es
Paquete R para análisis de poliadenilación alternativa (APA/IPA) en single-cell
(Seurat/Signac/PASTA): carga y anotación, QC, integración, ensamblaje del assay
polyA y análisis diferencial (RED, DPAU), más visualización. La meta es una
librería sólida que sea un pipeline integral de análisis single-cell polyA; evoluciona
según las necesidades del análisis de datos reales.

## Límites del proyecto
- Rama de trabajo: `test_PolyAIA` (siempre en evolución; es prueba y error).
  La memoria (`AGENTS.md`, `CLAUDE.md`, `docs/`) se versiona solo en esta rama.
- **Git: nunca hacer commit ni push.** Al cerrar un bloque, proponer los
  mensajes de commit (en inglés, estilo del historial: `WIP: <Función> ...`,
  `docs: ...`). La persona usuaria los ejecuta en su terminal y pide la
  presentación de commit/push cuando la necesita.
- **Nunca borrar archivos** sin indicación explícita. Solo crear o modificar.
- `../Mice_example/` (datos reales 3xTg/B6) es intocable: no se modifica, mueve ni borra.
- No editar un archivo mientras la persona usuaria lo está ejecutando; si solo
  pide una respuesta, contestar en texto sin editar.
- `man/*.Rd` y `NAMESPACE` los genera roxygen2 (`devtools::document()`): se
  edita el bloque `#'` en `R/` y se regenera; no editar los `.Rd` a mano.
- La validación es con datos reales en la máquina de la persona usuaria; en el
  entorno de Claude no hay R instalado, así que solo se revisa el código.
- Entorno: conda (no mamba).

## No romper
- Antes de proponer commit, revisar que el cambio funciona con las funciones que
  lo consumen (p. ej. filtro de `DEPsMatrix` -> `QCDensity`): un arreglo no debe
  necesitar otro commit.
- Revertir guardando: si se deshace una función añadida, el archivo modificado
  se copia antes a `backup/` (`<archivo>_<descripcion>_<YYYYMMDD_HHMM>.R`; `backup/`
  está en .gitignore) por si se quiere recuperar.
- Azimuth vive en `AzimuthAnnotation()`, llamada por `UploadSce(run_azimuth = TRUE)`;
  sube `future.globals.maxSize` (sin eso falla en silencio). Se intentó dejarla
  inline y se revirtió.
- `UploadSce` no combina muestras de referencias 10x distintas (hg19 + GRCh38): revertido.
- `PolyAPlot`: no reintroducir `coverage_label_hjust` (revertido). Mantener el
  estilo acordado: título "<célula> - <gen>" centrado; subtítulo
  "Intronic APA | REDi = x, FDR = y" (FDR = PASTA_padj primero, luego padj del
  test RED; sin "p -> d"); overlay con leyenda vertical entre eje y coverage y
  el control al fondo; sitios "p"/"d" (p1/d1… si hay más de dos). La señal de
  extremo 3' se descartó por ruidosa.
- No reintroducir APAdeg (`APAdegsMatrix`, `CompareDEGs`, `CompareDEGsPlot`): eliminado a propósito.
- `DEPsMatrix`: `pair_rule = "max"` por defecto (igual que MAAPER);
  `exclude_peaks = FALSE`; los warnings de glm no deben volver NA los p-valores del LRT.

## Registro obligatorio
Cada cambio se registra en `docs/CHANGELOG.md` con ID `CH-YYYYMMDD-NN` y sus
casos en `docs/VALIDACION.md`. `docs/ESTADO.md` se actualiza al terminar.
Nada se declara validado sin confirmación de la persona usuaria (que valida con
datos reales y avisa si algo falla: entonces se corrige o se revierte).
