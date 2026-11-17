# Registro de cambios — PolyAIA

Lo más reciente arriba. Estados: Implementado / Validado / Revertido.
El historial anterior al 2026-11-17 está en `git log` (rama `test_PolyAIA`).

## 2026-11-17 — SeuratPipeline multi-assay para WNN (CH-20261117-02)
- **Estado:** Implementado; pendiente de validación con datos reales (polyA, Harmony)
- **Razón:** correr `SeuratPipeline()` sobre el assay `polyA` sin sobrescribir las
  reducciones del RNA, como entrada de `FindMultiModalNeighbors()`.
- **Cambio:** argumentos `pca.name`, `pca.key`, `integration.key`, `keep.scale.data`;
  HVFs tomadas de `VariableFeatures(seurat_obj[[assay]])` (antes del DefaultAssay);
  claves derivadas del nombre y chequeo de claves duplicadas; `scale.data` se puede
  descartar tras la integración; `misc$n_dims_by_reduction` y
  `misc$cluster_prefix_by_assay`; `misc$n_dims`/`misc$cluster_prefix` solo se
  escriben para el DefaultAssay o si están vacíos.
- **Archivos:** `R/1_AnnotationSceData.R` (+ `man/SeuratPipeline.Rd` vía `devtools::document()`)
- **Validación:** probado con datos simulados (Seurat 5.0.1, `RPCAIntegration`); VAL-BASE-04

## 2026-11-17 — Línea base de la memoria del proyecto (CH-20261117-01)
- **Estado:** Implementado; pendiente de validación
- **Razón:** adoptar un sistema de memoria en archivos para que cualquier sesión retome el trabajo.
- **Cambio:** se crean `AGENTS.md`, `CLAUDE.md` y `docs/` (PLAN, ESTADO, CHANGELOG,
  VALIDACION) a partir del inventario del repo, las notas del proyecto y la entrevista.
  Se añaden `AGENTS.md`, `CLAUDE.md` y `docs/` a `.Rbuildignore` para que no entren
  en el paquete. Referencia: commit `984841d` (2026-11-13), con 3 archivos modificados sin commit.
- **Archivos:** `AGENTS.md`, `CLAUDE.md`, `docs/*.md`, `.Rbuildignore`
- **Validación:** VAL-BASE-01..12; VAL-CH-20261117-01-01
