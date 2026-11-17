# Estado actual — PolyAIA
_Última actualización: 2026-11-17._

## Dónde estamos
- Fase activa: Fase 1 — evolución continua en `test_PolyAIA`.
- Rama `test_PolyAIA` al día con origin; último commit `984841d` (2026-11-13, Rd FilterPeaks/PolyAPipeline).
- Según la persona usuaria, todo funciona como se espera en el uso con datos reales.
- Línea base de validación: ver `docs/VALIDACION.md` (pendiente de recorrer).

## Cambios sin commit (al adoptar la memoria)
- `R/1_AnnotationSceData.R`: `SeuratPipeline` multi-assay (`pca.name`, `pca.key`,
  `integration.key`, `keep.scale.data`; `misc$n_dims_by_reduction`,
  `misc$cluster_prefix_by_assay`). Por validar.
- `examples/PolyAIA_example.Rmd` (+ `.nb.html`): experimentos locales; usa
  `mock_SeuratPipeline`/`mock_ClusterTree` y `library(PolyAIA)` comentado.
  Revisar antes de commitearlo (el `.nb.html` es el escaparate en GitHub).

## En evolución
- `DEPsMatrixBulk()` marcado WIP en el historial.

## Decisiones técnicas vigentes
| Área | Decisión | Motivo |
| --- | --- | --- |
| Flujo | Editar el paquete directamente; GitHub es el respaldo | Los mocks quedaron obsoletos (oct 2026) |
| Git | Claude propone mensajes; la persona usuaria hace commit/push | Control del historial |
| Documentación | roxygen2 -> `man/`, `NAMESPACE` | Generados; no editar a mano |
| Azimuth | `AzimuthAnnotation()` separada, llamada por `UploadSce` | Reutilizable y diagnosticable |
| RED | `pair_rule = "max"` por defecto; tests LRT/QL/Fisher/t-test | Igual que MAAPER |
| DEGs | edgeR directo en pseudobulk con covariables | Libra no admite covariables |
| AD vs NS | Dos llamadas DEPsMatrix por sexo, fondo `NSfemale`/`NSmale` | Fondo emparejado por sexo |
| Organoides | Comparaciones siempre vs WT (AD, TBEV, HSV-1, AD+TBEV, AD+HSV-1) | Mismo set que MAAPER |

## Próximo paso acordado
- Recorrer la línea base de VALIDACION.md y decidir qué hacer con los cambios sin commit.
