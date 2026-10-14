# Genera man/figures/flujo.svg: el flujo de trabajo de ieiprio con los seis
# graficos de la galeria incrustados como miniaturas.
# Requiere Python 3 con Pillow. Se corre desde la raiz del paquete, despues de
# data-raw/figuras_readme.R:
#   python3 data-raw/flujo.py

import base64
import io
from pathlib import Path

from PIL import Image

FIG = Path("man/figures")

GRAFICOS = [
    ("galeria-1-cromatograma.png", "iei_cromatograma()", "Picos de cada base en la variante"),
    ("galeria-2-dotplot.png", "iei_dotplot()", "Lectura frente a la referencia"),
    ("galeria-6-clinvar.png", "iei_gen()", "Variantes de ClinVar en la proteína"),
    ("galeria-3-circos.png", "iei_circos()", "Candidatos sobre el genoma"),
    ("galeria-4-acordes.png", "plot(b, tipo = \"acordes\")", "Signos del caso y genes"),
    ("galeria-5-cariograma.png", "iei_cariograma()", "Candidatos en los cromosomas"),
]


if not (FIG / "galeria-6-clinvar.png").exists() and (FIG / "galeria-6-ranking.png").exists():
    GRAFICOS[2] = ("galeria-6-ranking.png", "plot(b)", "Genes candidatos por fenotipo")


def miniatura(nombre):
    im = Image.open(FIG / nombre).convert("RGB").resize((520, 390), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, "JPEG", quality=82, optimize=True)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def caja(x, y, w, h, titulo, sub, fondo, borde, tinta, sub_tinta="#344054", rx=8):
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fondo}"'
        + (f' stroke="{borde}"' if borde else "")
        + "/>"
        f'<text x="{x + w / 2}" y="{y + h / 2 - 4}" text-anchor="middle" font-weight="600" fill="{tinta}">{titulo}</text>'
        f'<text x="{x + w / 2}" y="{y + h / 2 + 13}" text-anchor="middle" font-size="11" fill="{sub_tinta}">{sub}</text>'
    )


partes = []
W, H = 960, 900
partes.append(
    f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
    f'viewBox="0 0 {W} {H}" font-family="Inter, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif" '
    f'role="img" aria-label="Flujo de trabajo de ieiprio, de la secuencia al reporte, con los gráficos de cada paso">'
)
partes.append(
    '<defs><marker id="f" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
    'orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="#475467"/></marker>'
)
for i in range(len(GRAFICOS)):
    partes.append(f'<clipPath id="c{i}"><rect width="288" height="216" rx="10"/></clipPath>')
partes.append("</defs>")
partes.append(f'<rect width="{W}" height="{H}" rx="14" fill="#ffffff"/>')

# Encabezados de columna
etiqueta = '<text x="{x}" y="34" font-size="11" font-weight="700" fill="#667085" letter-spacing="1.2">{t}</text>'
partes.append('<g font-size="12">')
partes.append(etiqueta.format(x=24, t="ENTRADA"))
partes.append(etiqueta.format(x=248, t="PROCESO"))
partes.append(etiqueta.format(x=794, t="SALIDA"))

# Entradas
for y, t, s in [(48, "VCF", "iei_leer_vcf()"), (116, "Sanger .ab1 .scf FASTA", "iei_sanger()"),
                (184, "PDF, Word o texto", "iei_variantes()")]:
    partes.append(caja(24, y, 170, 54, t, s, "#eef4fb", "#1d6fa5", "#0b2545"))

# Proceso
for x, t, s in [(248, "Filtrar", "panel, calidad, sexo"), (430, "Anotar", "VEP, ClinVar, gnomAD"),
                (612, "Priorizar", "puntaje y ACMG")]:
    partes.append(caja(x, 104, 148, 78, t, s, "#0b2545", None, "#ffffff", "#d0d5dd"))

# Salidas
for y, t, s in [(70, "Reporte HTML", "iei_reporte()"), (138, "Excel, VCF, JSON", "iei_exportar()")]:
    partes.append(caja(794, y, 142, 54, t, s, "#ecfdf3", "#12b76a", "#054f31"))

# Apoyos
for x, t, s in [(248, "Sexo y CNV", "iei_verificar_sexo()"), (430, "Frecuencia", "iei_gnomad()"),
                (612, "Fenotipo", "HPO, iei_buscar_fenotipo()")]:
    partes.append(caja(x, 250, 148, 50, t, s, "#f9fafb", "#d0d5dd", "#344054", "#475467"))
partes.append("</g>")

# Flechas
partes.append('<g stroke="#475467" stroke-width="1.4" fill="none" marker-end="url(#f)">')
for d in ["M194 75 C220 75 220 143 246 143", "M194 143 L246 143", "M194 211 C220 211 220 143 246 143",
          "M396 143 L428 143", "M578 143 L610 143", "M760 143 L792 97", "M760 143 L792 165"]:
    partes.append(f'<path d="{d}"/>')
partes.append("</g>")
partes.append('<g stroke="#98a2b3" stroke-width="1.2" stroke-dasharray="4 3" fill="none" marker-end="url(#f)">')
for x in (322, 504, 686):
    partes.append(f'<path d="M{x} 250 L{x} 186"/>')
partes.append("</g>")

# Graficos
partes.append('<line x1="24" y1="334" x2="936" y2="334" stroke="#e6edf6" stroke-width="2"/>')
partes.append('<text x="24" y="362" font-size="11" font-weight="700" fill="#667085" letter-spacing="1.2">GRÁFICOS</text>')
partes.append('<text x="104" y="362" font-size="11.5" fill="#667085">Cada paso tiene su vista. Arriba, la secuencia y el gen; abajo, el fenotipo sobre el genoma.</text>')

xs = [24, 336, 648]
ys = [378, 640]
for i, (archivo, funcion, texto) in enumerate(GRAFICOS):
    x = xs[i % 3]
    y = ys[i // 3]
    datos = miniatura(archivo)
    partes.append(
        f'<g transform="translate({x},{y})">'
        f'<rect x="-1" y="-1" width="290" height="218" rx="11" fill="#ffffff" stroke="#dce6f2"/>'
        f'<image width="288" height="216" clip-path="url(#c{i})" preserveAspectRatio="xMidYMid meet" '
        f'href="data:image/jpeg;base64,{datos}" xlink:href="data:image/jpeg;base64,{datos}"/>'
        f'<text x="2" y="236" font-size="12.5" font-weight="700" fill="#0b2545" '
        f'font-family="JetBrains Mono, SFMono-Regular, Menlo, Consolas, monospace">{funcion.replace(chr(34), "&quot;")}</text>'
        f'<text x="2" y="252" font-size="11.5" fill="#52514e">{texto}</text>'
        "</g>"
    )

partes.append("</svg>")
(FIG / "flujo.svg").write_text("".join(partes), encoding="utf-8")
print("flujo.svg:", round((FIG / "flujo.svg").stat().st_size / 1024), "KB")
