# Revisa que GitHub muestre bien las fórmulas de los archivos Markdown. Uso,
# desde la raíz (requiere `gh` con sesión iniciada):
#   julia src/revisar_formulas.jl                 # todos los .md del curso
#   julia src/revisar_formulas.jl archivo.md ...  # sólo esos archivos
# Pide a `gh api markdown` el HTML de cada archivo y comprueba que:
#   1. cada bloque `$$` de la fuente salga como fórmula;
#   2. cada fórmula en línea `$...$` salga idéntica a la fuente;
#   3. no queden `$` sueltos en el texto (fórmulas que no se convirtieron);
#   4. no haya Markdown (cursivas, negritas) dentro de una fórmula ni bloques
#      `<pre lang="math">` sin convertir.
# Errores frecuentes: una fórmula en línea partida entre dos renglones, una
# que termina en `)` seguida de otro `)` (`$f(x)$)`), `}_` (abre cursivas),
# `\%` y bloques `$$` de varias líneas dentro de un elemento de lista con
# fórmulas en línea.
# Termina con código 1 si encuentra algún problema.

const DIRECTORIOS = ["programa", "docs", "Class_notes", joinpath("Class_notes", "repaso")]

problemas = String[]
reportar(msg) = push!(problemas, msg)

"Decodifica las entidades HTML que usa GitHub (`&amp;gt;` aparece doblemente codificado)."
function decodificar(s)
    for _ in 1:2
        s = replace(s, "&lt;" => "<", "&gt;" => ">", "&quot;" => "\"", "&#39;" => "'",
                    "&#x27;" => "'", "&amp;" => "&")
    end
    s
end

normalizar(s) = replace(decodificar(s), r"\s+" => "")
recortar(s, n = 80) = (s = replace(s, r"\s+" => " "); length(s) > n ? first(s, n) * "…" : s)

"Quita los bloques de código y el código en línea de la fuente."
sin_codigo(src) = replace(replace(src, r"```.*?```"s => ""), r"`[^`\n]*`" => "")

archivos = isempty(ARGS) ?
    ["README.md"; [joinpath(d, f) for d in DIRECTORIOS for f in readdir(d) if endswith(f, ".md")]] :
    ARGS

total_formulas = 0
for ruta in archivos
    global total_formulas
    src = read(ruta, String)
    html = read(`gh api markdown -f mode=gfm -F text=@$ruta`, String)
    formulas = [m[1] for m in eachmatch(r"<math-renderer[^>]*>(.*?)</math-renderer>"s, html)]
    total_formulas += length(formulas)
    salidas = Set(normalizar(f) for f in formulas)

    # --- 1. Bloques $$ ----------------------------------------------------------
    fuente = sin_codigo(src)
    for m in eachmatch(r"\$\$.*?\$\$"s, fuente)
        normalizar(m.match) in salidas || reportar("$ruta: bloque sin convertir: $(recortar(m.match))")
    end

    # --- 2. Fórmulas en línea -------------------------------------------------
    for (i, linea) in enumerate(split(replace(fuente, r"\$\$.*?\$\$"s => ""), '\n'))
        for m in eachmatch(r"(?<![\\$])\$[^$\n]+\$", linea)
            normalizar(m.match) in salidas || reportar("$ruta: fórmula distinta de la fuente: $(recortar(m.match))")
        end
    end

    # --- 3. Signos $ sueltos en el texto -----------------------------------------
    texto = replace(html, r"<math-renderer.*?</math-renderer>"s => "")
    texto = replace(texto, r"<(pre|code)[^>]*>.*?</\1>"s => "")
    for linea in split(replace(texto, r"<[^>]+>" => ""), '\n')
        occursin('$', linea) && reportar("$ruta: \$ suelto en el texto: $(recortar(decodificar(strip(linea))))")
    end

    # --- 4. Markdown dentro de fórmulas y bloques sin convertir -----------------
    for f in formulas
        occursin(r"<(em|strong)>", f) && reportar("$ruta: Markdown dentro de una fórmula: $(recortar(decodificar(f)))")
    end
    n = count("<pre lang=\"math\"", html)
    n > 0 && reportar("$ruta: $n bloque(s) <pre lang=\"math\"> sin convertir")
end

if isempty(problemas)
    println("Sin problemas ($(length(archivos)) archivos, $total_formulas fórmulas).")
else
    foreach(println, problemas)
    println("\n$(length(problemas)) problema(s).")
    exit(1)
end
