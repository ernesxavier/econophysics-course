# Revisa la consistencia del material del curso. Uso, desde la raíz:
#   julia src/revisar_material.jl
# Comprueba que:
#   1. cada clave del índice del repaso exista y su etiqueta incluya esa sesión;
#   2. cada clave citada en una sesión (o en un taller) aparezca en las filas
#      del índice de esa sesión (o de las sesiones del módulo);
#   3. los enlaces relativos de las notas lleven a archivos y anclas existentes;
#   4. no queden menciones a sesiones inexistentes ni a "semanas" del curso.
# Termina con código 1 si encuentra algún problema.

const REPASO = joinpath("Class_notes", "repaso")
const AREAS = Dict('P' => "probabilidad.md", 'E' => "procesos_estocasticos.md",
                   'L' => "algebra_lineal.md", 'C' => "calculo.md",
                   'F' => "fisica_estadistica.md", 'S' => "estadistica.md")
const NSESIONES = 20

problemas = String[]
reportar(msg) = push!(problemas, msg)

"Sesiones de una etiqueta como \"Sesiones 3, 8, proyecto final\" (0 = proyecto final)."
function sesiones_de(etiqueta)
    startswith(etiqueta, "Todas") && return Set(0:NSESIONES)
    if occursin("en adelante", etiqueta)
        n = parse(Int, match(r"\d+", etiqueta).match)
        return Set([0; n:NSESIONES])
    end
    m = match(r"(\d+)–(\d+)", etiqueta)
    m === nothing || return Set(parse(Int, m[1]):parse(Int, m[2]))
    s = Set(parse(Int, x.match) for x in eachmatch(r"\d+", etiqueta))
    occursin("proyecto", etiqueta) && push!(s, 0)
    s
end

"Expande \"[E.B2–B3](...)\" en las claves E.B2 y E.B3."
function claves_en(texto)
    ks = String[]
    for m in eachmatch(r"\b([PELCFS])\.([BIA])(\d+)(?:–[BIA]?(\d+))?", texto)
        a = parse(Int, m[3])
        b = m[4] === nothing ? a : parse(Int, m[4])
        append!(ks, ["$(m[1]).$(m[2])$n" for n in a:b])
    end
    ks
end

# --- Etiquetas de cada concepto -------------------------------------------
etiquetas = Dict{String,Set{Int}}()
for archivo in values(AREAS)
    texto = read(joinpath(REPASO, archivo), String)
    for m in eachmatch(r"### ([PELCFS]\.[BIA]\d+)[^\n]*\n\*([^*\n]+)\*", texto)
        etiqueta = rstrip(m[2], '.')
        if !(startswith(etiqueta, "Sesi") || startswith(etiqueta, "Todas"))
            reportar("$(m[1]): etiqueta con formato desconocido: $(m[2])")
            continue
        end
        etiquetas[m[1]] = sesiones_de(etiqueta)
    end
end

# --- 1. Índice del repaso --------------------------------------------------
indice = Dict{Int,Set{String}}()
for linea in eachline(joinpath(REPASO, "README.md"))
    m = match(r"^\| (\d+|—) \|[^|]*\|(.*)\|$", linea)
    m === nothing && continue
    s = m[1] == "—" ? 0 : parse(Int, m[1])
    indice[s] = Set(claves_en(m[2]))
    for k in indice[s]
        if !haskey(etiquetas, k)
            reportar("Índice, sesión $s: la clave $k no existe")
        elseif !(s in etiquetas[k])
            reportar("Índice, sesión $s: la etiqueta de $k no incluye esa sesión")
        end
    end
end
sort(collect(keys(indice))) == 0:NSESIONES ||
    reportar("El índice del repaso no tiene una fila por sesión (más el proyecto final)")

# --- 2. Claves citadas en cada sesión ---------------------------------------
"Sesiones que cubre un archivo, según su título (\"# Sesión 3:\" o \"# Sesiones 1 y 2:\")."
function sesiones_del_archivo(ruta)
    titulo = first(filter(l -> startswith(l, "# Sesi"), readlines(ruta)))
    Set(parse(Int, x.match) for x in eachmatch(r"\d+", split(titulo, ':')[1]))
end

# Cada taller puede citar las claves de las sesiones de su módulo
const MODULOS = Dict(1 => Set(1:7), 2 => Set(8:14), 3 => Set(15:20))

for ruta in [readdir("src"; join = true); readdir("Class_notes"; join = true)]
    taller = match(r"/taller(\d)_[a-z_]+\.jl$", ruta)
    (taller !== nothing || occursin(r"/s\d\d_[a-z_]+\.(jl|md)$", ruta)) || continue
    cubre = taller === nothing ? sesiones_del_archivo(ruta) : MODULOS[parse(Int, taller[1])]
    # La última sesión también apunta al proyecto final (fila 0 del índice)
    alcance = NSESIONES in cubre ? union(cubre, Set([0])) : cubre
    permitidas = union((get(indice, s, Set{String}()) for s in alcance)...)
    for k in unique(claves_en(read(ruta, String)))
        haskey(etiquetas, k) || (reportar("$ruta: cita la clave inexistente $k"); continue)
        k in permitidas || reportar("$ruta: cita $k, que no está en el índice de las sesiones $(sort(collect(cubre)))")
    end
end

# --- 3. Enlaces relativos ---------------------------------------------------
anclas(ruta) = Set(m[1] for m in eachmatch(r"<a id=\"([^\"]+)\"></a>", read(ruta, String)))
for dir in ("Class_notes", REPASO, "docs", ".")
    for nombre in readdir(dir)
        endswith(nombre, ".md") || continue
        ruta = joinpath(dir, nombre)
        for m in eachmatch(r"\]\(([^)#\s]+\.md)(?:#([^)\s]+))?\)", read(ruta, String))
            destino = normpath(joinpath(dir, m[1]))
            if !isfile(destino)
                reportar("$ruta: enlace roto a $(m[1])")
            elseif m[2] !== nothing && !(m[2] in anclas(destino))
                reportar("$ruta: el ancla #$(m[2]) no existe en $(m[1])")
            end
        end
    end
end

# --- 4. Menciones a sesiones y semanas --------------------------------------
for dir in ("src", "Class_notes", REPASO)
    for nombre in readdir(dir)
        (endswith(nombre, ".jl") || endswith(nombre, ".md")) || continue
        ruta = joinpath(dir, nombre)
        for (i, linea) in enumerate(eachline(ruta))
            for m in eachmatch(r"\b[Ss]esi(?:ón|ones)\]?\s+(\d+)", linea)
                parse(Int, m[1]) <= NSESIONES || reportar("$ruta:$i: sesión inexistente: $(m.match)")
            end
            occursin(r"\b[Ss]emanas?\s+\d", linea) && reportar("$ruta:$i: menciona semanas del curso: $(strip(linea))")
        end
    end
end

if isempty(problemas)
    println("Sin problemas.")
else
    foreach(println, problemas)
    println("\n$(length(problemas)) problema(s).")
    exit(1)
end
