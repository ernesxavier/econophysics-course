# Convierte un archivo fuente con marcadores `#%%` (código visible) y
# `#%% hide` (código oculto) en un cuaderno de Pluto, lo ejecuta sin interfaz y
# guarda las salidas (imágenes y HTML) para revisarlas. Uso, desde la raíz:
#   julia --project=Code src/construir_cuaderno.jl src/s03_colas_pesadas.jl Code/s03_colas_pesadas.jl /tmp/salidas_s03
# El destino debe quedar en Code/ para que Pkg.activate(@__DIR__) encuentre el entorno.
using Pluto

fuente, destino, salidas = ARGS
mkpath(salidas)

bloques = split(read(fuente, String), r"^#%%"m)[2:end]
celdas = map(bloques) do b
    primera, resto = split(b, "\n"; limit = 2)
    c = Pluto.Cell(rstrip(resto))
    c.code_folded = occursin("hide", primera)
    c
end
nb = Pluto.Notebook(celdas, destino)
Pluto.save_notebook(nb)
println("Cuaderno con $(length(celdas)) celdas guardado en $destino")

sesion = Pluto.ServerSession()
nb = Pluto.SessionActions.open(sesion, destino; run_async = false)
errores = 0
for (i, c) in enumerate(nb.cells)
    if c.errored
        global errores += 1
        msg = c.output.body isa Dict ? get(c.output.body, :msg, c.output.body) : c.output.body
        println("ERROR en celda $i: ", first(c.code, 80), "\n   → ", first(string(msg), 400))
    end
    m = c.output.mime
    nombre = joinpath(salidas, lpad(i, 2, '0'))
    if m isa MIME"image/png"
        write(nombre * ".png", c.output.body)
    elseif m isa MIME"image/svg+xml"
        write(nombre * ".svg", c.output.body)
    elseif m isa MIME"text/html"
        write(nombre * ".html", c.output.body)
    end
end
println("Celdas con error: $errores de $(length(nb.cells))")
Pluto.SessionActions.shutdown(sesion, nb; async = false)
