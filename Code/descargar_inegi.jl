# Descarga la población por localidad del Censo de Población y Vivienda 2020
# (INEGI, Principales resultados por localidad, ITER) y la guarda en
# Data/localidades_2020.csv.
#
# Uso (desde la raíz del repositorio):
#     julia --project=Code Code/descargar_inegi.jl
#
# Se conservan las localidades con al menos 100 habitantes. Se excluyen las
# filas de totales (nacional, estatal y municipal) y los agrupados de
# localidades de una y dos viviendas (claves 9998 y 9999).

using CSV, DataFrames, Downloads, p7zip_jll

const URL = "https://www.inegi.org.mx/contenidos/programas/ccpv/2020/datosabiertos/iter/iter_00_cpv2020_csv.zip"
const MINIMO = 100

tmp = mktempdir()
zip = joinpath(tmp, "iter.zip")
println("Descargando ITER 2020 (≈ 37 MB)...")
Downloads.download(URL, zip; headers = ["User-Agent" => "Mozilla/5.0"])
run(pipeline(`$(p7zip()) x $zip -o$tmp -y`; stdout = devnull))
archivo = only(filter(endswith("CSV20.csv"),
                      readdir(joinpath(tmp, "iter_00_cpv2020", "conjunto_de_datos"); join = true)))

cols = [:ENTIDAD, :NOM_ENT, :MUN, :NOM_MUN, :LOC, :NOM_LOC, :POBTOT]
iter = CSV.read(archivo, DataFrame; select = cols, types = Dict(:POBTOT => String))
filter!(r -> r.ENTIDAD != 0 && r.MUN != 0 && r.LOC ∉ (0, 9998, 9999), iter)
iter.POBTOT = parse.(Int, iter.POBTOT)
println("Localidades: ", nrow(iter), ", población total: ", sum(iter.POBTOT))

filter!(:POBTOT => >=(MINIMO), iter)
sort!(iter, :POBTOT; rev = true)
rename!(iter, [:entidad, :nombre_entidad, :municipio, :nombre_municipio,
               :localidad, :nombre_localidad, :poblacion])

destino = joinpath(@__DIR__, "..", "Data", "localidades_2020.csv")
CSV.write(destino, iter)
println("Guardadas ", nrow(iter), " localidades con al menos $MINIMO habitantes en ",
        normpath(destino))
