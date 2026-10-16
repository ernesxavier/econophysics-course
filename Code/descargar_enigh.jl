# Descarga la Encuesta Nacional de Ingresos y Gastos de los Hogares (ENIGH)
# 2024 de INEGI y genera:
#
#   Data/enigh2024_hogares.csv   un renglón por hogar, con su factor de
#                                expansión y sus ingresos mensuales
#
# Uso (desde la raíz del repositorio):
#     julia --project=Code Code/descargar_enigh.jl
#
# Se usa la tabla `concentradohogar` (nueva serie). INEGI reporta los ingresos
# del trimestre; aquí se dividen entre 3 para tenerlos por mes, en pesos de
# 2024. El factor de expansión indica cuántos hogares del país representa cada
# hogar de la muestra: todas las estadísticas deben ponderarse con él.

using CSV, DataFrames, Downloads, p7zip_jll

const URL = "https://www.inegi.org.mx/contenidos/programas/enigh/nc/2024/datosabiertos/conjunto_de_datos_enigh2024_ns_csv.zip"
const TABLA = "conjunto_de_datos_concentradohogar_enigh2024_ns.csv"
const DATA = joinpath(@__DIR__, "..", "Data")

tmp = mktempdir()
zip = joinpath(tmp, "enigh.zip")
println("Descargando ", basename(URL), " (unos 100 MB)...")
Downloads.download(URL, zip; headers = ["User-Agent" => "Mozilla/5.0"])
run(pipeline(`$(p7zip()) e $zip -o$tmp $TABLA -r -y`; stdout = devnull))

t = CSV.read(joinpath(tmp, TABLA), DataFrame;
             select = [:ubica_geo, :tam_loc, :factor, :tot_integ, :ing_cor, :ingtrab, :transfer],
             types = Dict(:ubica_geo => String))
mensual(x) = round.(x ./ 3; digits = 2)
hogares = DataFrame(entidad = parse.(Int, first.(t.ubica_geo, 2)),
                    tam_loc = t.tam_loc,                 # 1: 100 mil hab. o más … 4: menos de 2,500
                    factor = t.factor,
                    integrantes = t.tot_integ,
                    ingreso = mensual(t.ing_cor),        # ingreso corriente
                    ingreso_trabajo = mensual(t.ingtrab),
                    transferencias = mensual(t.transfer))

CSV.write(joinpath(DATA, "enigh2024_hogares.csv"), hogares)
println("Guardados ", nrow(hogares), " hogares, que representan ",
        sum(hogares.factor), " hogares del país")
