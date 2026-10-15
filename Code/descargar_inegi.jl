# Descarga los Principales resultados por localidad (ITER) de los Censos de
# Población y Vivienda 2010 y 2020 (INEGI) y genera:
#
#   Data/localidades_2020.csv     población por localidad en 2020 (≥ 100 hab.)
#   Data/municipios_2010_2020.csv población por municipio en 2010 y 2020
#
# Uso (desde la raíz del repositorio):
#     julia --project=Code Code/descargar_inegi.jl
#
# Se excluyen las filas de totales nacional y estatal, y los agrupados de
# localidades de una y dos viviendas (claves 9998 y 9999). La fila con
# localidad 0 de cada municipio es el total municipal. En el archivo de
# municipios sólo quedan los que existían en ambos censos (13 municipios se
# crearon entre 2010 y 2020).

using CSV, DataFrames, Downloads, p7zip_jll

const URL_2020 = "https://www.inegi.org.mx/contenidos/programas/ccpv/2020/datosabiertos/iter/iter_00_cpv2020_csv.zip"
const URL_2010 = "https://www.inegi.org.mx/contenidos/programas/ccpv/2010/datosabiertos/iter_nal_2010_csv.zip"
const MINIMO = 100
const DATA = joinpath(@__DIR__, "..", "Data")

"Descarga y descomprime un ITER; devuelve la tabla con las columnas básicas."
function leer_iter(url)
    tmp = mktempdir()
    zip = joinpath(tmp, "iter.zip")
    println("Descargando ", basename(url), "...")
    Downloads.download(url, zip; headers = ["User-Agent" => "Mozilla/5.0"])
    run(pipeline(`$(p7zip()) x $zip -o$tmp -y`; stdout = devnull))
    csvs = [joinpath(d, f) for (d, _, fs) in walkdir(tmp) for f in fs
            if endswith(lowercase(f), ".csv") && filesize(joinpath(d, f)) > 10^7]
    t = CSV.read(only(csvs), DataFrame; normalizenames = true)
    rename!(lowercase, t)
    t = select(t, :entidad, :nom_ent, :mun, :nom_mun, :loc, :nom_loc, :pobtot)
    t.pobtot = [x isa Number ? Int(x) : tryparse(Int, string(x)) for x in t.pobtot]
    filter!(r -> r.entidad != 0 && r.mun != 0 && !isnothing(r.pobtot), t)
end

iter2020 = leer_iter(URL_2020)
iter2010 = leer_iter(URL_2010)

# Localidades 2020
loc = filter(r -> r.loc ∉ (0, 9998, 9999), iter2020)
println("Localidades 2020: ", nrow(loc), ", población total: ", sum(loc.pobtot))
filter!(:pobtot => >=(MINIMO), loc)
sort!(loc, :pobtot; rev = true)
loc = select(loc, :entidad => :entidad, :nom_ent => :nombre_entidad, :mun => :municipio,
             :nom_mun => :nombre_municipio, :loc => :localidad, :nom_loc => :nombre_localidad,
             :pobtot => :poblacion)
CSV.write(joinpath(DATA, "localidades_2020.csv"), loc)
println("Guardadas ", nrow(loc), " localidades con al menos $MINIMO habitantes")

# Municipios 2010 y 2020
mun(t) = select(filter(:loc => ==(0), t), :entidad, :mun, :nom_ent, :nom_mun, :pobtot)
m10, m20 = mun(iter2010), mun(iter2020)
println("Municipios 2010: ", nrow(m10), " (", sum(m10.pobtot), " hab.); 2020: ",
        nrow(m20), " (", sum(m20.pobtot), " hab.)")
m = innerjoin(select(m10, :entidad, :mun, :pobtot => :poblacion_2010),
              select(m20, :entidad, :mun, :nom_ent => :nombre_entidad,
                     :nom_mun => :nombre_municipio, :pobtot => :poblacion_2020),
              on = [:entidad, :mun])
m = select(m, :entidad, :nombre_entidad, :mun => :municipio, :nombre_municipio,
           :poblacion_2010, :poblacion_2020)
sort!(m, [:entidad, :municipio])
CSV.write(joinpath(DATA, "municipios_2010_2020.csv"), m)
println("Guardados ", nrow(m), " municipios presentes en ambos censos")
