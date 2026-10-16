# Descarga precios diarios de Yahoo Finance y los guarda en Data/precios.csv.
#
# Uso (desde la raíz del repositorio):
#     julia --project=Code Code/descargar_datos.jl
#
# Los notebooks leen el CSV, así que en clase no se depende de la API de Yahoo.
# Las fechas son fijas para que todos trabajen con los mismos datos.

using CSV, DataFrames, Dates, Downloads, JSON

const ACTIVOS = [
    ("^MXX",    "IPC (BMV)"),
    ("MXN=X",   "USD/MXN"),
    ("^GSPC",   "S&P 500"),
    ("BTC-USD", "Bitcoin"),
]
const INICIO = Date(2000, 1, 1)
const FIN    = Date(2026, 9, 30)

unix(d::Date) = round(Int, datetime2unix(DateTime(d)))

function yahoo(ticker; inicio = INICIO, fin = FIN)
    url = "https://query1.finance.yahoo.com/v8/finance/chart/" *
          replace(ticker, "^" => "%5E", "=" => "%3D") *
          "?period1=$(unix(inicio))&period2=$(unix(fin + Day(1)))&interval=1d"
    io = IOBuffer()
    Downloads.download(url, io; headers = ["User-Agent" => "Mozilla/5.0"])
    res = JSON.parse(String(take!(io)))["chart"]["result"][1]
    fechas = Date.(unix2datetime.(res["timestamp"]))
    cierre = res["indicators"]["adjclose"][1]["adjclose"]
    df = DataFrame(fecha = fechas, cierre = cierre)
    # Yahoo marca con `nothing` los días sin cotización
    filter!(r -> !isnothing(r.cierre) && inicio <= r.fecha <= fin, df)
    df.cierre = Float64.(df.cierre)
    unique!(df, :fecha)
end

# Al incluir este archivo desde otro script solo se define `yahoo`
if abspath(PROGRAM_FILE) == @__FILE__
    datos = mapreduce(vcat, ACTIVOS) do (ticker, nombre)
        df = yahoo(ticker)
        df.activo .= nombre
        println(rpad(nombre, 12), nrow(df), " días  (", first(df.fecha), " a ", last(df.fecha), ")")
        df
    end

    destino = joinpath(@__DIR__, "..", "Data", "precios.csv")
    mkpath(dirname(destino))
    CSV.write(destino, select(datos, :fecha, :activo, :cierre))
    println("Guardado en ", normpath(destino))
end
