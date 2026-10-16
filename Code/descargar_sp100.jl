# Descarga precios diarios de las acciones del S&P 100 y los guarda en
# Data/sp100_precios.csv (fecha, ticker, cierre) y Data/sp100_sectores.csv.
#
# Uso (desde la raíz del repositorio):
#     julia --project=Code Code/descargar_sp100.jl
#
# Composición del índice al 21 de septiembre de 2026 (Wikipedia, con datos del
# ETF OEF de iShares). Se omite GOOG porque es otra clase de acción de GOOGL.
# Sectores según la clasificación GICS.
#
# Solo se guardan las acciones que cotizan durante toda la ventana, para que la
# matriz de rendimientos quede completa. Usar la composición actual introduce
# sesgo de supervivencia: las empresas que salieron del índice no aparecen.

include(joinpath(@__DIR__, "descargar_datos.jl"))   # función `yahoo`

const ACCIONES = [
    ("AAPL",  "Apple Inc.",                              "Tecnología"),
    ("ABBV",  "AbbVie",                                  "Salud"),
    ("ABT",   "Abbott Laboratories",                     "Salud"),
    ("ACN",   "Accenture",                               "Tecnología"),
    ("ADBE",  "Adobe Inc.",                              "Tecnología"),
    ("AMAT",  "Applied Materials",                       "Tecnología"),
    ("AMD",   "Advanced Micro Devices",                  "Tecnología"),
    ("AMGN",  "Amgen",                                   "Salud"),
    ("AMT",   "American Tower",                          "Bienes raíces"),
    ("AMZN",  "Amazon",                                  "Consumo discrecional"),
    ("ANET",  "Arista Networks",                         "Tecnología"),
    ("AVGO",  "Broadcom",                                "Tecnología"),
    ("AXP",   "American Express",                        "Finanzas"),
    ("BA",    "Boeing",                                  "Industria"),
    ("BAC",   "Bank of America",                         "Finanzas"),
    ("BKNG",  "Booking Holdings",                        "Consumo discrecional"),
    ("BLK",   "BlackRock",                               "Finanzas"),
    ("BMY",   "Bristol Myers Squibb",                    "Salud"),
    ("BNY",   "BNY Mellon",                              "Finanzas"),
    ("BRK-B", "Berkshire Hathaway (Class B)",            "Finanzas"),
    ("C",     "Citigroup",                               "Finanzas"),
    ("CAT",   "Caterpillar Inc.",                        "Industria"),
    ("CMCSA", "Comcast",                                 "Comunicaciones"),
    ("COF",   "Capital One",                             "Finanzas"),
    ("COP",   "ConocoPhillips",                          "Energía"),
    ("COST",  "Costco",                                  "Consumo básico"),
    ("CRM",   "Salesforce",                              "Tecnología"),
    ("CSCO",  "Cisco",                                   "Tecnología"),
    ("CVS",   "CVS Health",                              "Salud"),
    ("CVX",   "Chevron Corporation",                     "Energía"),
    ("DE",    "Deere & Company",                         "Industria"),
    ("DELL",  "Dell Technologies",                       "Tecnología"),
    ("DHR",   "Danaher Corporation",                     "Salud"),
    ("DIS",   "Walt Disney Company (The)",               "Comunicaciones"),
    ("DUK",   "Duke Energy",                             "Servicios públicos"),
    ("EMR",   "Emerson Electric",                        "Industria"),
    ("FDX",   "FedEx",                                   "Industria"),
    ("GD",    "General Dynamics",                        "Industria"),
    ("GE",    "GE Aerospace",                            "Industria"),
    ("GEV",   "GE Vernova",                              "Industria"),
    ("GILD",  "Gilead Sciences",                         "Salud"),
    ("GM",    "General Motors",                          "Consumo discrecional"),
    ("GOOGL", "Alphabet Inc. (Class A)",                 "Comunicaciones"),
    ("GS",    "Goldman Sachs",                           "Finanzas"),
    ("HD",    "Home Depot",                              "Consumo discrecional"),
    ("IBM",   "IBM",                                     "Tecnología"),
    ("INTC",  "Intel",                                   "Tecnología"),
    ("INTU",  "Intuit",                                  "Tecnología"),
    ("ISRG",  "Intuitive Surgical",                      "Salud"),
    ("JNJ",   "Johnson & Johnson",                       "Salud"),
    ("JPM",   "JPMorgan Chase",                          "Finanzas"),
    ("KO",    "Coca-Cola Company (The)",                 "Consumo básico"),
    ("LIN",   "Linde plc",                               "Materiales"),
    ("LLY",   "Eli Lilly and Company",                   "Salud"),
    ("LMT",   "Lockheed Martin",                         "Industria"),
    ("LOW",   "Lowe's",                                  "Consumo discrecional"),
    ("LRCX",  "Lam Research",                            "Tecnología"),
    ("MA",    "Mastercard",                              "Finanzas"),
    ("MCD",   "McDonald's",                              "Consumo discrecional"),
    ("MDLZ",  "Mondelēz International",                  "Consumo básico"),
    ("MDT",   "Medtronic",                               "Salud"),
    ("META",  "Meta Platforms",                          "Comunicaciones"),
    ("MMM",   "3M",                                      "Industria"),
    ("MO",    "Altria",                                  "Consumo básico"),
    ("MRK",   "Merck & Co.",                             "Salud"),
    ("MS",    "Morgan Stanley",                          "Finanzas"),
    ("MSFT",  "Microsoft",                               "Tecnología"),
    ("MU",    "Micron Technology",                       "Tecnología"),
    ("NEE",   "NextEra Energy",                          "Servicios públicos"),
    ("NFLX",  "Netflix, Inc.",                           "Comunicaciones"),
    ("NOW",   "ServiceNow",                              "Tecnología"),
    ("NVDA",  "Nvidia",                                  "Tecnología"),
    ("ORCL",  "Oracle Corporation",                      "Tecnología"),
    ("PANW",  "Palo Alto Networks",                      "Tecnología"),
    ("PEP",   "PepsiCo",                                 "Consumo básico"),
    ("PFE",   "Pfizer",                                  "Salud"),
    ("PG",    "Procter & Gamble",                        "Consumo básico"),
    ("PLTR",  "Palantir Technologies",                   "Tecnología"),
    ("PM",    "Philip Morris International",             "Consumo básico"),
    ("QCOM",  "Qualcomm",                                "Tecnología"),
    ("RTX",   "RTX Corporation",                         "Industria"),
    ("SBUX",  "Starbucks",                               "Consumo discrecional"),
    ("SCHW",  "Charles Schwab Corporation",              "Finanzas"),
    ("SNDK",  "Sandisk",                                 "Tecnología"),
    ("SO",    "Southern Company",                        "Servicios públicos"),
    ("T",     "AT&T",                                    "Comunicaciones"),
    ("TMO",   "Thermo Fisher Scientific",                "Salud"),
    ("TMUS",  "T-Mobile US",                             "Comunicaciones"),
    ("TSLA",  "Tesla, Inc.",                             "Consumo discrecional"),
    ("TXN",   "Texas Instruments",                       "Tecnología"),
    ("UBER",  "Uber",                                    "Industria"),
    ("UNH",   "UnitedHealth Group",                      "Salud"),
    ("UNP",   "Union Pacific Corporation",               "Industria"),
    ("UPS",   "United Parcel Service",                   "Industria"),
    ("USB",   "U.S. Bancorp",                            "Finanzas"),
    ("V",     "Visa Inc.",                               "Finanzas"),
    ("VZ",    "Verizon",                                 "Comunicaciones"),
    ("WFC",   "Wells Fargo",                             "Finanzas"),
    ("WMT",   "Walmart",                                 "Consumo básico"),
    ("XOM",   "ExxonMobil",                              "Energía"),
]
const INICIO_SP100 = Date(2015, 1, 1)

# Errores de ajuste de Yahoo: en la escisión de Fortive el precio ajustado de
# Danaher sube 61% en un día. Se reescalan los precios anteriores para que el
# rendimiento de ese día sea cero (el real fue de alrededor de ±1%).
const CORRECCIONES = [("DHR", Date(2016, 7, 5))]

function corregir!(df, ticker)
    for (t, dia) in CORRECCIONES
        t == ticker || continue
        i = findfirst(==(dia), df.fecha)
        df.cierre[1:i-1] .*= df.cierre[i] / df.cierre[i-1]
    end
    df
end

precios = DataFrame[]
omitidas = String[]
for (ticker, nombre, sector) in ACCIONES
    df = yahoo(ticker; inicio = INICIO_SP100)
    # Cubre la ventana si cotiza desde la primera semana
    if first(df.fecha) > INICIO_SP100 + Day(7)
        println(rpad(ticker, 7), "omitida: empieza el ", first(df.fecha))
        push!(omitidas, ticker)
        continue
    end
    corregir!(df, ticker)
    df.ticker .= ticker
    push!(precios, df)
    sleep(0.3)   # no saturar la API
end
datos = reduce(vcat, precios)

# Solo días en que cotizan todas las acciones
dias = combine(groupby(datos, :fecha), nrow => :n)
completos = Set(dias.fecha[dias.n .== length(precios)])
println("Días con todas las acciones: ", length(completos), " de ", nrow(dias))
filter!(r -> r.fecha in completos, datos)

sectores = DataFrame([(ticker = t, nombre = n, sector = s)
                      for (t, n, s) in ACCIONES if !(t in omitidas)])

destino = joinpath(@__DIR__, "..", "Data")
CSV.write(joinpath(destino, "sp100_precios.csv"), select(datos, :fecha, :ticker, :cierre))
CSV.write(joinpath(destino, "sp100_sectores.csv"), sectores)
println(nrow(sectores), " acciones guardadas; omitidas: ", join(omitidas, ", "))
