function parseSummary(fname)
    file = open(fname)
    lines = readlines(file)
    p = Dict{String,Float64}()
    for line in lines
        els = split(line, " ")
        if els[1] == "parameterSummary"
            k = els[2]
            if startswith(k, "exrep1.") == false
                str = strip(SubString(els[3], 6), ';')
                if str != ""
                    v = parse(Float64, str)
                    p[k] = v
                end
            end
        end
    end
    return p
end
