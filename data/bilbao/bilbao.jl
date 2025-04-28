using HTTP, Gumbo, JSON

"""
    encode_form(data::Dict{String,String})
    # Arguments
    - `data::Dict{String,String}`: A dictionary containing form data to be encoded.
    # Returns
    - `String`: The encoded form data as a URL-encoded string.
"""

function encode_form(data::Dict{String,String})
    return join(["$(HTTP.escape(s))=$(HTTP.escape(v))" for (s, v) in data], "&")
end

"""
    This function retrieves the HTML document containing the table of generators for a given space group number from Bilbao Crystallographic Server.
    get_gen_doc(gnum::Integer)
    # Arguments
    - `gnum::Integer`: The space group number for which to retrieve the generators.
    # Returns
    - `Gumbo.HTMLDocument`: The parsed HTML document containing the generators.
"""
function get_gen_doc(gnum::Integer)::HTMLDocument
    url = "https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-getgen"

    form_data = Dict(
        "list" => "Standard/Default Setting",
        "what" => "gen",
        "gnum" => string(gnum)
    )

    encoded_body = encode_form(form_data)

    headers = [
        "Content-Type" => "application/x-www-form-urlencoded"
    ]

    response = HTTP.post(url, headers, encoded_body)

    if response.status == 200
        html_string = String(response.body)
        return parsehtml(html_string)
    else
        error("Request failed with status: $(response.status)")
    end
end


function parse_gen_doc(doc::HTMLDocument)::Dict{String,Any}
    d=Dict{String,Any}()
    d["name"]=parse_group_name(doc) # The name of the space group
    # The rows of the table
    rows=doc.root[2][6][1][1].children
    numgens = length(rows)-2
    generators=Vector{Dict{String,String}}(undef, numgens)
    d["generators"]=generators
    for i=1:numgens
        gen=Dict{String,String}()
        # The columns of the table
        cols=rows[i+2]
        gen["xyz"]=strip(cols[2][1].text) # The generator in XYZ format
        m=cols[3][1][1][1][2][1][1].text # The generator in 3x4 matrix format
        r=split(m, '\n') # rows of the 3x4 matrix
        io_a=IOBuffer()
        io_b=IOBuffer()
        write(io_a, '[')
        write(io_b, '[')
        for j=1:3
            c=split(r[j]) # split the row into columns
            join(io_a, c[1:3], " ") # join first three columns
            write(io_b, replace(c[4], "/" => "//")) # write the last column
            if j<3
                write(io_a, ";")
                write(io_b, ",")
            end
        end
        write(io_a, ']')
        write(io_b, ']')
        gen["a"]=String(take!(io_a))
        gen["b"]=String(take!(io_b))
        generators[i]=gen
    end
    d
end



function parse_group_name(doc::HTMLDocument)::String
    sub_chars = Dict('0'=>'₀', '1'=>'₁', '2'=>'₂', '3'=>'₃', '4'=>'₄', '5'=>'₅', '6'=>'₆', '7'=>'₇', '8'=>'₈', '9'=>'₉')
    header=doc.root[2][4].children
    io = IOBuffer()
    for e in header[2:end] # Skip the preamble
        if e isa HTMLText
            write(io, e.text)
        elseif e isa HTMLElement{:i}
            write(io, e[1].text)
        elseif e isa HTMLElement{:sub} 
            for c in e[1].text
                if haskey(sub_chars, c)
                    write(io, sub_chars[c])
                else
                    write(io, c)
                end
            end
        end
    end
    rstrip(String(take!(io)))
end

function get_generators3D()
    generators=Vector{Any}()
    for gnum in 1:230
        doc = get_gen_doc(gnum)
        d = parse_gen_doc(doc)
        println(d["name"])
        push!(generators, d)
        #sleep(1) # Sleep for 1 second to avoid overwhelming the server
    end
    # Save the generators to a JSON file
    script_dir = @__DIR__
    json_file = joinpath(script_dir, "generators3D.json")
    open(json_file, "w") do io
        JSON.print(io, generators, 4)
    end
end