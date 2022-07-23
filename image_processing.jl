using Statistics
using DataFrames
using GLMakie
using CSV
using Gtk
using Tables
using Polynomials
using StatsPlots
using HypothesisTests
using Distributions
using Symbolics
using StatsBase
using StatsModels
using AbstractPlotting
using ImageView
using RandomNumbers
using TikzPictures
using CurveFit
using TexTables
using Images

##

function picking_files()

    dataframe_path = raw"C:\Users\olive\Desktop\Uni\Summer Programs\DESY\Photos"

    # dataframe_path = raw"C:\Users\liebeoli\Desktop\Functional Cellulose-lignin-coating on Porous Materials\Photos";

    files = open_dialog("Chose a file", GtkNullContainer(), String["*.jpg"], select_multiple=true)

    entries = readdir(dataframe_path, join=true)

    for i in range(1, length(entries))
        entries[i] = split(entries[i], raw"\\")[7][1:end-4]
    end

    solvent = String[]

    producer_material = String[]

    iterate_files = []

    pulses = []

    for file_pick in files
        push!(solvent, split(file_pick, raw"_")[1][end-3:end])
        push!(producer_material, split(file_pick, raw"_")[2])
        push!(pulses, parse(Int64, split(file_pick, raw"_")[3][1:end-1]))
        push!(iterate_files, split(file_pick, raw"_")[1])

    end

    sort(pulses, rev = true)
    
    images = []

    for file in files
        push!(images, load(file))
    end

    down_dimension = 7

    # imshow(images[1])

    colours_init = zeros(RGB{Float64}, down_dimension, 250)

    for image in images
        global k = 1
        for i in (1:1:size(image)[1])
            for l in (1:1:size(image)[2])
                if float(red(image[i, l])) < 0.5 && float(green(image[i, l])) < 0.5 && float(blue(image[i, l])) > 0.7
                    # println(Int(255*float(red(image[i, l]))), raw", " , Int(255*float(green(image[i, l]))), raw", " , Int(255*float(blue(image[i, l]))))
                    for t in range(1, down_dimension)
                        colours_init[t, k] = image[i, l]
                    end
                    global k += 1
                end
            end
        end
    end

    sch = RGB{Float64}(0.0, 0.0, 0.0)

    indexArray = findall(x -> x == sch, colours_init)

    for i in indexArray
        println("There are ", i[2], " blue pixels in the photo")
        global g = i[2]
        break
    end

    colours = zeros(RGB{Float64}, down_dimension, g - 1)

    for i in range(1, g - 1)
        for t in range(1, down_dimension)
            colours[t, i] = colours_init[1, i]
        end
    end

    imshow(colours)

    return solvent, producer_material, pulses, entries
end


##

picking_files();