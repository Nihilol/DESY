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

    dataframe_path = raw"C:\Users\liebeoli\Desktop\Functional Cellulose-lignin-coating on Porous Materials\Photos";

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

    # imshow(images[1])

    for image in images
        for i in (1:500:size(image)[1])
            for l in (1:1:size(image)[2])
                println(float(red(image[i, l])), raw", " , float(green(image[i, l])), raw", " , float(blue(image[i, l])))
            end
        end
    end

    return solvent, producer_material, pulses, entries
end

##

picking_files();