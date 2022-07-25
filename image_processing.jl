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

    # dataframe_path = raw"C:\Users\olive\Desktop\Uni\Summer Programs\DESY\Photos"

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

    down_dimension = 10

    length_of_colourscheme= 16000

    colours_init_blue = zeros(RGB{Float64}, down_dimension, length_of_colourscheme)

    colours_init_red = zeros(RGB{Float64}, down_dimension, length_of_colourscheme)

    colours_init_green = zeros(RGB{Float64}, down_dimension, length_of_colourscheme)

    

    for image in images
        global b = 1
        global r = 1
        global g = 1
        for i in (1:1:size(image)[1])
            for l in (1:1:size(image)[2])
                if float(red(image[i, l])) < 0.5 && float(green(image[i, l])) < 0.5 && float(blue(image[i, l])) > 0.6
                    for t in range(1, 10)
                        colours_init_blue[t, b] = image[i, l]
                    end
                    global b += 1
                end
                if float(red(image[i, l])) > 0.85 && float(green(image[i, l])) < 0.5 && float(blue(image[i, l])) < 0.5
                    for t in range(1, 10)
                        colours_init_red[t, r] = image[i, l]
                    end
                    global r += 1
                end
                if float(red(image[i, l])) < 0.5 && float(green(image[i, l])) > 0.8 && float(blue(image[i, l])) < 0.5
                    for t in range(1, 10)
                        colours_init_green[t, g] = image[i, l]
                    end
                    global g += 1
                end
            end 
        end
        global b = b
        global r = r
        global g = g
    end

    println(b - 1, " Is the amount of blue. ", r - 1, " Is the amount of red. ",  g - 1, " Is the amount of green")

    # sch = RGB{Float64}(0.0, 0.0, 0.0)

    # indexArray = findall(x -> x == sch, colours_init_blue)

    # for i in indexArray
    #     println("There are ", i[2], " blue pixels in the photo")
    #     global g = i[2]
    #     break
    # end

    colours_blue = zeros(RGB{Float64}, down_dimension, b - 1)
    colours_red = zeros(RGB{Float64}, down_dimension, r - 1)
    colours_green = zeros(RGB{Float64}, down_dimension, g - 1)

    for i in range(1, b - 1)
        for t in range(1, down_dimension)
            colours_blue[t, i] = colours_init_blue[1, i]
        end
    end
    for i in range(1, r - 1)
        for t in range(1, down_dimension)
            colours_red[t, i] = colours_init_red[1, i]
        end
    end
    for i in range(1, g - 1)
        for t in range(1, down_dimension)
            colours_green[t, i] = colours_init_green[1, i]
        end
    end

    if colours_blue != empty
        imshow(colours_blue)
    else
        println("There are no blue pixels in this photo")
    end
    if colours_red != empty
        imshow(colours_red)
    else
        println("There are no green pixels in this photo")
    end
    if colours_green != empty
        imshow(colours_green)
    else
        println("There are no green pixels in this photo")
    end
    return solvent, producer_material, pulses, entries
end


##

picking_files();