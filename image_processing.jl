# using Statistics
# using DataFrames
# using CSV
using Gtk
# using Tables
# using Polynomials
# using StatsPlots
# using HypothesisTests
# using Distributions
# using Symbolics
# using StatsBase
# using StatsModels
using Images
using ImageView
# using RandomNumbers
# using TikzPictures
# using CurveFit
# using TexTables
using GLMakie

##

function picking_files()

    seperation_of_copies = false

    show_plots = false

    solvent = String[]

    zoom = String[]

    producer_material = String[]

    iterate_files = []

    pulses = []

    down_dimension = 10

    length_of_colourscheme= 16000

    colours_init_blue = zeros(RGB{Float64}, down_dimension, length_of_colourscheme)

    colours_init_red = zeros(RGB{Float64}, down_dimension, length_of_colourscheme)

    colours_init_green = zeros(RGB{Float64}, down_dimension, length_of_colourscheme)

    # dataframe_path = raw"C:\Users\olive\Desktop\Uni\Summer Programs\DESY\Photos"

    dataframe_path = raw"C:\Users\liebeoli\Desktop\Functional Cellulose-lignin-coating on Porous Materials\Photos";

    files = open_dialog("Chose a file", GtkNullContainer(), String["*.jpg"], select_multiple=true)

    entries = readdir(dataframe_path, join=true)

    for i in range(1, length(entries))
        entries[i] = split(entries[i], raw"\\")[7][1:end-4]
    end

    for file_pick in files
        push!(solvent, split(file_pick, raw"_")[1][end-3:end])
        push!(producer_material, split(file_pick, raw"_")[2])
        push!(pulses, parse(Int64, split(file_pick, raw"_")[3][1:end-1]))
        push!(iterate_files, split(file_pick, raw"_")[1])
        # println(split(file_pick, raw"_")[5][1:end - 5])
        push!(zoom, split(file_pick, raw"_")[5][1:end - 5])
    end

    sort(pulses, rev = true)
    
    images = []

    for file in files
        push!(images, load(file))
    end

    global b = 1
    global r = 1
    global g = 1
    

    for image in images
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

    println(b - 1, " Is the amount of blue pixels, ", r - 1, " Is the amount of red pixels, ",  g - 1, " Is the amount of green pixels, with solvent: ",
                             solvent[1], ", at ", pulses[1], " number of pulses, with ", zoom[1], " zoom.")

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

    if seperation_of_copies == true
        if colours_blue != empty
            for i in range(1, b - 1)
                for l in range(i + 1, b - 1)
                    if float(green(colours_blue[1, i])) - 0.00001 < float(green(colours_blue[1, l])) < float(green(colours_blue[1, i])) + 0.00001 && 
                        float(red(colours_blue[1, i])) - 0.00001 < float(red(colours_blue[1, l])) < float(red(colours_blue[1, i])) + 0.00001 &&
                        float(blue(colours_blue[1, i])) - 0.00001 < float(blue(colours_blue[1, l])) < float(blue(colours_blue[1, i])) + 0.00001
                        for t in range(1, 10)
                            colours_blue[t, l] = RGB{Float64}(0.0, 0.0, 0.0)
                        end
                    end
                end
            end
            search = RGB{Float64}(0.0, 0.0, 0.0)
            indexArray_blue = findall(x -> x != search, colours_blue[1,:])
            global blue_seperator = length(indexArray_blue)
            println(blue_seperator, " Is the corrected colour-scheme for blue.")
            colour_seperated_blue = zeros(RGB{Float64}, Int(round(blue_seperator/10)), blue_seperator)
            w = 1
            for i in range(1, b - 1)
                if colours_blue[1, i] != RGB{Float64}(0.0, 0.0, 0.0)
                    for p in range(1, Int(round(blue_seperator/10)))
                        colour_seperated_blue[p, w] = colours_blue[1, i]
                    end
                    w += 1
                end
            end
            if show_plots == true
                imshow(colour_seperated_blue)
            end
        else
            println("There are no blue pixels in this photo")
        end
        if colours_red != empty
            for i in range(1, r - 1)
                for l in range(i + 1, r - 1)
                    if float(green(colours_red[1, i])) - 0.00001 < float(green(colours_red[1, l])) < float(green(colours_red[1, i])) + 0.00001 && 
                        float(red(colours_red[1, i])) - 0.00001 < float(red(colours_red[1, l])) < float(red(colours_red[1, i])) + 0.00001 &&
                        float(blue(colours_red[1, i])) - 0.00001 < float(blue(colours_red[1, l])) < float(blue(colours_red[1, i])) + 0.00001
                        for t in range(1, 10)
                            colours_red[t, l] = RGB{Float64}(0.0, 0.0, 0.0)
                        end
                    end
                end
            end
            search = RGB{Float64}(0.0, 0.0, 0.0)
            indexArray_red = findall(x -> x != search, colours_red[1,:])
            global red_seperator = length(indexArray_red)
            println(red_seperator, " Is the corrected colour-scheme for red.")
            colour_seperated_red = zeros(RGB{Float64}, Int(round(red_seperator/10)), red_seperator)
            q = 1
            for i in range(1, r - 1)
                if colours_red[1, i] != RGB{Float64}(0.0, 0.0, 0.0)
                    for p in range(1, Int(round(red_seperator/10)))
                        colour_seperated_red[p, q] = colours_red[1, i]
                    end
                    q += 1
                end
            end
            if show_plots == true
                imshow(colour_seperated_red)
            end
        else
            println("There are no red pixels in this photo")
        end
        if colours_green != empty
            for i in range(1, g - 1)
                for l in range(i + 1, g - 1)
                    if float(green(colours_green[1, i])) - 0.00001 < float(green(colours_green[1, l])) < float(green(colours_green[1, i])) + 0.00001 && 
                        float(red(colours_green[1, i])) - 0.00001 < float(red(colours_green[1, l])) < float(red(colours_green[1, i])) + 0.00001 &&
                        float(blue(colours_green[1, i])) - 0.00001 < float(blue(colours_green[1, l])) < float(blue(colours_green[1, i])) + 0.00001
                        for t in range(1, 10)
                            colours_green[t, l] = RGB{Float64}(0.0, 0.0, 0.0)
                        end
                    end
                end
            end
            search = RGB{Float64}(0.0, 0.0, 0.0)
            indexArray_green = findall(x -> x != search, colours_green[1,:])
            global green_seperator = length(indexArray_green)
            println(green_seperator, " Is the corrected colour-scheme for green.")
            colour_seperated_green = zeros(RGB{Float64}, Int(round(green_seperator/10)), green_seperator)
            e = 1
            for i in range(1, g - 1)
                if colours_green[1, i] != RGB{Float64}(0.0, 0.0, 0.0)
                    for p in range(1, Int(round(green_seperator/10)))
                        colour_seperated_green[p, e] = colours_green[1, i]
                    end
                    e += 1
                end
            end
            if show_plots == true
                imshow(colour_seperated_green)
            end
        else
            println("There are no green pixels in this photo.")
        end
    else
        if colours_green != empty && show_plots == true
            imshow(colours_green)
        end
        if colours_red != empty && show_plots == true
            imshow(colours_red)
        end
        if colours_blue != empty && show_plots == true
            imshow(colours_blue)
        end
    end
    return solvent, producer_material, pulses, entries, b, r, g
end


##


function picking_files_and_plotting()

    amount_of_photos = 4

    solvent_list = String[]

    pulse_list = Float32[]

    blue_pixels = Float32[]

    red_pixels = Float32[]

    green_pixels = Float32[]

    for i in range(1, amount_of_photos)
        solvent, producer_material, pulses, entries, blue, red, green = picking_files();
        push!(solvent_list, solvent[1])
        push!(pulse_list, pulses[1])
        push!(blue_pixels, blue[1])
        push!(red_pixels, red[1])
        push!(green_pixels, green[1])
    end


    fig1 = GLMakie.Figure()
    ax1 = GLMakie.Axis(fig1[1,1], xlabel = raw"Pulse number", ylabel = raw"Amount of Coloured Pixels in The Photo")
    GLMakie.scatter!(ax1, pulse_list, blue_pixels, markersize = 10, color = :blue, label = "Blue pixels")
    GLMakie.scatter!(ax1, pulse_list, red_pixels, markersize = 10, color = :red, label = "Red pixels")
    GLMakie.scatter!(ax1, pulse_list, green_pixels, markersize = 10, color = :green, label = "Green pixels")
    axislegend(ax1, position=:lt, labelsize = 30)

    display(fig1)
    return 0
end


picking_files_and_plotting();