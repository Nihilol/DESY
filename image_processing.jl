# using Statistics
using DataFrames
using CSV
using Gtk
using Tables
# using Polynomials
# using StatsPlots
# using HypothesisTests
# using Distributions
using Symbolics
# using StatsBase
# using StatsModels
using Images
using ImageView
# using RandomNumbers
# using TikzPictures
# using CurveFit
# using TexTables
using GLMakie
using Flux
using HCubature

##

pick_images = false

density_calculations = false

pick_csvs = true
## 


# dataframe_path = raw"C:\Users\olive\Desktop\Uni\Summer Programs\DESY\Photos"

dataframe_path = raw"C:\Users\liebeoli\Desktop\Functional Cellulose-lignin-coating on Porous Materials\Photos"


function picking_images()

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
    return solvent, producer_material, pulses, entries, b, r, g, zoom
end


##


function picking_images_and_plotting()

    amount_of_photos = 3

    solvent_list = String[]

    pulse_list = Float32[]

    blue_pixels = Float32[]

    red_pixels = Float32[]

    green_pixels = Float32[]

    zoom_list = String[]

    for i in range(1, amount_of_photos)
        solvent, producer_material, pulses, entries, blue, red, green, zoom = picking_images();
        push!(solvent_list, solvent[1])
        push!(pulse_list, pulses[1])
        push!(blue_pixels, blue[1])
        push!(red_pixels, red[1])
        push!(green_pixels, green[1])
        push!(zoom_list, zoom[1])
    end


    fig1 = GLMakie.Figure()
    ax1 = GLMakie.Axis(fig1[1,1], xlabel = raw"Pulse number", ylabel = raw"Amount of Coloured Pixels in The Photo")
    GLMakie.scatter!(ax1, pulse_list, blue_pixels, markersize = 10, color = :blue, label = "Blue pixels")
    GLMakie.scatter!(ax1, pulse_list, red_pixels, markersize = 10, color = :red, label = "Red pixels")
    GLMakie.scatter!(ax1, pulse_list, green_pixels, markersize = 10, color = :green, label = "Green pixels")
    axislegend(ax1, position=:lt, labelsize = 30)

    title_of_plot = string(dataframe_path, raw"\Graphs\\", solvent_list[1], raw"_", zoom_list[1], raw"s", raw".png")

    println(title_of_plot)

    Makie.save(title_of_plot, fig1)

    display(fig1)
    return 0
end


##

if pick_images == true
    picking_images_and_plotting();
end

if density_calculations == true
    function calc_dens(n, σ_x, σ_y, A, x_0, y_0)
        density = 1
        mass = 600
        function peaks()
            aa = LinRange(-15, 15, n)
            bb = LinRange(-15, 15, n)
            g = [A * exp.(-(((a.-x_0).^2)/(2*σ_x^2) + ((b.-y_0).^2)/(2*σ_y^2))) for a in aa, b in bb]
            return (aa, bb, g)
        end
        f(x, y) = A*exp(-((x-x_0)^2/(2*σ_x^2) + (y-y_0)^2/(2*σ_y^2)))
        f(v) = f(v...)
        a0, b0 = -10, 10
        a1, b1 = -10, 10
        area = hcubature(f, (a0, a1), (b0, b1))[1]
        if mass/density - 1 < area < mass/density + 1
            a, b, c = peaks()
            fig1 = Figure()
            ax1 = GLMakie.Axis3(fig1[1, 1])
            hm = GLMakie.surface!(ax1, a, b, c)
            Colorbar(fig1[1, 2], hm, height=Relative(0.5))
            display(fig1)
        end
        return σ_x, σ_y, A
    end
    for i in (0.1:2)
        for j in (1:100)
            calc_dens(49, i, i, j, 0, 0);
        end
    end
end

function picking_csv_and_plotting()
    data_frame_path = raw"C:\Users\liebeoli\Desktop\Functional Cellulose-lignin-coating on Porous Materials\SE CSV Files"

    files = open_dialog("Chose a file", GtkNullContainer(), String["*.csv"], select_multiple=true)

    entries = readdir(data_frame_path, join = true)

    pulses = []

    dates = []

    coating = []

    total_files = []

    place = []

    for file_pick in files
        push!(pulses, split(split(file_pick, raw"_")[1], raw"\\")[end])
        push!(dates, string(split(split(file_pick, raw"_")[4], raw"\\")[1], raw"_",  split(split(file_pick, raw"_")[5], raw"\\")[1][1:end-4]))
        push!(coating, split(split(file_pick, raw"_")[2], raw"\\")[1])
    end
    for file in entries, i in eachindex(pulses)
        if occursin(pulses[i], file) && occursin(dates[i], file) && occursin(coating[i], file)
            push!(total_files, file)
        end
    end

    global figure = Figure()

    global axis = GLMakie.Axis(figure[1, 1])

    for file in total_files

        Dataframe = DataFrame(CSV.File(file))

        df = Dataframe[:, 1]

        Wavelength = Float32[]

        Intensity = Float32[]

        place = split(split(file, raw"_")[3], raw"\\")[1]

        pulse = split(split(file, raw"_")[1], raw"\\")[end]

        for i in eachindex(df)
            if i == 1
                i = 2
            end
            push!(Wavelength, parse(Float64, df[i][1:10]))
            push!(Intensity, parse(Float64, df[i][12:18]))
        end

        colours = [:crimson, :dodgerblue, :slateblue1, :sienna1, :orchid1]

        label = string(place, raw" and ", pulse)

        GLMakie.scatterlines!(axis, Wavelength, Intensity, markersize = 3, label = label)

    end

    axislegend(axis, position=:rt, labelsize = 10)

    display(figure)
end

if pick_csvs == true
    picking_csv_and_plotting();
end


