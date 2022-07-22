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

function picking_files()

    dataframe_path = raw"C:\Users\olive\Desktop\Zipped"

    files = open_dialog("Chose a file", GtkNullContainer(), String["*.csv"], select_multiple=true)

    if occursin("nooil", files[1]) || occursin("No Oil", files[1])
        println("It does contain oil")
        dataframe_path = raw"C:\Users\olive\Desktop\Zipped\No Oil"
    end

    entries = readdir(dataframe_path, join=true)

    conc_files = []

    concentrations = String[]

    time_vector = Vector{Float64}[]

    for file_pick in files
        push!(conc_files, split(file_pick, raw"_")[1])
        push!(concentrations, last(split(split.(file_pick, raw"_")[1], raw"\\")))

    end

    full_mean = Vector{Float64}[]

    full_std = Vector{Float64}[]

    for check_file in conc_files
        mean_data_simple = Float64[]
        time_vector_simple = Float64[]
        std_data_simple = Float64[]
        for file_pick2 in entries
            i = 1
            if occursin(check_file, file_pick2) &! occursin("radius", file_pick2)
                if occursin("Oil", dataframe_path)
                    file_path = string(raw"C:\Users\olive\Desktop\Zipped\\", concentrations[i])
                    start_radius = Tables.matrix(CSV.File(string(file_path, "radiusZipped.csv")))[1]
                else
                    start_radius= Tables.matrix(CSV.File(string(check_file, "radiusZipped.csv")))[1]
                end
                # println(start_radius)
                matrix_of_data = Tables.matrix(DataFrame(CSV.File(file_pick2, select=[1])))
                push!(mean_data_simple, mean(matrix_of_data/start_radius))
                minute = split(file_pick2, raw"_")[3][7:end-4]
                push!(time_vector_simple, parse(Float64, minute))
                push!(std_data_simple, std(matrix_of_data/start_radius))
            end
        end
        push!(full_std, std_data_simple)
        push!(full_mean, mean_data_simple)
        push!(time_vector, time_vector_simple)
    end


    println("You have chosen ", length(concentrations), " different concentrations")
    println("The concentrations are", concentrations)
    return full_mean, time_vector, concentrations, full_std
    empty!(full_mean)
    empty!(time_vector)
    println(time_vector)
    empty!(full_std)
    empty!(concentrations)
end

##

function find_radius()
    path_to_zipped = raw"C:\Users\olive\Desktop\Zipped"
    file = open_dialog("Chose a file", GtkNullContainer(), String["*.csv"])
    data = Tables.matrix(DataFrame(CSV.File(file)))
    x1 = data[4]
    x2 = data[5]
    x3 = data[6]
    y1 = data[7]
    y2 = data[8]
    y3 = data[9]
    temp = x2^2 + y2^2
    bc = (x1^2 + y1^2 - temp) / 2
    cd = (temp - x3^2 - y3^2) / 2
    det = (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2)
    if (abs(det) < 1.0e-10)
        return None
    end
    # Center of circle
    cx = (bc*(y2 - y3) - cd*(y1 - y2)) / det
    cy = ((x1 - x2) * cd - (x2 - x3) * bc) / det
    radius = [((cx - x1)^2 + (cy - y1)^2)^(.5)]
    full_file = string(path_to_zipped, "\\", split(file, raw"\\")[7][1:end-4], "Zipped", ".csv")
    CSV.write(full_file, DataFrame([radius], [:radius]))
    return radius
end

##


function plot_data(full_mean::Vector, full_std::Vector, time_vector::Vector)
    colours = [:crimson, :dodgerblue, :slateblue1, :sienna1, :orchid1]
    fig1 = Figure()
    fig2 = Figure()
    title_of_saved_figure = ""
    conc = Float64[305, 255, 195, 98, 19.8]
    slopes = Float64[]
    ax1 = Axis(fig1[1, 1], xlabel = raw"Time in minutes",
     ylabel = raw"Normalised Length of interfacial region", xlabelsize = 30, ylabelsize = 27)#, limits = (minimum(time_vector[2]) - 1.0
     #, maximum(time_vector[1]) + 3.0, 0, maximum(full_mean[2]) + 1/2*maximum(full_mean[2])))
    ax2 = Axis(fig2[1,1], xlabel = raw"Concentration in mM",
     ylabel = raw"Normalised Linear slope of concentration growth", xlabelsize = 30, ylabelsize = 21)
    for i in eachindex(full_mean)
        if occursin("19comma7mM", concentrations[i])
            concentrations[i] = "20mM"
        end
        linfit = Polynomials.fit(time_vector[i], full_mean[i], 1)
        push!(slopes, (tryparse(Float64, split(split((String(string(linfit))), raw" ")[3], raw"*")[1])))
        error_on_time = fill(0.3, length(time_vector[i]))
        error_on_mean = (full_std[i]/sqrt(length(full_mean[i])))
        errorbars!(ax1, time_vector[i], full_mean[i], error_on_time, color = colours[i], whiskerwidth = 4,
        direction = :x)
        errorbars!(ax1, time_vector[i], full_mean[i], error_on_mean, color = colours[i], whiskerwidth = 4,
        direction = :y)
        GLMakie.scatter!(ax1, time_vector[i], full_mean[i], markersize = 4, color=colours[i],
            label=concentrations[i])
        title_of_saved_figure = string(title_of_saved_figure, concentrations[i])
        empty!(error_on_time)
        empty!(error_on_mean)
    end
    GLMakie.scatter!(ax2, conc, slopes, markersize = 10, color = :black)
    println(slopes)
    # placement = raw"C:\Users\olive\Desktop\Uni\Bachelor\6. Semester\Bachelor\\"
    # full_title = string(placement, title_of_saved_figure, raw".png")
    # vlines!(ax1, [41,53], ymin = [0.1, 0.1], ymax = [0.7, 0.7], color = :black, linestyle = :dash)
    # hlines!(ax1, [50,40], xmin = [41, 41.0], xmax = [53, 53.0], color = :black, linestyle = :dash)
    # GLMakie.lines!(ax1, [41, 53], [0.0028, 0.0028], linestyle = :dash, color = :black)
    # GLMakie.lines!(ax1, [41, 53], [0.019, 0.019], linestyle = :dash, color = :black)
    axislegend(ax1, position=:lt, labelsize = 30, fontsize = 50)
    #Makie.save(full_title, fig1)
    #Makie.save(raw"C:\Users\olive\Desktop\Uni\Bachelor\6. Semester\Bachelor\slopes.png", fig2)
    display(fig2)
    # display(fig1)
end

##

full_mean, time_vector, concentrations, full_std = picking_files();

##

plot_data(full_mean, full_std, time_vector)

function zip_files()
    dataframe_path = raw"C:\Users\olive\Desktop\Non-zipped"

    files = open_dialog("Chose a file", GtkNullContainer(), String["*.csv"], select_multiple=true)

    if occursin("No Oil", files[1])
        println("It does contain oil")
        dataframe_path = raw"C:\Users\olive\Desktop\Non-zipped\No Oil"
    end

    entries = readdir(dataframe_path, join=true)

    conc_files = String[]

    for file_split in files
        splitting = split(file_split, raw"dot")
        push!(conc_files, splitting[1])
    end

    data_vector_simple = Float64[]

    for check_file in conc_files
        data_vector_simple = Float64[]
        for file in entries
            if occursin(check_file, file)
                loaded_data = vec(Tables.matrix(DataFrame(CSV.File(file, select=[2]))))
                data_vector_simple = [data_vector_simple; loaded_data]   # Might need to add the absolute value here
                path_to_leftover_file = string(raw"C:\Users\olive\Desktop\Non-zipped\Leftovers",
                 raw"\\", split(file, raw"\\")[6])
                cp(file, path_to_leftover_file, force = true)
                rm(file)
            end
        end
        zip_path = raw"C:\Users\olive\Desktop\Zipped\\"
        full_file = string(zip_path, split(check_file, raw"\\")[6], ".csv")
        if occursin("No Oil", dataframe_path)
            zip_path = raw"C:\Users\olive\Desktop\Zipped\No Oil\\"
            full_file = string(zip_path, split(check_file, raw"\\")[7], ".csv")
        end
        println(full_file)
        println("This file has been added to Zipped files: ", full_file)
        CSV.write(full_file, DataFrame([data_vector_simple], [:lengths]))
    end
    empty!(entries)
    empty!(conc_files)
    empty!(files)
end


##

zip_files();

##


find_radius()


function angles_distribution()
    file = open_dialog("Chose a file", GtkNullContainer(), String["*.csv"])
    vector_of_data = Tables.matrix(DataFrame(CSV.File(file, select=[2])))
    title_of_plot = string(split(file, raw"\\")[7][1:5], " - Parallel Flow Orientation")
    angles = Float64[]
    ranges = 1:2:length(vector_of_data) - 1
    for i in eachindex(vector_of_data)
        if vector_of_data[i] < 0
            vector_of_data[i] = 180 - abs(vector_of_data[i]);
        end
    end
    for i in ranges 
        if vector_of_data[i] > vector_of_data[i + 1]
            θ = vector_of_data[i]
            ϕ = vector_of_data[i + 1]
        else
            θ = vector_of_data[i + 1]
            ϕ = vector_of_data[i]
        end
        angle = abs(θ - ϕ)
        if angle > 90
            angle = 180 - angle
        end
        push!(angles, angle)
    end
    sample_size = length(angles)
    number_of_bins = Int64(ceil(log2(sample_size) + 1))
    xlabel_fig = string("Angle between fiber and micro-groove \n Sample size: ", sample_size)
    fig_hist = Figure()
    axhist = Axis(fig_hist[1,1], xlabel = xlabel_fig, ylabel = "Normalised value", title = title_of_plot,
        xlabelsize = 30, ylabelsize = 40)
    hist!(axhist, angles, bins = number_of_bins, normalization = :pdf, bar_labels = :values, color = :values,
     strokewidth = 0.5, stroke_color = (:black, 0.5), label_size = 10)
    # title_of_saved_figure = string(raw"C:\Users\olive\Desktop\Uni\Bachelor\6. Semester\Bachelor\\", title_of_plot, ".png")
    # Makie.save(title_of_saved_figure, fig_hist)
    display(fig_hist)
    # uniform_fit = StatsPlots.fit(Uniform, angles)
    # StatsPlots.qqplot(angles, uniform_fit)
    # ExactOneSampleKSTest(angles, uniform_fit)
end

angles_distribution() 