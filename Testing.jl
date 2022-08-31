using GLMakie: Figure
using GLMakie


f = Figure()

GLMakie.Axis(f[1, 1])

labels = [raw"Hello", raw"Yo", raw"What", raw"Up"]

markersizes = [5, 10, 15, 20]
colors = [:red, :green, :blue, :orange]

for ms in markersizes, color in colors
    scatter!(randn(5, 2), markersize = ms, color = color)
end

group_size = [MarkerElement(marker = :circle, color = colour,
    strokecolor = :transparent,
    markersize = 40) for colour in colors]

legends = [Legend(f,
    [group_size],
    [string.(labels)],
    ["Greeting"])]

f[1, 2] = legends[1]


# for l in legends[4:6]
#     l.orientation = :horizontal
#     l.tellheight = true
#     l.tellwidth = false
# end

# legends[2].titleposition = :left
# legends[5].titleposition = :left

# legends[3].nbanks = 2
# legends[5].nbanks = 2
# legends[6].nbanks = 2

f